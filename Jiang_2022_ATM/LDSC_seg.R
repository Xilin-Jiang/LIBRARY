##################################################################################
# LDSC_seg  for cell-type enrichment (Finucane et al. 2018 Nature Genetics)
##################################################################################

#################################################################################
# step 1 -------- install LDSC from https://github.com/bulik/ldsc
#################################################################################

#################################################################################
# step 3 -------- preparing summary statistics data
#################################################################################
# need to run analysis with both -assoc for plink and -- logistic; the latter doesn't have A1 A2 column, which is required for sum stats.
library(dplyr)
simple_assoc <- read.table("./LDSC_seg/250.2_topic7.assoc", header = T)
logistic_assoc <- read.table("./LDSC_seg/250.2_topic7.assoc.logistic", header = T)

logistic_a1_a2 <- logistic_assoc %>%
  left_join(select(simple_assoc, SNP, A1, A2), by = c("SNP", "A1"))
write.table(logistic_a1_a2, paste0("./LDSC_seg/250.2_topic7.logistic.a2"), sep="\t", row.names = FALSE ,quote = F)
# ----- now use the function in .sh file to map it

# both stratified LDSC and simple LDSC is performed in the .sh script

#################################################################################
# step 4 ---------- extract heritability estimation from ldsc results
#################################################################################
DIR <- "/users/mcvean/xilin/xilin/Multimorbidity-biobank/LDSC_SEG/"
common_disease_within_topics <- read.table("/users/mcvean/xilin/xilin/Multimorbidity-biobank/BOLT_LMM_subtype_list.txt")
# need to add all disease incidence case for comparison
common_disease_within_topics <- common_disease_within_topics

# everything we want to save
pheh2g <- matrix(NA, nrow = dim(common_disease_within_topics)[1], ncol = 1)
seh2g <- matrix(NA, nrow = dim(common_disease_within_topics)[1], ncol = 1)
for(i in 1:dim(common_disease_within_topics)[1]){
  ds_id <-  common_disease_within_topics$V1[i]
  topic_id <- common_disease_within_topics$V2[i]
  print(ds_id)
  try({
    # extract heritability for first subtype
    h2g.1 <- readLines(paste(DIR, ds_id,"_topic", topic_id,"/", ds_id, "_topic", topic_id,"_imputed.h2g.log",sep=""))
    pheh2g[i,1] <- as.numeric(str_split(h2g.1[length(h2g.1) - 6], "\\s+")[[1]][5])
    se_h2g <- str_split(h2g.1[length(h2g.1) - 6], "\\s+")[[1]][6]
    seh2g[i,1] <- as.numeric(substr(se_h2g, 2,nchar(se_h2g)-1))
  })
}
common_disease_within_topics$h2g <- pheh2g
common_disease_within_topics$seh2g <- seh2g
common_disease_within_topics <- common_disease_within_topics %>%
  rename(disease = V1, topic = V2, age = V3, N = V4) %>%
  mutate(z_score = pheh2g/seh2g)
write.csv(common_disease_within_topics ,paste(DIR, "h2g_imputed.csv",sep=""), row.names = FALSE )

#################################################################################
# step 5 -------- run LDSC seg with saved jackknife coefficients
#################################################################################
##########################################################
# use the jackknife samples to compute difference p-value
##########################################################
CTS_dir <- "~/Desktop/comorbidity/Multi-morbidity_biobank/CTS_results/CTS_results/"
disease_ldsc_seg <- read.table("BOLT_LMM_subtype_list.txt", header =  F) %>%
  arrange(V1)
names(disease_ldsc_seg) <- c("disease", "topic", "age",  "subtype_size")

cts_names <- read.table(paste0(CTS_dir, "153.2_topic4_imputed.cell_type_results.txt"), header = T) %>%
  pull(Name)

# topic_id <- 7
# disease_id <- 250.2

strong_genetic_subtype <- read.csv("Association_analysis/subtypes_Fst_matched_topic.csv") %>%
  filter(p_fst < 0.011) %>%
  pull(disease)

all_ldsc_seg <- read.csv("h2g_imputed.csv", header =  T)
names(all_ldsc_seg) <- c("disease", "topic", "age", "subtype_size", "h2g", "h2g_se", "h2g_zscore")

ds_high_h2g <- all_ldsc_seg %>%
  filter(disease %in% strong_genetic_subtype, topic != "all", h2g_zscore > 5) %>%
  left_join(select(phe_phecode, phecode, phenotype), by = c("disease" = "phecode") )


# ds_high_h2g <- all_ldsc_seg %>%
#   filter(h2g_zscore > 5, topic != "all") %>%
#   arrange(disease) %>%
#   left_join(select(phe_phecode, phecode, phenotype), by = c("disease" = "phecode") )

LDSC_seg_subtp_diff <- list()
for(ds_tp_id in 1:dim(ds_high_h2g)[1]){
  topic_id <- ds_high_h2g$topic[ds_tp_id]
  disease_id <- ds_high_h2g$disease[ds_tp_id]
  print(disease_id)

  cts_mean <- c()
  cts_se <- c()
  cts_p <- c()
  for(num_cts in 1:length(cts_names)){
    cts_id <- cts_names[num_cts]
    subtp_cts <- read.table(paste0(CTS_dir, disease_id, "_topic", topic_id,"_imputed.", cts_id, ".part_delete")) %>%
      pull(V54)

    all_cts <- read.table(paste0(CTS_dir, disease_id, "_topicall_imputed.", cts_id, ".part_delete")) %>%
      pull(V54)
    Jacknife_samples <- (subtp_cts - all_cts)
    cts_mean[num_cts] <- mean(Jacknife_samples)
    cts_se[num_cts] <- sqrt(var(Jacknife_samples) * (length(Jacknife_samples) - 1)^2/ length(Jacknife_samples) )
    cts_p[num_cts] <- (1-pnorm(cts_mean[num_cts]/cts_se[num_cts]))
  }
  LDSC_seg_subtp_diff[[ds_tp_id]] <- data.frame(Name = cts_names, diff_mean = cts_mean, diff_se = cts_se, P = cts_p) %>%
    mutate(topic = topic_id, disease = disease_id)

}
LDSC_seg_subtp_diff <- bind_rows(LDSC_seg_subtp_diff)
LDSC_seg_subtp_diff$FDR <- p.adjust(LDSC_seg_subtp_diff$P , method = "fdr")
LDSC_seg_subtp_diff %>%
  arrange(P)


plot_ds <- LDSC_seg_subtp_diff
plt_qq <- data.frame(y = sort(-log10(plot_ds$P)), x = sort(-log10(runif(length(plot_ds$P)))))

ggplot(plt_qq) +
  geom_point(aes(x = x, y = y), size = 0.5, alpha = 0.5) +
  geom_abline(slope = 1, linetype = "dashed", color = red) +
  theme(legend.position = "none",panel.background=element_blank()) +
  xlab(expression("Expected" * -log[10](P)) ) + ylab(expression("Observed" * -log[10](P))) +
  ggtitle(paste0("LDSC: Subtypes v.s. all disease"))


