# genetic tools
##################################################
# MR analysis
##################################################
ibrary(MendelianRandomization)
eggerMR <- function(loadings, disease, genetics, method = "ivw") {
  d_pvalue <- rep(NA, dim(genetics)[2])
  d_effects <- rep(NA, dim(genetics)[2])
  d_se <- rep(NA, dim(genetics)[2])
  t_pvalue <- rep(NA, dim(genetics)[2])
  t_effects <- rep(NA, dim(genetics)[2])
  t_se <- rep(NA, dim(genetics)[2])
  for(g_id in 1:dim(genetics)[2]){
    g2T <- lm(loadings ~ genetics[,g_id])
    t_pvalue[g_id] <- summary(g2T)$coefficients[2,4]
    t_effects[g_id] <- summary(g2T)$coefficients[2,1]
    t_se[g_id] <- summary(g2T)$coefficients[2,2]

    g2D <- glm(disease ~ genetics[,g_id],
               family = binomial)
    d_pvalue[g_id] <- summary(g2D)$coefficients[2,4]
    d_effects[g_id] <- summary(g2D)$coefficients[2,1]
    d_se[g_id] <- summary(g2D)$coefficients[2,2]
  }

  IV_topic <- which(p.adjust(t_pvalue,method = "fdr") < 0.05)
  IV_disease <- which(p.adjust(d_pvalue,method = "fdr") < 0.05)

  if(length(IV_topic) > 0){
    t2d_MRdata <- mr_input(bx = t_effects[IV_topic], bxse = t_se[IV_topic],
                           by = d_effects[IV_topic], byse = d_se[IV_topic])
    if(method == "egger"){
      topic2diseaseMR <- mr_egger(t2d_MRdata)
    }else if(method == "ivw"){
      topic2diseaseMR <- mr_ivw(t2d_MRdata)
    }

  }else{
    topic2diseaseMR <- NULL
  }

  if(length(IV_disease) > 0){
    d2t_MRdata <- mr_input(bx = d_effects[IV_disease], bxse = d_se[IV_disease],
                           by = t_effects[IV_disease], byse = t_se[IV_disease])
    if(method == "egger"){
      disease2topicMR <- mr_egger(d2t_MRdata)
    }else if(method == "ivw"){
      disease2topicMR <- mr_ivw(d2t_MRdata)
    }

  }else{
    disease2topicMR <- NULL
  }
  MR_results <- list(topic2diseaseMR, disease2topicMR)
  return(MR_results)
}
# small example

