##########################################
# tissue specific analysis using LDSC
##########################################
# step 1: making the sum states file: starting from 250.2_topic7.logistic.a2 from the .R data
cd LDSC_seg/ldsc
size=11750 # sample size of the GWAS
source activate ldsc
cd ..
./ldsc/munge_sumstats.py \
        --sumstats 250.2_topic7.assoc \
        --N ${size} \
        --out 250.2_topic7 \
        --merge-alleles w_hm3.snplist

# step 2: download tissue specific data,also the baselineLD data (to adjust for the annotation baseline)
# following https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses
cts_name=Multi_tissue_gene_expr
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/${cts_name}_1000Gv3_ldscores.tgz
# tar -xvzf ${cts_name}_1000Gv3_ldscores.tgz
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
# tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
# tar -xvzf weights_hm3_no_hla.tgz

# step 3: LDSC_SEG
mkdir CTS_results
./ldsc/ldsc.py  \
	--h2-cts  250.2_topic7.sumstats.gz \
	--ref-ld-chr-cts Multi_tissue_gene_expr.ldcts \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--w-ld-chr weights_hm3_no_hla/weights. \
	--out CTS_results/250.2_topic7

# step 3': simple LDSC for heritability estimation
./ldsc/ldsc.py  \
        --h2  250.2_topic7.sumstats.gz \
        --ref-ld-chr /users/mcvean/xilin/xilin/UK_biobank/eur_w_ld_chr/ \
        --w-ld-chr  /users/mcvean/xilin/xilin/UK_biobank/eur_w_ld_chr/  \
        --out 250.2_topic7.h2g

# module load R/3.6.2-foss-2019b
# Rscript /users/mcvean/xilin/xilin/Multimorbidity-biobank/extract_h2_imputed.R

#########################################################################################################
#  Step 5 -- perform LDSC for each cell type while saving the JackKnife results
#########################################################################################################
# note here --print-delete-vals will save the jackknife coefficients; it could not be done with the wrapper --h2-cts
while read ct files;
do
/apps/well/ldsc/20180517/ldsc.py  \
        --h2  ${ds_id}_topic${topic_id}/${ds_id}_topic${topic_id}.bgen.sumstats.gz \
        --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline.,${files} \
        --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
        --w-ld-chr weights_hm3_no_hla/weights. \
        --overlap-annot --print-cov --print-coefficients --print-delete-vals \
        --out CTS_results/${ds_id}_topic${topic_id}_imputed.${ct}
done < Multi_tissue_gene_expr.ldcts










