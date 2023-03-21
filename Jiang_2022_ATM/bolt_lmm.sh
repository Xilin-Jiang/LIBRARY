#################################################
# step 1: create a imputed file with less SNPs
#################################################
# from LDSC_seg/run_create_noX.sh # create a bfile list with MAF > 0.001 and chromosome from 1-23
/apps/well/plink/1.90b2n/plink \
                  --bfile /well/mcvean/xilin/UK_biobank/all_chromosome \
                  --maf 0.001 \
                  --chr 1-23 \
                  --make-bed --out /well/mcvean/xilin/UK_biobank/all_chromosome_maf0.001_chr1_23

###################################################################################
# step 2: create a samples (needs to be smaller and case > 10%); in the .R file
###################################################################################
topic_id=$(awk -v var="${1}" 'FNR == var { print $2 }' ../BOLT_LMM_subtype_list.txt )
ds_id=$(awk -v var="${1}" 'FNR == var { print $1}' ../BOLT_LMM_subtype_list.txt )
/users/mcvean/xilin/xilin/softwares/BOLT-LMM_v2.4/bolt \
                  --bfile=/well/mcvean/xilin/UK_biobank/all_chromosome_maf0.001_chr1_23 \
                  --remove=bolt.in_plink_but_not_imputed.FID_IID.968.txt \
                  --phenoFile=${ds_id}_topic${topic_id}/${ds_id}_topic${topic_id}_keep.txt \
                  --phenoCol=phenotype \
                  --lmm \
                  --LDscoresFile=/apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
                  --numThreads=8 \
                  --statsFile=${ds_id}_topic${topic_id}/${ds_id}_topic${topic_id}.lmm.stats \
                  --bgenFile=/well/ukbb-wtchg/v3/imputation/ukb_imp_chr{1:22}_v3.bgen \
                  --bgenMinMAF=1e-3 \
                  --bgenMinINFO=0.5 \
                  --sampleFile=/users/mcvean/xilin/xilin/UK_biobank/ukb1062_imp_chr1_v2_s487398.sample \
                  --statsFileBgenSnps=${ds_id}_topic${topic_id}/${ds_id}_topic${topic_id}.bolt.bgen.stat \
                  --verboseStats \
                  2>&1 | tee ${ds_id}_topic${topic_id}/${ds_id}_topic_${topic_id}.lmm.log # log output written to stdout and stderr

# basic args:
# --bfile: prefix of PLINK genotype files (bed/bim/fam)
# --remove: list of individuals to remove (FID IID)
# --exclude: list of SNPs to exclude (rs###)
# --phenoFile: phenotype file
# --phenoCol: column of phenotype file containing phenotypes
# --covarFile: covariate file
# --covarCol: column(s) containing categorical covariate (multiple ok)
# --qCovarCol: column(s) containing quantitative covariates (array format)
# --modelSnps: subset of SNPs to use in GRMs
# --lmm: flag to perform default BOLT-LMM mixed model association
# --LDscoresFile: reference LD Scores (data table in separate download)
# --numThreads: multi-threaded execution

# additional args for association testing on imputed SNPs:
# --statsFile: output file for association statistics at PLINK-format SNPs
# --dosageFile: file(s) containing additional dosage-format SNPs (multiple ok)
# --dosageFidIidFile: file containing FIDs and IIDs for dosage-format SNPs
# --statsFileDosageSnps: output file for assoc stats at dosage-format SNPs
# --impute2FileList: file listing chroms and IMPUTE2-format additional SNPs
# --impute2FidIidFile: file containing FIDs and IIDs for IMPUTE2-format SNPs
# --impute2CallThresh: minimum pAA+pAB+pBB for calling IMPUTE2-format SNPs
# --statsFileImpute2Snps: output file for assoc stats at IMPUTE2-format SNPs
# --dosage2FileList: file listing map and 2-dosage format additional SNPs
# --statsFileDosage2Snps: output file for assoc stats at 2-dosage format SNPs
