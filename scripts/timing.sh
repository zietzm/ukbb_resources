#!/bin/bash

set -e

PHENOTYPE_ROOT=/data1/home/mnz2108/git/maxgcp-analysis/data/phenotypes
PT=$PHENOTYPE_ROOT/timing_pheno.tsv
CV=$PHENOTYPE_ROOT/covar_regenie.tsv
GT_ROOT=/data1/home/mnz2108/data_resources/ukbiobank/hapmap3_genotypes
GT=$GT_ROOT/hapmap3_variants_white_british
GTSUB=$GT_ROOT/hapmap3_variants_white_british_500k_mac100

SAIGE=/data1/home/mnz2108/git/SAIGE/extdata

TP=/data1/home/mnz2108/git/maxgcp-analysis/data/timing

################################################################################
# Plink 1 timing
################################################################################

/usr/bin/time -p -o $TP/plink1_1pt.txt plink \
	--bfile $GT \
	--pheno $PT \
	--pheno-name q_30700_0 \
	--prune \
	--covar $CV \
	--vif 100000 \
	--memory 25000 \
	--linear hide-covar \
	--out $TP/gwas/plink1_1pt

################################################################################
# Plink 2 timing
################################################################################

/usr/bin/time -p -o $TP/plink2_1pt.txt plink2 \
	--pfile $GT \
	--pheno $PT \
	--pheno-name q_30700_0 \
	--covar $CV \
	--require-covar age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
	--vif 100000 \
	--memory 25000 \
	--glm hide-covar log10 \
	--out $TP/gwas/plink2_1pt

/usr/bin/time -p -o $TP/plink2_5pt.txt plink2 \
	--pfile $GT \
	--pheno $PT \
	--pheno-name q_30700_0,q_30710_0,q_30720_0,q_30730_0,q_30740_0 \
	--covar $CV \
	--require-covar age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
	--vif 100000 \
	--memory 25000 \
	--glm hide-covar log10 \
	--out $TP/gwas/plink2_5pt

/usr/bin/time -p -o $TP/plink2_10pt.txt plink2 \
	--pfile $GT \
	--pheno $PT \
	--pheno-name q_30700_0,q_30710_0,q_30720_0,q_30730_0,q_30740_0,q_30750_0,q_30760_0,q_30770_0,q_30780_0,q_30790_0 \
	--covar $CV \
	--require-covar age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
	--vif 100000 \
	--memory 25000 \
	--glm hide-covar log10 \
	--out $TP/gwas/plink2_10pt

################################################################################
# Regenie timing (step 1)
################################################################################

/usr/bin/time -p -o $TP/regenie_1pt_step1.txt regenie \
	--step 1 \
	--pgen $GTSUB \
	--ref-first \
	--covarFile $CV \
	--covarColList "age,sex,PC{1:10}" \
	--phenoFile $PT \
	--phenoColList q_30700_0 \
	--strict \
	--force-qt \
	--bsize 1000 \
	--out $TP/gwas/regenie_1pt

/usr/bin/time -p -o $TP/regenie_5pt_step1.txt regenie \
	--step 1 \
	--pgen $GTSUB \
	--ref-first \
	--covarFile $CV \
	--covarColList "age,sex,PC{1:10}" \
	--phenoFile $PT \
	--phenoColList q_30700_0,q_30710_0,q_30720_0,q_30730_0,q_30740_0 \
	--strict \
	--force-qt \
	--bsize 1000 \
	--out $TP/gwas/regenie_5pt

/usr/bin/time -p -o $TP/regenie_10pt_step1.txt regenie \
	--step 1 \
	--pgen $GTSUB \
	--ref-first \
	--covarFile $CV \
	--covarColList "age,sex,PC{1:10}" \
	--phenoFile $PT \
	--phenoColList q_30700_0,q_30710_0,q_30720_0,q_30730_0,q_30740_0,q_30750_0,q_30760_0,q_30770_0,q_30780_0,q_30790_0 \
	--strict \
	--force-qt \
	--bsize 3500 \
	--out $TP/gwas/regenie_10pt

################################################################################
# Regenie timing (step 2)
################################################################################

/usr/bin/time -p -o $TP/regenie_1pt_step2.txt regenie \
	--step 2 \
	--pgen $GT \
	--ref-first \
	--covarFile $CV \
	--covarColList "age,sex,PC{1:10}" \
	--phenoFile $PT \
	--phenoColList q_30700_0 \
	--bsize 1000 \
	--pred $TP/gwas/regenie_1pt_pred.list \
	--out $TP/gwas/regenie_1pt

/usr/bin/time -p -o $TP/regenie_5pt_step2.txt regenie \
	--step 2 \
	--pgen $GT \
	--ref-first \
	--covarFile $CV \
	--covarColList "age,sex,PC{1:10}" \
	--phenoFile $PT \
	--phenoColList q_30700_0,q_30710_0,q_30720_0,q_30730_0,q_30740_0 \
	--bsize 1000 \
	--pred $TP/gwas/regenie_5pt_pred.list \
	--out $TP/gwas/regenie_5pt

/usr/bin/time -p -o $TP/regenie_10pt_step2.txt regenie \
	--step 2 \
	--pgen $GT \
	--ref-first \
	--covarFile $CV \
	--covarColList "age,sex,PC{1:10}" \
	--phenoFile $PT \
	--phenoColList q_30700_0,q_30710_0,q_30720_0,q_30730_0,q_30740_0,q_30750_0,q_30760_0,q_30770_0,q_30780_0,q_30790_0 \
	--bsize 1000 \
	--pred $TP/gwas/regenie_10pt_pred.list \
	--out $TP/gwas/regenie_10pt

################################################################################
# SAIGE timing (sparse GRM)
################################################################################

xsv select -d $'\t' FID,IID $PT |
	xsv fmt -t $'\t' -o $PHENOTYPE_ROOT/white_british_samples.txt

plink2 \
	--pfile $GTSUB \
	--keep $PHENOTYPE_ROOT/white_british_samples.txt \
	--make-bed \
	--out $GTSUB

/usr/bin/time -p -o $TP/saige_sparse_grm.txt Rscript $SAIGE/createSparseGRM.R \
	--plinkFile=$GTSUB \
	--nThreads=40 \
	--outputPrefix="${GTSUB}_sparse_kinship" \
	--numRandomMarkerforSparseKin=2000 \
	--relatednessCutoff=0.125

################################################################################
# SAIGE timing (step 1)
################################################################################

# Put the covariates and phenotypes together
xsv join -d $'\t' 1 $PT 1 $CV |
	xsv select '!FID[1]-IID[1]' |
	xsv fmt -t $'\t' \
		>$PHENOTYPE_ROOT/timing_pheno_covar.tsv

# 1000 random markers
shuf -n 1000 $GT.bim | awk '{print $2}' >${GT}.1000Markers.list

# Make a new plink file containing those 1000 markers
plink2 \
	--bfile $GTSUB \
	--extract ${GT}.1000Markers.list \
	--make-bed \
	--out ${GTSUB}_random1000

# Fit the null model and estimate a variance ratio
# q_30700_0,q_30710_0,q_30720_0,q_30730_0,q_30740_0,q_30750_0,q_30760_0,q_30770_0,q_30780_0,q_30790_0
/usr/bin/time -p -o $TP/saige_1pt_step1.txt Rscript $SAIGE/step1_fitNULLGLMM.R \
	--sparseGRMFile="${GTSUB}_sparse_kinship_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx" \
	--sparseGRMSampleIDFile="${GTSUB}_sparse_kinship_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
	--useSparseGRMtoFitNULL=TRUE \
	--plinkFile=${GT}_random1000 \
	--phenoFile=$PHENOTYPE_ROOT/timing_pheno_covar.tsv \
	--skipVarianceRatioEstimation=FALSE \
	--phenoCol=q_30700_0 \
	--covarColList=age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
	--qCovarColList=sex \
	--sampleIDColinphenoFile=IID \
	--traitType=quantitative \
	--invNormalize=TRUE \
	--nThreads=35 \
	--outputPrefix=$TP/gwas/saige_1pt \
	--IsOverwriteVarianceRatioFile=TRUE

################################################################################
# SAIGE timing (step 2)
################################################################################

/usr/bin/time -p -o $TP/saige_1pt_step2.txt Rscript $SAIGE/step2_SPAtests.R \
	--bedFile=${GT}.bed \
	--bimFile=${GT}.bim \
	--famFile=${GT}.fam \
	--AlleleOrder=ref-first \
	--SAIGEOutputFile=$TP/gwas/saige_1pt.result.txt \
	--minMAF=0 \
	--minMAC=20 \
	--GMMATmodelFile=$TP/gwas/saige_1pt.rda \
	--varianceRatioFile=$TP/gwas/saige_1pt.varianceRatio.txt \
	--LOCO=FALSE \
	--sparseGRMFile="${GTSUB}_sparse_kinship_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx" \
	--sparseGRMSampleIDFile="${GTSUB}_sparse_kinship_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
	--is_fastTest=TRUE
