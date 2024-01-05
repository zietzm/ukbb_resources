#!/bin/bash

# This script performs genetic analysis on the UK Biobank data.
# We start by running GWAS on the UK Biobank data, using Plink
# and Regenie.

# Stop if any command fails
set -e

SERVER=mimir

if [ "$SERVER" = "eir" ]; then
	# IID_FILE=/data2/michael/data_resources/ukbiobank/ukb_white_british_ids.tsv

	# Path to genotypes pgen/pvar/psam files
	GENOTYPES=/data2/michael/data_resources/ukbiobank/hapmap3_genotypes/hapmap3_variants_white_british

	# Paths to write phenotype and covariate files
	PHENOTYPE_ROOT=/data2/michael/maxgcp/phenotypes
	mkdir -p $PHENOTYPE_ROOT
	PHENOTYPES=$PHENOTYPE_ROOT/pheno.tsv
	COVARIATES=$PHENOTYPE_ROOT/covar.tsv

	# Paths to write GWAS results
	GWAS_PATH=/data2/michael/maxgcp/gwas_results
	mkdir -p $GWAS_PATH
elif [ "$SERVER" = "mimir" ]; then
	# IID_FILE=/data1/home/mnz2108/data_resources/ukbiobank/ukb_white_british_ids.tsv

	# Path to genotypes pgen/pvar/psam files
	GT_ROOT=/data1/home/mnz2108/data_resources/ukbiobank/hapmap3_genotypes
	GENOTYPES=$GT_ROOT/hapmap3_variants_white_british

	# Paths to write phenotype and covariate files
	PHENOTYPE_ROOT=/data1/home/mnz2108/git/maxgcp-analysis/data/phenotypes
	mkdir -p $PHENOTYPE_ROOT
	PHENOTYPES=$PHENOTYPE_ROOT/pheno.tsv
	COVARIATES=$PHENOTYPE_ROOT/covar.tsv

	# Paths to write GWAS results
	GWAS_PATH=/data1/home/mnz2108/git/maxgcp-analysis/data/gwas_results
	mkdir -p $GWAS_PATH

	# Path to write LD score regression results
	H2_PATH=/data1/home/mnz2108/git/maxgcp-analysis/data/h2_results
	mkdir -p $H2_PATH

	# Path to write genetic correlation results
	RG_PATH=/data1/home/mnz2108/git/maxgcp-analysis/data/rg_results
	mkdir -p $RG_PATH
else
	echo "Unknown server $SERVER"
	exit 1
fi

# cp /huggin/data2/mnz2108/igwas/pt/01_07_full_table.tsv $PHENOTYPE_ROOT/01_07_full_table.tsv
#
# xsv select eid,eid,age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
#   $PHENOTYPE_ROOT/01_07_full_table.tsv \
#   | sed 's/eid,eid/#FID,IID/' \
#   | xsv join 1 /dev/stdin 1 $IID_FILE \
#   | xsv select '#FID-PC10' \
#   | xsv fmt -t $'\t' \
#   > $COVARIATES
#
# xsv select 'eid,eid,b_A01-' $PHENOTYPE_ROOT/01_07_full_table.tsv \
#   | sed 's/eid,eid/#FID,IID/' \
#   | xsv join 1 /dev/stdin 1 $IID_FILE \
#   | xsv select '!#FID[1]-IID[1]' \
#   | xsv fmt -t $'\t' \
#   > $PHENOTYPES

# xsv select 'FID,IID,q_30700_0,q_30710_0,q_30720_0,q_30730_0,q_30740_0,q_30750_0,q_30760_0,q_30770_0,q_30780_0,q_30790_0' \
# 	$PHENOTYPE_ROOT/pheno_regenie.tsv |
# 	xsv fmt -t $'\t' \
# 		>$PHENOTYPE_ROOT/timing_pheno.tsv

plink2 \
	--pfile $GENOTYPES \
	--make-bed \
	--out $GENOTYPES

# plink2 \
# 	--pfile $GENOTYPES \
# 	--pheno $PHENOTYPES \
# 	--covar $COVARIATES \
# 	--require-covar age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
# 	--vif 100000 \
# 	--memory2500 25000 \
# 	--glm hide-covar log10 \
# 	--out $GWAS_PATH/plink_white_british
#
# plink2 \
# 	--pfile $GENOTYPES \
# 	--pheno $PHENOTYPE_ROOT/quant_pheno_final.tsv \
# 	--covar $COVARIATES \
# 	--require-covar age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
# 	--vif 100000 \
# 	--memory 25000 \
# 	--glm hide-covar log10 \
# 	--out $GWAS_PATH/plink_white_british

##############################################################################
########################### Regenie GWAS #####################################
##############################################################################

# Prepare data for regenie step 1
# 500k variants
# MAC filter 100

# plink2 \
# 	--pfile $GENOTYPES \
# 	--mac 100 \
# 	--thin-count 500000 \
# 	--make-pgen \
# 	--out "$GT_ROOT"/hapmap3_variants_white_british_500k_mac100

# Edit the phenotype and covariate files. Regenie wants "FID IID", not "#FID IID"
# head -n 1 $PHENOTYPES | sed 's/#//' > $PHENOTYPE_ROOT/pheno_regenie.tsv
# tail -n +2 $PHENOTYPES >> $PHENOTYPE_ROOT/pheno_regenie.tsv

# head -n 1 $COVARIATES | sed 's/#//' > $PHENOTYPE_ROOT/covar_regenie.tsv
# tail -n +2 $COVARIATES >>$PHENOTYPE_ROOT/covar_regenie.tsv

# function do_regenie {
# 	CHUNK=$1
# 	CHUNKSTART=$((10 * CHUNK + 3))
# 	CHUNKEND=$((10 * CHUNK + 12))
#
# 	# 2638 phenotypes
# 	if [ "$CHUNK" -eq 263 ]; then
# 		CHUNKEND=2640
# 	fi
#
# 	xsv select "1-2,${CHUNKSTART}-${CHUNKEND}" $PHENOTYPE_ROOT/pheno_regenie.tsv |
# 		xsv fmt -t $'\t' \
# 			>"${PHENOTYPE_ROOT}/pheno_regenie_filtered_${CHUNK}.tsv"
#
# 	regenie \
# 		--step 1 \
# 		--pgen "$GT_ROOT/hapmap3_variants_white_british_500k_mac100" \
# 		--ref-first \
# 		--covarFile $PHENOTYPE_ROOT/covar_regenie.tsv \
# 		--covarColList "age,sex,PC{1:10}" \
# 		--phenoFile "$PHENOTYPE_ROOT/pheno_regenie_filtered_$CHUNK.tsv" \
# 		--strict \
# 		--force-qt \
# 		--bsize 1000 \
# 		--lowmem \
# 		--out "$GWAS_PATH/regenie_white_british_chunk_$CHUNK"
# }
#
# for CHUNK in {0..263}; do
# 	do_regenie "$CHUNK"
# done

# regenie \
# 	--step 1 \
# 	--pgen "$GT_ROOT"/hapmap3_variants_white_british_500k_mac100 \
# 	--ref-first \
# 	--covarFile $PHENOTYPE_ROOT/covar_regenie.tsv \
# 	--covarColList "age,sex,PC{1:10}" \
# 	--phenoFile $PHENOTYPE_ROOT/pheno_regenie_filtered.tsv \
# 	--strict \
# 	--force-qt \
# 	--bsize 20000 \
# 	--lowmem \
# 	--out $GWAS_PATH/regenie_white_british

# regenie step 2 on the 100 random phenotypes
# regenie \
# 	--step 2 \
# 	--pgen $GENOTYPES \
# 	--ref-first \
# 	--covarFile $COVARIATES \
# 	--covarColList "age,sex,PC{1:10}" \
# 	--phenoFile $PHENOTYPES \
# 	--bsize 1000 \
# 	--pred $GWAS_PATH/regenie_white_british_pred.list \
# 	--out $GWAS_PATH/regenie_white_british
