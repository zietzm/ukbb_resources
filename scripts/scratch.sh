#!/bin/bash

set -e

full_gwas() {
	UKB_DIR=/data1/deep_storage/ukbiobank
	BGEN_DIR=$UKB_DIR/imp_bgen_files

	# Paths to write phenotype and covariate files
	PHENOTYPE_ROOT=/data1/home/mnz2108/git/maxgcp-analysis/data/phenotypes
	PHENOTYPES=$PHENOTYPE_ROOT/pheno.tsv
	COVARIATES=$PHENOTYPE_ROOT/covar.tsv

	# Paths to write GWAS results
	GWAS_PATH=/data1/home/mnz2108/git/maxgcp-analysis/data/gwas_results

	for chr in {1..22}; do
		echo "Processing chromosome $chr"
		plink2 \
			--bgen $BGEN_DIR/ukb_imp_chr"${chr}"_v3.bgen ref-first \
			--sample $UKB_DIR/genotypes/ukb41019_imp_chr"${chr}"_v3_s487282.sample \
			--pheno $PHENOTYPES \
			--pheno-name b_A01 \
			--covar $COVARIATES \
			--require-covar age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
			--vif 100000 \
			--memory 25000 \
			--glm hide-covar log10 \
			--out $GWAS_PATH/timing_chr"${chr}"
	done
}

time full_gwas
