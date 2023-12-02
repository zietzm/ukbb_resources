#!/bin/bash

# This script performs genetic analysis on the UK Biobank data.
# We start by running GWAS on the UK Biobank data, using Plink
# and Regenie.

IID_FILE=/data2/michael/data_resources/ukbiobank/ukb_white_british_ids.tsv

# Path to genotypes pgen/pvar/psam files
GENOTYPES=/data2/michael/data_resources/ukbiobank/hapmap3_genotypes/hapmap3_variants_white_british

PHENOTYPE_ROOT=/data2/michael/maxgcp/phenotypes
PHENOTYPES=$PHENOTYPE_ROOT/pheno.tsv
COVARIATES=$PHENOTYPE_ROOT/covar.tsv

GWAS_PATH=/data2/michael/maxgcp/gwas_results

cp /huggin/data2/mnz2108/igwas/pt/01_07_full_table.tsv $PHENOTYPE_ROOT/01_07_full_table.tsv

xsv select eid,eid,age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  $PHENOTYPE_ROOT/01_07_full_table.tsv \
  | sed 's/eid,eid/#FID,IID/' \
  | xsv join 1 /dev/stdin 1 $IID_FILE \
  | xsv fmt -t $'\t' \
  > $COVARIATES

xsv select 'eid,b_A01-' $PHENOTYPE_ROOT/01_07_full_table.tsv \
  | xsv fmt -t $'\t' \
  > $PHENOTYPES

plink2 \
  --pfile $GENOTYPES \
  --pheno $PHENOTYPES \
  --covar $COVARIATES \
  --require-covar age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --vif 100000 \
  --memory 35000 \
  --glm hide-covar log10 \
  --out $GWAS_PATH/plink_white_british

regenie \
  --step 1 \
  --pgen $GENOTYPES \
  --ref-first \
  --covarFile $COVARIATES \
  --covarColList "age,sex,PC{1:10}" \
  --phenoFile $PHENOTYPES \
  --bsize 10000 \
  --out $GWAS_PATH/regenie_white_british

# regenie step 2 on the 100 random phenotypes
regenie \
  --step 2 \
  --pgen $GENOTYPES \
  --ref-first \
  --covarFile $COVARIATES \
  --covarColList "age,sex,PC{1:10}" \
  --phenoFile $PHENOTYPES \
  --bsize 1000 \
  --pred $GWAS_PATH/regenie_white_british_pred.list \
  --out $GWAS_PATH/regenie_white_british
