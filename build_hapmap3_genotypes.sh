#!/bin/bash

# This script builds the data_resources/ukbiobank/hapmap3_genotypes/ directory
# for imputed SNP data from the UK Biobank. It downloads the HapMap3 variant list
# from LDAK, and then uses plink to extract the SNPs from the imputed data.
# In addition to GNU coreutils, this assumes that Plink, Bgenix, and XSV are
# installed in the PATH. This is supposed to be run on eir, a Tatonetti Lab server.

mkdir -p /data2/michael/data_resources/ukbiobank/hapmap3_genotypes/
cd /data2/michael/data_resources/ukbiobank/hapmap3_genotypes/

# The HapMap3 tagging files for LDAK is available from the LDAK website
wget https://genetics.ghpc.au.dk/doug/bld.ldak.hapmap.gbr.tagging.gz
gunzip bld.ldak.hapmap.gbr.tagging.gz

wget https://genetics.ghpc.au.dk/doug/ldak.thin.hapmap.gbr.tagging.gz
gunzip ldak.thin.hapmap.gbr.tagging.gz

# Get the variant list from the LDAK tagging file. Note, variants are formatted as
# CHR:POS (e.g. 1:10000) from Chr37/hg19. This file has the variant names in the
# first column, so we use cut to extract the first column. We have to add a flag
# for the delimiter, which is a space. We also remove the first row, which is a
# header.
cut -d' ' -f1 ldak.thin.hapmap.gbr.tagging \
  | tail -n +2 \
  | sort -k1,1 \
  > hapmap3_variants.txt

# The data we have from the UK Biobank is also in the Chr37/hg19 build. To prove
# this, take a look at /huggin/data2/ukb_data/22418_genotyping_plink/ukb_snp_qc.txt.
# The first line after the header gives: rs28659788 ... 1 723307 C G ...
# This is the variant rs28659788, which is on chromosome 1 at position 723307 in
# Chr37/hg19 (see https://www.ncbi.nlm.nih.gov/snp/?term=rs28659788).

# Imputed variants are in .bgen files, so we will use bgenix and the UK Biobank
# provided index files to extract the variants we want. Unfortunately, we have
# identifiers for HapMap3 variants only in CHR:BP format, so we'll first map
# these to rsIDs using the UK Biobank provided variant list, then query the
# bgen files using rsIDs, then convert the bgen files to plink format, then
# map the rsIDs back to CHR:BP format for easy use with LDAK.

BGEN_DIR=/huggin/data2/ukb_data/22828_imputed_genotype_bgen
BGENIX_DIR=/huggin/data2/ukb_data/100319_imputed_genotype_bgi

process_chromosome() {
  chr=$1

  BGEN_DIR=/huggin/data2/ukb_data/22828_imputed_genotype_bgen
  BGENIX_DIR=/huggin/data2/ukb_data/100319_imputed_genotype_bgi
  OUT_DIR=/data2/michael/data_resources/ukbiobank/hapmap3_genotypes/

  echo "Processing chromosome $chr"

  # Create a map between rsid and CHR:BP for the HapMap3 variants,
  # using only the variants available in the imputed data, on the current
  # chromosome. These files have two columns: CHR:BP and rsID.
  bgenix \
    -g $BGEN_DIR/ukb22828_imp_chr${chr}_v3.bgen \
    -i $BGENIX_DIR/ukb100319_imp_chr${chr}_v3.bgen.bgi \
    -list \
    > $OUT_DIR/chr${chr}_variants_bgenix.txt

  # Since the bgen files are on huggin, and since we have to transfer data
  # between eir and huggin, we'll do all this processing locally instead of
  # in the pipe.
  cat $OUT_DIR/chr${chr}_variants_bgenix.txt \
    | tail -n +3 \
    | awk '{print $3":"$4"\t"$2}' - \
    | sed 's/^0//g' \
    | sort -k1,1 \
    | join -1 1 -2 1 $OUT_DIR/hapmap3_variants.txt - \
    > $OUT_DIR/chr${chr}_variants_map.txt

  # Extract the rsIDs from the map file (to query the bgen files).
  cut -d' ' -f2 $OUT_DIR/chr${chr}_variants_map.txt \
    > $OUT_DIR/chr${chr}_rsids.txt

  # Extract the variants from the bgen files.
  bgenix \
    -g $BGEN_DIR/ukb22828_imp_chr${chr}_v3.bgen \
    -i $BGENIX_DIR/ukb100319_imp_chr${chr}_v3.bgen.bgi \
    -incl-rsids $OUT_DIR/chr${chr}_rsids.txt \
    > $OUT_DIR/chr${chr}_hapmap3_variants.bgen

  # Convert the bgen files to plink2 format.
  plink2 \
    --bgen $OUT_DIR/chr${chr}_hapmap3_variants.bgen \
    --sample $BGEN_DIR/ukb22828_c${chr}_b0_v3_s487280.sample \
    --rm-dup exclude-mismatch \
    --make-pgen \
    --out $OUT_DIR/chr${chr}_hapmap3_variants

  # Map the rsIDs back to CHR:BP format.
  plink2 \
    --pfile $OUT_DIR/chr${chr}_hapmap3_variants \
    --update-name $OUT_DIR/chr${chr}_variants_map.txt 1 2 \
    --make-pgen \
    --out $OUT_DIR/chr${chr}_hapmap3_variants_mapped

  echo "$OUT_DIR/chr${chr}_hapmap3_variants_mapped" >> $OUT_DIR/chr_merge_list.txt
}

export -f process_chromosome

# Process each chromosome in parallel.
parallel --jobs 22 process_chromosome ::: {1..22}

# Merge the plink files.
plink2 \
  --pmerge-list chr_merge_list.txt pfile \
  --make-pgen \
  --out hapmap3_variants

# Create a file with eids for White British individuals.
# Note field 21000 (ethnic background) https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=21000
# Note coding https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=1001
# Note value 1001 = White British
xsv search -s '"21000-0.0"' "1001" /huggin/data2/ukb_data/ukb23674_data_pull/100065_ethnicity.csv \
  | xsv search -s '"21000-1.0"' "(^1001$|^$)" \
  | xsv search -s '"21000-2.0"' "(^1001$|^$)" \
  | xsv select eid,eid \
  | xsv sort \
  | tail -n +2 \
  | cat <(echo "#FID,IID") - \
  | xsv fmt -t '\t' \
  > /data2/michael/data_resources/ukbiobank/ukb_white_british_ids.tsv

# Filter genotypes for White British individuals.
plink2 \
  --pfile hapmap3_variants \
  --keep /data2/michael/data_resources/ukbiobank/ukb_white_british_ids.txt \
  --make-pgen \
  --out hapmap3_variants_white_british
