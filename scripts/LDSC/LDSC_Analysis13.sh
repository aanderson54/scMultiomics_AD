#!/bin/bash

export MKL_NUM_THREADS=8

BED=$1
#BIMFILE=$2
INPUT_BED=$2
OUTPUT_DIR=$3
#BIM_DIR=$4
samp=$4
#CHR=$7
GWAS=$5


study=`basename ${GWAS} | cut -d "." -f1`
RESULTS_DIR="/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/LDSC_all_links_results/"

python /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/scripts/ldsc.py --h2 ${GWAS} --w-ld-chr /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/LDscore_GRCh38/weights/weights.hm3_noMHC. --ref-ld-chr ${OUTPUT_DIR}/${samp}.,/gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/LDscore_GRCh38/baselineLD_v2.2/baselineLD. --overlap-annot --frqfile-chr /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/1000G_Phase3_frq/1000G.EUR.QC. --out ${RESULTS_DIR}/${samp}_${study}.GWAS --print-coefficients &>/gpfs/gpfs1/home/irodriguez/Myers_Lab/tmp/${samp}.baseline.ldsclog
