#!/bin/bash

export MKL_NUM_THREADS=8

BED=$1
#BIMFILE=$2
INPUT_BED=$2
OUTPUT_DIR=$3
BIM_DIR=$4
samp=$5
#CHR=$7
GWAS=$6

study=`basename ${GWAS} | cut -d "." -f1`
RESULTS_DIR="/cluster/home/irodriguez/Myers_Lab/20220430_RNAseq_Alasooetal-2018_iPSCs_Macrophages_challenged/LDSC_ATACseq/LDSC_analysis/LDSC_results/"

for CHR in {1..22};do
BIMFILE=${BIM_DIR}/1000G.EUR.hg38.${CHR}.bim

python /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/scripts/make_annot.py --bed-file ${BED} --bimfile ${BIMFILE} --annot-file ${OUTPUT_DIR}/${samp}.${CHR}.annot.gz &>/gpfs/gpfs1/home/irodriguez/Myers_Lab/tmp/${samp}.${CHR}.annotlog



python /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/scripts/ldsc.py --l2 --bfile ${BIM_DIR}/1000G.EUR.hg38.${CHR} --ld-wind-cm 1 --annot ${OUTPUT_DIR}/${samp}.${CHR}.annot.gz --thin-annot --out ${OUTPUT_DIR}/${samp}.${CHR} \
--print-snps /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/1000G_EUR_Phase3_baseline/print_snps.txt &>/gpfs/gpfs1/home/irodriguez/Myers_Lab/tmp/${samp}.${CHR}.ldsclog
done


python /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/scripts/ldsc.py --h2 ${GWAS} --w-ld-chr /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/LDscore_GRCh38/weights/weights.hm3_noMHC. --ref-ld-chr ${OUTPUT_DIR}/${samp}.,/gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/LDscore_GRCh38/baselineLD_v2.2/baselineLD. --overlap-annot --frqfile-chr /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/1000G_Phase3_frq/1000G.EUR.QC. --out ${RESULTS_DIR}/${samp}_${study}.GWAS --print-coefficients &>/gpfs/gpfs1/home/irodriguez/Myers_Lab/tmp/${samp}.baseline.ldsclog
