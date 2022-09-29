#!/bin/bash


LIB=$1

BAMS=/cluster/home/aanderson/myers/multiomics_cellranger/211206/${LIB}/outs/gex_possorted_bam.bam
OUT_DIR=/cluster/home/aanderson/myers/cellsnp/211206/$LIB/OUT_${LIB}_1000G
BARCODES=/cluster/home/aanderson/myers/cellsnp/211206/barcodes/barcodes_${LIB}.tsv
REF2=/cluster/home/aanderson/software/gatk_resource/retry/1000G_phase1.snps.high_confidence.hg38.vcf.gz

cellsnp-lite -s $BAMS -b $BARCODES -O $OUT_DIR -p 30  --gzip  -R $REF2

