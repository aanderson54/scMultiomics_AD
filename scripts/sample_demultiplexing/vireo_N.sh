#!/bin/bash


LIB=$1

CELL_DATA=/cluster/home/aanderson/myers/cellsnp/211206/$LIB/OUT_${LIB}_1000G
OUT_DIR=/cluster/home/aanderson/myers/cellsnp/211206/$LIB/vireo_${LIB}_1000G
DONOR_GT_FILE=REF2=/cluster/projects/ADFTD/BrainTF-ADsupp/croo/ADsuppGenomes_dbSNP-154_wClinvar.vcf.gz


vireo -c $CELL_DATA  -o $OUT_DIR -N 2

