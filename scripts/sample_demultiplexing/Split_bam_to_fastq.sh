#!/bin/bash


BAM_FILE=$1
SAMPLE=$2
export PATH=$PATH:~/software/cellranger-arc-2.0.0/


#samtools view -H $BAM_FILE > ${SAMPLE}_2_SAM_header

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
#samtools view $BAM_FILE | LC_ALL=C grep -F -f ${SAMPLE}.txt > ${SAMPLE}_2_filtered_SAM_body

# Combine header and body
#cat ${SAMPLE}_2_SAM_header ${SAMPLE}_2_filtered_SAM_body > ${SAMPLE}_2_filtered.sam

# Convert filtered.sam to BAM format
#samtools view -b ${SAMPLE}_2_filtered.sam > atac_${SAMPLE}_2.bam
samtools collate -o  ${SAMPLE}_tmp.bam atac_${SAMPLE}_2.bam 
samtools fixmate -m ${SAMPLE}_tmp.bam ${SAMPLE}_tmp2.bam
samtools sort -o atac_${SAMPLE}_2.bam ${SAMPLE}_tmp2.bam
samtools markdup -r  atac_${SAMPLE}_2.bam  atac_${SAMPLE}.bam

#mkdir ${SAMPLE}_2_atac

/cluster/home/aanderson/software/cellranger-arc-2.0.1/lib/bin/bamtofastq --nthreads=8 atac_${SAMPLE}.bam ${SAMPLE}_atac
