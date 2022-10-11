##LDSC analysis pipeline AD multiomics
Author: Ivan Rodriguez-Nunez


1) Formatting GWAS sumstats (For AD sumstats I removed HLA and APOE regions)
gzcat /Users/ivanrodriguez/Documents/Ivan/Myers_Lab/LDSC-SEG/GWAS_sumstats/PMID_33589840/GCST90012877/GCST90012877_buildGRCh37.tsv.gz | awk '!($3 == 19 && $4 >= 45000000 && $4 <= 45800000)' | gzip > /Users/ivanrodriguez/Documents/Ivan/Myers_Lab/LDSC-SEG/GWAS_sumstats/PMID_33589840/GCST90012877/AD_Schwartzentruber_not_APOE.txt.gz
gzcat /Users/ivanrodriguez/Documents/Ivan/Myers_Lab/LDSC-SEG/GWAS_sumstats/PMID_33589840/GCST90012877/AD_Schwartzentruber_not_APOE.txt.gz  | awk '!($3 == 6 && $4 >= 28477797 && $4 <= 33448354)' | gzip > /Users/ivanrodriguez/Documents/Ivan/Myers_Lab/LDSC-SEG/GWAS_sumstats/PMID_33589840/GCST90012877/AD_Schwartzentruber_not_APOE_HLA.txt.gz 			  

python ./munge_sumstats.py \
--sumstats ./AD_Schwartzentruber_not_APOE_HLA.txt.gz \
--ignore variant_id,GWAS_BETA,GWAS_SE,GWAS_P,GWAX_UKBB_BETA,GWAX_UKBB_SEmGWAX_UKBB_P,DIRECT,HET_P \
--snp SNP_ID \
--a1 effect_allele \
--a2 other_allele \
--info INFO \
--N-cas 75024 \
--N-con 397844 \
--out ./AD_Schwartzentruber_not_APOE_HLA \
--merge-alleles ./w_hm3.snplist \

2) Script to create annotation files and run LDSC

/cluster/home/irodriguez/Myers_Lab/LDscore/scripts/LDSC_Analysis12.sh

bsub -n 12 -o ./testlog.txt /cluster/home/irodriguez/Myers_Lab/LDscore/scripts/bsub_LDSC12.sh \
-d /cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/bed_all_links_filtered_files/ \
-o /cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/annotation_files \
-b /cluster/home/irodriguez/Myers_Lab/LDscore/LDscore_GRCh38/plink_files \
-g /cluster/home/irodriguez/Myers_Lab/LDscore/GWAS_summary/AD_Kunkle_Stage1_not_APOE_HLA.sumstats.gz 

3) Script to run LDSC with annotation files already created 

/cluster/home/irodriguez/Myers_Lab/LDscore/scripts/LDSC_Analysis13.sh

bsub -n 12 -o ./testlog.txt /cluster/home/irodriguez/Myers_Lab/LDscore/scripts/bsub_LDSC13.sh \
-d /cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/bed_all_links_filtered_files/ \
-o /cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/annotation_files \
-g /gpfs/gpfs1/home/irodriguez/Myers_Lab/LDscore/GWAS_summary/AD_Schwartzentruber_not_APOE_HLA.sumstats.gz

4) Script for multitest correction of LDSC results
/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/Partitioned_heritability_files_20220511.R

5) Script to make heatmap
/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/LDSC_heatmap_all_links_filtered_BH_20220729.R