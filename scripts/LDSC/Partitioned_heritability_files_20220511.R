# open R 4.1.0
bsub -R "rusage[mem=10000]" -n 10 -Is bash
module load g/gsl/1.15   
module load g/htslib/1.9-devel
module load g/build_tools/llvm/8.0
module load g/util/gcc/8.2
module load g/R/4.1.0

R

library(scales)
library(dplyr)
library(tidyr)
library(stringr)


# Setting work_dir and Rvar file
work_dir = "/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/Partitioned_heritability_all_links_files/"
setwd(work_dir)

filelist = list.files(pattern = "/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/LDSC_all_links_results/*_AD_Kunkle_Stage1_not_APOE_HLA.GWAS.results")

dir_results_data = "/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/LDSC_all_links_results"
all_partitioned_heritability_file_list = list.files(path=dir_results_data,pattern="*_AD_Kunkle_Stage1_not_APOE_HLA.GWAS.results",full.names=T)

#assuming tab separated values with a header    
datalist = lapply(all_partitioned_heritability_file_list, function(x)read.table(x, header=T))

for (i in 1:length(datalist)){datalist[[i]]<-cbind(datalist[[i]],all_partitioned_heritability_file_list[i])}

for (i in 1:length(datalist)){datalist[[i]] <- datalist[[i]] %>%
  mutate(Coefficient_p = pnorm(`Coefficient_z.score`, lower.tail = FALSE)) %>%
  mutate(
    Coefficient_holm = p.adjust(Coefficient_p, method = "holm"),
    Enrichment_holm = p.adjust(Enrichment_p, method = "holm"),
    Coefficient_holm_cutoff =
      max(Coefficient_p[Coefficient_holm < 0.05], na.rm = TRUE),
    Enrichment_holm_cutoff =
      max(Enrichment_p[Enrichment_holm < 0.05], na.rm = TRUE)) %>%
  mutate(sig_coef = Coefficient_holm < 0.05)}

for (i in 1:length(datalist)){datalist[[i]] <- datalist[[i]] %>%
  mutate(
    Coefficient_BH = p.adjust(Coefficient_p, method = "BH"),
    Enrichment_BH = p.adjust(Enrichment_p, method = "BH"),
    Coefficient_BH_cutoff =
      max(Coefficient_p[Coefficient_BH < 0.05], na.rm = TRUE),
    Enrichment_holm_cutoff =
      max(Enrichment_p[Enrichment_BH < 0.05], na.rm = TRUE)) %>%
  mutate(sig_coef = Coefficient_BH < 0.05)}

result <- lapply(datalist, function(x) {x[1,]})

#assuming the same header/columns for all files
datafr = do.call("rbind", result)

colnames(datafr)[11] <- "Functional_category"

datafr$Functional_category <- as.factor(datafr$Functional_category)

levels(datafr$Functional_category) <- gsub("/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/LDSC_all_links_results/", "", levels(datafr$Functional_category), fixed=TRUE)
levels(datafr$Functional_category) <- gsub("_AD_Kunkle_Stage1_not_APOE_HLA.GWAS.results", "", levels(datafr$Functional_category), fixed=TRUE)
levels(datafr$Functional_category) <- gsub("scAD_all_links_", "", levels(datafr$Functional_category), fixed=TRUE)

datafr <- datafr[order(datafr$Enrichment_BH),]

write.csv(datafr, file ="/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/RINdrop_filtered_20220511/LDSC_analysis/Partitioned_heritability_all_links_files/scAD_Partitioned_Heritability_all_links_AD_Kunkle_Stage1_not_APOE_HLA.csv")
