# R script to generate heatmaps for the LDSC results with cell-type
# specific feature links from the single cell Multiomics study with controls and AD cases

#Load packages
library(ComplexHeatmap)
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(circlize)

# Set working directory
work_dir = "/Users/ivanrodriguez/Documents/Ivan/Myers_Lab/AD_supp/20220729_LDSC/LDSC_all_links_results"
setwd(work_dir)

# List of LDSC_sc Multiomics results to represent in a heatmap
filelist = list.files(pattern = "/Users/ivanrodriguez/Documents/Ivan/Myers_Lab/AD_supp/20220511_LDSC/LDSC_all_links_results/*.csv")

dir_results_data = "/Users/ivanrodriguez/Documents/Ivan/Myers_Lab/AD_supp/20220511_LDSC/LDSC_all_links_results"
all_partitioned_heritability_file_list = list.files(path=dir_results_data,pattern="*.csv",full.names=T)

#assuming tab separated values with a header    
datalist = lapply(all_partitioned_heritability_file_list, function(x)read.csv(x))

# add a column with the file name
for (i in 1:length(datalist)){datalist[[i]]<-cbind(datalist[[i]],all_partitioned_heritability_file_list[i])}

#assuming the same header/columns for all files
datafr = do.call("rbind", datalist)

# Formatting the dataframe
colnames(datafr)[22] <- "GWAS"

datafr$Functional_category <- as.factor(datafr$Functional_category)
datafr$GWAS <- as.factor(datafr$GWAS)

levels(datafr$GWAS) <- gsub("/Users/ivanrodriguez/Documents/Ivan/Myers_Lab/AD_supp/20220511_LDSC/LDSC_all_links_results/scAD_Partitioned_Heritability_all_links_", "", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub(".csv", "", levels(datafr$GWAS))

levels(datafr$GWAS) <- gsub("AD_Bellenguez_GRCh38_not_APOE_HLA", "AD, Bellenguez et al.(2022)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("AD_Jansen_not_APOE_HLA", "AD, Jansen et al.(2019)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("AD_Kunkle_Stage1_not_APOE_HLA", "AD, Kunkle et al.(2019)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("AD_Schwartzentruber_not_APOE_HLA", "AD, Schwartzentruber et al.(2021)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("AD_Wightman_not_APOE_HLA", "AD, Wightman et al.(2021)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("ALS_Rheenen", "ALS, van Rheenen et al.(2021)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("ASD", "ASD, Grove et al.(2019)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("^BD", "BD, Ruderfer et al.(2018)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("BMI", "BMI, Locke et al.(2015)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("Height", "Height, Wood et al.(2014)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("IBD", "IBD, Liu et al.(2015)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("MDD", "MDD, Wray et al.(2018)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("MS", "MS, Sawcer et al.(2011)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("RA", "RA, Okada et al.(2014)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("SCZ_Pardinas", "SCZ, Pardiñas et al.(2018)", levels(datafr$GWAS))
levels(datafr$GWAS) <- gsub("SLE", "SLE, Julià et al.(2018)", levels(datafr$GWAS))


# Copy functional category column
datafr$Functional_category2 <- datafr$Functional_category

# Split original functional category column
datafr <- datafr %>% tidyr::separate(Functional_category, c("cell.type", "group"), sep='_', extra='drop')

# Replace values in group column
datafr$group <- datafr$group %>% tidyr::replace_na('All')

datafr$group <- factor(datafr$group, levels = c("All", "common", "AD", "Ctrl"))
levels(datafr$group) <- gsub("common", "Common", levels(datafr$group))
levels(datafr$group) <- gsub("Ctrl", "Control", levels(datafr$group))


# Create new columns with logical value for multitest correction significance
datafr$sig_holm <- datafr$Enrichment_holm < 0.05 & datafr$Coefficient_holm < 0.05
datafr$sig_BH <- datafr$Enrichment_BH < 0.05 & datafr$Coefficient_BH < 0.05


# Order dataframe to create matrix later
datafr <- datafr[with(datafr,
                      order(factor(datafr$GWAS, 
                                   levels = c("AD, Kunkle et al.(2019)",
                                              "AD, Jansen et al.(2019)",
                                              "AD, Schwartzentruber et al.(2021)",
                                              "AD, Wightman et al.(2021)",
                                              "AD, Bellenguez et al.(2022)",
                                              "MS, Sawcer et al.(2011)",
                                              "ALS, van Rheenen et al.(2021)",
                                              "ASD, Grove et al.(2019)",
                                              "MDD, Wray et al.(2018)",
                                              "BD, Ruderfer et al.(2018)",
                                              "SCZ, Pardiñas et al.(2018)",
                                              "IBD, Liu et al.(2015)",
                                              "RA, Okada et al.(2014)",
                                              "SLE, Julià et al.(2018)",
                                              "BMI, Locke et al.(2015)",
                                              "Height, Wood et al.(2014)")),
                            datafr$cell.type,
                            datafr$group)),]

# Create new dataframes with coefficient z-scores and enrichment to generate matrices
Z_score <- datafr[,c(12:13,11,15,20,23:24)]

Enrichment <- datafr[,c(12:13,6,16,21,23:24)]

# Z-score dataframes in wide format for matrices
Z_score  <- dcast(setDT(Z_score), Functional_category2 ~ GWAS, value.var = c("cell.type", 
                                                                             "group", 
                                                                             "Coefficient_z.score",
                                                                             "Coefficient_holm",
                                                                             "Coefficient_BH"))


Z_score <- Z_score[,c(1,4,3,5:6,2,14,7:8,13,9,16,12,15,17,10:11,
                      20,19,21:22,18,30,23:24,29,25,32,28,31,33,26:27,
                      36,35,37:38,34,46,39:40,45,41,48,44,47,49,42:43,
                      52,51,53:54,50,62,55:56,61,57,64,60,63,65,58:59,
                      68,67,69:70,66,78,71:72,77,73,80,76,79,81,74:75)]

# Enrichment dataframes in wide format for matrices
Enrichment <- dcast(setDT(Enrichment), Functional_category2 ~ GWAS, value.var = c("cell.type",
                                                                                  "group",
                                                                                  "Enrichment",
                                                                                  "Enrichment_holm",
                                                                                  "Enrichment_BH"))

Enrichment <- Enrichment[,c(1,4,3,5:6,2,14,7:8,13,9,16,12,15,17,10:11,
                            20,19,21:22,18,30,23:24,29,25,32,28,31,33,26:27,
                            36,35,37:38,34,46,39:40,45,41,48,44,47,49,42:43,
                            52,51,53:54,50,62,55:56,61,57,64,60,63,65,58:59,
                            68,67,69:70,66,78,71:72,77,73,80,76,79,81,74:75)]

# Rename coefficient z scores columns
names(Z_score)[names(Z_score) == 'Coefficient_z.score_AD, Bellenguez et al.(2022)'] <- 'AD, Bellenguez et al.(2022)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_AD, Jansen et al.(2019)'] <- 'AD, Jansen et al.(2019)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_AD, Kunkle et al.(2019)'] <- 'AD, Kunkle et al.(2019)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_AD, Schwartzentruber et al.(2021)'] <- 'AD, Schwartzentruber et al.(2021)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_AD, Wightman et al.(2021)'] <- 'AD, Wightman et al.(2021)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_ALS, van Rheenen et al.(2021)'] <- 'ALS, van Rheenen et al.(2021)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_ASD, Grove et al.(2019)'] <- 'ASD, Grove et al.(2019)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_BD, Ruderfer et al.(2018)'] <- 'BD, Ruderfer et al.(2018)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_BMI, Locke et al.(2015)'] <- 'BMI, Locke et al.(2015)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_Height, Wood et al.(2014)'] <- 'Height, Wood et al.(2014)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_IBD, Liu et al.(2015)'] <- 'IBD, Liu et al.(2015)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_MDD, Wray et al.(2018)'] <- 'MDD, Wray et al.(2018)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_MS, Sawcer et al.(2011)'] <- 'MS, Sawcer et al.(2011)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_RA, Okada et al.(2014)'] <- 'RA, Okada et al.(2014)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_SCZ, Pardiñas et al.(2018)'] <- 'SCZ, Pardiñas et al.(2018)'
names(Z_score)[names(Z_score) == 'Coefficient_z.score_SLE, Julià et al.(2018)'] <- 'SLE, Julià et al.(2018)'

# Rename Enrichmnet columns
names(Enrichment)[names(Enrichment) == 'Enrichment_AD, Bellenguez et al.(2022)'] <- 'AD, Bellenguez et al.(2022)'
names(Enrichment)[names(Enrichment) == 'Enrichment_AD, Jansen et al.(2019)'] <- 'AD, Jansen et al.(2019)'
names(Enrichment)[names(Enrichment) == 'Enrichment_AD, Kunkle et al.(2019)'] <- 'AD, Kunkle et al.(2019)'
names(Enrichment)[names(Enrichment) == 'Enrichment_AD, Schwartzentruber et al.(2021)'] <- 'AD, Schwartzentruber et al.(2021)'
names(Enrichment)[names(Enrichment) == 'Enrichment_AD, Wightman et al.(2021)'] <- 'AD, Wightman et al.(2021)'
names(Enrichment)[names(Enrichment) == 'Enrichment_ALS, van Rheenen et al.(2021)'] <- 'ALS, van Rheenen et al.(2021)'
names(Enrichment)[names(Enrichment) == 'Enrichment_ASD, Grove et al.(2019)'] <- 'ASD, Grove et al.(2019)'
names(Enrichment)[names(Enrichment) == 'Enrichment_BD, Ruderfer et al.(2018)'] <- 'BD, Ruderfer et al.(2018)'
names(Enrichment)[names(Enrichment) == 'Enrichment_BMI, Locke et al.(2015)'] <- 'BMI, Locke et al.(2015)'
names(Enrichment)[names(Enrichment) == 'Enrichment_Height, Wood et al.(2014)'] <- 'Height, Wood et al.(2014)'
names(Enrichment)[names(Enrichment) == 'Enrichment_IBD, Liu et al.(2015)'] <- 'IBD, Liu et al.(2015)'
names(Enrichment)[names(Enrichment) == 'Enrichment_MDD, Wray et al.(2018)'] <- 'MDD, Wray et al.(2018)'
names(Enrichment)[names(Enrichment) == 'Enrichment_MS, Sawcer et al.(2011)'] <- 'MS, Sawcer et al.(2011)'
names(Enrichment)[names(Enrichment) == 'Enrichment_RA, Okada et al.(2014)'] <- 'RA, Okada et al.(2014)'
names(Enrichment)[names(Enrichment) == 'Enrichment_SCZ, Pardiñas et al.(2018)'] <- 'SCZ, Pardiñas et al.(2018)'
names(Enrichment)[names(Enrichment) == 'Enrichment_SLE, Julià et al.(2018)'] <- 'SLE, Julià et al.(2018)'

# Copy datafr to order heatmap rows
datafr2 <- datafr

# Order heatmap rows
datafr2_AD <- subset(datafr2, datafr2$GWAS=="AD, Jansen et al.(2019)")

# Match Z-score and enrichment row orders with the datafr2_AD for heatmap 
Z_score <- Z_score[match(datafr2_AD$Functional_category2, Z_score$Functional_category2),]

Enrichment <- Enrichment[match(datafr2_AD$Functional_category2, Enrichment$Functional_category2),]

# Convert Z-score and Enrichment DF to matrix for heatmap
Z_score_matrix <- data.matrix(Z_score[,34:49])
rownames(Z_score_matrix) = Z_score$`cell.type_AD, Kunkle et al.(2019)`
colnames(Z_score_matrix) =  colnames(Z_score[1,34:49])

Enrichment_matrix <- data.matrix(Enrichment[,34:49])
rownames(Enrichment_matrix) = Z_score$`cell.type_AD, Kunkle et al.(2019)`
colnames(Enrichment_matrix) =  colnames(Z_score[1,34:49])

# Create matrices to annotate heatmaps with significance

# Match functional category order of the all, distal and promoter DFs for matrix sig 
matrix_sig_holm <- matrix(datafr2$sig_holm, 24, 16)
matrix_sig_BH <- matrix(datafr2$sig_BH, 24, 16)

matrix_holm_sig <- matrix(datafr2$Enrichment_holm, 24, 16)
matrix_BH_sig <- matrix(datafr2$Enrichment_BH, 24, 16)

rownames(matrix_sig_holm) = Z_score$`cell.type_AD, Jansen et al.(2019)`

rownames(matrix_sig_BH) = Z_score$`cell.type_AD, Jansen et al.(2019)`

# Make annotations to characterize the heatmap 
block = rowAnnotation(foo = anno_block(gp = gpar(fill = c("Ast"="darkgoldenrod1",
                                                          "Exc"="cornflowerblue",
                                                          "Inh"="seagreen3",
                                                          "Mic"="mediumorchid3",
                                                          "Oli"="coral3",
                                                          "OPC"="firebrick")),
                                       labels = c("Ast","Exc","Inh", 
                                                  "Mic","Oli",
                                                  "OPC"), 
                                       labels_gp = gpar(col = "white", fontsize = 9)),
                      group = Z_score$`group_AD, Jansen et al.(2019)`,
                      gp = gpar(col = "black"), show_annotation_name = FALSE,
                      col = list(group = c("All" = "black",
                                           "Common" = "grey64",
                                           "AD" = "red",
                                           "Control" = "blue")))

top_block = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("black", "grey", "red4","blue4")),
                                               labels = c("Alzheimer", "Brain-related", 
                                                          "IMID", "Other"), 
                                               labels_gp = gpar(col = "white", fontsize = 9)))


# Size feature links in bp
size <- read.table("/Users/ivanrodriguez/Documents/Ivan/Myers_Lab/AD_supp/20220511_LDSC/scAD_all_links_npeaks_size.txt")
colnames(size) <- c("N_links","Functional_category2")

size$Functional_category2 <- as.factor(size$Functional_category2)

levels(size$Functional_category2) <- gsub("/cluster/home/irodriguez/Myers_Lab/AD_supp/LDSC/bed_all_links_filtered_files/", "", levels(size$Functional_category2))
levels(size$Functional_category2) <- gsub("scAD_all_links_", "", levels(size$Functional_category2))
levels(size$Functional_category2) <- gsub(".bed", "", levels(size$Functional_category2))

size <- subset(size, size$N_links> 0)
size <- subset(size, size$Functional_category2 !=  "total")

size <- droplevels(size)

size$N_links_log2 <- log2(size$N_links)

size <- size[match(datafr2_AD$Functional_category2, size$Functional_category2),]

bar = rowAnnotation(n.links.log2 = anno_barplot(size$N_links_log2, 
                                               fontsize = 5),
                    annotation_name_rot = 90)


# Split heatmap by cell type and group
split<- Z_score$`cell.type_AD, Jansen et al.(2019)`



# Split heatmap by GWAS
GWAS <- as.data.frame(unique(datafr$GWAS))

GWAS <- GWAS %>%
  mutate(GWAS_type = case_when(str_detect(unique(datafr$GWAS), "AD," ) ~ "Alzheimer",
                               str_detect(unique(datafr$GWAS), "^BD," ) ~ "Brain-related",
                               str_detect(unique(datafr$GWAS), "SCZ," ) ~ "Brain-related",
                               str_detect(unique(datafr$GWAS), "MDD," ) ~ "Brain-related",
                               str_detect(unique(datafr$GWAS), "ASD," ) ~ "Brain-related",
                               str_detect(unique(datafr$GWAS), "ALS," ) ~ "Brain-related",
                               str_detect(unique(datafr$GWAS), "MS," ) ~ "Brain-related",
                               str_detect(unique(datafr$GWAS), "RA" ) ~ "IMID",
                               str_detect(unique(datafr$GWAS), "IBD" ) ~ "IMID",
                               str_detect(unique(datafr$GWAS), "SLE" ) ~ "IMID",
                               TRUE ~ "Other"))

split_GWAS <- GWAS$GWAS_type 

# Heatmap for coefficient Z score annotates with enrichment value if it is significant
# with BH correction
# All peaks
heatmap <- Heatmap(Z_score_matrix, name = "Z_score",
                   col = colorRamp2(breaks = c(-3.25, 0, 3.25), c("blue", "white", "red"), transparency = .2),
                   heatmap_legend_param = list(
                     title = "Z_score", at = c(-6,-4,-2, 0, 2, 4, 6)),
                   cluster_rows = FALSE,
                   show_column_dend = FALSE,
                   show_row_names = F,
                   row_order = rownames(Z_score$`cell.type_AD, Kunkle et al.(2019)`),
                   column_order = colnames(Z_score_matrix),
                   column_title = "Feature links",
                   column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                   row_names_gp = gpar(fontsize = 11),
                   show_column_names = T,
                   height = unit(7.5, "cm"),
                   width = unit(7, "cm"),
                   border = TRUE,
                   rect_gp = gpar(col = "black", lwd = 1),
                   row_split =split,
                   row_title = NULL,
                   row_title_rot = 0,
                   column_split = split_GWAS,
                   left_annotation = block,
                   top_annotation = top_block,
                   right_annotation = bar,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if(matrix_sig_BH[i, j] == TRUE)
                       grid.text(sprintf("%.0f", Enrichment_matrix[i, j]), x, y, 
                                 gp = gpar(fontsize = 6, col = "white", fontface = "bold"))
                   })


draw(heatmap)

