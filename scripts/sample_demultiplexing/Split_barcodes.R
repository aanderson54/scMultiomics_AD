
library(Seurat)
library(sctransform)
library(dplyr)
library(scDblFinder)
library(biomaRt) 
library(scater)
library(SingleCellExperiment)
library(SeuratDisk)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(scales)
library(future) 
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2)



###################################################################################

gtf<-rtracklayer::import("~/myers/scAnalysis/genes.gtf")
chrY<-gtf[seqnames(gtf) =="chrY",]
genes<-unique(chrY$gene_name)
#####################





###################################################################################
data.data <- Read10X(data.dir = "/cluster/home/aanderson/myers/multiomics_cellranger/211206/2/outs/filtered_feature_bc_matrix/")
rna_counts <- data.data$`Gene Expression`
atac_counts<-data.data$Peaks
data2 <- CreateSeuratObject(counts = rna_counts,project="chrY")
data2<-NormalizeData(data2)
tmp_assay<-subset(data2[["RNA"]],features =genes)
Key(tmp_assay)<-"y_"
data2[["Y"]]<-tmp_assay
#####################
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% "chrY"
atac_counts <- atac_counts[as.vector(grange.use), ]
Y_assay<-CreateAssayObject(atac_counts)
data2[["ATAC_Y"]] <- Y_assay
#####################
md2<-data2@meta.data
md2$ratio<-md2$nCount_Y/md2$nCount_RNA
png("~/myers/scAnalysis/211206/sam2_atacY_rnaY.png",width=4,height=4,units="in",res=500,type="cairo")
ggplot(md2, aes(x=nCount_Y, y=nCount_ATAC_Y, color=nFeature_RNA))+geom_jitter(size=0.3)+ylim(0,40)+xlim(0,300)
dev.off()
md2$ratio<-md2$nFeature_Y/md2$nFeature_RNA
md2$sex2<-ifelse(md2$nFeature_RNA<500 & md2$nCount_Y<5, "unassigned",
                 ifelse(md2$ratio>0.0005,"M","F")) #this ratio doesn't stand
png("~/myers/scAnalysis/211206/sam2_rnaY_cellsnp.png",width=4,height=4,units="in",res=500,type="cairo")
ggplot(md2, aes(y=nCount_Y, x=nCount_RNA, color=donor_id))+geom_point(size=0.1)+xlim(0,300000)+ylim(0,300)
dev.off()



vireo<-read.table("/cluster/home/aanderson/myers/cellsnp/211206/2/vireo_2_1000G/donor_ids.tsv")
colnames(vireo)<-vireo[1,]
vireo<-vireo[-1,]
rownames(vireo)<-vireo$cell
data2<-AddMetaData(data2, vireo)

aggregate(nCount_Y~donor_id, data2@meta.data, quantile) 
#F: BNT1261 Ctrl
#M: BEB18034 PD
barcodes<-gsub("1","2",colnames(data2))
f_barcodes<-barcodes[data2$donor_id=="donor1"]
m_barcodes<-barcodes[data2$donor_id=="donor0"]

write.csv(f_barcodes, "~/myers/scAnalysis/211206/sample2_barcodes_BNT1261.csv")
write.csv(m_barcodes, "~/myers/scAnalysis/211206/sample2_barcodes_BEB18034.csv")

###################################################################################

###################################################################################
data.data <- Read10X(data.dir = "/cluster/home/aanderson/myers/multiomics_cellranger/211206/3/outs/filtered_feature_bc_matrix/")
rna_counts <- data.data$`Gene Expression`
atac_counts<-data.data$Peaks
data3 <- CreateSeuratObject(counts = rna_counts,project="chrY")
data3<-NormalizeData(data3)
tmp_assay<-subset(data3[["RNA"]],features =genes)
Key(tmp_assay)<-"y_"
data3[["Y"]]<-tmp_assay
#####################
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% "chrY"
atac_counts <- atac_counts[as.vector(grange.use), ]
Y_assay<-CreateAssayObject(atac_counts)
data3[["ATAC_Y"]] <- Y_assay
#####################
md3<-data3@meta.data
md3$ratio<-md3$nCount_Y/md3$nCount_RNA
png("~/myers/scAnalysis/211206/sam3_atacY_rnaY.png",width=4,height=4,units="in",res=500,type="cairo")
ggplot(md3, aes(x=nCount_Y, y=nCount_ATAC_Y, color=nFeature_RNA))+geom_jitter(size=0.3)
dev.off()
md3$ratio<-md3$nFeature_Y/md3$nFeature_RNA
md3$sex2<-ifelse(md3$nFeature_RNA<500 & md3$nCount_Y<5, "unassigned",
                 ifelse(md3$ratio>0.0005,"M","F")) #this ratio doesn't stand
png("~/myers/scAnalysis/211206/sam3_rnaY_cellsnp.png",width=4,height=4,units="in",res=500,type="cairo")
ggplot(md3, aes(y=nCount_Y, x=nCount_RNA, color=donor_id))+geom_point(size=0.1)
dev.off()





vireo<-read.table("/cluster/home/aanderson/myers/cellsnp/211206/3/vireo_3_1000G/donor_ids.tsv")
colnames(vireo)<-vireo[1,]
vireo<-vireo[-1,]

rownames(vireo)<-vireo$cell
data3<-AddMetaData(data3, vireo)

data3$ID<-ifelse(data3$prob_max>0.9, data3$best_singlet, data3$donor_id)

aggregate(nCount_Y~ID, data3@meta.data, quantile) 



#F: BEB19157  PD
#M: 1230    ctrl
barcodes<-gsub("1","3",colnames(data3))

barcodes_BEB19157<-barcodes[data3$donor_id=="donor1"]  #female
barcodes_1230<-barcodes[data3$donor_id=="donor0"] #male

write.csv(barcodes_BEB19157, "~/myers/scAnalysis/211206/sample3_barcodes_BEB19157.csv")
write.csv(barcodes_1230, "~/myers/scAnalysis/211206/sample3_barcodes_1230.csv")

###################################################################################


###################################################################################
data.data <- Read10X(data.dir = "/cluster/home/aanderson/myers/multiomics_cellranger/211206/4/outs/filtered_feature_bc_matrix/")
rna_counts <- data.data$`Gene Expression`
atac_counts<-data.data$Peaks
data4 <- CreateSeuratObject(counts = rna_counts,project="chrY")
data4<-NormalizeData(data4)
genes<-pos[which(pos$chromosome_name=="Y"),]$hgnc_symbol
tmp_assay<-subset(data4[["RNA"]],features =genes)
Key(tmp_assay)<-"y_"
data4[["Y"]]<-tmp_assay
#####################
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% "chrY"
atac_counts <- atac_counts[as.vector(grange.use), ]
Y_assay<-CreateAssayObject(atac_counts)
data4[["ATAC_Y"]] <- Y_assay
data4.2<-subset(data4, subset = nFeature_RNA < 10000 & nFeature_RNA > 200)
#####################
md4<-data4@meta.data
md4$sex<-ifelse(md4$nCount_Y>1,"M","F")
md4$ratio<-md4$nFeature_Y/md4$nFeature_RNA
md4$sex2<-ifelse(md4$nFeature_RNA<500 & md4$nCount_Y<5, "unassigned",
                 ifelse(md4$ratio>0.0005,"M","F"))
png("~/myers/scAnalysis/211206/sam4_rnaY.png",width=4,height=4,units="in",res=500,type="cairo")
ggplot(md4, aes(y=nCount_Y, x=nCount_RNA, color=sex2))+geom_point(size=0.1)
dev.off()





vireo<-read.table("/cluster/home/aanderson/myers/cellsnp/211206/4/vireo_4_1000G/donor_ids.tsv")
colnames(vireo)<-vireo[1,]
vireo<-vireo[-1,]
rownames(vireo)<-vireo$cell
data4<-AddMetaData(data4, vireo)

aggregate(nCount_Y~donor_id, data4@meta.data, quantile) 
#F: 3329  AD
#M: HBQS  ALS
barcodes<-gsub("1","4",colnames(data4))
f_barcodes<-barcodes[data4$donor_id=="donor1"]
m_barcodes<-barcodes[data4$donor_id=="donor0"]

write.csv(f_barcodes, "~/myers/scAnalysis/211206/sample4_barcodes_3329.csv")
write.csv(m_barcodes, "~/myers/scAnalysis/211206/sample4_barcodes_HBQS.csv")

###################################################################################


###################################################################################
data.data <- Read10X(data.dir = "/cluster/home/aanderson/myers/multiomics_cellranger/211206/5/outs/filtered_feature_bc_matrix/")
rna_counts <- data.data$`Gene Expression`
atac_counts<- data.data$Peaks
data5 <- CreateSeuratObject(counts = rna_counts,project="chrY")
data5<-NormalizeData(data5)
tmp_assay<-subset(data5[["RNA"]],features =genes)
Key(tmp_assay)<-"y_"
data5[["Y"]]<-tmp_assay
#####################
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% "chrY"
atac_counts <- atac_counts[as.vector(grange.use), ]
Y_assay<-CreateAssayObject(atac_counts)
data5[["ATAC_Y"]] <- Y_assay
#####################
md5<-data5@meta.data
md5$ratio<-md5$nCount_Y/md5$nCount_RNA
md5$sex2<-ifelse(md5$nFeature_RNA<500 & md5$nCount_Y<5, "unassigned",
                 ifelse(md5$ratio>0.001,"M","F"))

png("~/myers/scAnalysis/211206/sam5_rnaY_cellSNP.png",width=4,height=4,units="in",res=600,type="cairo")
ggplot(md5, aes(y=nCount_Y, x=nCount_RNA, color=donor_id))+geom_point(size=0.3)
dev.off()




vireo<-read.table("/cluster/home/aanderson/myers/cellsnp/211206/5/vireo_5_1000G/donor_ids.tsv")
colnames(vireo)<-vireo[1,]
vireo<-vireo[-1,]
rownames(vireo)<-vireo$cell
data5<-AddMetaData(data5, vireo)

aggregate(nCount_Y~donor_id, data5@meta.data, quantile) 
#F: NT1271    Ctrl 
#M: BEB18062  ALS
barcodes<-gsub("1","5",colnames(data5))
f_barcodes<-barcodes[data5$donor_id=="donor0"]
m_barcodes<-barcodes[data5$donor_id=="donor1"]

write.csv(f_barcodes, "~/myers/scAnalysis/211206/sample5_barcodes_NT1271.csv")
write.csv(m_barcodes, "~/myers/scAnalysis/211206/sample5_barcodes_BEB18062.csv")

###################################################################################



###################################################################################
data.data <- Read10X(data.dir = "/cluster/home/aanderson/myers/multiomics_cellranger/211206/9/outs/filtered_feature_bc_matrix/")
rna_counts <- data.data$`Gene Expression`
atac_counts<-data.data$Peaks
data9 <- CreateSeuratObject(counts = rna_counts,project="chrY")
data9<-NormalizeData(data9)
tmp_assay<-subset(data9[["RNA"]],features =genes)
Key(tmp_assay)<-"y_"
data9[["Y"]]<-tmp_assay
#####################
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% "chrY"
atac_counts <- atac_counts[as.vector(grange.use), ]
Y_assay<-CreateAssayObject(atac_counts)
data9[["ATAC_Y"]] <- Y_assay
#####################
md9<-data9@meta.data
md9$ratio<-md9$nCount_Y/md9$nCount_RNA
md9$sex2<-ifelse(md9$nFeature_RNA<500 & md9$nCount_Y<5, "unassigned",
                 ifelse(md9$ratio>0.0005,"M","F"))
png("~/myers/scAnalysis/211206/sam9_rnaY_cellsnp.png",width=4,height=4,units="in",res=500,type="cairo")
ggplot(md9, aes(y=nCount_Y, x=nCount_RNA, color=donor_id))+geom_jitter(size=0.3)
dev.off()





vireo<-read.table("/cluster/home/aanderson/myers/cellsnp/211206/9/vireo_9_1000G/donor_ids.tsv")
colnames(vireo)<-vireo[1,]
vireo<-vireo[-1,]
rownames(vireo)<-vireo$cell
data9<-AddMetaData(data9, vireo)

aggregate(nCount_Y~donor_id, data9@meta.data, quantile) 
#F: 4313 AD
#M: 4482   AD
barcodes<-gsub("1","9",colnames(data9))
f_barcodes<-barcodes[which(data9$donor_id=="donor1")]
m_barcodes<-barcodes[which(data9$donor_id=="donor0")]

write.csv(f_barcodes, "~/myers/scAnalysis/211206/sample9_barcodes_4313.csv")
write.csv(m_barcodes, "~/myers/scAnalysis/211206/sample9_barcodes_4482.csv")

###################################################################################

###################################################################################
data.data <- Read10X(data.dir = "/cluster/home/aanderson/myers/multiomics_cellranger/211206/10/outs/filtered_feature_bc_matrix/")
rna_counts <- data.data$`Gene Expression`
atac_counts<-data.data$Peaks
data10 <- CreateSeuratObject(counts = rna_counts,project="chrY")
data10<-NormalizeData(data10)
tmp_assay<-subset(data10[["RNA"]],features =genes)
Key(tmp_assay)<-"y_"
data10[["Y"]]<-tmp_assay
#####################
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% "chrY"
atac_counts <- atac_counts[as.vector(grange.use), ]
Y_assay<-CreateAssayObject(atac_counts)
data10[["ATAC_Y"]] <- Y_assay
#####################
md10<-data10@meta.data
md10$ratio<-md10$nCount_Y/md10$nCount_RNA
md10$sex2<-ifelse(md10$nFeature_RNA<500 & md10$nCount_Y<5, "unassigned",
                  ifelse(md10$ratio>0.0005,"M","F"))
png("~/myers/scAnalysis/211206/sam10_rnaY_cellsnp.png",width=4,height=4,units="in",res=500,type="cairo")
ggplot(md10, aes(y=nCount_Y, x=nCount_RNA, color=donor_id))+geom_jitter(size=0.3)
dev.off()



vireo<-read.table("/cluster/home/aanderson/myers/cellsnp/211206/10/vireo_10_1000G/donor_ids.tsv")
colnames(vireo)<-vireo[1,]
vireo<-vireo[-1,]
rownames(vireo)<-vireo$cell
data10<-AddMetaData(data10, vireo)

aggregate(nCount_Y~donor_id, data10@meta.data, quantile) 
#F: HCT17HEX Ctrl
#M: HCTZZT  Ctrl
barcodes<-gsub("1","10",colnames(data10))
f_barcodes<-barcodes[which(data10$donor_id=="donor0")]
m_barcodes<-barcodes[which(data10$donor_id=="donor1")]

write.csv(f_barcodes, "~/myers/scAnalysis/211206/sample10_barcodes_HCT17HEX.csv")
write.csv(m_barcodes, "~/myers/scAnalysis/211206/sample10_barcodes_HCTZZT.csv")

###################################################################################

###################################################################################
data.data <- Read10X(data.dir = "/cluster/home/aanderson/myers/multiomics_cellranger/211206/11/outs/filtered_feature_bc_matrix/")
rna_counts <- data.data$`Gene Expression`
atac_counts <-data.data$Peaks
data11 <- CreateSeuratObject(counts = rna_counts,project="chrY")
data11<-NormalizeData(data11)
tmp_assay<-subset(data11[["RNA"]],features =genes)
Key(tmp_assay)<-"y_"
data11[["Y"]]<-tmp_assay
#####################
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% "chrY"
atac_counts <- atac_counts[as.vector(grange.use), ]
Y_assay<-CreateAssayObject(atac_counts)
data11[["ATAC_Y"]] <- Y_assay
#####################
md11<-data11@meta.data
md11$ratio<-md11$nCount_Y/md11$nCount_RNA
md11$sex2<-ifelse(md11$nFeature_RNA<500 & md11$nCount_Y<5, "unassigned",
                  ifelse(md11$ratio>0.0005,"M","F"))
md11$donor_id[1]<-"unassigned"
png("~/myers/scAnalysis/211206/sam11_rnaY_cellSNP.png",width=4,height=4,units="in",res=500,type="cairo")
ggplot(md11, aes(y=nCount_Y, x=nCount_RNA, color=donor_id))+geom_jitter(size=0.3)
dev.off()




vireo<-read.table("/cluster/home/aanderson/myers/cellsnp/211206/11/vireo_11_1000G/donor_ids.tsv")
colnames(vireo)<-vireo[1,]
vireo<-vireo[-1,]

vireo$ID<-ifelse(vireo$prob_max>0.75, vireo$best_singlet, "unassigned")
rownames(vireo)<-vireo$cell
data11<-AddMetaData(data11, vireo)

aggregate(nCount_Y~donor_id, data11@meta.data, quantile) 
#F: 4443 AD
#M: 4305 AD
barcodes<-gsub("1","11",colnames(data11))
f_barcodes<-barcodes[data11$donor_id=="donor0"]
m_barcodes<-barcodes[data11$donor_id=="donor1"]

write.csv(f_barcodes, "~/myers/scAnalysis/211206/sample11_barcodes_4443.csv")
write.csv(m_barcodes, "~/myers/scAnalysis/211206/sample11_barcodes_4305.csv")

###################################################################################
