### GWAS of different cell types from EXTEND data set
setwd('/mnt/data1/GWAS_bloodcelltypes/EXTEND/')
library(methylumi)
library(wateRmelon)
require(gdata)
library(minfi)
library(ggplot2)
require(gridExtra)
library(plyr)
require(IlluminaHumanMethylationEPICmanifest)
library(dplyr)
library(tidyr)

####################
#### SNP QC
####################
#This was conducted by gemma and the script is located "/mnt/data1/EXTEND/GWAS_bloodcelltypes/BDR_QC.r"
# HOWEVER NEED TO CHECK THROUGH
# library(snpStats)
# pathM <- paste('/mnt/data1/EXTEND/Genotypes/Imputed/EXTEND_Unrelated_EUR_QCd', c(".bed", ".bim", ".fam"), sep = "")
# SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
# Take one bim map (all 3 maps are based on the same ordered set of SNPs)
# map <- SNP_M$map
# fam <- SNP_M$fam
# SNP_M$genotypes
# colnames(map) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")


####################
#### Methylation QC
####################
#This was conducted by gemma and script is located "/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/merging_EXTEND_1_2_passed_samples_withGeno.Rmd"
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EXTEND_batch1_2_genoQCd_Normalised.rdat")
# No RGSET was found
# ! CREATE RGSET FROM ORGININAl --- HOWEVER SAMPLES ONLY CONTAIN PASSED! Run This chink once
# Batch1
# batch1 <- list.files('/mnt/data1/EXTEND/Methylation/idats/batch1')
# batch1 <- as.data.frame(batch1)
# colnames(batch1)[colnames(batch1) =='batch1'] <- 'Basename'
# batch1$Basename <- as.character(batch1$Basename)
# batch1 <- separate(data = batch1, col = Basename,
                      # into = c("Chip", "Chip_Position", "colour"), sep="_")
# batch1$Basename <- paste(batch1$Chip, batch1$Chip_Position, sep="_")
# batch <- unique(batch1$Basename)
# batch <- batch[batch %in% pheno$Basename] #these only reads in the idats that passed QC
# batch1 <- as.data.frame(batch)
# colnames(batch1)[colnames(batch1) =='batch'] <- 'Basename'
# batch1$Basename <- as.character(batch1$Basename)
# idatPath<-c('/mnt/data1/EXTEND/Methylation/idats/batch1')
# RGset_batch1 <- read.metharray.exp(base = idatPath, targets = batch1, force = TRUE)
# save(RGset_batch1, file="EXTEND_batch1_RGset.rdat")
# 
# #Batch2
# batch2 <- list.files('/mnt/data1/EXTEND/Methylation/idats/batch2')
# batch2 <- as.data.frame(batch2)
# colnames(batch2)[colnames(batch2) =='batch2'] <- 'Basename'
# batch2$Basename <- as.character(batch2$Basename)
# batch2 <- separate(data = batch2, col = Basename,
                   # into = c("Chip", "Chip_Position", "colour"), sep="_")
# batch2$Basename <- paste(batch2$Chip, batch2$Chip_Position, sep="_")
# batch <- unique(batch2$Basename)
# batch <- batch[batch %in% pheno$Basename] #these only reads in the idats that passed QC
# batch2 <- as.data.frame(batch)
# colnames(batch2)[colnames(batch2) =='batch'] <- 'Basename'
# batch2$Basename <- as.character(batch2$Basename)
# idatPath<-c('/mnt/data1/EXTEND/Methylation/idats/batch2')
# RGset_batch2 <- read.metharray.exp(base = idatPath, targets = batch2, force = TRUE)
# save(RGset_batch2, file="EXTEND_batch2_RGset.rdat")
#to check that these idats match the pheno, batch1(776) + batch2(246) = 1022(pheno)

#load the rgsets
# load("EXTEND_batch2_RGset.rdat")
# load("EXTEND_batch1_RGset.rdat")
# RGset <- combineArrays(RGset_batch1, RGset_batch2)
# save(RGset, file="EXTEND_RGset.rdat")
# rm(RGset_batch1, RGset_batch2)

#load mset
# load("EXTEND_batch2_mset.rdat")
# load("EXTEND_batch1_mset.rdat")
# mset <- combo(mset_batch1, mset_batch2)
# save(mset, file="EXTEND_mset.rdat")

load("EXTEND_CellCounts.rdat") #this was done in the 'resolving_est_celltype.R'
pcs <- read.table("/mnt/data1/EXTEND/GWAS_bloodcelltypes/EXTEND_Unrelated_EUR_QCd.pca.eigenvec", stringsAsFactors = F, header = T)

Cellcounts <- cbind(pheno$Age, pheno$Sex)
cell_counts <- as.data.frame(cell_counts)
Cellcounts <- cbind(Cellcounts, cell_counts)

#give pcs basename rownames
IID <- pheno$IID[match(pheno$Basename, rownames(Cellcounts))]
#rownames(Cellcounts) <- IID
Cellcounts <- cbind(Cellcounts, IID)
Cellcounts <- as.data.frame(Cellcounts, stringsAsFactors = F)
##Reformat the cov file to include FID, IID, covariates
library(snpStats)
pathM <- paste('/mnt/data1/EXTEND/Genotypes/Imputed/EXTEND_Unrelated_EUR_QCd', c(".bed", ".bim", ".fam"), sep = "")
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
fam <- SNP_M$fam

cov2 <- cbind(Cellcounts, pcs[match(Cellcounts$IID, pcs$X10K06515PH),3:12])
cov2 <- cbind(fam[match(cov2$IID, fam$pedigree),1:2], cov2)
cov2 <- cov2[ , -which(names(cov2) %in% c("IID"))]
colnames(cov2)<-c("FID","IID", "Age", "Sex", "CD8T","CD4T", "NK",
                 "Bcell","Mono","Gran",
                 "Eos", "Neu", paste("PC", 1:10, sep = ""))
cov2$Age <- as.numeric(as.character(cov2$Age))
cov2$Sex <- as.character(cov2$Sex)
cov2$Sex[cov2$Sex == 'Female'] <- 2
cov2$Sex[cov2$Sex == 'Male'] <- 1
cov2$Sex <- as.numeric(cov2$Sex)-1
cov2 <- cov2[!(is.na(cov2$FID)),]
write.table(cov2, 'covariates.txt',sep = "\t", quote = FALSE, row.names = F)
ex_sex <- cov2[,c(1,2,4)]
ex_sex$Sex <- as.numeric(ex_sex$Sex)+1
ex_sex$Sex[is.na(ex_sex$Sex)] <- 0
write.table(ex_sex, 'EXTEND_sex_update.txt',sep = "\t", quote = FALSE, row.names = F)

#pheno file
cell_phn <- cov2[,c("FID","IID", "CD8T","CD4T", "NK",
                      "Bcell","Mono","Gran",
                      "Eos", "Neu")]
write.table(cell_phn, 'cell_phenotype.txt',sep = "\t", quote = FALSE, row.names = F)

#covariates, age,sex,pcs
cell_cov <- cov2[,-which(names(cov2) %in% c("CD8T","CD4T", "NK",
                    "Bcell","Mono","Gran",
                    "Eos", "Neu"))]
write.table(cell_cov, 'cov_agesexpcs.txt',sep = "\t", quote = FALSE, row.names = F)


####### PCA visualisation ######
extend_eur <- read.table("/mnt/data1/EXTEND/GWAS_bloodcelltypes/EXTEND_Unrelated_EUR_QCd.pca.eigenvec", header = T)
library(factoextra)
library(bigpca)
pcs <-extend_eur[,3:22]

pdf("pca.pdf")
plot(pcs[,1], pcs[,2], xlab = "PC1", ylab = "PC2")
dev.off()



