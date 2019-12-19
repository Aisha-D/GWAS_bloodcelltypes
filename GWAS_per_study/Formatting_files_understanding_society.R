## This script takes in the genotype and methylation data to estimate cell types (phenotype files)
## And the covariates (age,sex,pcs) for PLINK GWAS
setwd('/mnt/data1/GWAS_bloodcelltypes/Scripts/')
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
library(snpStats)
library(factoextra)
load("/mnt/data1/EPICQC/UnderstandingSociety/US_Betas_Pheno.rda")


#############       Phenotype txt file for PLINK          ##################
#Remove old cell type calculations
pheno <- pheno[,-which(names(pheno) %in% c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"))]

#Cell counts
source("/mnt/data1/reference_files/BloodCellPropCalc/cellpropfunctions.r") #loads functions
load('/mnt/data1/reference_files/BloodCellPropCalc/Bloodcoefs_withNeu.rdat') #loads the bloodcoefs to calculate cell proportion
library(genefilter)
library(quadprog)
library(matrixStats)

betas <- dat #or if it is a betas matrix 
#here we only use the probes required to estiamte cell counts
coefs2 <- coefs[rownames(coefs) %in% rownames(betas),]
betas2 <- betas[rownames(coefs2),]
counts <- projectCellType(Y = betas2, coefCellType = coefs2)
rownames(counts) <- colnames(betas2)
save(counts, file = '/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/cellcounts.rdat')
boxplot(counts, las = 2)


#############       Covariates txt file for PLINK         ###################
pathM <- paste("/mnt/data1/EPICQC/UnderstandingSociety/Genotypes/Imputed/data_filtered_1", c(".bed", ".bim", ".fam"), sep = "")
##This contains European samples only^
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
fam <- SNP_M$fam #Take the fam file to add the IID and FID to covaraites table and make sure the methylation and genotype files match

dim(fam) #1111 samples genotype info
dim(pheno) #1175 samples methylation info

#Drop some pheno files which exist in the fam file
fam2 <- separate(data = fam, col = member,
                 into = c("ID", "ID2"), sep="_") #the ID is repeated so remove that
pheno <- pheno[which(pheno$Essex.Reference %in% fam2$ID),] #this drops the number of samples to 1111
rownames(pheno) <- pheno$barcode
betas <- betas[,which(colnames(betas) %in% rownames(pheno))]
rm(fam2)

#read in the PCA eigenvectors and visualise any samples that are extreme
pca <- read.table("/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/data_filtered_1.pca.eigenvec", header = T)
library(bigpca)
pcs <-pca[,3:22]
pdf("/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/pca.pdf")
plot(pcs[,1], pcs[,2], xlab = "PC1", ylab = "PC2")
dev.off()

#cbind the pheno,cellcounts,pca
#rename the sample IDs in the first two columns of pca file
pca <- separate(data = pca, col = X312_312,
                 into = c("ID", "ID2"), sep="_")
#match the number of samples from pheno and cellcoutns
counts <- counts[which(rownames(counts) %in% rownames(pheno)),]
cov <- cbind(pheno, counts[match(rownames(counts), rownames(pheno)),])
cov <- cbind(cov, pca[match(cov$Essex.Reference, pca$ID), 4:13])
#remove columns which are not needed
cov <- cov[ , -which(names(cov) %in% c("Unique.Tube.Barcode","Rack.Barcode" ,         
                                          "Destination.Well","MethArray_ChipPosition",
                                          "MethArray_Chip","barcode","bloodprocess"))]
#MAKE the IID and FID columns match the fam style with _
cov$IID <- paste(cov$Essex.Reference, cov$Essex.Reference, sep = "_")
cov$FID <- paste(cov$Essex.Reference, cov$Essex.Reference, sep = "_")
cov <- cov[,-1]
cov <- cov[,c("FID","IID","confage","nsex","CD8T","CD4T","NK" ,          
                "Bcell","Mono","Gran","Eos","Neu" ,         
                "X.0.0312902","X0.0133535","X.0.0100221","X0.0341959","X.0.0570412",  
                "X.0.00131213","X0.00949711","X0.0419407","X.0.00691399","X.0.000153751"
                )]
colnames(cov)<-c("FID","IID", "Age", "Sex", "CD8T","CD4T", "NK",
                  "Bcell","Mono","Gran",
                  "Eos", "Neu", paste("PC", 1:10, sep = ""))
cov2 <- cov[,c("FID","IID","Age","Sex", paste("PC", 1:10, sep = ""))]
#covariates file
cov2$Age <- as.numeric(as.character(cov2$Age))
cov2$Sex <- as.numeric(cov2$Sex)
write.table(cov2, '/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/cov_agesexpcs.txt',sep = "\t", quote = FALSE, row.names = F)

#pheno file for cells
cellphen <- cov[,c("FID","IID", "CD8T","CD4T", "NK",
               "Bcell","Mono","Gran",
               "Eos", "Neu")]
write.table(cellphen, '/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/cell_phenotype.txt',sep = "\t", quote = FALSE, row.names = F)

#sex update column
sexupd <- cov[,c("FID","IID","Sex")]
sexupd$Sex[is.na(sexupd$Sex)] <- 0
write.table(sexupd, '/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/US_sex_update.txt',sep = "\t", quote = FALSE, row.names = F)
