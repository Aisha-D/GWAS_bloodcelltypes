## Compare age and sex influence on cell distribution
setwd("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis/")


#Extend
load("/mnt/data1/GWAS_bloodcelltypes/EXTEND/EXTEND_CellCounts.rdat")
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EUR_unlrelated_QCd/EXTEND_batch1_2_genoQCd_Normalised.rdat")

pheno2 = pheno[,c(1,5,6)]
rownames(pheno2) <- pheno2$Basename
pheno2 = pheno2[,-1]

pheno3 = merge(pheno2, cell_counts, by = "row.names")
rownames(pheno3) = pheno3$Basename
pheno3 = pheno3[,3:12]

pdf("EXTEND_AGE.pdf")
par(mfrow = c(2,2)) 
plot(pheno3$Age, pheno3$CD8T)
plot(pheno3$Age, pheno3$CD4T)
plot(pheno3$Age, pheno3$NK)
plot(pheno3$Age, pheno3$Bcell)
plot(pheno3$Age, pheno3$Mono)
plot(pheno3$Age, pheno3$Gran)
plot(pheno3$Age, pheno3$Eos)
plot(pheno3$Age, pheno3$Neu)
dev.off()

scatterplot(pheno3$CD8T ~ pheno3$Age, data = pheno3)

pdf("EXTEND_AGE.pdf")
par(mfrow = c(2,2)) 
plot(pheno3$Age, pheno3$CD8T)
abline(lm(pheno3$CD8T ~ pheno3$Age, data = pheno3), col="red") # regression line (y~x)
plot(pheno3$Age, pheno3$CD4T)
abline(lm(pheno3$CD4T ~ pheno3$Age, data = pheno3), col="red")
plot(pheno3$Age, pheno3$NK)
abline(lm(pheno3$NK ~ pheno3$Age, data = pheno3), col="red")
plot(pheno3$Age, pheno3$Bcell)
abline(lm(pheno3$Bcell ~ pheno3$Age, data = pheno3), col="red")
plot(pheno3$Age, pheno3$Mono)
abline(lm(pheno3$Mono ~ pheno3$Age, data = pheno3), col="red")
plot(pheno3$Age, pheno3$Gran)
abline(lm(pheno3$Gran ~ pheno3$Age, data = pheno3), col="red")
plot(pheno3$Age, pheno3$Eos)
abline(lm(pheno3$Eos ~ pheno3$Age, data = pheno3), col="red")
plot(pheno3$Age, pheno3$Neu)
abline(lm(pheno3$Neu ~ pheno3$Age, data = pheno3), col="red")
dev.off()

library(ggpubr)
ggscatter(pheno3, x = 'CD8T', y = 'Age',
          add = 'reg.line',
          add.params = list(color = 'blue', fill = 'lightgray'),
          conf.int = T) +
  stat_cor(method = 'pearson') +
  labs(x = 'CD8T', y = 'Age')


pdf("EXTEND_SEX.pdf")
par(mfrow = c(2,2)) 
boxplot(pheno3$CD8T ~ pheno3$Sex)
abline(lm(pheno3$CD8T ~ pheno3$Sex, data = pheno3), col="red") # regression line (y~x)
boxboxplot(pheno3$CD4T ~ pheno3$Sex)
abline(lm(pheno3$CD4T ~ pheno3$Sex, data = pheno3), col="red")
boxplot(pheno3$NK ~ pheno3$Sex)
abline(lm(pheno3$NK ~ pheno3$Sex, data = pheno3), col="red")
boxplot(pheno3$Bcell ~ pheno3$Sex)
abline(lm(pheno3$Bcell ~ pheno3$Sex, data = pheno3), col="red")
boxplot(pheno3$Mono ~ pheno3$Sex)
abline(lm(pheno3$Mono ~ pheno3$Sex, data = pheno3), col="red")
boxplot(pheno3$Gran ~ pheno3$Sex)
abline(lm(pheno3$Gran ~ pheno3$Sex, data = pheno3), col="red")
boxplot(pheno3$Eos ~ pheno3$Sex)
abline(lm(pheno3$Eos ~ pheno3$Sex, data = pheno3), col="red")
boxplot(pheno3$Neu ~ pheno3$Sex)
abline(lm(pheno3$Neu ~ pheno3$Sex, data = pheno3), col="red")
dev.off()

library(ggpubr)
ggscatter(pheno3, x = 'CD8T', y = 'Sex',
          add = 'reg.line',fill = 'Sex',
          add.params = list(color = 'blue', fill = 'lightgray'),
          conf.int = T) +
  stat_cor(method = 'pearson') +
  labs(x = 'CD8T', y = 'Sex')

library(ForestPMPlot)
install.packages("ForestPMPlot")
tmp =  matrix(c(1,1,1,1,1,1,
                1,1,0,1,1,1,
                0,0,0,1,0,0,
                0,1,0,0,1,1), 6, 4)
tmp

library(nonpar)
cochrans.q(extend$)
