## Summary of cell type distribution in each cohort 
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(patchwork)
library(gtable)

###########################################
#####     Load Cell Count data    #########
###########################################
## Extend
load("/mnt/data1/GWAS_bloodcelltypes/EXTEND/cell_counts.rdat")
extend_cellcounts = as.data.frame(counts)
rm(counts)
## EuGEI
load('/mnt/data1/GWAS_bloodcelltypes/EuGEI/cellcounts.rdat')
eugei_cellcounts = as.data.frame(counts)
rm(counts)
load('/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/cellcounts.rdat')
us_cellcounts = as.data.frame(counts)
rm(counts)

#############################################
####     Plot Distributions Per Study    ####
#############################################
celltypes = colnames(extend_cellcounts)

pdf("Cell distribution per study.pdf")
g = lapply(1:length(celltypes), function(i) {
  j = celltypes[i]
  ggplot(extend_cellcounts, aes(extend_cellcounts[,j])) +
    geom_histogram() +
    theme(panel.background = element_blank()) +
    labs(title =paste('Extend',j , sep = '_'), x = "", y = "")
})
grid.arrange(grobs = g, ncol = 4)

g = lapply(1:length(celltypes), function(i) {
  j = celltypes[i]
  ggplot(eugei_cellcounts, aes(eugei_cellcounts[,j])) +
    geom_histogram() +
    theme(panel.background = element_blank()) +
    labs(title =paste('EuGEI',j , sep = '_'), x = "", y = "")
})
grid.arrange(grobs = g, ncol = 4)

g = lapply(1:length(celltypes), function(i) {
  j = celltypes[i]
  ggplot(us_cellcounts, aes(us_cellcounts[,j])) +
    geom_histogram() +
    theme(panel.background = element_blank()) +
    labs(title =paste('UnderSoc',j , sep = '_'), x = "", y = "")
})
grid.arrange(grobs = g, ncol = 4)
dev.off()

#############################################
####     Plot CellDist per celltype     #####
#############################################
CellDist = function(celltypes = celltypes, dat = dat,
                    studies = studies) {
  #  x = vector of cell types that are in all dataframes
  #  dat = a list of dataframes you want to use
  #  studies =  a vector of names of studies which is listed in the same order
  #             as the dfs
  s = length(studies)
  for ( i in 1:length(celltypes)) {
    j = celltypes[i]
    print(j)
    par(mfrow = c(2, 2))
    for(k in 1:length(studies)){
      plot(density(dat[[k]][,j]), main = paste(studies[k], j, sep= "_"), xlab = "")
    }
  }
}

pdf("Cell Distribution per celltype (density).pdf")
CellDist(celltypes = colnames(extend_cellcounts), 
         dat = list(extend_cellcounts, eugei_cellcounts, us_cellcounts),
         studies = c("Extend", "EuGEI", "UnderSoc"))
dev.off()



#############################################
#       Plot CellDist per celltype          #
#       Overlay the studies on one          #
#       graph                               #
#############################################

##write a function to plot 
#This function makes dataframe for one cell type
CellDens = function(dat = dat, celltype = celltype, study = study){
  
  ## dat = list of dataframes
  ## celltypes = should be one cell type of interest
  ## study = name of studies(same order as list of dataframe)
  
  #create the merged dataframe of one cell type
  celldf = c()
  for(i in 1:length(dat)) {
    df = as.data.frame(dat[i], stringsAsFactors = FALSE)
    studyname = study[i]
    df1 = as.data.frame(cbind(df[,celltype], 
                                  rep(studyname, nrow(df))),
                            stringsAsFactors = FALSE)
    colnames(df1) <- c("CellType","Study")
    celldf <- rbind(celldf, df1)
  }
  celldf$CellType <- as.numeric(celldf$CellType)

  p = ggplot(celldf, aes(x = CellType, fill = Study)) +
    geom_density(alpha = 0.5) +
    labs(x = "", y = "Density",
         title = paste(celltype, " Distribution")) +
           theme(panel.background = element_blank())
  print(p)
  return(p)
}

CellDens(dat = list(extend_cellcounts,eugei_cellcounts,us_cellcounts),
         celltype = 'CD4T',
         study = c("Ex", "Eu","US"))


#lets loop through to plot all types of celltypes
pdf('CellDistribution.pdf', onefile = T)
for(i in 1:ncol(extend_cellcounts)){
  cell = colnames(extend_cellcounts)[i]
  CellDens(dat = list(extend_cellcounts,eugei_cellcounts,us_cellcounts),
           celltype = cell,
           study = c("Ex", "Eu","US"))
}
dev.off()
