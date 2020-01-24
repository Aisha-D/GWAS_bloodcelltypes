## Script that generates the coeffs table to create counts table

source("/mnt/data1/reference_files/BloodCellPropCalc/cellpropfunctions.r") #loads functions
load('/mnt/data1/reference_files/BloodCellPropCalc/BloodCoefs.rdat') #loads the bloodcoefs to calculate cell proportion
library(genefilter)
library(quadprog)
library(matrixStats)

############ Create coeffs from sorted blood data
load("/mnt/data1/reference_files/BloodCellPropCalc/450k_reference.rdat") #contains beta data
library(FlowSorted.Blood.450k)
ref_metadata <- as.data.frame(FlowSorted.Blood.450k@colData)
ref_beta <- as.data.frame(getBeta(FlowSorted.Blood.450k))
celltypes = c("Gran","CD4T","CD8T","Bcell","Mono","NK","Neu", "Eos")

samplenames <- ref_metadata[which(ref_metadata$CellType %in% celltypes),'Sample_Name']
# CD14 is marker for monocyte, CD19 is marker for B cell, CD56 is a marker NK cells
ref_beta <- ref_beta[,samplenames]
ref_metadata <- ref_metadata[which(ref_metadata$CellType %in% celltypes),]



# Create coefs 
compositeCellType = "Blood"
processMethod = "auto"
probeSelect = "auto"
cellTypes = c("Gran","CD4T","CD8T","Bcell","Mono","NK","Neu", "Eos")
referencePlatform = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k")
returnAll = FALSE
meanPlot = FALSE
verbose = TRUE
cellInd = ref_metadata
rawbetas= ref_beta

## If using Neu will have to save change p value in the function
compData <- pickCompProbes(rawbetas= ref_beta, cellInd= ref_metadata, 
                           cellTypes = cellTypes, numProbes = 50, 
                           probeSelect = probeSelect, p.value = 1e-5)

coefs <- compData$coefEsts

# save(coefs, file = /mnt/data1/reference_files/BloodCellPropCalc/Bloodcoefs_withNeu.rdat)
