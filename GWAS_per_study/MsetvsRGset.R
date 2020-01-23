### Correlation plot comparing the RGset and Mset of each cell type from each study
library(dplyr)
library(tidyr)
library(genefilter)
library(quadprog)
library(matrixStats)
library(snpStats)
library(gridExtra)
library(ggpubr)

##################### Estimate Cell Count Function ###################
# Internal functions -----------------------------------------------------------

pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50,
                           compositeCellType = compositeCellType,
                           probeSelect = probeSelect) {
  .isMatrixBackedOrStop(mSet)
  splitit <- function(x) {
    split(seq_along(x), x)
  }
  
  p <- getBeta(mSet)
  pd <- as.data.frame(colData(mSet))
  if (!is.null(cellTypes)) {
    if (!all(cellTypes %in% pd$CellType))
      stop("elements of argument 'cellTypes' is not part of ",
           "'mSet$CellType'")
    keep <- which(pd$CellType %in% cellTypes)
    pd <- pd[keep,]
    p <- p[,keep]
  }
  # NOTE: Make cell type a factor
  pd$CellType <- factor(pd$CellType, levels = cellTypes)
  ffComp <- rowFtests(p, pd$CellType)
  prof <- vapply(
    X = splitit(pd$CellType),
    FUN = function(j) rowMeans2(p, cols = j),
    FUN.VALUE = numeric(nrow(p)))
  r <- rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2, -1, 0) + ncol(compTable)] <-
    c("low", "high", "range")
  tIndexes <- splitit(pd$CellType)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  
  if (probeSelect == "any") {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < 1e-5, ]
      yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
      c(rownames(yAny)[seq(numProbes * 2)])
    })
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < 1e-5, ]
      yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
      yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
      c(rownames(yUp)[seq_len(numProbes)],
        rownames(yDown)[seq_len(numProbes)])
    })
  }
  
  trainingProbes <- unique(unlist(probeList))
  p <- p[trainingProbes,]
  
  pMeans <- colMeans2(p)
  names(pMeans) <- pd$CellType
  
  form <- as.formula(
    sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+")))
  phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
  colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
  if (ncol(phenoDF) == 2) {
    # Two group solution
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
  } else {
    # > 2 groups solution
    tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts
  }
  
  list(
    coefEsts = coefEsts,
    compTable = compTable,
    sampleMeans = pMeans)
}

projectCellType <- function(Y, coefCellType, contrastCellType = NULL,
                            nonnegative = TRUE, lessThanOne = FALSE) {
  if (is.null(contrastCellType)) {
    Xmat <- coefCellType
  } else {
    Xmat <- tcrossprod(coefCellType, contrastCellType)
  }
  
  nCol <- dim(Xmat)[2]
  if (nCol == 2) {
    Dmat <- crossprod(Xmat)
    mixCoef <- t(
      apply(Y, 2, function(x) solve(Dmat, crossprod(Xmat, x))))
    colnames(mixCoef) <- colnames(Xmat)
    return(mixCoef)
  } else {
    nSubj <- dim(Y)[2]
    
    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)
    
    if (nonnegative) {
      if (lessThanOne) {
        Amat <- cbind(rep(-1, nCol), diag(nCol))
        b0vec <- c(-1, rep(0, nCol))
      } else {
        Amat <- diag(nCol)
        b0vec <- rep(0, nCol)
      }
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve.QP(
          Dmat = Dmat,
          dvec = crossprod(Xmat[obs,], Y[obs,i]),
          Amat = Amat,
          bvec = b0vec)$sol
      }
    } else {
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
      }
    }
    mixCoef
  }
}

validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]
  
  if (is.null(L.forFstat)) {
    # NOTE: All non-intercept coefficients
    L.forFstat <- diag(sizeModel)[-1,]
    colnames(L.forFstat) <- colnames(xTest)
    rownames(L.forFstat) <- colnames(xTest)[-1]
  }
  
  # Initialize various containers
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()
  
  if (verbose) cat("[validationCellType] ")
  # Loop over each CpG
  for (j in seq_len(M)) {
    # Remove missing methylation values
    ii <- !is.na(Y[j, ])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j,]
    
    if (j %% round(M / 10) == 0 && verbose) cat(".") # Report progress
    
    # Try to fit a mixed model to adjust for plate
    try({
      if (!is.null(modelBatch)) {
        fit <- try(
          lme(modelFix, random = modelBatch, data = pheno[ii, ]))
        # NOTE: If LME can't be fit, just use OLS
        OLS <- inherits(fit, "try-error")
      } else {
        OLS <- TRUE
      }
      
      if (OLS) {
        fit <- lm(modelFix, data = pheno[ii, ])
        fitCoef <- fit$coef
        sigmaResid[j] <- summary(fit)$sigma
        sigmaIcept[j] <- 0
        nClusters[j] <- 0
      } else {
        fitCoef <- fit$coef$fixed
        sigmaResid[j] <- fit$sigma
        sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
        nClusters[j] <- length(fit$coef$random[[1]])
      }
      coefEsts[j,] <- fitCoef
      coefVcovs[[j]] <- vcov(fit)
      
      useCoef <- L.forFstat %*% fitCoef
      useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
    })
  }
  if (verbose) cat(" done\n")
  
  # Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1
  
  # Get P values corresponding to F statistics
  Pval <- 1 - pf(Fstat, sizeModel, degFree)
  
  list(
    coefEsts = coefEsts,
    coefVcovs = coefVcovs,
    modelFix = modelFix,
    modelBatch = modelBatch,
    sigmaIcept = sigmaIcept,
    sigmaResid = sigmaResid,
    L.forFstat = L.forFstat,
    Pval = Pval,
    orderFstat = order(-Fstat),
    Fstat = Fstat,
    nClusters = nClusters,
    nObserved = nObserved,
    degFree = degFree)
}

# Exported functions -----------------------------------------------------------

estimateCellCounts <- function(rgSet, compositeCellType = "Blood",
                               processMethod = "auto", probeSelect = "auto",
                               cellTypes = c("CD8T", "CD4T", "NK", "Bcell",
                                             "Mono", "Gran"),
                               referencePlatform = c(
                                 "IlluminaHumanMethylation450k",
                                 "IlluminaHumanMethylationEPIC",
                                 "IlluminaHumanMethylation27k"),
                               returnAll = FALSE, meanPlot = FALSE,
                               verbose = TRUE, ...) {
  
  # Check inputs
  .isMatrixBackedOrStop(rgSet, "estimateCellCounts")
  .isRGOrStop(rgSet)
  rgSet <- as(rgSet, "RGChannelSet")
  referencePlatform <- match.arg(referencePlatform)
  rgPlatform <- sub(
    "IlluminaHumanMethylation",
    "",
    annotation(rgSet)[which(names(annotation(rgSet)) == "array")])
  platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
  if ((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes)) {
    message("[estimateCellCounts] Consider including 'nRBC' in argument 'cellTypes' for cord blood estimation.\n")
  }
  referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
  subverbose <- max(as.integer(verbose) - 1L, 0L)
  if (!require(referencePkg, character.only = TRUE)) {
    stop(sprintf("Could not find reference data package for compositeCellType '%s' and referencePlatform '%s' (inferred package name is '%s')",
                 compositeCellType, platform, referencePkg))
  }
  data(list = referencePkg)
  referenceRGset <- get(referencePkg)
  if (rgPlatform != platform) {
    rgSet <- convertArray(
      object = rgSet,
      outType = referencePlatform,
      verbose = subverbose)
  }
  if (!"CellType" %in% names(colData(referenceRGset))) {
    stop(sprintf("the reference sorted dataset (in this case '%s') needs to have a phenoData column called 'CellType'"),
         names(referencePkg))
  }
  if (sum(colnames(rgSet) %in% colnames(referenceRGset)) > 0) {
    stop("the sample/column names in the user set must not be in the ",
         "reference data ")
  }
  if (!all(cellTypes %in% referenceRGset$CellType)) {
    stop(sprintf("all elements of argument 'cellTypes' needs to be part of the reference phenoData columns 'CellType' (containg the following elements: '%s')",
                 paste(unique(referenceRGset$cellType), collapse = "', '")))
  }
  if (length(unique(cellTypes)) < 2) {
    stop("At least 2 cell types must be provided.")
  }
  if ((processMethod == "auto") &&
      (compositeCellType %in% c("Blood", "DLPFC"))) {
    processMethod <- "preprocessQuantile"
  }
  if ((processMethod == "auto") &&
      (!compositeCellType %in% c("Blood", "DLPFC"))) {
    processMethod <- "preprocessNoob"
  }
  processMethod <- get(processMethod)
  if ((probeSelect == "auto") && (compositeCellType == "CordBlood")) {
    probeSelect <- "any"
  }
  if ((probeSelect == "auto") && (compositeCellType != "CordBlood")) {
    probeSelect <- "both"
  }
  
  if (verbose) {
    message("[estimateCellCounts] Combining user data with reference ",
            "(flow sorted) data.\n")
  }
  newpd <- DataFrame(
    sampleNames = c(colnames(rgSet), colnames(referenceRGset)),
    studyIndex = rep(
      x = c("user", "reference"),
      times = c(ncol(rgSet), ncol(referenceRGset))),
    stringsAsFactors = FALSE)
  referencePd <- colData(referenceRGset)
  combinedRGset <- combineArrays(
    object1 = rgSet,
    object2 = referenceRGset,
    outType = "IlluminaHumanMethylation450k")
  colData(combinedRGset) <- newpd
  colnames(combinedRGset) <- newpd$sampleNames
  rm(referenceRGset)
  
  if (verbose) {
    message("[estimateCellCounts] Processing user and reference data ",
            "together.\n")
  }
  if (compositeCellType == "CordBlood") {
    # NOTE: Here Shan wants to discard probes that they have decided
    #       shouldn't be used, for example multi-mapping probes. This is
    #       done by only using probes with names in the comptable.
    #       This is kind of ugly, and dataset dependent.
    combinedMset <- processMethod(combinedRGset, verbose = subverbose)
    compTable <- get(paste0(referencePkg, ".compTable"))
    combinedMset <- combinedMset[
      which(rownames(combinedMset) %in% rownames(compTable)),]
  } else {
    combinedMset <- processMethod(combinedRGset)
  }
  rm(combinedRGset)
  
  # Extract normalized reference data
  referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
  colData(referenceMset) <- as(referencePd, "DataFrame")
  mSet <- combinedMset[, combinedMset$studyIndex == "user"]
  colData(mSet) <- as(colData(rgSet), "DataFrame")
  rm(combinedMset)
  
  if (verbose) {
    message("[estimateCellCounts] Picking probes for composition ",
            "estimation.\n")
  }
  compData <- pickCompProbes(
    mSet = referenceMset,
    cellTypes = cellTypes,
    compositeCellType = compositeCellType,
    probeSelect = probeSelect)
  coefs <- compData$coefEsts
  write.table(coefs, 'Extend_coeffs.csv')
  # TODO: Shouldn't be necessary to rm() anything
  rm(referenceMset)
  
  if (verbose) message("[estimateCellCounts] Estimating composition.\n")
  counts <- projectCellType(getBeta(mSet)[rownames(coefs), ], coefs)
  rownames(counts) <- colnames(rgSet)
  
  if (meanPlot) {
    smeans <- compData$sampleMeans
    smeans <- smeans[order(names(smeans))]
    sampleMeans <- c(
      colMeans2(
        x = getBeta(mSet),
        rows = match(rownames(coefs), rownames(mSet))),
      smeans)
    sampleColors <- c(
      rep(1, ncol(mSet)),
      1 + as.numeric(factor(names(smeans))))
    plot(sampleMeans, pch = 21, bg = sampleColors)
    legend("bottomleft",
           c("blood", levels(factor(names(smeans)))),
           col = 1:7,
           pch = 15)
  }
  if (returnAll) {
    return(list(
      counts = counts,
      compTable = compData$compTable,
      
      normalizedData = mSet))
  } else {
    counts
  }
}

# Global variables -------------------------------------------------------------
utils::globalVariables(c("channel"))
# Other hardcoded variables ----------------------------------------------------

.default.27k.annotation  <- "ilmn12.hg19"
.default.450k.annotation <- "ilmn12.hg19"
.default.epic.annotation <- "ilm10b4.hg19"
.metharray.types <- c("IlluminaHumanMethylation450k",
                      "IlluminaHumanMethylationEPIC",
                      "IlluminaHumanMethylation27k")
.seqnames.order.all <- c(paste0("chr", c(1:22, "X", "Y")), "multi", "unmapped")
.seqnames.order <- paste0("chr", c(1:22, "X", "Y"))

# Internal functions -----------------------------------------------------------

logit2 <- function(x) log2(x) - log2(1 - x)

ilogit2 <- function(x) 2^x / (1 + 2^x)

.show.annotation <- function(annotation, indent = "  ") {
  cat("Annotation\n")
  if (length(annotation) == 1) {
    cat(sprintf("%sarray: %s\n", indent, annotation))
  } else {
    sapply(seq(along = annotation), function(ii) {
      cat(sprintf("%s%s: %s\n",
                  indent,
                  names(annotation)[ii],
                  annotation[ii]))
    })
  }
}

.show.preprocessMethod <- function(preprocessMethod) {
  if (length(preprocessMethod) == 3 && is.null(names(preprocessMethod))) {
    names(preprocessMethod) <- c("rg.norm", "minfi", "manifest")
  }
  if (length(preprocessMethod) == 0) {
    preprocessMethod <- c(
      rg.norm = "unknown",
      minfi = "unknown",
      manifest = "unknown")
  }
  cat("Preprocessing\n")
  cat(sprintf("  Method: %s\n  minfi version: %s\n  Manifest version: %s\n",
              preprocessMethod["rg.norm"],
              preprocessMethod["minfi"],
              preprocessMethod["manifest"]))
}

.getManifestString <- function(annotation) {
  if (length(annotation) == 1) {
    if (annotation == "Unknown") {
      stop("Cannot get Manifest object for an 'Unknown' array")
    }
    return(paste0(annotation, "manifest"))
  }
  if ("array" %in% names(annotation)) {
    if (annotation["array"] == "Unknown") {
      stop("Cannot get Manifest object for an 'Unknown' array")
    }
    return(paste0(annotation["array"], "manifest"))
  }
  stop("unable to get the manifest string for this object")
}

.getAnnotationString <- function(annotation) {
  if (length(annotation) == 1) {
    if (annotation == "Unknown") {
      stop("Cannot get Annotation object for an 'Unknown' array")
    }
    return(sprintf("%sanno", annotation))
  }
  if (all(c("array", "annotation") %in% names(annotation))) {
    if (annotation["array"] == "Unknown") {
      stop("Cannot get Annotation object for an 'Unknown' array")
    }
    return(
      sprintf("%sanno.%s", annotation["array"], annotation["annotation"]))
  }
  stop("unable to get the annotation string for this object")
}

.betaFromMethUnmeth <- function(Meth, Unmeth, object, offset = 0,
                                betaThreshold = 0, minZero = TRUE) {
  stopifnot(offset >= 0)
  stopifnot(betaThreshold >= 0 & betaThreshold <= 0.5)
  if (minZero) {
    Meth <- pmax2(Meth, 0)
    Unmeth <- pmax2(Unmeth, 0)
  }
  beta <- Meth / (Meth + Unmeth + offset)
  if (betaThreshold > 0) {
    beta <- pmin2(pmax2(beta, betaThreshold), 1 - betaThreshold)
  }
  beta
}

.checkAssayNames <- function(object, names) {
  nms <- names(assays(object, withDimnames = FALSE))
  if (!all(names %in% nms)) {
    return(sprintf(
      "object of class '%s' needs to have assay slots with names '%s'",
      class(object),
      paste0(names, collapse = ", ")))
  } else {
    NULL
  }
}

.digestMatrix <- function(mat, digits = 6) {
  content <- sprintf(paste0("%.", digits, "f"), mat)
  # NOTE: Handling signed zero as per IEEE specs
  zero <- paste(c("0.", rep("0", digits)), collapse = "")
  content[content == paste0("-", zero)] <- zero
  digest::digest(c(content, rownames(mat), colnames(mat)))
}

.digestVector <- function(vec, digits = 6) {
  content <- sprintf(paste0("%.", digits, "f"), vec)
  # NOTE: Handling signed zero as per IEEE specs
  zero <- paste(c("0.", rep("0", digits)), collapse = "")
  content[content == paste0("-", zero)] <- zero
  digest::digest(content)
}

.isGenomicOrStop <- function(object) {
  if (!is(object, "GenomicMethylSet") && !is(object, "GenomicRatioSet")) {
    stop("object is of class '", class(object), "', but needs to be of ",
         "class 'GenomicMethylSet' or 'GenomicRatioSet'")
  }
}

.isMethylOrStop <- function(object) {
  if (!is(object, "MethylSet") && !is(object, "GenomicMethylSet")) {
    stop("object is of class '", class(object), "', but needs to be of ",
         "class 'MethylSet' or 'GenomicMethylSet'")
  }
}

.isMethylOrRatio <- function(object) {
  is(object, "MethylSet") ||
    is(object, "GenomicMethylSet") ||
    is(object, "RatioSet") ||
    is(object, "GenomicRatioSet")
}

.isRGOrStop <- function(object) {
  if (!is(object, "RGChannelSet")) {
    stop("object is of class '", class(object), "', but needs to be of ",
         "class 'RGChannelSet' or 'RGChannelSetExtended'")
  }
}

.isMatrixBacked <- function(object) {
  stopifnot(is(object, "SummarizedExperiment"))
  all(vapply(assays(object), is.matrix, logical(1L)))
}

.isDelayedArrayBacked <- function(object) {
  stopifnot(is(object, "SummarizedExperiment"))
  all(vapply(assays(object), is, logical(1L), "DelayedArray"))
}

.isMatrixBackedOrWarning <- function(object, FUN) {
  if (.isDelayedArrayBacked(object)) {
    warning("Memory usage may be high because '", FUN, "()' is not yet ",
            "optimized for use with DelayedArray-backed minfi objects.",
            call. = FALSE,
            immediate. = TRUE)
  } else if (!.isMatrixBacked(object)) {
    warning("Memory usage may be high because '", FUN, "()' is not yet ",
            "optimized for use with non-matrix-backed minfi objects.",
            call. = FALSE,
            immediate. = TRUE)
  }
}

.isMatrixBackedOrStop <- function(object, FUN) {
  if (!.isMatrixBacked(object)) {
    stop("'", FUN, "()' only supports matrix-backed minfi objects.",
         call. = FALSE)
  }
}

.is27k <- function(object) {
  annotation(object)["array"] == "IlluminaHumanMethylation27k"
}

.is450k <- function(object) {
  annotation(object)["array"] == "IlluminaHumanMethylation450k"
}

.isEPIC <- function(object) {
  annotation(object)["array"] == "IlluminaHumanMethylationEPIC"
}

.harmonizeSex <- function(vector) {
  # TODO: This function is not yet implemented
  stop("function not done")
  validMale <- c("M", "MALE")
  validFemale <- c("F", "FEMALE")
  # validUnknown <- c("U", "Unknown")
  if (is.factor(vector)) vector <- as.character(vector)
  if (!is.character(vector)) {
    stop("[.harmonizeSet] argument 'vector' needs to be either a ",
         "character or a factor")
  }
  vector <- toupper(vector)
  vector[vector %in% validMale] <- "M"
  vector[vector %in% validFemale] <- "F"
  if (any(!vector %in% c("M", "F"))) {
    stop("[.harmonizeSet] could not harmonize the vector argument to be ",
         "either 'M' or 'F'")
  }
  vector
}

.harmonizeDataFrames <- function(x, y) {
  stopifnot(is(x, "DataFrame"))
  stopifnot(is(y, "DataFrame"))
  x.only <- setdiff(names(x), names(y))
  y.only <- setdiff(names(y), names(x))
  if (length(x.only) > 0) {
    df.add <- x[1, x.only, drop = FALSE]
    is.na(df.add[1, ]) <- TRUE
    y <- cbind(y, df.add)
  }
  if (length(y.only) > 0) {
    df.add <- y[1, y.only, drop = FALSE]
    is.na(df.add[1, ]) <- TRUE
    x <- cbind(x, df.add)
  }
  list(x = x, y = y[, names(x)])
}

.pDataAdd <- function(object, df) {
  stopifnot(is(df, "data.frame") || is(df, "DataFrame"))
  pD <- colData(object)
  if (any(names(df) %in% names(pD))) {
    alreadyPresent <- intersect(names(df), names(pD))
    warning(sprintf(
      "replacing the following columns in colData(object): %s",
      paste(alreadyPresent, collapse = ", ")))
    pD[, alreadyPresent] <- df[, alreadyPresent]
    df <- df[, !names(df) %in% alreadyPresent]
  }
  if (ncol(df) > 0) {
    # NOTE: Work around for bug in cbind(DataFrame, DataFrame)
    rownam <- rownames(pD)
    pD <- cbind(pD, df)
    rownames(pD) <- rownam
  }
  colData(object) <- pD
  object
}

.pDataFix <- function(df) {
  characterColumns <- c(
    "Slide", "Array", "Sample_Name", "Basename", "SampleID")
  for (col in characterColumns) {
    if (col %in% names(df))
      df[[col]] <- as.character(df[[col]])
  }
  df
}

.NA_type <- function(type) {
  c(vector(type), NA)
}

# Exported functions -----------------------------------------------------------
getMethSignal <- function(object, what = c("Beta", "M"), ...) {
  what <- match.arg(what)
  switch(what,
         "Beta" = getBeta(object, ...),
         "M" = getM(object, ...)
  )
}



############## Load data ########################

### EuGEI ##
load("/mnt/data1/EuGEI/QC/GeorginasQC/All_Plates_Blood_WithRepeats/JustEuGEIresults/EuGEIBloodSamples_Normalised.rdat")
keep = read.table('/mnt/data1/GWAS_bloodcelltypes/EuGEI/keep.txt', header = T, stringsAsFactors = F)
#934 samples
pheno_fam <- pheno[which(pheno$Geno.CHIP.Location %in% keep$FID),]
pheno_eu <- pheno_fam
betas_eu <- betas
rm(pheno,betas, pheno_fam,keep)

## Extend ##
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EUR_unlrelated_QCd/EXTEND_batch1_2_genoQCd_Normalised.rdat")
#1022 samples

pathM <- paste('/mnt/data1/EXTEND/Genotypes/Imputed/EXTEND_Unrelated_EUR_QCd', c(".bed", ".bim", ".fam"), sep = "")
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
fam <- SNP_M$fam

pheno_ex <- pheno[which(pheno$IID %in% fam$pedigree),]
betas_ex <- betas[,which(colnames(betas) %in% pheno_eu$Basename)]
rm(pheno, betas,fam, SNP_M, pathM)

## UnderSoc ##
load("/mnt/data1/EPICQC/UnderstandingSociety/US_Betas_Pheno.rda")
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
betas <- dat[,which(colnames(dat) %in% rownames(pheno))]

pheno_us <- pheno
betas_us <- betas
rm(fam2, fam, SNP_M, dat, betas, pheno, pathM)

###############################################################



############# Run RGset ################
library(methylumi)
library(minfi)
library(wateRmelon)
require(gdata)
library(minfi)
library(ggplot2)
require(gridExtra)
library(plyr)
require(IlluminaHumanMethylationEPICmanifest)
library(dplyr)
library(tidyr)

## EuGEI - no RGset ... create one
# idatPath <- c('/mnt/data1/EuGEI/idats/')
# Basename <- as.data.frame(pheno_eu$Basename)
# colnames(Basename) <- 'Basename'
# Basename$Basename <- as.character(Basename$Basename)
# RGset_eugei <- read.metharray.exp(base = idatPath, targets = Basename, force = TRUE)
# save(RGset_eugei, file = "/mnt/data1/GWAS_bloodcelltypes/EuGEI/RGset_eugei.rdat")

# load("/mnt/data1/GWAS_bloodcelltypes/EuGEI/RGset_eugei.rdat")
# eugei_counts = estimateCellCounts(RGset_eugei, compositeCellType = "Blood", 
#                    processMethod = "auto", probeSelect = "auto", 
#                    cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "Neu", "Eos"))
# 
# save(eugei_counts, file =  "/mnt/data1/GWAS_bloodcelltypes/EuGEI/RGset_cell_count.rdat")


## Extend
load("/mnt/data1/GWAS_bloodcelltypes/EXTEND/EXTEND_RGset.rdat")
RGset_extend <- RGset
rm(RGset)
extend_counts = estimateCellCounts(RGset_extend, compositeCellType = "Blood", 
                                  processMethod = "auto", probeSelect = "auto", 
                                  cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "Neu", "Eos"))
#filter for passed samples as rgset contains all studies


## UnderSoc
#Make RGset but idat files are missing?? Create list of sampeles to move from MDrive to knight
#Make Basename df again for US
Basename <- as.data.frame(pheno_us$barcode)
colnames(Basename) <- 'Basename'
Basename$Basename <- as.character(Basename$Basename)
# write.csv(Basename, "US_Samples.csv")

idatPath = c('/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/idats')
RGset_undersoc <- read.metharray.exp(base = idatPath, targets = Basename, force = TRUE)
save(RGset_undersoc, file = "/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/RGset_undersoc.rdat")

load('/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/RGset_undersoc.rdat')
undersoc_counts = estimateCellCounts(RGset_undersoc, compositeCellType = "Blood",
                   processMethod = "auto", probeSelect = "auto",
                   cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "Neu", "Eos"))

save(undersoc_counts, file =  "/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/RGset_cell_count.rdat")

########### Statistical Analysis #############

###########
## EuGEI
###########
load("/mnt/data1/GWAS_bloodcelltypes/EuGEI/RGset_cell_count.rdat")
RGset_eugei <- eugei_counts
load("/mnt/data1/GWAS_bloodcelltypes/EuGEI/cellcounts.rdat")
Mset_eugei <- counts
rm(eugei_counts, counts)
#Remove extra samples from Mset_eugei
Mset_eugei <- Mset_eugei[which(rownames(Mset_eugei) %in% rownames(RGset_eugei)),]

#Match order of samples
RGset_eugei <- RGset_eugei[order(rownames(RGset_eugei)),]
Mset_eugei <- Mset_eugei[order(rownames(Mset_eugei)),]
identical(rownames(Mset_eugei), rownames(RGset_eugei))

#order 
RGset_eugei <- RGset_eugei[,colnames(Mset_eugei)]
results_eugei <- c()

for (i in 1:ncol(Mset_eugei)){
  print(i)
  y = Mset_eugei[,i]
  x = RGset_eugei[,i]
  temp <- t.test(y , x)
  corrtest <- cor.test(x, y, method="pearson")
  results_eugei <- rbind(results_eugei, cbind(temp$p.value, corrtest$p.value, corrtest$estimate))

}
rownames(results_eugei) <- colnames(Mset_eugei)
colnames(results_eugei) <- c('T.Test P.value', "Pearson Correlation P.value", 
                              "Pearson Correlation R")


estcell_eugei <- as.data.frame(cbind(Mset_eugei[,1], RGset_eugei[,1]))
colnames(estcell_eugei) <- c("Mset", "RGset")

library(gridExtra)
library(ggpubr)

pdf('EuGEI Correlations.pdf', onefile = TRUE)
for(i in 1:ncol(Mset_eugei)){
  estcell_eugei <- as.data.frame(cbind(Mset_eugei[,i], RGset_eugei[,i]))
  colnames(estcell_eugei) <- c("Mset", "RGset")
  print(ggscatter(data = estcell_eugei, x = 'Mset', y = "RGset",
                  add = 'reg.line',
                  add.params = list(color = 'blue', fill = 'lightgray'),
                  conf.int = T) +
          stat_cor(method = 'pearson') +
          labs(title = paste("eugei", colnames(Mset_eugei)[i], sep = "_"),
               x = 'Mset', y = 'RGset') +
          theme(plot.title = element_text(hjust = 0.5)))
}
dev.off()
###########
## Extend
###########
load("/mnt/data1/GWAS_bloodcelltypes/EXTEND/cell_counts_using_rgset.rdat")
RGset_extend <- cell_counts
load("/mnt/data1/GWAS_bloodcelltypes/EXTEND/cell_counts.rdat")
Mset_extend <- counts
rm(cell_counts, counts)

#Remove extra samples from Mset_extend
Mset_extend <- Mset_extend[which(rownames(Mset_extend) %in% rownames(RGset_extend)),]

#Match order of samples
RGset_extend <- RGset_extend[order(rownames(RGset_extend)),]
Mset_extend <- Mset_extend[order(rownames(Mset_extend)),]
identical(rownames(Mset_extend), rownames(RGset_extend))

#order 
RGset_extend <- RGset_extend[,colnames(Mset_extend)]

results_extend <- c()
for (i in 1:ncol(Mset_extend)){
  y = Mset_extend[,i]
  x = RGset_extend[,i]
  temp <- t.test(y,x)
  corrtest <- cor.test(x, y, method="pearson")
  results_extend <- rbind(results_extend, cbind(temp$p.value, corrtest$p.value, corrtest$estimate))
}

rownames(results_extend) <- colnames(Mset_extend)
colnames(results_extend) <- c('T.Test P.value', "Pearson Correlation P.value", 
                                "Pearson Correlation R")

pdf('extend Correlations.pdf', onefile = TRUE)
for(i in 1:ncol(Mset_extend)){
  estcell_extend <- as.data.frame(cbind(Mset_extend[,i], RGset_extend[,i]))
  colnames(estcell_extend) <- c("Mset", "RGset")
  print(ggscatter(data = estcell_extend, x = 'Mset', y = "RGset",
                  add = 'reg.line',
                  add.params = list(color = 'blue', fill = 'lightgray'),
                  conf.int = T) +
          stat_cor(method = 'pearson') +
          labs(title = paste("extend", colnames(Mset_extend)[i], sep = "_"),
               x = 'Mset', y = 'RGset') +
          theme(plot.title = element_text(hjust = 0.5)))
}
dev.off()
###########
#### Understanding Society
############
load('/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/RGset_cell_count.rdat')
RGset_undersoc <- undersoc_counts
load("/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/cellcounts.rdat")
Mset_undersoc <- counts
rm(counts, undersoc_counts)

#Remove extra samples from Mset_undersoc
Mset_undersoc <- Mset_undersoc[which(rownames(Mset_undersoc) %in% rownames(RGset_undersoc)),]

#Match order of samples
RGset_undersoc <- RGset_undersoc[order(rownames(RGset_undersoc)),]
Mset_undersoc <- Mset_undersoc[order(rownames(Mset_undersoc)),]
identical(rownames(Mset_undersoc), rownames(RGset_undersoc))

#order 
RGset_undersoc <- RGset_undersoc[,colnames(Mset_undersoc)]

#t test loop
results_undersoc <- c()
for (i in 1:ncol(Mset_undersoc)){
  print(i)
  y = Mset_undersoc[,i]
  x = RGset_undersoc[,i]
  temp <- t.test(y , x)
  corrtest <- cor.test(x, y, method="pearson")
  results_undersoc <- rbind(results_undersoc, cbind(temp$p.value, corrtest$p.value, corrtest$estimate))
}

rownames(results_undersoc) <- colnames(Mset_undersoc)
colnames(results_undersoc) <- c('T.Test P.value', "Pearson Correlation P.value", 
                                "Pearson Correlation R")

grid.table(results_undersoc)

pdf('UnderSoc Correlations.pdf', onefile = TRUE)
for(i in 1:ncol(Mset_undersoc)){
  estcell_undersoc <- as.data.frame(cbind(Mset_undersoc[,i], RGset_undersoc[,i]))
  colnames(estcell_undersoc) <- c("Mset", "RGset")
  print(ggscatter(data = estcell_undersoc, x = 'Mset', y = "RGset",
            add = 'reg.line',
            add.params = list(color = 'blue', fill = 'lightgray'),
            conf.int = T) +
    stat_cor(method = 'pearson') +
    labs(title = paste("undersoc", colnames(Mset_undersoc)[i], sep = "_"),
         x = 'Mset', y = 'RGset') +
    theme(plot.title = element_text(hjust = 0.5)))
}
dev.off()

########## Tables & Histrogram summaries ##############
pdf('Correlation Undersoc.pdf', onefile = T)
grid.table(results_undersoc)
dev.off()
pdf('Correlation Extend.pdf', onefile = T)
grid.table(results_extend)
dev.off()
pdf('Correlation EuGEI.pdf', onefile = T)
grid.table(results_eugei)
dev.off()

corretbl <- cbind(results_eugei[,3], results_extend[,3], results_undersoc[,3])
corretbl <- as.data.frame(corretbl)
colnames(corretbl) <- c('Eugei', 'Extend', 'UnderSoc')

corretbl$celltype <- rownames(corretbl)
corretbl2 <- melt(corretbl, id.vars = 'celltype')

pdf('correlation summary.pdf')
ggplot(corretbl2, aes(x = factor(corretbl2$celltype), y = corretbl2$value,
                      fill = factor(corretbl2$variable))) +
  geom_bar(stat = 'identity',position=position_dodge()) +
  xlab("Cell Type") + ylab("Pearson Correlation (R)") +
  scale_fill_hue(name="Study") 
dev.off()

### Does Gender distribution differ between studies? ####
Gen_us <- as.data.frame(pheno_us$nsex)
Gen_us$Study <- rep('UnderSoc', nrow(Gen_us))
colnames(Gen_us) <- c('Sex', 'Study')
Gen_us$Sex <- as.numeric(as.character(Gen_us$Sex))
Gen_us$Sex[Gen_us$Sex == 2] <- 'Female'
Gen_us$Sex[Gen_us$Sex == 1] <- 'Male'
Gen_us$Sex[Gen_us$Sex == 0] <- NA


Gen_eu <- as.data.frame(pheno_eu$Sex)
Gen_eu$Study <- rep('EuGEI', nrow(Gen_eu))
colnames(Gen_eu) <- c('Sex', 'Study')
Gen_eu$Sex <- as.character(Gen_eu$Sex)

Gen_ex <- as.data.frame(pheno_ex$Sex)
Gen_ex$Study <- rep('Extend', nrow(Gen_ex))
colnames(Gen_ex) <- c('Sex', 'Study')
Gen_ex$Sex <- as.character(Gen_ex$Sex)

Gen <- rbind(Gen_eu, Gen_ex, Gen_us)
Gen <- Gen[which(!is.na(Gen$Sex)),]
pdf('Sex distribution.pdf')
ggplot(Gen, aes(fill= Gen$Sex, x =  Gen$Study)) +
  geom_bar(position=position_dodge()) +
  xlab("Study") +
  scale_fill_hue(name="Sex")
dev.off()

### Does Age distribution differ between studies?
Age_us <- as.data.frame(pheno_us$confage)
Age_us$Study <- rep('UnderSoc', nrow(Age_us))
colnames(Age_us) <- c('Age', 'Study')
Age_us$Age <- as.numeric(as.character(Age_us$Age))



Age_eu <- as.data.frame(pheno_eu$Age)
Age_eu$Study <- rep('EuGEI', nrow(Age_eu))
colnames(Age_eu) <- c('Age', 'Study')
Age_eu$Age <- as.numeric(as.character(Age_eu$Age))

Age_ex <- as.data.frame(pheno_ex$Age)
Age_ex$Study <- rep('Extend', nrow(Age_ex))
colnames(Age_ex) <- c('Age', 'Study')
Age_ex$Age <- as.numeric(as.character(Age_ex$Age))

allAges <- rbind(Age_eu, Age_ex, Age_us)
Age <- as.data.frame(Age)

pdf('Age Distribution.pdf')
ggplot(allAges, aes(allAges$Age, fill = allAges$Study)) +
  geom_density(alpha = 0.5) +
  xlab("Age") +
  scale_fill_hue(name="Study")
dev.off()

### Boxplot of each cell type 
pdf('Cell Type Distribution per study.pdf', onefile = T)
par(mfrow = c(1,2))
boxplot(Mset_eugei, las = 2, main = 'Mset_eugei')
boxplot(RGset_eugei, las = 2, main = 'RGset_eugei')

boxplot(Mset_extend, las = 2, main = 'Mset_extend')
boxplot(RGset_extend, las = 2, main = 'RGset_extend')

boxplot(Mset_undersoc, las = 2, main = 'Mset_undersoc')
boxplot(RGset_undersoc, las = 2, main = 'RGset_undersoc')
dev.off()

### Do the probes in the RGSet change each time? And why?
