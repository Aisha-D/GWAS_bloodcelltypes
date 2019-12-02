##This script looks into the issue with estiamtecellcounts from minfi. 
##Neutrophils do not have enough probes to carry the estimation on. We will compare using rgset and mset
##Then look at the other issue of what p-value is suitable for each cell type

setwd('/mnt/data1/aisha/GWAS_bloodcelltypes/')

#######################################
###        PART1:Using RGSET        ###
#######################################
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

load("/mnt/data1/EXTEND/GWAS_bloodcelltypes/EXTEND_RGset.rdat")
# https://github.com/hansenlab/minfi/issues/167 describes well what the issue is
############ Subset rgset ###########
# keep <- RGset@colData$Basename[1:10]
# eset2 = RGset[, sampleNames(RGset) %in% keep]

#break the estimate cell comp function
# https://github.com/hansenlab/minfi/issues/167 describes well what the issue is
#functions that are preset outside estimate cell functions
################## DO NOT EDIT THESE FUNCTIONS  ###################

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

######  functions estiamte cell relaies on ######### Can edit these
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
      y <- x[x[, "p.value"] < 1e-5, ] #########! ADJUST HERE TO DECREASE NA ##########
      yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
      c(rownames(yAny)[seq(numProbes * 2)])
    })
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < 1e-5, ] #########! ADJUST HERE TO DECREASE NA ##########
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

cell_counts <- estimateCellCounts(RGset,compositeCellType = "Blood",
                                  processMethod = "auto", probeSelect = 'auto',
                                  cellTypes = c("CD8T", "CD4T", "NK", "Bcell",
                                                "Mono", "Gran","Eos", "Neu"))

save(cell_counts, file = 'cell_counts_using_rgset.rdat')


############################################
###              PART2:MSET              ###
############################################
rm(list=setdiff(ls(), 'RGset'))
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EXTEND_batch1_2_genoQCd_Normalised.rdat")

## need to use custom functions
library(genefilter)
library(quadprog)
validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]
  
  if(is.null(L.forFstat)) {
    L.forFstat <- diag(sizeModel)[-1,] # All non-intercept coefficients
    colnames(L.forFstat) <- colnames(xTest) 
    rownames(L.forFstat) <- colnames(xTest)[-1] 
  }
  
  ## Initialize various containers
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()
  
  if(verbose)
    cat("[validationCellType] ")
  for(j in 1:M) { # For each CpG
    ## Remove missing methylation values
    ii <- !is.na(Y[j,])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j,]
    
    if(j%%round(M/10)==0 && verbose)
      cat(".") # Report progress
    
    try({ # Try to fit a mixed model to adjust for plate
      if(!is.null(modelBatch)) {
        fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
        OLS <- inherits(fit,"try-error") # If LME can't be fit, just use OLS
      } else
        OLS <- TRUE
      
      if(OLS) {
        fit <- lm(modelFix, data=pheno[ii,])
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
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
    })
  }
  if(verbose)
    cat(" done\n")
  ## Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1
  
  ## Get P values corresponding to F statistics
  Pval <- 1-pf(Fstat, sizeModel, degFree)
  
  out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
              sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval,
              orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, nObserved=nObserved,
              degFree=degFree)
  
  out
}

pickCompProbes <- function(rawbetas, cellInd, cellTypes = NULL, numProbes = 50, probeSelect = probeSelect) {
  ## p is matrix of beta values
  ## cellInd is vector denoting cell type 
  splitit <- function(x) {
    split(seq(along=x), x)
  }
  
  #   p <- getBeta(mSet)
  #   pd <- as.data.frame(colData(mSet))
  if(!is.null(cellTypes)) {
    if(!all(cellTypes %in% as.character(cellInd)))
      stop("elements of argument 'cellTypes' is not part of 'cellInd'")
    keep <- which(as.character(cellInd) %in% cellTypes)
    rawbetas <- rawbetas[,keep]
    cellInd<-cellInd[keep]
  }
  ## make cell type a factor 
  cellInd <- factor(cellInd)
  ffComp <- rowFtests(rawbetas, cellInd) #each row has ftests for each cell type factor
  prof <- sapply(splitit(cellInd), function(i) rowMeans(rawbetas[]))
  r <- matrixStats::rowRanges(rawbetas)
  compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range") 
  tIndexes <- splitit(cellInd)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(rawbetas))
    x[i] <- 1
    return(rowttests(rawbetas, factor(x)))
  })
  
  if (probeSelect == "any"){
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[,"p.value"] < 1e-5,]
      yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]      
      c(rownames(yAny)[1:(numProbes*2)])
    })
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[,"p.value"] < 1e-5,]
      yUp <- y[order(y[,"dm"], decreasing=TRUE),]
      yDown <- y[order(y[,"dm"], decreasing=FALSE),]
      c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
    })
  }
  
  trainingProbes <- unique(unlist(probeList))
  rawbetas <- rawbetas[trainingProbes,]
  
  pMeans <- colMeans(rawbetas)
  names(pMeans) <- cellInd
  
  form <- as.formula(sprintf("y ~ %s - 1", paste(levels(cellInd), collapse="+")))
  phenoDF <- as.data.frame(model.matrix(~cellInd-1))
  colnames(phenoDF) <- sub("cellInd", "", colnames(phenoDF))
  if(ncol(phenoDF) == 2) { # two group solution
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(rawbetas))
  } else { # > 2 group solution
    tmp <- validationCellType(Y = rawbetas, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts
  }
  
  out <- list(coefEsts = coefEsts, compTable = compTable,
              sampleMeans = pMeans)
  return(out)
}

projectCellType <- function(Y, coefCellType, contrastCellType=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
  if(is.null(contrastCellType))
    Xmat <- coefCellType
  else
    Xmat <- tcrossprod(coefCellType, contrastCellType) 
  
  nCol <- dim(Xmat)[2]
  if(nCol == 2) {
    Dmat <- crossprod(Xmat)
    mixCoef <- t(apply(Y, 2, function(x) { solve(Dmat, crossprod(Xmat, x)) }))
    colnames(mixCoef) <- colnames(Xmat)
    return(mixCoef)
  } else {
    nSubj <- dim(Y)[2]
    
    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)
    
    if(nonnegative){
      if(lessThanOne) {
        Amat <- cbind(rep(-1, nCol), diag(nCol))
        b0vec <- c(-1, rep(0, nCol))
      } else {
        Amat <- diag(nCol)
        b0vec <- rep(0, nCol)
      }
      for(i in 1:nSubj) {
        obs <- which(!is.na(Y[,i])) 
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), Amat, b0vec)$sol
      }
    } else {
      for(i in 1:nSubj) {
        obs <- which(!is.na(Y[,i])) 
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
      }
    }
    return(mixCoef)
  }
}

## need to select from mset
sampleSheet <- pheno
rawbetas<-betas
sampleSheet<-sampleSheet[match(colnames(rawbetas), sampleSheet$Basename),]

compositeCellType = "Blood"
processMethod = "auto"
probeSelect = "auto"
cellTypes = c("CD8T","CD4T", "NK", "Bcell","Mono","Gran","Eos","Neu")
cellInd <- cellTypes
referencePlatform = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k")
returnAll = FALSE
meanPlot = FALSE
verbose = TRUE

if(verbose) cat("[estimateCellCounts] Picking probes for composition estimation.\n")
compData <- pickCompProbes(rawbetas=rawbetas, cellInd=cellInd, cellTypes = cellTypes, numProbes = 50, probeSelect = probeSelect)
coefs <- compData$coefEsts

if(verbose) cat("[estimateCellCounts] Estimating composition.\n")
counts <- projectCellType(rawbetas[rownames(coefs), ], coefs)
rownames(counts) <- colnames(rawbetas)

save(counts, file = 'cell_counts_using_mset.rdat')

############################################
###              PART3:P.VAL             ###
############################################
#we are interested in looking at the pvalue of each cell type where noise decreases

load("/mnt/data1/aisha/GWAS_bloodcelltypes/cell_counts_using_mset.rdat")
load("/mnt/data1/aisha/GWAS_bloodcelltypes/cell_counts_using_rgset.rdat")

cell_counts <- as.data.frame(cell_counts)
counts <- as.data.frame(counts)

##first we compare the boxplot of mset and rgset - are the variation the same?
pdf("EstCellCounts rgset vs mset.pdf", width = 10, height = 8)
par(mfrow = c(1,2))
boxplot(cell_counts, las = 2,  main="RGset")
boxplot(counts, las = 2, main="Mset")
dev.off()

#Noise model explained well by - https://stats.stackexchange.com/questions/351707/noise-in-regression-data
#now we model noise - as mset is easier we will work with using mset than rgset


