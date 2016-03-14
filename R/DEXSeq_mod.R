## MOdified 2 functions from the DEXSeq package; so that it's possible to read
## our Cigene RefSeq Gene IDs with DEXseq.

## The initial problem is that DEXSeq tries to seperate exon number and geneID
## by ':'. This causes problems since we have already a ':' in our cigene
## RefSeq IDs. So the solution is to gsub the ':' from the geneID into a '-'
## This is done within the 'DEXSeqDataSet.MOD' function. 

DEXSeqDataSet.MOD <- 
    function (countData, sampleData, design = ~sample + exon + condition:exon, 
    featureID, groupID, featureRanges = NULL, transcripts = NULL) 
{
    stopifnot(class(countData) %in% c("matrix", "data.frame"))
    countData <- as.matrix(countData)
    stopifnot(class(featureID) %in% c("character", "factor"))
    stopifnot(class(groupID) %in% c("character", "factor"))
    stopifnot(class(sampleData) %in% c("data.frame"))
    stopifnot(length(groupID) == nrow(countData))
    stopifnot(length(featureID) == length(groupID))
    stopifnot(nrow(sampleData) == ncol(countData))
    modelFrame <- cbind(sample = rownames(sampleData), sampleData)
    modelFrame <- rbind(cbind(modelFrame, exon = "this"), cbind(modelFrame, 
        exon = "others"))
    rownames(modelFrame) <- NULL
    colData <- DataFrame(modelFrame)
    if (!"exon" %in% all.vars(design)) {
        stop("The formula does not specify a contrast with the variable 'exon'")
    }
    allVars <- all.vars(design)
    if (any(!allVars %in% colnames(colData))) {
        notPresent <- allVars[!allVars %in% colnames(colData)]
        notPresent <- paste(notPresent, collapse = ",")
        stop(sprintf("the variables '%s' of the parameter 'design' are not specified in the columns of the sampleData", 
            notPresent))
    }
    ## Modified
    if (any(grepl(" |:", groupID) | grepl(" |:", featureID))) {
        warning("empty spaces or ':' characters were found either in your groupIDs or in your featureIDs, these will be removed from the identifiers")
        groupID <- gsub(" ", "", groupID)
        featureID <- gsub(" |:", "-", featureID)
    }
    rownames(countData) <- paste(groupID, featureID, sep = ":")
    forCycle <- split(1:nrow(countData), as.character(groupID))
    others <- lapply(forCycle, function(i) {
        sct <- countData[i, , drop = FALSE]
        rs <- t(sapply(1:nrow(sct), function(r) colSums(sct[-r, 
            , drop = FALSE])))
        rownames(rs) <- rownames(sct)
        rs
    })
    others <- do.call(rbind, others)
    stopifnot(all(rownames(countData) %in% rownames(others)))
    others <- others[rownames(countData), ]
    nCountData <- cbind(countData, others)
    if (!is.null(featureRanges)) {
        stopifnot(class(featureRanges) %in% c("GRanges", "GRangesList"))
        se <- SummarizedExperiment(nCountData, colData = colData, 
            rowRanges = featureRanges)
    }
    else {
        se <- SummarizedExperiment(nCountData, colData = colData)
    }
    names(assays(se))[1] = "counts"
    mcols(se)$featureID <- featureID
    mcols(se)$groupID <- groupID
    mcols(se)$exonBaseMean <- rowMeans(countData)
    mcols(se)$exonBaseVar <- rowVars(countData)
    if (!is.null(transcripts)) {
        mcols(se)$transcripts <- transcripts
    }
    rownames(se) <- paste(groupID, featureID, sep = ":")
    dds <- DESeqDataSet(se, design, ignoreRank = TRUE)
    maxGene <- names(which.max(table(groupID)))
    rows <- mcols(dds)$groupID %in% maxGene
    numExons <- sum(rows)
    exonCol <- rep(factor(featureID[rows]), nrow(sampleData))
    modelFrame <- data.frame(sample = rep(rownames(sampleData), 
        each = numExons), exon = exonCol)
    varNames <- colnames(sampleData)
    for (i in varNames) {
        modelFrame[[i]] <- rep(sampleData[[i]], each = numExons)
    }
    modelFrame$dispersion <- NA
    modelFrame$sizeFactor <- NA
    modelFrame$count <- NA
    dxd <- new("DEXSeqDataSet", dds, modelFrameBM = modelFrame)
    return(dxd)
}


DEXSeqDataSetFromHTSeq.MOD <-
    function (countfiles, sampleData,
                                    design = ~sample + exon + condition:exon, 
                                    flattenedfile = NULL) 
{
    if (!all(sapply(countfiles, class) == "character")) {
        stop("The countfiles parameter must be a character vector")
    }
    lf <- lapply(countfiles, function(x) read.table(x, header = FALSE, 
        stringsAsFactors = FALSE))
    if (!all(sapply(lf[-1], function(x) all(x$V1 == lf[1]$V1)))) 
        stop("Count files have differing gene ID column.")
    dcounts <- sapply(lf, `[[`, "V2")
    rownames(dcounts) <- lf[[1]][, 1]
    dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", 
        ]
    ##rownames(dcounts) <- sub(":", ":E", rownames(dcounts))
    ## Modified
    fix.ids <- function(x){
        len = nchar(x)
        ex = substr(x, len-2, len)
        gene = substr(x, 1, len-4)
        return(paste(gene, ex, sep =":E"))
    }
    rownames(dcounts) <- sapply(rownames(dcounts), fix.ids, USE.NAMES = FALSE )
    colnames(dcounts) <- countfiles
    splitted <- strsplit(rownames(dcounts), ":E")
    exons <- paste('E',sapply(splitted, "[[", 2), sep ='')
    genesrle <- sapply(splitted, "[[", 1)
    if (!is.null(flattenedfile)) {
        aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, 
            header = FALSE)
        colnames(aggregates) <- c("chr", "source", "class", "start", 
            "end", "ex", "strand", "ex2", "attr")
        aggregates$strand <- gsub("\\.", "*", aggregates$strand)
        aggregates <- aggregates[which(aggregates$class == "exonic_part"), 
            ]
        aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
        aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", 
            aggregates$attr)
        transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", 
            aggregates$attr)
        transcripts <- strsplit(transcripts, "\\+")
        exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", 
            aggregates$attr)
        exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, 
            end = aggregates$end), strand = aggregates$strand)
        names(exoninfo) <- paste(aggregates$gene_id, exonids, 
            sep = ":E")
        names(transcripts) <- rownames(exoninfo)
        if (!all(rownames(dcounts) %in% names(exoninfo))) {
            stop("Count files do not correspond to the flattened annotation file")
        }
        matching <- match(rownames(dcounts), names(exoninfo))
        stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
        stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
        dxd <- DEXSeqDataSet.MOD(dcounts, sampleData, design, exons, 
            genesrle, exoninfo[matching], transcripts[matching])
        return(dxd)
    }
    else {
        dxd <- DEXSeqDataSet.MOD(dcounts, sampleData, design, exons, 
            genesrle)
        return(dxd)
    }
}
