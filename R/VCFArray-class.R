### -------------------------
### classes
### -------------------------

setClassUnion("VcfFile_OR_VcfStack", c("VcfFile", "VcfStack"))

setClass("VCFArraySeed",
         contains = "Array",
         slots = c(vcffile = c("VcfFile_OR_VcfStack"),
                   name = "character",
                   dim = "integer",
                   dimnames = "list",
                   gr = "GRanges"))

### -------------------------
### VCFArraySeed methods
### -------------------------

setMethod("dim", "VCFArraySeed", function(x) x@dim)
setGeneric("vcffile", function(x) standardGeneric("vcffile"))
setMethod("vcffile", "VCFArraySeed", function(x) x@vcffile)
setMethod("rowRanges", "VCFArraySeed", function(x) x@gr)

setMethod("show", "VCFArraySeed", function(object)
{
    vcf <- vcffile(object)
    if (is(vcf, "VcfFile")) {
        cat("VCFArraySeed\n",
            "VCF file path: ", path(vcf), "\n",
            "VCF index path: ", index(vcf), "\n",
            "array data: ", object@name, "\n",
            "dim: ", paste(dim(object), collapse=" x "), "\n",
            sep="")
    } else if (is(vcf, "VcfStack")) {
        vcffiles <- files(vcf)
        cat("VCFArraySeed\n",
            "VcfStack object with ", nrow(vcf), " files and ", ncol(vcf), " samples", "\n", 
            "VCF file path: \n", paste(unname(sapply(vcffiles, path))), "\n",
            ## "VCF index path: ", index(vcf), "\n",
            "array data: ", object@name, "\n",
            "dim: ", paste(dim(object), collapse=" x "), "\n",
            sep="")
    }
})

#' @import GenomicRanges
.extract_array_from_VCFArray <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    tp <- ifelse(x@name == "GT", "character", "integer")  ## FIXME
    if (any(ans_dim == 0L)){
        ans <- eval(parse(text = tp))(0)  ## return integer(0) / character(0)
        dim(ans) <- ans_dim
    } else {
        vcf <- vcffile(x)
        ridx <- index[[1]]
        if(is.null(ridx))
            ridx <- seq_len(nrow(x))
        cidx <- index[[2]]
        if(is.null(cidx))
            cidx <- seq_len(ncol(x))
        gr <- x@gr
        param <- ScanVcfParam(which = gr[gr$pos %in% ridx],
                              samples = colnames(x)[cidx])
        res <- readGeno(vcf, x@name, param = param)
        ans <- res
    }
    ans
}

setMethod("extract_array", "VCFArraySeed", .extract_array_from_VCFArray)

### ---------------------------
### VCFArraySeed constructor
### ---------------------------

#' @import VariantAnnotation
## #' @importFrom Rsamtools countTabix

VCFArraySeed <- function(file = character(), index = character(), name = character())
{
    ## if (!isSingleString(file))
    ##     stop(wmsg(
    ##         "'file' must be a single string specifying the path to ",
    ##         "the vcf file where the assay data is located."))
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name of ",
                  "the assay data corresponding to the vcf 'FORMAT' field."))

    if (is(file, "VcfFile")) {
        vcf <- file
        if (!is.na(index(vcf)) && length(index)) {
            stop("'index' cannot be used when 'VcfFile' ",
                 "input already has the index file.")
        } else if (is.na(index(vcf))) {
            if (length(index)) {
                index(vcf) <- index
            } else {
                ## index(vcf) <- paste(path(vcf), "tbi", sep = ".")
                vcf <- indexVcf(vcf)
            }
        }
    } else if (is(file, "RangedVcfStack")) {
        vcf <- file
    } else if(isSingleString(file)) {
        if(file.exists(file)) file <- normalizePath(file)  ## in base R
        if (!length(index)) index = paste(file, "tbi", sep = ".")
        vcf <- VcfFile(file, index = index)
    }
    ## read the header info
    if (is(vcf, "VcfStack")) {
        header <- scanVcfHeader(files(vcf)[[1]])
    } else {
        header <- scanVcfHeader(vcf)
    }
    stopifnot(name %in% rownames(geno(header)))
    nsamps <- length(samples(header))
    
    ## lightweight filter. Only return REF, rowRanges
    if (is(vcf, "RangedVcfStack")) {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA, which = rowRanges(vcf))
        readvcf <- readVcfStack(vcf, param = param)
    } else {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA, samples = )
        readvcf <- readVcf(vcf, genome = "hg19", param = param)
    }
    gr <- granges(rowRanges(readvcf)) 
    gr$pos <- seq_along(gr)
    
    nvars <- length(gr)
    
    new("VCFArraySeed", vcffile = vcf, name = name,
        dim = c(nvars, nsamps),
        dimnames = list(names(gr), samples(header)),
        gr = gr)
}

### -------------------
### VCFArray class
### -------------------
## ' @importClassesFrom DelayedArray DelayedArray DelayedMatrix
setClass("VCFArray", contains = "DelayedArray")
setClass("VCFMatrix", contains=c("DelayedMatrix", "VCFArray"))
setMethod("matrixClass", "VCFArray", function(x) "VCFMatrix")
setAs("VCFArray", "VCFMatrix", function(from) new("VCFMatrix", from))
setAs("VCFMatrix", "VCFArray", function(from) from)
setAs(
    "ANY", "VCFMatrix",
    function(from) as(as(from, "VCFArray"), "VCFMatrix"))


### -----------------
### Validity check
### -----------------

.validate_VCFArray <- function(x)
{
    if (!is(x@seed, "VCFArraySeed"))
        return(wmsg("'x@seed' must be a VCFArraySeed object"))
    TRUE
}

setValidity2("VCFArray", .validate_VCFArray)

### --------------
### VCFArray constructor
### --------------

setMethod(
    "DelayedArray", "VCFArraySeed",
    function(seed) new_DelayedArray(seed, Class="VCFArray")  ## need "extract_array" to work.
    )

VCFArray <- function(file = character(), index = character(), name=NA)
{
    if (is(file, "VCFArraySeed")) {
        if (!missing(name))
            stop(wmsg(
                "VCFArray() must be called with a single argument ",
                "when passed an VCFArraySeed object"))
        seed <- file
    } ## else if (is(file, "VcfStack")) {
    ##     NULL
    ## }
    else {
        if (is.na(name)) {
            header <- scanVcfHeader(file) 
            geno <- rownames(geno(header))
            if (length(geno) == 1) {
                name <- geno
            } else {
                message('The Available values for "name" argument are: ',
                               paste(geno, collapse=" "), "\n",
                               "Please specify, otherwise, ",
                               'The default value of "GT" will be returned.', "\n")
                name <- "GT"
            }
        }
        seed <- VCFArraySeed(file, index = index, name = name)
    }
    DelayedArray(seed)   ## does the automatic coercion to VCFMatrix if 2-dim.
}

## setMethod("seed", "VCFArray", function(x) x@seed)
### -------------
### example 
### -------------

## setMethod("example", "VCFArray", function (topic = "ANY", package = "VCFArray")
## {
##     fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
##     seed <- VCFArraySeed(fl, "GT")
##     ## seed <- VCFArraySeed(fl, "DS")
##     ## seed <- VCFArraySeed(fl, "GL")
##     VCFArray(seed)
## })
