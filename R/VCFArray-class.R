### -------------------------
### classes
### -------------------------

#' @import GenomicFiles
setClassUnion("VcfFile_OR_RangedVcfStack", c("VcfFile", "RangedVcfStack"))

setClass("VCFArraySeed",
         contains = "Array",
         slots = c(vcffile = c("VcfFile_OR_RangedVcfStack"),
                   vcfheader = "VCFHeader",
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
            "VcfStack object with ", nrow(vcf), " files and ",
            ncol(vcf), " samples", "\n", 
            "VCF file path: ",
            "'", path(vcffiles[[1]]), "' and ", nrow(vcf)-1, " more...", "\n",
            "VCF index path: ",
            "'", index(vcffiles[[1]]), "' and ", nrow(vcf)-1, " more...", "\n",
            "array data: ", object@name, "\n",
            "dim: ", paste(dim(object), collapse=" x "), "\n",
            sep="")
    }
})

.get_VCFArraySeed_type <- function(seed, pfix, name)
{
    tp <- eval(parse(text = pfix))(seed@vcfheader)[name, "Type"]  ## FIXME: geno/info/fixed
    map <- c(Integer = "integer", String = "character", Float = "numeric")
    tp <- map[tp]
    tp
}

.get_VCFArraySeed_basic_param <- function(seed, pfix, name)
{
    if (pfix == "geno") {
        param <- ScanVcfParam(fixed = NA, info = NA)
    } else if (pfix == "fixed") {
        param <- ScanVcfParam(fixed = name, info = NA, geno = NA)
    } else if (pfix == "info") {
        param <- ScanVcfParam(fixed = NA, info = name, geno = NA)
    }
    param
}
.readVcf_for_class <- function(vcf, param, pfix, name)
{
    if(is(vcf, "VcfFile")) {
        res <- readVcf(vcf, genome = "hg19", param = param)
    } else if (is(vcf, "RangedVcfStack")) {
        res <- readVcfStack(vcf, param = param)
    }
    res <- eval(parse(text = pfix))(res)[[name]]
    if (is(res, "DataFrame")){   ## FIXME: for some "IntegerList"
                                 ## columns from info(), fails for
                                 ## dim.
        res <- res[[1]]
    }
    res
}

#' @import GenomicRanges
.extract_array_from_VCFArray <- function(x, index)
{
    ## browser()
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    pfix <- sub("/.*", "", x@name)
    name <- sub(".*/", "", x@name)

    if (any(ans_dim == 0L)){
        tp <- .get_VCFArraySeed_type(x, pfix, name)
        ans <- eval(parse(text = tp))(0)  ## return integer(0) / character(0)
        dim(ans) <- ans_dim
    } else {
        vcf <- vcffile(x)
        for(i in seq_along(index)) {
            if(is.null(index[[i]]))
                index[[i]] <- seq_len(ans_dim[i])
        }
        gr <- x@gr

        ## set basic params
        param <- .get_VCFArraySeed_basic_param(x, pfix, name)
        vcfWhich(param) <- gr[gr$pos %in% index[[1]] ]
        if (pfix == "geno" && length(ans_dim) > 1) {
            vcfSamples(param) <- colnames(x)[ index[[2]] ]
        }

        ## read array data from VcfFile/RangedVCfStack object.
        ans <- .readVcf_for_class(vcf, param, pfix, name)

        ## final touch to return array (1D) / to subset >2D arrays.
        if (length(ans_dim) == 1) {
            dim(ans) <- ans_dim
        } else if (length(ans_dim) > 2) {
            index.ext <- c(vector("list", 2), index[-c(1:2)])
            ans <- extract_array(ans, index.ext)
        }
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
    if (is(file, "VcfFile")) {
        vcf <- file
        if (!is.na(index(vcf)) && length(index)) {
            stop("'index' cannot be used when 'VcfFile' ",
                 "input already has the index file.")
        } else if (is.na(index(vcf))) {
            if (length(index)) {
                index(vcf) <- index
            } else {
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
        header <- scanVcfHeader(vcf)   ## FIXME: add the "scanVcfHeader,VcfStack".
    }
    geno <- rownames(geno(header))
    fixed <- names(fixed(header))
    info <- rownames(info(header))
    msg <- paste('The Available values for "name" argument are: \n',
               "fixed(", length(fixed), "): ", paste(fixed, collapse = " "), "\n",
               "info(", length(info), "): ", paste(info, collapse = " "), "\n",
               "geno(", length(geno), "): ", paste(geno, collapse = " "), "\n",
               sep = "")
                  ## paste(geno, collapse=" "), "\n")

    ## check "name" argument
    if (missing(name) || !name %in% c(fixed, info, geno))
        stop(msg, "Please specify corectly!")
    
    ## lightweight filter. Only return REF, rowRanges
    if (is(vcf, "RangedVcfStack")) {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA, which = rowRanges(vcf))
        readvcf <- readVcfStack(vcf, param = param)
    } else {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA)
        readvcf <- readVcf(vcf, genome = "hg19", param = param)
    }
    gr <- granges(rowRanges(readvcf)) 
    gr$pos <- seq_along(gr)

    ## check the category of geno/info/fixed
    pfix <- ifelse(name %in% geno, "geno",
            ifelse(name %in% fixed, "fixed",
            ifelse(name %in% info, "info", NULL)))
    
    ## dims
    nvars <- length(gr)
    nsamps <- length(samples(header))
    dims <- nvars
    dimnames <- list(names(gr))

    if (pfix == "geno") {
        dims[2] <- nsamps
        dimnames[[2]] <- samples(header)
        extradim <- as.integer(geno(header)[name, "Number"]) ## FIXME:
                                                             ## geno()/info()/fixed()
        if (!is.na(extradim) && extradim != 1) {
            dims <- c(dims, extradim)
            dimnames <- c(dimnames, list(as.character(seq_len(extradim))))
        }
    }
    
    new("VCFArraySeed", vcffile = vcf, vcfheader = header,
        name = paste(pfix, name, sep = "/"),
        dim = dims, dimnames = dimnames, 
        gr = gr)
}

### -------------------
### VCFArray class
### -------------------
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
    }
    else {
        seed <- VCFArraySeed(file, index = index, name = name)
    }
    DelayedArray(seed)   ## does the automatic coercion to VCFMatrix if 2-dim.
}

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
