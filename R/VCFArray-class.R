### -------------------------
### classes
### -------------------------

setClass("VCFArray", contains = "DelayedArray")
setClass("VCFArraySeed",
         contains = "Array",
         slots = c(path = "character",
                   name = "character",
                   dim = "integer",
                   dimnames = "list"))

### -------------------------
### VCFArraySeed methods
### -------------------------

setMethod("dim", "VCFArraySeed", function(x) x@dim)
setMethod("path", "VCFArraySeed", function(object) object@path)

setMethod("show", "VCFArraySeed", function(object)
{
    cat("VCFArraySeed\n",
        "file path: ", path(object), "\n",
        "array data: ", object@name, "\n",
        "dim: ", paste(dim(object), collapse=" x "), "\n",
        sep="")
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
        vcf <- VcfFile(path(x))
        header <- scanVcfHeader(vcf)
        ridx <- index[[1]]
        ## if(is.null(ridx))
        ##     cidx <- seq_along(samples(scanVcfHeader(vcf)))
        cidx <- index[[2]]
        if(is.null(cidx))
            cidx <- seq_along(samples(scanVcfHeader(vcf)))

        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA)  ## lightweight
                                                                 ## filter. Only
                                                                 ## return
                                                                 ## REF,
                                                                 ## rowRanges
        gr <- granges(rowRanges(readVcf(vcf, genome = "hg19", param = param))) 
        ## gr <- granges(rowRanges(readVcf(vcf)))  ## FIXME, readVcf() heavy? 
        gr$pos <- seq_along(gr)
        param <- ScanVcfParam(fixed = NA, info = NA, geno = x@name,
                              which = gr[gr$pos %in% ridx],
                              samples = samples(header)[cidx])
        res <- readVcf(vcf, genome = "hg19", param = param)
        ans <- geno(res)[[1]]
        ## res <- scanVcf(vcf, param = gr[gr$pos %in% ridx])
        ## ans <- lapply(seq_len(length(res)), function(i) res[[i]]$GENO[[x@name]])
        ## ans <- do.call(rbind, ans)
    }
    ans
}

setMethod("extract_array", "VCFArraySeed", .extract_array_from_VCFArray)

### ---------------------------
### VCFArraySeed constructor
### ---------------------------

#' @import VariantAnnotation
#' @importFrom Rsamtools countTabix
VCFArraySeed <- function(path = character(), name = character())
{
    if (!isSingleString(path))
        stop(wmsg(
            "'path' must be a single string specifying the path to ",
            "the vcf file where the assay data is located."))
    if (!isSingleString(name))
        stop("'name' must be a single string.")
    if(file.exists(path)) path <- normalizePath(path)  ## in base R
    
    vcf <- VcfFile(path, index = paste(path, "tbi", sep = "."))
    header <- scanVcfHeader(vcf)

    stopifnot(name %in% rownames(geno(header)))
    nsamps <- length(samples(header))
    nvars <- countTabix(vcf)[[1]]
    
    new("VCFArraySeed", path = path, name = name,
        dim = c(nvars, nsamps),
        dimnames = list(seq_len(nvars), samples(header)))
}

### --------------
### VCFArray 
### --------------

setMethod(
    "DelayedArray", "VCFArraySeed",
    function(seed) new_DelayedArray(seed, Class="VCFArray")  ## need "extract_array" to work.
    )

VCFArray <- function(path, name=NA)
{
    if (is(path, "VCFArraySeed")) {
        if (!missing(name))
            stop(wmsg(
                "VCFArray() must be called with a single argument ",
                "when passed an VCFArraySeed object"))
        seed <- path
    } else {
        if (is.na(name)) {
            vcf <- VcfFile(path, index = paste(path, "tbi", sep="."))
            header <- scanVcfHeader(vcf)
            geno <- rownames(geno(header))
            if (length(geno) == 1) {
                name <- geno
            } else {
                message(paste0('The Available values for "name" argument are: ',
                               paste(geno, collapse=" "), "\n",
                               "Please specify, otherwise, ",
                               'The default value of "GT" will be returned.', "\n"))
                name <- "GT"
            }
        }
        seed <- VCFArraySeed(path, name)
    }
    DelayedArray(seed)   ## does the automatic coercion to VCFMatrix if 2-dim.
}

### -------------------
### other classes
### -------------------

## setClass("VCFMatrix", contains=c("DelayedMatrix", "VCFArray"))
## setMethod("matrixClass", "VCFArray", function(x) "VCFMatrix")
## setAs("VCFArray", "VCFMatrix", function(from) new("VCFMatrix", from))
## setAs("VCFMatrix", "VCFArray", function(from) from)
## setAs(
##     "ANY", "VCFMatrix",
##     function(from) as(as(from, "VCFArray"), "VCFMatrix"))


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
