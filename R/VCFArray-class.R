### -------------------------
### classes
### -------------------------
#' @import methods DelayedArray 
setClassUnion("VcfFile_OR_RangedVcfStack", c("VcfFile", "RangedVcfStack"))

setClass("VCFArraySeed",
         contains = "Array",
         slots = c(vcffile = c("VcfFile_OR_RangedVcfStack"),
                   ## vcfheader = "VCFHeader",
                   name = "character",
                   dim = "integer",
                   dimnames = "list",
                   gr = "GRanges"))

### -------------------------
### VCFArraySeed methods
### -------------------------

#' @export
setMethod("dim", "VCFArraySeed", function(x) x@dim)
#' @export
setMethod("dimnames", "VCFArraySeed", function(x) x@dimnames)
#' @export
setGeneric("vcffile", function(x) standardGeneric("vcffile"))
setMethod("vcffile", "VCFArraySeed", function(x) x@vcffile)
#' @export
setMethod("rowRanges", "VCFArraySeed", function(x) x@gr)

.header <- function(file)
{
    if (is(file, "RangedVcfStack")) {
        header <- scanVcfHeader(files(file)[[1]])
    } else {
        header <- scanVcfHeader(file)   ## FIXME: add the "scanVcfHeader,VcfStack".
    }
    header
}
    

#' @export
availableNames <- function(file)
{
    header <- .header(file)
    geno <- rownames(geno(header))
    ## fixed <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
    fixed <- c("REF", "ALT", "QUAL", "FILTER")
    info <- rownames(info(header))
    list(fixed = fixed, info = info, geno = geno)
}

.availableNames_msg <- function(file)
{
    avail <- availableNames(file)
    msg <- paste('The Available values for "name" argument are: \n',
                 "fixed(", length(avail$fixed), "): ", paste(avail$fixed, collapse = " "), "\n",
                 "info(", length(avail$info), "): ", paste(avail$info, collapse = " "), "\n",
                 "geno(", length(avail$geno), "): ", paste(avail$geno, collapse = " "), "\n",
                 sep = "")
    msg
}

#' @export
#' @import VariantAnnotation GenomicFiles
#' @importFrom Rsamtools index index<-
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

## For generating an R null object of type ... in VCF. 
.get_VCFArraySeed_type <- function(seed, pfix, name)
{
    hdr <- .header(vcffile(seed))
    if (pfix %in% c("info", "geno")) {
        tp <- eval(parse(text = pfix))(hdr)[name, "Type"] 
    } else if (name %in% c("REF", "ALT", "FILTER")) {
        tp <- "Character"
    } else if (name == "QUAL") {
        tp <- "Float"
    }
    map <- c(Integer = "integer",  Float = "numeric", Flag = "character",
             String = "character", Character = "character")
    tp <- map[tp]
    tp
}

.get_VCFArraySeed_basic_param <- function(seed, pfix, name)
{
    if (pfix == "geno") {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = name)
    } else if (pfix == "info") {
        param <- ScanVcfParam(fixed = NA, info = name, geno = NA)
        ##    } else if (pfix == "fixed" && name %in% c("CHROM", "POS", "ID", "REF")) {
    } else if (pfix == "fixed" && name == "REF") {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA)
    } else if (pfix == "fixed") {
        param <- ScanVcfParam(fixed = name, info = NA, geno = NA)
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
    if(is(res, "XStringSetList")) {
        res <- array(res@unlistData)
    }else if (is(res, "list_OR_List")) {
        res <- array(res)
    }
    res
}

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
        ## if (name %in% c("CHROM", "POS", "ID", "REF")) {
        ## if (name == "REF") {
        ##     ans <- mcols(gr)[[name]][ index[[1]] ]
        ## } else {
        ## set basic params
        param <- .get_VCFArraySeed_basic_param(x, pfix, name)
        vcfWhich(param) <- gr[gr$pos %in% index[[1]] ]
        if (pfix == "geno" && length(ans_dim) > 1) {
            vcfSamples(param) <- colnames(x)[ index[[2]] ]
        }        
        ## read array data from VcfFile/RangedVCfStack object.
        ans <- .readVcf_for_class(vcf, param, pfix, name)
        ##}
        
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

#' @export
setMethod("extract_array", "VCFArraySeed", .extract_array_from_VCFArray)

### ---------------------------
### VCFArraySeed constructor
### ---------------------------

#' @import GenomicRanges S4Vectors
#' 
VCFArraySeed <- function(file, index = character(), name = character())
{
    ## browser()
    if(isSingleString(file)) {
        file <- VcfFile(file)
    } else if (is(file, "VcfFile")) {
        if (!is.na(index(file)) && length(index)) {
            stop("'index' cannot be used when ",
                 "input already has the index file.")
        } else if (is.na(index(file))) {
            if (length(index)) {
                index(file) <- index
            } else {
                if (file.exists(path(file))) {
                    file <- indexVcf(file)
                } else {
                    stop("Please specify the \"index\" file for the remote VCF file.")
                }
            }
        }
    }

    avail <- availableNames(file)
    ## check "name" argument (case sensitive)
    if (missing(name) || !name %in% unname(unlist(avail)))
        stop(.availableNames_msg(file), "Please specify corectly!")

    ## lightweight filter. Only return REF, rowRanges
    if (is(file, "RangedVcfStack")) {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA, which = rowRanges(file))
        readvcf <- readVcfStack(file, param = param)
    } else {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA)
        readvcf <- readVcf(file, genome = "hg19", param = param)
    }
    gr <- granges(rowRanges(readvcf))
    gr$pos <- seq_along(gr)
    ## gr <- rowRanges(readvcf)
    ## mcols(gr) <- DataFrame(REF = mcols(gr)$REF, pos = seq_along(gr))
    
    ## check the category of geno/info/fixed
    pfix <- ifelse(name %in% avail$geno, "geno",
            ifelse(name %in% avail$fixed, "fixed",
            ifelse(name %in% avail$info, "info", NULL)))
    
    ## header
    header <- .header(file)
    
    ## dims
    nvars <- length(gr)
    nsamps <- length(samples(header))
    dims <- nvars
    dimnames <- list(names(gr))

    if (pfix == "geno") {
        dims[2] <- nsamps
        dimnames[[2]] <- samples(header)

        extradim <- as.integer(geno(header)[name, "Number"]) 
        if (!is.na(extradim) && extradim != 1) {
            dims <- c(dims, extradim)
            dimnames <- c(dimnames, list(as.character(seq_len(extradim))))
        }
    }
    
    new("VCFArraySeed",
        vcffile = file,
        ## vcfheader = header,
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

### -------------------
### VCFArray methods
### -------------------
setMethod("vcffile", "VCFArray", function(x) vcffile(seed(x)))

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

#' @export
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
