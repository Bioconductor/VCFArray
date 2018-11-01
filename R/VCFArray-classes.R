### -------------------------
### classes
### -------------------------
#' @import methods
#' @import BiocGenerics
#' @import DelayedArray
#' @importFrom tools file_path_as_absolute
#' 
setClassUnion("VcfFile_OR_RangedVcfStack", c("VcfFile", "RangedVcfStack"))

setClass("VCFArraySeed",
         contains = "Array",
         slots = c(vcffile = c("VcfFile_OR_RangedVcfStack"),
                   pfix = "character",
                   name = "character",
                   dim = "integer",
                   dimnames = "list",
                   gr = "GRanges",
                   pos = "integer"))

.extract_array_from_VCFArraySeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    pfix <- x@pfix
    name <- x@name

    if (any(ans_dim == 0L)){
        tp <- .get_VCFArraySeed_type(x, pfix, name)
        ans <- get(tp)(0)  ## return integer(0) / character(0) for 0
                           ## dim.
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
        vcfWhich(param) <- gr[x@pos %in% index[[1]] ]
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

#' VCFArray constructor and coercion methods.
#'
#' @name extract_array
#' @export
#' @description \code{extract_array}: the function to extract data
#'     from a \code{VCF} file, by taking \code{VCFArraySeed} as
#'     input. This function is required by the \code{DelayedArray} for
#'     the seed contract.
#' @param x the VCFArraySeed object
#' @param index in \code{extract_array()}, an unnamed list of
#'     subscripts as positive integer vectors, one vector per
#'     dimension in \code{x}. Empty and missing subscripts
#'     (represented by \code{integer(0)} and \code{NULL} list
#'     elements, respectively) are allowed. The subscripts can contain
#'     duplicated indices. They cannot contain NAs or non-positive
#'     values.
#' @aliases extract_array,VCFArraySeed-method
#' @rdname VCFArray-classes

setMethod("extract_array", "VCFArraySeed",
          .extract_array_from_VCFArraySeed)

### ---------------------------
### VCFArraySeed constructor
### ---------------------------

#' @import GenomicRanges
#' @import S4Vectors
#' 
VCFArraySeed <- function(file, vindex = character(),
                         name = character())
{
    ## check "file" argument
    if (!(isSingleString(file) || is(file, "VcfFile_OR_RangedVcfStack"))) {
        stop("The \"file\" argument must be either character string ",
             "(indicating the VCF file path), ",
             "or \"VcfFile\" object, or \"RangedVcfStack\" object. ")
    }
    
    ## check "name" argument (case sensitive)
    avail <- vcfFields(file)
    if (missing(name) || !name %in% unname(unlist(avail)))
        stop(.availableNames_msg(file), "Please specify corectly!")

    ## check "vindex" argument
    if(isSingleString(file)) {
        file <- VcfFile(file)
    } else if (is(file, "VcfFile")) {
        if (!is.na(index(file)) && length(vindex)) {
            stop("'vindex' cannot be used when ",
                 "input already has the index file.")
        } else if (is.na(index(file))) {
            if (length(vindex)) {
                index(file) <- vindex
            } else {
                if (file.exists(path(file))) {
                    file <- indexVcf(file)
                } else {
                    stop("Please specify the ",
                         "\"vindex\" file for the remote VCF file.")
                }
            }
        }
    }

    ## lightweight filter. Only return REF, rowRanges
    if (is(file, "RangedVcfStack")) {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA,
                              which = rowRanges(file))
        readvcf <- readVcfStack(file, param = param)
    } else {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA)
        readvcf <- readVcf(file, genome = "hg19", param = param)
    }
    gr <- granges(rowRanges(readvcf))

    pos <- seq_along(gr)
    
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
        ## convert into integer, "G/A/R/." will be NA. 
        if (!is.na(extradim) && extradim != 1) {
            dims <- c(dims, extradim)
            dimnames <- c(dimnames,
                          list(as.character(seq_len(extradim))))
        }
    }
    
    new("VCFArraySeed",
        vcffile = file,
        pfix = pfix,
        name = name,
        dim = dims, dimnames = dimnames, 
        gr = gr,
        pos = pos
        )
}

### --------------------------------
### VCFArray and VCFMatrix objects
### --------------------------------

### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user
### will see and manipulate VCFArray and VCFMatrix objects instead of
### DelayedArray and DelayedMatrix objects.

#' @exportClass VCFArray
#' @rdname VCFArray-classes
#' @aliases VCFArray-class matrixClass,VCFArray-method
#' @param file takes values for charater string (specifying the VCF
#'     file path), \code{VcfFile} object, and \code{RangedVcfStack}
#'     object.
#' @param vindex in \code{VCFArray()}, the character string specifying
#'     the index file path. This argument is required if an remote VCF
#'     file is used for the \code{file} argument.
#' @param name the data entry from VCF file to be read into
#'     VCFArraySeed / VCFArray. For \code{VCFArray}. This argument
#'     should always be specified.
#' @return \code{VCFArray} class object.

setClass("VCFArray", contains = "DelayedArray")

#' @name VCFMatrix
#' @exportClass VCFMatrix
#' @aliases VCFMatrix-class
#' @rdname VCFArray-classes

setClass("VCFMatrix", contains=c("DelayedMatrix", "VCFArray"))

## for internal use only.
setMethod("matrixClass", "VCFArray", function(x) "VCFMatrix")

### Automatic coercion method from VCFArray to VCFMatrix (muted for
### higher dimensions) this function works only when VCFArray is
### 2-dimensional, otherwise, it fails.

#' @name coerce
#' @export
#' @aliases coerce,VCFArray,VCFMatrix-method
#'     coerce,VCFMatrix,VCFArray-method coerce,ANY,VCFMatrix-method
#' @rdname VCFArray-classes
#' 
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
    function(seed) new_DelayedArray(seed, Class="VCFArray")
    ## need "extract_array" to work.
    )

#' @description \code{VCFArray}: The function to convert data entries
#'     inside VCF file into the \code{VCFArray} instance.
#' @export
#' @aliases VCFArray-method
#' @rdname VCFArray-classes
#' @examples
#' fl <- system.file("extdata", "chr22.vcf.gz",
#'                   package="VariantAnnotation")
#' va <- VCFArray(fl, name = "GT")
#' va
#' vcf <- VariantAnnotation::VcfFile(fl)
#' va1 <- VCFArray(vcf, name = "GT")
#' va1
#' all.equal(va, va1)
#'
#' \dontrun{
#' ## RangedVcfStack class
#' library(GenomicFiles)
#' extdata <- system.file(package="GenomicFiles", "extdata")
#' files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)
#' names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))
#' seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")
#' stack <- VcfStack(files, seqinfo)
#' gr <- as(seqinfo(stack)[rownames(stack)], "GRanges")
#' rgstack <- RangedVcfStack(stack, rowRanges = gr)  ## RangedVcfStack
#' va2 <- VCFArray(rgstack, name = "SB")
#' va2
#' }
#' ## coercion
#' as(va[1:10, ], "array")


VCFArray <- function(file, vindex = character(),
                     name=NA)
{
    if (is(file, "VCFArraySeed")) {
        if (!missing(name))
            stop(wmsg(
                "VCFArray() must be called with a single argument ",
                "when passed an VCFArraySeed object"))
        seed <- file
    }
    else {
        seed <- VCFArraySeed(file, vindex = vindex, name = name)
    }
    DelayedArray(seed)   ## does the automatic coercion to VCFMatrix
                         ## if 2-dim.
}

