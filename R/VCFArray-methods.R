### -------------------------
### VCFArraySeed methods
### -------------------------

#' VCFArraySeed or VCFArray related methods, slot getters and setters.
#' @description \code{dim}, \code{dimnames}: dimension and dimnames of
#'     object contained in the VCF file.
#' @param x the \code{VCFArray} or \code{VCFArraySeed} objects.
#' @return \code{dim}: the integer vector of dimensions for
#'     \code{VCFArray} or \code{VCFArraySeed} objects.
#' @rdname VCFArray-methods
#' @exportMethod dim
#' @examples
#' fl <- system.file("extdata", "chr22.vcf.gz",
#'                   package="VariantAnnotation")
#' va <- VCFArray(fl, name = "GT")
#' dim(va)
#' dimnames(va)
#' vcffile(va)
#' seed(va)
#' dim(seed(va))

setMethod("dim", "VCFArraySeed", function(x) x@dim)

#' @rdname VCFArray-methods
#' @export
#' @return \code{dimnames}: the unnamed list of dimension names for
#'     \code{VCFArray} and \code{VCFArraySeed} objects.
 
setMethod("dimnames", "VCFArraySeed", function(x) x@dimnames)

#' @rdname VCFArray-methods
#' @export
setGeneric("vcffile", function(x) standardGeneric("vcffile"))

#' @description \code{vcffile}: extract the \code{VcfFile} object
#'     corresponding to the backend VCF file.
#' @aliases vcffile vcffile,VCFArraySeed vcffile,VCFArray
#' @return \code{vcffile}: the \code{VcfFile} object corresponding to
#'     the backend VCF file.
#' @rdname VCFArray-methods

setMethod("vcffile", "VCFArraySeed", function(x) x@vcffile)

#' @description \code{rowRanges}: extract the \code{rowRanges}
#'     information from the backend VCF file.
#' @rdname VCFArray-methods
#' @export

setMethod("rowRanges", "VCFArraySeed", function(x) x@gr) 

#' @export
#' @import VariantAnnotation
#' @import GenomicFiles
#' @importFrom Rsamtools index index<-
#' @param object the \code{VCFArraySeed} object.
#' @rdname VCFArray-methods

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

### -------------------
### VCFArray methods
### -------------------

#' @export
#' @rdname VCFArray-methods
setMethod("vcffile", "VCFArray", function(x) vcffile(seed(x)))

#' @export
#' @rdname VCFArray-methods
setMethod("rowRanges", "VCFArray", function(x) seed(x)@gr) 

