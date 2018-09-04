## #' availableNames
## #' @name availableNames
## #' @description this function takes the VCF file path or `VcfFile`
## #'     object as input, and returns all available data entries that
## #'     could be converted into `VCFArray` instances. The returned data
## #'     entry names can be used as input for the `name` argument in
## #'     `VCFArray` constructor. 
## #' @param file takes values for charater string (specifying the VCF
## #'     file path), \code{VcfFile} object, and \code{RangedVcfStack}
## #'     object.
## #' @return An named list of available data entries from the backend
## #'     VCF file that could be converted into a \code{VCFArray}
## #'     instance. The names are \code{geno}, \code{fixed}, and
## #'     \code{info}, which represents the category each data entry
## #'     belongs to.
## #' @export
## #' @examples
## #' fl <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
## #' avail <- availableNames(fl)
## #' avail
## #' VCFArray(fl, name = avail$info[1])

## availableNames <- function(file)
## {
##     header <- .header(file)
##     geno <- rownames(geno(header))
##     fixed <- c("REF", "ALT", "QUAL", "FILTER")
##     info <- rownames(info(header))
##     list(fixed = fixed, info = info, geno = geno)
## }
