.header <- function(file)
{
    if (is(file, "RangedVcfStack")) {
        header <- scanVcfHeader(files(file)[[1]])
    } else {
        header <- scanVcfHeader(file)
    }
    header
}

.availableNames_msg <- function(file)
{
    avail <- vcfFields(file)
    msg <- paste(
        'The available values for "name" argument are: \n',
        "fixed(", length(avail$fixed), "): ",
        paste(avail$fixed, collapse = " "), "\n",
        "info(", length(avail$info), "): ",
        paste(avail$info, collapse = " "), "\n",
        "geno(", length(avail$geno), "): ",
        paste(avail$geno, collapse = " "), "\n",
        sep = "")
    msg
}

.pfixFun <- function(x)
{
    get(x, envir = getNamespace("VariantAnnotation"))
}

## For generating an R null object of type ... in VCF. 
.get_VCFArraySeed_type <- function(seed, pfix, name)
{
    hdr <- .header(vcffile(seed))
    if (pfix %in% c("info", "geno")) {
        tp <- .pfixFun(pfix)(hdr)[name, "Type"]
    } else if (name %in% c("REF", "ALT", "FILTER")) {
        tp <- "Character"
    } else if (name == "QUAL") {
        tp <- "Float"
    } else {
        return(NULL)
    }
    map <- c(Integer = "integer",  Float = "numeric", Flag = "character",
             String = "character", Character = "character")
    tp <- map[tp]
    tp
}

.get_VCFArraySeed_basic_param <- function(seed, pfix, name)
{
    if (pfix == "geno") {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = name, which = seed@gr)
    } else if (pfix == "info") {
        param <- ScanVcfParam(fixed = NA, info = name, geno = NA, which = seed@gr)
    } else if (pfix == "fixed" && name == "REF") {
        param <- ScanVcfParam(fixed = NA, info = NA, geno = NA, which = seed@gr)
    } else if (pfix == "fixed") {
        param <- ScanVcfParam(fixed = name, info = NA, geno = NA, which = seed@gr)
    } else {
        return(NULL)
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
    res <- .pfixFun(pfix)(res)[[name]]
    if(is(res, "XStringSetList")) {
        res <- array(res@unlistData)
    }else if (is(res, "list_OR_List")) {
        res <- array(res)
    }
    res
}
