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

.extract_array_from_VCFArray <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    
    tp <- ifelse(x@name == "GT", "character", "integer")  ## FIXME
    if (any(ans_dim == 0L)){
        ans <- eval(parse(text = tp))(0)  ## return integer(0) / character(0)
        dim(ans) <- ans_dim
    } else {
        vcf <- VcfFile(path(x))
        ridx <- index[[1]]  ## FIXME: here ignoring the column index. 
        gr <- granges(rowRanges(readVcf(vcf)))  ## FIXME, readVcf() heavy? 
        gr$pos <- seq_along(gr)
        res <- scanVcf(vcf, param = gr[gr$pos %in% ridx])
        ans <- lapply(seq_len(length(res)), function(i) res[[i]]$GENO[[x@name]])
        ans <- do.call(rbind, ans)
        cidx <- index[[2]]
        ans <- ans[, cidx]
    }
    ans
}

setMethod("extract_array", "VCFArraySeed", .extract_array_from_VCFArray)

### ---------------------------
### VCFArraySeed constructor
### ---------------------------

VCFArraySeed <- function(path = character(), name = character())
{
    vcf <- VcfFile(path)
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
            vcf <- VcfFile(path)
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_VCFArray <- function(x)
{
    if (!is(x@seed, "VCFArraySeed"))
        return(wmsg("'x@seed' must be a VCFArraySeed object"))
    TRUE
}

setValidity2("VCFArray", .validate_VCFArray)
