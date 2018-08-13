readVcfStack <- function (x, i, j = colnames(x), param = ScanVcfParam()) 
{
    stopifnot(is(x, "VcfStack"))
    if ((!missing(i) || !missing(j)) && !missing(param)) 
        stop("'i' and 'j' cannot be used with 'param'")
    gr <- NULL
    if (missing(param) && missing(i) && is(x, "RangedVcfStack")) {
        gr <- rowRanges(x)
        i = intersect(names(files(x)), as.character(seqnames(gr)))
    }
    else if (missing(param) && missing(i)) {
        gr <- GRanges()
        i = names(files(x))
    }
    else if (missing(param) && is(i, "GRanges")) {
        gr <- i
        i = unique(seqnames(i))
    }
    else if (missing(param)) {
        if (is.numeric(i)) 
            i = names(files(x))[i]
        gr <- GRanges()
    }
    else {
        gr <- GRanges(vcfWhich(param))
        i = intersect(names(files(x)), as.character(seqnames(gr)))
    }
    x = x[i]
    if (is.numeric(j)) {
        j <- colnames(x)[j]
    }
    else if (!missing(param)) {
        j <- vcfSamples(param)
    }
    genome <- genome(x)
    vcfSamples(param) <- j
    vcfWhich(param) <- gr
    vcf <- lapply(names(files(x)), function(i, files, genome, 
        param) {
        file <- files[[i]]
        if (length(vcfWhich(param))) 
            vcfWhich(param) <- vcfWhich(param)[i]
        readVcf(file, genome, param)
    }, files(x), genome, param)
    do.call(rbind, vcf)
}

