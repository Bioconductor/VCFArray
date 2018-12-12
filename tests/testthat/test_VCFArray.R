test_that("VCFArraySeed arguments check works", {

    ## singleString input
    fl <- system.file("extdata", "chr22.vcf.gz",
                      package="VariantAnnotation")
    vcf <- VcfFile(fl)
    seed <- VCFArraySeed(fl, name = "GT")
    expect_s4_class(seed, "VCFArraySeed")

    ## "index": input is indexed "VcfFile", and the "index" argument
    ## is not NULL.
    expect_error(VCFArraySeed(vcf, index=index(vcf), name="DS"))
    
    index(vcf) <- NA
    seed <- VCFArraySeed(vcf, name="DS")
    expect_true(validObject(seed))
    
    index(vcf) <- NA
    index <- paste(path(vcf), "tbi", sep=".")
    seed <- VCFArraySeed(vcf, vindex = index, name="DS")
    expect_true(validObject(seed))
    expect_equal(index(vcffile(seed)), index)

    ## "name"
    expect_error(VCFArraySeed(fl, name = "notValidName"))
    expect_error(VCFArraySeed(fl))

}) 

test_that("VCFArraySeed and VCFArray constructor works", {

    fl <- system.file("extdata", "chr22.vcf.gz",
                      package="VariantAnnotation")
    ## geno()
    seed <- VCFArraySeed(fl, name = "GT")
    expect_true(validObject(seed))
    expect_identical(dim(seed), c(10376L, 5L))

    ## info()
    seed <- VCFArraySeed(fl, name = "LDAF")
    expect_equal(dim(seed), 10376L)  
    va <- VCFArray(seed)
    expect_true(validObject(va))
    expect_identical(dim(seed), dim(va))
    
    ## fixed()
    seed <- VCFArraySeed(fl, name = "REF")
    expect_equal(dim(seed), 10376L)
    va <- VCFArray(seed)
    va1 <- VCFArray(fl, name = "REF")
    expect_equal(va, va1)

    ## XStringSetList, IntegerList...
    va <- VCFArray(fl, name = "CIEND")
    expect_true(validObject(va))
    ## DelayedArray::type(va), "IntegerList"
    
    ##----------------
    ## RangedVcfStack
    ##----------------

    extdata <- system.file(package = "GenomicFiles", "extdata")
    files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)[1:2]
    names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))
    seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")
    stack <- GenomicFiles::VcfStack(files, seqinfo)
    gr <- as(GenomicFiles::seqinfo(stack)[rownames(stack)], "GRanges")
    rgstack <- GenomicFiles::RangedVcfStack(stack, rowRanges = gr)  
    
    ## geno()
    seed <- VCFArraySeed(rgstack, name = "GT")
    expect_identical(dim(seed), c(318L, 3L))

    ## fixed()
    seed <- VCFArraySeed(rgstack, name = "FILTER")
    expect_identical(dim(seed), 318L)
    va <- VCFArray(seed)
    expect_s4_class(va, "VCFArray")

    ## info()
    seed <- VCFArraySeed(rgstack, name = "set")
    expect_identical(dim(seed), 318L)

    ## 3-dim array
    seed <- VCFArraySeed(rgstack, name = "SB")
    va <- VCFArray(seed)
    expect_identical(dim(va), c(318L, 3L, 4L))
    expect_s4_class(va, "VCFArray")
    
    va1 <- VCFArray(rgstack, name = "SB")
    expect_identical(va, va1)
})

