test_that("VCFArraySeed constructor works", {

    ## singleString input
    fl <- system.file("extdata", "chr22.vcf.gz",
                      package="VariantAnnotation")
    vcf <- VcfFile(fl)
    seed <- VCFArraySeed(fl, name = "GT")
    expect_s4_class(seed, "VCFArraySeed")

    ##----------------
    ## VcfFile / character
    ##----------------

    ## input is "VcfFile" with index file, and the "index" argument is
    ## not NULL.
    expect_error(VCFArraySeed(vcf, index=index(vcf), name="DS"))
    
    index(vcf) <- NA
    seed <- VCFArraySeed(vcf, name="DS")
    expect_true(validObject(seed))
    
    index(vcf) <- NA
    index <- paste(path(vcf), "tbi", sep=".")
    seed <- VCFArraySeed(vcf, index = index, name="DS")
    expect_true(validObject(seed))
    expect_equal(index(vcffile(seed)), index)

    ## geno()
    seed <- VCFArraySeed(fl, name = "GT")
    expect_true(validObject(seed))
    expect_identical(dim(seed), c(10376L, 5L))

    ## seed1 <- VCFArraySeed(fl, name = "gt")
    ## expect_equal(seed, seed1)

    expect_error(VCFArraySeed(fl, name = "any"))
    expect_error(seed <- VCFArraySeed(fl))
    
    ## info()
    seed <- VCFArraySeed(fl, name = "LDAF")
    expect_equal(dim(seed), 10376L)  
        
    ## fixed()
    seed <- VCFArraySeed(fl, name = "REF")
    expect_equal(dim(seed), 10376L)
    
    ##----------------
    ## RangedVcfStack
    ##----------------

    rgstackFile <- system.file("extdata", "rgstack.rda", package = "VCFArray")
    rgstack <- readRDS(rgstackFile)

    ## geno()
    seed <- VCFArraySeed(rgstack, name = "GT")  ## success
    expect_identical(dim(seed), c(1000L, 3L))

    seed <- VCFArraySeed(rgstack, name = "AD")  ## warning... 

    hdr <- scanVcfHeader(files(rgstack)[[1]])

    ## fixed: REF, ALT, FILTER, QUAL,
    ## info: AC, all
    ## geno: all
    ## infos <- rownames(info(hdr))
    ## for (i in seq_along(infos)) {
    ##     seed <- VCFArraySeed(rgstack, name = infos[i])
    ##     va <- VCFArray(seed)
    ##     print(infos[i])
    ##     print(va)
    ## }

    ## genos <- rownames(geno(hdr))
    ## for (i in seq_along(genos)) {
    ##     seed <- VCFArraySeed(rgstack, name = genos[i])
    ##     va <- VCFArray(seed)
    ##     print(genos[i])
    ##     print(va)
    ## }

    ## 3-dim array
    seed <- VCFArraySeed(rgstack, name = "SB")
    va <- VCFArray(seed)
    expect_identical(dim(va), c(1000L, 3L, 4L))
    
    va1 <- VCFArray(rgstack, name = "SB")
    expect_identical(va, va1)
})

test_that("VCFArray constructor works", {
    fl <- system.file("extdata", "chr22.vcf.gz",
                      package="VariantAnnotation")
    seed <- VCFArraySeed(fl, name = "GT")
    va <- VCFArray(seed)
    expect_s4_class(va, "VCFMatrix")
    vasubset <- va[1:12, ]  ## simple operation degrades "VCFMatrix"
                            ## into "DelayedMatrix".
    expect_s4_class(vasubset, "DelayedMatrix")

    va <- VCFArray(fl, name = "LDAF")
    expect_s4_class(va, "VCFArray")
    expect_equal(dim(va), 10376L)

    va <- VCFArray(fl, name = "CIEND")
    expect_true(validObject(va))

})
