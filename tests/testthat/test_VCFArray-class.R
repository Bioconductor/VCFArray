test_that("VCFArraySeed constructor works", {

    ## singleString input
    fl <- system.file("extdata", "chr22.vcf.gz",
                      package="VariantAnnotation")
    seed <- VCFArraySeed(fl, name = "GT")
    expect_s4_class(seed, "VCFArraySeed")
    
    ## "VcfFile" input
    vcf <- VcfFile(fl)
    seed <- VCFArraySeed(vcf, name="DS")
    expect_true(validObject(seed))
    ## expect_equal(VCFArraySeed(fl, name="DS"), seed)
    
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

    ## "RangedVcfStack" input
    rgstackFile <- system.file("extdata", "rgstack.rda", package = "VCFArray")
    rgstack <- readRDS(rgstackFile)
    
    ## geno()
    seed <- VCFArraySeed(rgstack, name = "GT")  ## success
    expect_true(validObject(seed))
    expect_identical(dim(seed), c(1000L, 3L))

    seed1 <- VCFArraySeed(rgstack, name = "gt")
    expect_identical(seed, seed1)

    expect_error(VCFArraySeed(rgstack, name = "any"))
    expect_error(seed <- VCFArraySeed(rgstack))
    
    
    ## info()
    seed <- VCFArraySeed(fl, name = "LDAF")
    expect_equal(dim(seed), 10376L)  
        
    ## fixed()
    seed <- VCFArraySeed(fl, name = "REF")
    expect_equal(dim(seed), 10376L)
    
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
    expect_s4_class(VAsubset, "DelayedMatrix")

    va <- VCFArray(fl, name = "LDAF")
    expect_s4_class(va, "VCFArray")
    expect_equal(dim(va), 10376L)

    va <- VCFArray(fl, name = "CIEND")
    expect_true(validObject(va))

})
