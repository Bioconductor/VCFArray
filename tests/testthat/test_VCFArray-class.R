test_that("VCFArraySeed constructor works", {

    ## singString input
    fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
    seed <- VCFArraySeed(fl, name = "GT")
    expect_s4_class(seed, "VCFArraySeed")
    
    ## "VcfFile" input
    vcf <- VcfFile(fl)
    seed <- VCFArraySeed(vcf, name="DS")
    expect_true(validObject(seed))
    expect_equal(VCFArraySeed(fl, name="DS"), seed)
    
    ## input is "VcfFile" with index file, and the "index" argument is not NULL. 
    expect_error(VCFArraySeed(vcf, index=index(vcf), name="DS"))
    
    index(vcf) <- NA
    seed <- VCFArraySeed(vcf, name="DS")
    expect_true(validObject(seed))
    
    index(vcf) <- NA
    index <- paste(path(vcf), "tbi", sep=".")
    seed <- VCFArraySeed(vcf, index = index, name="DS")
    expect_true(validObject(seed))
    expect_equal(index(vcffile(seed)), index)
})

test_that("VCFArray constructor works", {
    fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
    seed <- VCFArraySeed(fl, name = "GT")
    VA <- VCFArray(seed)
    expect_s4_class(VA, "VCFMatrix")
    VAsubset <- VA[1:12, ]  ## simple operation degrades "VCFMatrix" into "DelayedMatrix". 
    expect_s4_class(VAsubset, "DelayedMatrix")


})
