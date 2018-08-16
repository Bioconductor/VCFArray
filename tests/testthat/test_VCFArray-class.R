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
    extdata <- system.file(package="GenomicFiles", "extdata")
    files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)
    names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))
    seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")
    stack <- VcfStack(files, seqinfo)
    gr <- as(seqinfo(stack)[rownames(stack)], "GRanges")
    rgstack <- RangedVcfStack(stack, rowRanges = gr)  ## RangedVcfStack
                                                      ## object
                                                      ## (rowRanges()
                                                      ## available)
    
    seed <- VCFArraySeed(rgstack, name = "GT")  ## success
    expect_true(validObject(seed))
    expect_identical(dim(seed), c(1000L, 3L))

    expect_error(VCFArraySeed(rgstack, name = "any"))
    expect_error(seed <- VCFArraySeed(rgstack))

    ## Fixed / info
    seed <- VCFArraySeed(fl, name = "LDAF")
    expect_equal(dim(seed), 10376L)  

    seed <- VCFArraySeed(fl, name = "AVGPOST")

    vcfinfo <- as.data.frame(info(seed@vcfheader)[,1:2])
    infos <- rownames(vcfinfo)
    for (i in seq_along(infos)) {
        seed <- VCFArraySeed(fl, name = infos[i])
        va <- VCFArray(seed)
        print(va)
    }

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
    VA <- VCFArray(seed)
    expect_s4_class(VA, "VCFMatrix")
    VAsubset <- VA[1:12, ]  ## simple operation degrades "VCFMatrix"
                            ## into "DelayedMatrix".
    expect_s4_class(VAsubset, "DelayedMatrix")

    va <- VCFArray(fl, name = "LDAF")
    expect_s4_class(va, "VCFArray")
    expect_equal(dim(va), 10376L)

})
