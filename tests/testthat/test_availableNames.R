test_that("availableNames works", {
    fl <- system.file("extdata", "chr22.vcf.gz",
                      package="VariantAnnotation")
    avail <- availableNames(fl)
    expect_true(is(avail, "list"))
    expect_equal(names(avail), c("fixed", "info", "geno"))
    expect_equal(unname(lengths(avail)), c(4L, 22L, 3L))

    vcf <- VcfFile(fl)
    avail1 <- availableNames(vcf)
    expect_identical(avail, avail1)

    rgstackFile <- system.file("extdata", "rgstack.rds", package = "VCFArray")
    rgstack <- readRDS(rgstackFile)
    avail2 <- availableNames(rgstack)

    expect_identical(names(avail2), names(avail))
    expect_identical(unname(lengths(avail2)), c(4L, 26L, 9L))
})
