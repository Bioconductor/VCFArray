fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
seed <- VCFArraySeed(fl, "GT")
VCFArray(seed)
seed <- VCFArraySeed(fl, "DS")
VCFArray(seed)
seed <- VCFArraySeed(fl, "GL")
VCFArray(seed)


bg <- system.file("extdata", "CEU_Exon.vcf.bgz", package="SeqArray")
VAseed <- VCFArraySeed(bg, "GT")
VAseed
dim(VAseed)
path(VAseed)
DelayedArray(VAseed)
identical(VCFArray(VAseed), VCFArray(bg, "GT"))

VAseed1 <- VCFArraySeed(bg, "DP")
VAseed1
dim(VAseed1)
path(VAseed1)
DelayedArray(VAseed1)

fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
VcfFile(fl)
seed <- VCFArraySeed(fl, "GT")
VCFArray(seed)

chr22url <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
chr22url.tbi <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"

seed <- VCFArraySeed(chr22url, "GT")
va <- VCFArray(seed)

vcf <- VcfFile(chr22url, index=chr22url.tbi, yieldSize = 10000)
header <- scanVcfHeader(vcf)

## figure out how to select 
param <- ScanVcfParam(info = sample(rownames(info(header)), 2), geno = "GT", samples= samples(header)[1:10])
param <- ScanVcfParam(info = NA, fixed = NA, geno = "GT")
a <- scanVcf(vcf, param = param)

param <- ScanVcfParam(fixed = NA, geno=NA, info=NA)
a <- readVcf(vcf, param = param)  ## return only rowRanges, REF
## how to extract the a$rowRanges? Could be used in matching the row indexes in "extract_array".
## 
rowRanges(a)

param <- ScanVcfParam(fixed = NA, info = NA, geno = "GT")
a <- readVcf(vcf, param = param)
geno(a)

a$REF
a$ALT
a$QUAlL
a$FILTER
a$INFO
dim(a$<NA>$GENO$GT)
