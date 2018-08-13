library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- VcfFile(fl)
seed <- VCFArraySeed(fl, name = "GT")
VCFArray(seed)
VCFArray(seed)[1:12, ]  ## simple operation degrades "VCFMatrix" into "DelayedMatrix". 

table(isIndel(readVcf(vcffile(seed))))
## FALSE  TRUE 
##  9970   406 
 table(isInsertion(readVcf(vcffile(seed))))
## FALSE  TRUE 
## 10202   174 
table(isSNV(readVcf(vcffile(seed))))
## FALSE  TRUE 
##   407  9969 
