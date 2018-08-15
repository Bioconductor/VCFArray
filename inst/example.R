library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- VcfFile(fl)
seed <- VCFArraySeed(fl, name = "GT")
va <- VCFArray(seed)

options(error=recover)
trace(extract_array, tracer = quote(str(index)), exit = quote(print(returnValue())))
va


VCFArray(seed)[1:12, ]  ## simple operation degrades "VCFMatrix" into "DelayedMatrix". 

seed <- VCFArraySeed(fl, name = "DS")
VCFArray(seed)
seed <- VCFArraySeed(fl, name = "GL")
VCFArray(seed)

bg <- system.file("extdata", "CEU_Exon.vcf.bgz", package="SeqArray")
VAseed <- VCFArraySeed(bg, name = "GT")
VAseed
dim(VAseed)
vcffile(VAseed)
DelayedArray(VAseed)
VCFArray(VAseed)
all.equal(VCFArray(VAseed), VCFArray(bg, name = "GT"))

VAseed1 <- VCFArraySeed(bg, name = "DP")
DelayedArray(VAseed1)

chr22url <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
chr22url.tbi <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"

## VcfFile(, index = ), will check for file.exists(), so must manually specify the "index" argument.
## 

seed <- VCFArraySeed(chr22url, name = "GT")
va <- VCFArray(seed)

va[1:5, 1:5]  ## upper-left
va[1:5, seq(to=ncol(va), length.out=5)]  ## upper-right, FAILED. 
va[seq(to=nrow(va), length.out=5), 1:5]  ## lower-left
va[seq(to=nrow(va), length.out=5), seq(to=ncol(va), length.out=5)]  ## lower-right
va
## <1105538 x 2504> DelayedMatrix object of type "character":
##             HG00096 HG00097 HG00099 ... NA21143 NA21144
## rs559462325 "0|0"   "0|0"   "0|0"   .   "0|0"   "0|0"  
## rs181691356 "0|0"   "0|0"   "0|0"   .   "0|0"   "0|0"  
## rs548263598 "0|0"   "0|0"   "0|0"   .   "0|0"   "0|0"  
## rs561987868 "0|0"   "0|0"   "0|0"   .   "0|0"   "0|0"  
## rs531010746 "0|0"   "0|0"   "0|0"   .   "0|0"   "0|0"  
## ...         .       .       .       .   .       .      
## rs374001814 "0|0"   "1|0"   "0|0"   .   "0|0"   "0|0"  
## rs149048580 "0|0"   "0|0"   "0|1"   .   "0|0"   "1|0"  
## rs555877612 "0|0"   "0|0"   "0|0"   .   "0|0"   "0|0"  
## rs574115117 "0|0"   "0|0"   "0|0"   .   "0|0"   "0|0"  
## rs541185110 "0|0"   "0|0"   "0|0"   .   "0|0"   "0|0"  


vcf <- VcfFile(chr22url, index=chr22url.tbi, yieldSize = 10000)
header <- scanVcfHeader(vcf)
scanVcfHeader(chr22url)  ## equivalent

### -------------------------------
## figure out how to select
### -------------------------------

library(GenomicFiles)
extdata <- system.file(package="GenomicFiles", "extdata")
files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)
names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))
stack <- VcfStack(files)
files(stack)  ## VcfFileList of length 7

files(stack)[[1]]

gr <- GRanges(c("7:1-159138000", "X:1-155270560"))
## stack[i, j], i pass into ScanVcfParam(which), j pass into ScanVcfParam(samples). 
param <- ScanVcfParam(geno="GT", which = gr, samples = colnames(stack)[c(1,3)])
rvs <- readVcfStack(stack, param = param)
geno(rvs)
gt <- geno(rvs)[["GT"]]

headers <- lapply(files(stack), scanVcfHeader)


### -------------------------------
## VcfStack as input
### -------------------------------

extdata <- system.file(package="GenomicFiles", "extdata")
files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)
names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))

## input data.frame describing the length of each sequence, coerce to
## 'Seqinfo' object
seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")
stack <- VcfStack(files, seqinfo)
gr <- as(seqinfo(stack)[rownames(stack)], "GRanges")
rgstack <- RangedVcfStack(stack, rowRanges = gr)  ## RangedVcfStack object (rowRanges() available)
param <- ScanVcfParam(fixed = NA, info = NA, geno = NA, which = rowRanges(rgstack))
readvcf <- readVcfStack(stack, param = param)

seed <- VCFArraySeed(rgstack, name = "GT")  ## success
va <- VCFArray(seed)

### -------------------------------
## figure out how to select
### -------------------------------

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
