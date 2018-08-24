## generate an RangedVcfStack object as example for parsing raw data into VCFArray. 

extdata <- system.file(package="GenomicFiles", "extdata")
files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)
names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))
seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")
stack <- VcfStack(files, seqinfo)
gr <- as(seqinfo(stack)[rownames(stack)], "GRanges")
rgstack <- RangedVcfStack(stack, rowRanges = gr)  ## RangedVcfStack
saveRDS(rgstack, file = "inst/extdata/rgstack.rds")
