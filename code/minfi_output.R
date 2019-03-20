library(minfi)

baseDir <- "input/"
targets <- readRDS("output/sample_key.RDS")
targets$Basename <- targets$basename
RGS_file <- "output/RGCS.RDS"
if(!file.exists(RGS_file)){
    met_data <- read.metharray.exp(targets = targets)
    saveRDS(met_data, RGS_file)
} else {
    met_data <- readRDS(RGS_file)
}

outputObjs <- list(dims = dim(met_data),
    samp_info = head(pData(met_data)),
    probe_info = head(getProbeInfo(met_data)))
saveRDS(outputObjs, file = "output/minfo_objs.RDS")

