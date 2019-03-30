## Sample key production

library(recount)


filenames <- list.files('input/')
filenames <- filenames[grep(filenames, pattern="GSM")]

basenames <- list.files("input/", full.names=TRUE)
basenames <- basenames[grep(basenames, pattern="GSM")]
basenames <- sub(basenames, pattern="_....idat.gz", replacement="")
basenames <- basenames[!duplicated(basenames)]


gsms <- sub(filenames, pattern="(GSM.*?)_.*", replacement="\\1")
gsms <- gsms[!duplicated(gsms)]

sample_list_file <- "output/sample_list.RDS"
if(!file.exists(sample_list_file)){
	sample_list <- lapply(gsms, function(x) geo_characteristics(geo_info(x)))
	saveRDS(sample_list, file = sample_list_file)
} else {
	sample_list <- readRDS(sample_list_file)
}

sample_key <- Reduce(sample_list, f = rbind)
sample_key$gsm <- gsms
sample_key$basename <- basenames
saveRDS(sample_key, "output/sample_key.RDS")
