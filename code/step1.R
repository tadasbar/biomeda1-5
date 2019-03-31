library(minfi)

# 1. Annotation prep

sample_key = readRDS("output/sample_key.RDS")

age_range = range(as.numeric(sample_key$age))
sexes =  unique(sample_key$sex)
# which columns have the same value in all rows
non_redundand = !(sapply(sample_key, function(x) length(unique(x))) == 1)
sample_key <- sample_key[, non_redundand] 

colnames(sample_key)[colnames(sample_key) == "basename"] = "Basename"

# 2. RGCS prep

met_data = read.metharray.exp(target = sample_key)
saveRDS(met_data, file = "output/RGCS.RDS")

# 3. detection Ps

met_detectionP = detectionP(met_data)
met_detectionP_bad = met_detectionP > 0.01
good_samples = colMeans(met_detectionP_bad) < 0.01

met_data = met_data[, good_samples]
met_detectionP = met_detectionP[,good_samples]

# 4. Normalization

met_data = preprocessFunnorm(met_data)

# 5.
