library(minfi)
library(annmatrix)

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

# 5. Removal of probes that have > 0.01 

# since met_data is smaller now
met_detectionP <- met_detectionP[
	match(rownames(met_data), rownames(met_detectionP)),
	match(colnames(met_data), colnames(met_detectionP))
]

good_probes <- rowMeans(met_detectionP > 0.01) < 0.01

met_data <- met_data[good_probes, ]

# 6. Remove CH probes and SNP probes

met_data <- dropMethylationLoci(met_data)

# 7. Extract the 3 tables from the prepared object file

met_annmat <- annmatrix(x = getBeta(met_data),
	rann = getAnnotation(met_data),
	cann = pData(met_data))

saveRDS(met_annmat, "output/step1_annmat.RDS")

# 8. Implement the IAC outlier removal method

#removal of samples with IAC which deviate from the mean by 2sd

outlier_no = 100
sd_cutoff = -3
while(outlier_no > 0){
	cor_mat <- cor(met_annmat)
	meancorrs <- scale(colMeans(cor_mat))
	stopifnot(sum(rownames(meancorrs) != colnames(met_annmat)) == 0)
	met_annmat <- met_annmat[,meancorrs > sd_cutoff]
	outlier_no <- sum(meancorrs < sd_cutoff)
}

# 9. Add a short quality control step

#qc by relation to island
cpg_means <- rowMeans(met_annmat)
names(cpg_means) <- met_annmat@Relation_to_Island

saveRDS(met_annmat, "output/final_annmat.RDS")
