# Sample data cannot do this part since it is only a small subset of 450k cpg.

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(RColorBrewer)
library(minfiData)
library(ChAMP)

# read CpGs position file
hg19 <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations

# read Sample info (Sample names have to be provided in the dataframe)
# In the example, age, gender and disease stage are provided as well. 
# Patients are divided into 2 groups according to their age.
targets <- readRDS("Data/Clinical_LUAD.rds")
targets$NC <- factor(targets$NC, ordered = T, levels = rev(levels(factor(targets$NC))))
targets$age_qtile <- "over_65"
targets$age_qtile[which(targets$age < 65)] <- "under_65"
targets$stage <- gsub("a", "", targets$stage)
targets$stage <- gsub("b", "", targets$stage)
targets <- subset(targets, stage == "ii" | stage == "i" | NC == "Control")

# read Beta-value
# Initially, the columns are Sample name and the rows are cg*******.
# After the data being normalized and DMPs being picked, it will be transposed when doing following analysis.
# rds file was recommended for the storage of the large matrix, it can shorten the reading time.
Beta <- readRDS("Data/LUAD.rds")
# Check the format of the beta-matrix, if your matrix have already been processed, skip the following 5 lines.
rownames(Beta) <- as.character(Beta[,1])
Beta <- Beta[,-1]
hg19 <- subset(hg19, chr %in% paste("chr", c(1:22), sep = ""))
Beta <- Beta[which(rownames(Beta) %in% rownames(hg19)),]
Beta <- as.matrix(Beta)


# Here we remove CpGs that contain NA values. It's the requirement of the normalization procedure. 
keep <- apply(Beta, 1, function(x){length(which(is.na(x))) <= 0*length(x)})
Beta <- Beta[keep,]

# Visualise what the data looks like before and after normalization
par(mfrow=c(1,2))
densityPlot(Beta[,which(targets$NC == "Control")], main="Beta (Control)", legend=FALSE)
densityPlot(Beta[,which(targets$NC == "Case")], main="Beta (Case)", legend=FALSE)
# Attention!! this requires high memory limit, if you get a error message saying that memory insufficient, use the memory.limit(32000) change limit to 32G
Beta_Norm <- champ.norm(beta=Beta, arraytype="450k", cores=4)
densityPlot(Beta[,which(targets$NC == "Control")], main="Beta (Control)", legend=FALSE)
densityPlot(Beta[,which(targets$NC == "Case")], main="Beta (Case)", legend=FALSE)

# save normalized Beta-values
saveRDS(Beta_Norm, "Data/LUAD_norm.rds")
