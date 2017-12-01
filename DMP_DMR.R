library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

#######################################################
# If you are with a normalized data, start from here! #
#######################################################
################# Preliminary Analysis ################
targets <- readRDS("Data/Clinical_LUAD.rds")
targets$NC <- factor(targets$NC, ordered = T, levels = rev(levels(factor(targets$NC))))
targets$age_qtile <- "over_65"
targets$age_qtile[which(targets$age < 65)] <- "under_65"
targets$stage <- gsub("a", "", targets$stage)
targets$stage <- gsub("b", "", targets$stage)
targets <- subset(targets, stage == "ii" | stage == "i" | NC == "Control")
Beta <- readRDS("Data/LUAD_norm.rds")
# make sure that your sample info and beta value have the same order
Beta <- Beta[,targets$Sample]
# change Beta-value to M-value, the M-value is better for data analysis, and beta is better for interpretation.
M <- log(Beta/(1-Beta))

# Some samples with strange data distribution, remove them
densityPlot(M, sampGroups=targets$NC, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$NC)), 
       text.col=brewer.pal(8,"Dark2"))
donkeep <- which(apply(Beta, 2, function(x){length(which(x >= 0.999999 | x <= 0.000001))}) > nrow(Beta)*0.001)
if(length(donkeep) != 0){
  Beta <- Beta[,-donkeep]
  M <- M[,-donkeep]
  targets <- targets[-donkeep,]
}
densityPlot(M, sampGroups=targets$NC, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$NC)), 
       text.col=brewer.pal(8,"Dark2"))


pal <- brewer.pal(8,"Dark2")
# MDS plots to look at largest sources of variation
# option dim.plot = c(1,2) means to see the first 2 priciple components
par(mfrow=c(1,2))
plotMDS(M, top=1000, gene.selection="common", 
        col=pal[factor(targets$NC)], pch = 10, dim.plot = c(1,2))
legend("topleft", legend=levels(factor(targets$NC)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$gender)], pch = 10)
legend("topleft", legend=levels(factor(targets$gender)), text.col=pal,
       bg="white", cex=0.85)

plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$age_qtile)], pch = 10)
legend("topleft", legend=levels(factor(targets$age_qtile)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$stage)], pch = 10)
legend("topleft", legend=levels(factor(targets$stage)), text.col=pal,
       bg="white", cex=0.85)


# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(M, top=1000, gene.selection="common", 
        col=pal[factor(targets$NC)], pch = 10, dim.plot = c(1,3))
legend("topleft", legend=levels(factor(targets$NC)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common", 
        col=pal[factor(targets$NC)], pch = 10, dim.plot = c(2,3))
legend("topleft", legend=levels(factor(targets$NC)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common", 
        col=pal[factor(targets$NC)], pch = 10, dim.plot = c(3,4))
legend("topleft", legend=levels(factor(targets$NC)), text.col=pal,
       bg="white", cex=0.85)

plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$gender)], pch = 10, dim.plot = c(1,3))
legend("topleft", legend=levels(factor(targets$gender)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$gender)], pch = 10, dim.plot = c(2,3))
legend("topleft", legend=levels(factor(targets$gender)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$gender)], pch = 10, dim.plot = c(3,4))
legend("topleft", legend=levels(factor(targets$gender)), text.col=pal,
       bg="white", cex=0.85)

plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$age_qtile)], pch = 10, dim.plot = c(1,3))
legend("topleft", legend=levels(factor(targets$age_qtile)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$age_qtile)], pch = 10, dim.plot = c(2,3))
legend("topleft", legend=levels(factor(targets$age_qtile)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$age_qtile)], pch = 10, dim.plot = c(3,4))
legend("topleft", legend=levels(factor(targets$age_qtile)), text.col=pal,
       bg="white", cex=0.85)

plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$stage)], pch = 10, dim.plot = c(1,3))
legend("topleft", legend=levels(factor(targets$stage)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$stage)], pch = 10, dim.plot = c(2,3))
legend("topleft", legend=levels(factor(targets$stage)), text.col=pal,
       bg="white", cex=0.85)
plotMDS(M, top=1000, gene.selection="common",  
        col=pal[factor(targets$stage)], pch = 10, dim.plot = c(3,4))
legend("topleft", legend=levels(factor(targets$stage)), text.col=pal,
       bg="white", cex=0.85)



###################### DMPs part ######################
# In fact, following line is the same as cpg.annotate()

# this is the factor of interest (here is Sample and Control)
cellType <- factor(targets$NC)
# this is the individual effect that we need to account for (Age group)
individual <- factor(targets$age_qtile) 
# use the above to create a design matrix
design <- model.matrix(~0+cellType+individual, data=targets)
colnames(design) <- c(levels(cellType),levels(individual)[-1])

# fit the linear model 
fit <- lmFit(M, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(Case-Control,
                            levels=design)
# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, coef = "Case - Control", number = nrow(M))

# fdr in betafit are not used, we only extract the beta-fold-change for each CpG
betafit <- lmFit(ilogit2(M), design)
betafit <- contrasts.fit(betafit, contMatrix)
betafit <- eBayes(betafit)
betatt <- topTable(betafit, coef = "Case - Control", number = nrow(M))
# use fdr calculated from M-values and betafc calculated from beta-value
m <- match(rownames(tt), rownames(betatt))
tt$betafc <- betatt$logFC[m]
m <- match(rownames(M), rownames(tt))
tt <- tt[m, ]

# Annotate the DMPs position
grset <- makeGenomicRatioSetFromMatrix(M, 
                                       array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                       mergeManifest = TRUE, what = "M")
anno <- getAnnotation(grset)
stat <- tt$t
annotated <- data.frame(ID = rownames(tt), stat = stat, 
                        CHR = anno$chr, pos = anno$pos, betafc = tt$betafc, 
                        indfdr = tt$adj.P.Val, is.sig = tt$adj.P.Val < 0.05)
class(annotated) <- "annot"
myAnnotation <- annotated
saveRDS(myAnnotation, "Results/myAnnotation_LUAD.rds")

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# get the table of results for the first contrast (naive - rTreg)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(M),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)
saveRDS(DMPs, "Results/DMPs_LUAD.rds")
# Attention! betafc in DMPs is not correct(go back to myAnnotation and get the right numbers)! 
# it is just for seeing the annotation of each CpGs!


# plot the top 100 most significantly differentially methylated CpGs 
par(mfrow=c(4,5))
sapply(rownames(DMPs)[1:100], function(cpg){
  plotCpg(Beta, cpg=cpg, pheno=targets$NC, ylab = "Beta values")
})

###################### DMRs part ######################
# recomended lambda and C for 450k are 1000 and 2.
# lambda is the maxinum distance between 2 DMPs.
# min.cpgs is the minumun number of DMPs in a DMR.
# there will be a new fdr calculated by DMR, based on the fdr of DMPs in a DMR.
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, min.cpgs = 2)
head(DMRs$results)
nrow(DMRs$results)
# Total numbers of DMPs covered in each DMR
sum(DMRs$results$no.cpgs)
# Generate the bed for DMRs
starts <- as.integer(unlist(strsplit(unlist(strsplit(DMRs$results$coord, split = ":"))[seq(2,(2*nrow(DMRs$results)),2)], 
                                     split = "-"))[seq(1,(nrow(DMRs$results)*2-1),2)])
ends <- as.integer(unlist(strsplit(unlist(strsplit(DMRs$results$coord, split = ":"))[seq(2,(2*nrow(DMRs$results)),2)], 
                                   split = "-"))[seq(2,(nrow(DMRs$results)*2),2)])
chr <- unlist(strsplit(DMRs$results$coord, split = ":"))[seq(1,(2*nrow(DMRs$results)-1),2)]
# total length of DMRs
sum(ends-starts+1)
bed_LUAD <- data.frame(chr = chr, start = starts, end = (ends+1))
write.table(bed_LUAD, "Results/bed_DMR_LUAD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# If you only want to keep DMRs with abs(log(mean-betafc)) over 0.1
# DMRs$results <- DMRs$results[which(abs(DMRs$results$meanbetafc)>0.1),]
saveRDS(DMRs, "Results/DMRs_LUAD.rds")

# convert the regions to annotated genomic ranges
data(dmrcatedata)
results.ranges <- extractRanges(DMRs, genome = "hg19")
saveRDS(results.ranges, "Results/results.ranges_LUAD.rds")


# Save the first 1000 DMPs for variable selection
Beta_new <- Beta[rownames(DMPs)[1:1000],]
M_new <- M[rownames(DMPs)[1:1000],]
saveRDS(Beta_new, "Results/Beta_LUAD_top1000.rds")
saveRDS(M_new, "Results/M_LUAD_top1000.rds")

