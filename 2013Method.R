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
library(plyr)

########### define DMRs function ###########
# the origional function only count the number of DMPs covered by DMRs
# here we also count the number of Hyper and Hypo-methylation in DMRs
dmrcate_redefine <- function (object, lambda = 1000, C = NULL, p.adjust.method = "BH", 
                              pcutoff = "fdr", consec = FALSE, conseclambda = 10, betacutoff = NULL, 
                              min.cpgs = 2, mc.cores = 1) 
{
  stopifnot(is(object, "annot"))
  stopifnot(lambda >= 1)
  stopifnot(pcutoff == "fdr" | (0 <= pcutoff & pcutoff <= 1))
  stopifnot(C >= 0.2)
  if (consec & is.null(conseclambda)) {
    stop("Consecutive CpG bandwidth must be specified")
  }
  object <- data.frame(ID = object$ID, weights = abs(object$stat), 
                       CHR = as.character(object$CHR), pos = object$pos, betafc = object$betafc, 
                       indfdr = object$indfdr, is.sig = object$is.sig, is.hyper = object$is.hyper, is.hypo = object$is.hypo)
  object <- object[order(object$CHR, object$pos), ]
  if (is.null(C) & !consec) {
    if (nrow(object) < 9e+05) {
      C = 2
    }
    else {
      C = 50
    }
  }
  if (consec) {
    lambda = conseclambda
    message(paste("Consecutive mode specified, lambda is now set at", 
                  conseclambda, "consecutive CpGs."))
    if (is.null(C)) {
      stop("Error: argument C must be specified (in CpG sites) for consecutive mode.")
    }
    object$realcoordforconsec <- object$pos
    object$pos <- unlist(sapply(as.numeric(table(object$CHR)), 
                                function(x) 1:x))
  }
  lag = lambda
  chr.unique <- unique(c(as.character(object$CHR)))
  fitted <- mclapply(chr.unique, fitParallel, object = object, 
                     consec = consec, conseclambda = conseclambda, lambda = lambda, 
                     C = C, mc.cores = mc.cores)
  object <- rbind.fill(fitted)
  object$fdr <- p.adjust(object$raw, method = p.adjust.method)
  if (pcutoff == "fdr") {
    nsig <- sum(object$is.sig)
    if (nsig == 0) {
      txt <- "The FDR you specified in cpg.annotate() returned no significant CpGs, hence there are no DMRs.\n    Try specifying a value of 'pcutoff' in dmrcate() and/or increasing 'fdr' in cpg.annotate()."
      stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
    }
    pcutoff <- sort(object$fdr)[nsig]
  }
  object$sig <- (object$fdr <= pcutoff)
  if (nrow(object) == 0) {
    txt <- "No signficant regions found. Try increasing the value of\n    'pcutoff' in dmrcate() and/or 'fdr' in cpg.annotate()."
    stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
  }
  message("Demarcating regions...")
  chr.N <- as.character(object$CHR)
  pos.N <- object$pos
  sig.N <- object$sig
  N <- length(sig.N)
  n.K <- which(sig.N)
  K <- length(n.K)
  stopifnot(K >= 2)
  pos.K <- pos.N[n.K]
  chr.K <- chr.N[n.K]
  jump_chr.k <- (chr.K[-1] != chr.K[-K])
  jump_pos.k <- (diff(pos.K) > lag)
  jump.k <- (jump_chr.k | jump_pos.k)
  ksegments.A2 <- Segment(jump.k)
  A <- nrow(ksegments.A2)
  kstart.A <- ksegments.A2[, "start"]
  kend.A <- ksegments.A2[, "end"]
  realpos.K <- pos.K
  if (consec) {
    realpos.N <- object$realcoordforconsec
    realpos.K <- realpos.N[n.K]
  }
  start.A <- realpos.K[kstart.A]
  end.A <- realpos.K[kend.A]
  chr.A <- chr.K[kstart.A]
  stopifnot(all(chr.K[kend.A] == chr.A))
  fmt <- "%s:%1d-%1d"
  coord.A <- sprintf(fmt, chr.A, start.A, end.A)
  nstart.A <- n.K[kstart.A]
  nend.A <- n.K[kend.A]
  width.A <- nend.A + 1 - nstart.A
  a.Z <- rep(seq(A), width.A)
  fn <- function(a) seq(from = nstart.A[a], to = nend.A[a])
  l.listA <- lapply(seq(A), fn)
  n.Z <- unlist(l.listA)
  region.N <- rep(NA_integer_, N)
  region.N[n.Z] <- a.Z
  levels <- seq(A)
  region.N <- factor(region.N, levels = levels)
  no_cpg.A <- c(table(region.N))
  region.hyper <- region.N[which(object$is.hyper == T)]
  region.hypo <- region.N[which(object$is.hypo == T)]
  no_cpg.hyper <- c(table(region.hyper))
  no_cpg.hypo <- c(table(region.hypo))
  REGIONSTAT <- function(field, fn) {
    x.N <- object[[field]]
    x.R <- tapply(x.N, region.N, fn)
    c(x.R)
  }
  fn_Stouffer <- function(x) pnorm(sum(qnorm(x))/sqrt(length(x)))
  fn_max <- function(x) x[which.max(abs(x))]
  results <- data.frame(coord = coord.A, no.cpgs = no_cpg.A, no.hyper =  no_cpg.hyper, no.hypo =  no_cpg.hypo,
                        minfdr = REGIONSTAT("fdr", min), Stouffer = REGIONSTAT("indfdr", 
                                                                               fn_Stouffer), maxbetafc = REGIONSTAT("betafc", fn_max), 
                        meanbetafc = REGIONSTAT("betafc", mean), row.names = seq(A), 
                        stringsAsFactors = FALSE)
  keep <- (results$no.cpgs >= min.cpgs)
  results <- results[keep, ]
  if (!is.null(betacutoff)) {
    keep <- (abs(results$meanbetafc) >= betacutoff)
    results <- results[keep, ]
  }
  o <- order(results$Stouffer, -results$no.cpgs)
  results <- results[o, ]
  message("Done!")
  output <- NULL
  output$input <- object
  output$results <- results
  output$cutoff <- pcutoff
  class(output) <- "dmrcate.output"
  output
}


########### LUAD part ###########
targets <- readRDS("../450k/TCGA/LUAD/Clinical_LUAD.rds")
targets$NC <- factor(targets$NC, ordered = T, levels = rev(levels(factor(targets$NC))))
targets$age_qtile <- "over_65"
targets$age_qtile[which(targets$age < 65)] <- "under_65"
targets$stage <- gsub("a", "", targets$stage)
targets$stage <- gsub("b", "", targets$stage)
targets <- subset(targets, stage == "ii" | stage == "i" | NC == "Control")

Beta <- readRDS("../450k/TCGA/LUAD/LUAD_norm.rds")
Beta <- Beta[,targets$Sample]
M <- log(Beta/(1-Beta))

donkeep <- which(apply(Beta, 2, function(x){length(which(x >= 0.999999 | x <= 0.000001))}) > 500)
Beta <- Beta[,-donkeep]
M <- M[,-donkeep]
targets <- targets[-donkeep,]


cellType <- factor(targets$NC)
individual <- factor(targets$age_qtile) 
design <- model.matrix(~0+cellType+individual, data=targets)
colnames(design) <- c(levels(cellType),levels(individual)[-1])
fit <- lmFit(M, design)
contMatrix <- makeContrasts(Case-Control,
                            levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, coef = "Case - Control", number = nrow(M))
betafit <- lmFit(ilogit2(M), design)
betafit <- contrasts.fit(betafit, contMatrix)
betafit <- eBayes(betafit)
betatt <- topTable(betafit, coef = "Case - Control", number = nrow(M))
m <- match(rownames(tt), rownames(betatt))
tt$betafc <- betatt$logFC[m]
m <- match(rownames(M), rownames(tt))
tt <- tt[m, ]
grset <- makeGenomicRatioSetFromMatrix(M, 
                                       array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                       mergeManifest = TRUE, what = "M")
anno <- getAnnotation(grset)
stat <- tt$t
annotated <- data.frame(ID = rownames(tt), stat = stat, 
                        CHR = anno$chr, pos = anno$pos, betafc = tt$betafc, 
                        indfdr = tt$adj.P.Val, is.sig = tt$adj.P.Val < 0.001)
N_index <- which(targets$NC == "Control")
T_index <- which(targets$NC == "Case")
meanbeta_N <- apply(Beta[,N_index], 1, mean)
meanbeta_T <- apply(Beta[,T_index], 1, mean)
SDbeta_N <- apply(Beta[,N_index], 1, sd)
SDbeta_T <- apply(Beta[,T_index], 1, sd)
meanbeta_diff <- meanbeta_T - meanbeta_N
percent_ToverN <- NA
for(i in 1:length(meanbeta_N)){
  percent_ToverN[i] <- length(which(Beta[i,T_index] >= (meanbeta_N[i] + 2*SDbeta_N[i])))/length(T_index)
}
percent_NoverT <- NA
for(i in 1:length(meanbeta_T)){
  percent_NoverT[i] <- length(which(Beta[i,N_index] >= (meanbeta_T[i] + 2*SDbeta_T[i])))/length(N_index)
}
annotated$meanbeta_diff <- meanbeta_diff
annotated$is.hyper <- (meanbeta_N < 0.25 & percent_ToverN >= 0.7 & meanbeta_diff >= 0.2)
# 8462
annotated$is.hypo <- (percent_NoverT >= 0.7 & meanbeta_diff <= -0.2)
# 472
myAnnotation <- subset(annotated, is.sig == TRUE & abs(meanbeta_diff)>=0.2)
ID_LUAD <- myAnnotation$ID
LUAD_hypo <- annotated$ID[which(annotated$is.hypo == T)]
LUAD_hyper <- annotated$ID[which(annotated$is.hyper == T)]
class(myAnnotation) <- "annot"


DMRs <- dmrcate_redefine(myAnnotation, lambda=250, C=2, min.cpgs = 2)
length(myAnnotation$meanbeta_diff)
# 22110
nrow(DMRs$results)
# 1193
sum(DMRs$results$no.cpg)
# 2529
sum(DMRs$results$no.hyper)
# 953
sum(DMRs$results$no.hypo)
# 50
saveRDS(DMRs, "../Results/updated/2013method_LUAD_DMR.rds")

ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(M),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)

annotated$is.hyper <- (meanbeta_diff >= 0.2)
# 8462
annotated$is.hypo <- (meanbeta_diff <= -0.2)
# 472
myAnnotation <- subset(annotated, is.sig == TRUE & abs(meanbeta_diff)>=0.2)
class(myAnnotation) <- "annot"
DMRs_test <- dmrcate_redefine(myAnnotation, lambda=250, C=2, min.cpgs = 2)
sum(DMRs_test$results$no.hyper)
# 1591
sum(DMRs_test$results$no.hypo)
# 938


########### LUSC part ###########

targets <- readRDS("../450k/TCGA/LUSC/Clinical_LUSC.rds")
targets$NC <- factor(targets$NC, ordered = T, levels = rev(levels(factor(targets$NC))))
targets$age_qtile <- "over_65"
targets$age_qtile[which(targets$age < 65)] <- "under_65"
targets$stage <- gsub("stage ", "", targets$stage)
targets$stage <- gsub("a", "", targets$stage)
targets$stage <- gsub("b", "", targets$stage)
targets <- subset(targets, stage == "ii" | stage == "i" | NC == "Control")

Beta <- readRDS("../450k/TCGA/LUSC/LUSC_norm.rds")
Beta <- Beta[,targets$Sample]
M <- log(Beta/(1-Beta))

donkeep <- which(apply(Beta, 2, function(x){length(which(x >= 0.999999 | x <= 0.000001))}) > 500)
if(length(donkeep) != 0){
  Beta <- Beta[,-donkeep]
  M <- M[,-donkeep]
  targets <- targets[-donkeep,]
}

cellType <- factor(targets$NC)
individual <- factor(targets$age_qtile) 
design <- model.matrix(~0+cellType+individual, data=targets)
colnames(design) <- c(levels(cellType),levels(individual)[-1])
fit <- lmFit(M, design)
contMatrix <- makeContrasts(Case-Control,
                            levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, coef = "Case - Control", number = nrow(M))
betafit <- lmFit(ilogit2(M), design)
betafit <- contrasts.fit(betafit, contMatrix)
betafit <- eBayes(betafit)
betatt <- topTable(betafit, coef = "Case - Control", number = nrow(M))
m <- match(rownames(tt), rownames(betatt))
tt$betafc <- betatt$logFC[m]
m <- match(rownames(M), rownames(tt))
tt <- tt[m, ]
grset <- makeGenomicRatioSetFromMatrix(M, 
                                       array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                       mergeManifest = TRUE, what = "M")
anno <- getAnnotation(grset)
stat <- tt$t
annotated <- data.frame(ID = rownames(tt), stat = stat, 
                        CHR = anno$chr, pos = anno$pos, betafc = tt$betafc, 
                        indfdr = tt$adj.P.Val, is.sig = tt$adj.P.Val < 0.001)
N_index <- which(targets$NC == "Control")
T_index <- which(targets$NC == "Case")
meanbeta_N <- apply(Beta[,N_index], 1, mean)
meanbeta_T <- apply(Beta[,T_index], 1, mean)
SDbeta_N <- apply(Beta[,N_index], 1, sd)
SDbeta_T <- apply(Beta[,T_index], 1, sd)
meanbeta_diff <- meanbeta_T - meanbeta_N
percent_ToverN <- NA
for(i in 1:length(meanbeta_N)){
  percent_ToverN[i] <- length(which(Beta[i,T_index] >= (meanbeta_N[i] + 2*SDbeta_N[i])))/length(T_index)
}
percent_NoverT <- NA
for(i in 1:length(meanbeta_T)){
  percent_NoverT[i] <- length(which(Beta[i,N_index] >= (meanbeta_T[i] + 2*SDbeta_T[i])))/length(N_index)
}
annotated$meanbeta_diff <- meanbeta_diff
annotated$is.hyper <- (meanbeta_N < 0.25 & percent_ToverN >= 0.7 & meanbeta_diff >= 0.2)
# 8839
annotated$is.hypo <- (percent_NoverT >= 0.7 & meanbeta_diff <= -0.2)
# 2148
myAnnotation <- subset(annotated, is.sig == TRUE & abs(meanbeta_diff)>=0.2)
ID_LUSC <- myAnnotation$ID
LUSC_hypo <- annotated$ID[which(annotated$is.hypo == T)]
LUSC_hyper <- annotated$ID[which(annotated$is.hyper == T)]
class(myAnnotation) <- "annot"

DMRs <- dmrcate_redefine(myAnnotation, lambda=250, C=2, min.cpgs = 2)
length(myAnnotation$meanbeta_diff)
# 41551
nrow(DMRs$results)
# 3657
sum(DMRs$results$no.cpg)
# 8123
sum(DMRs$results$no.hyper)
# 1728
sum(DMRs$results$no.hypo)
# 395
saveRDS(DMRs, "../Results/updated/2013method_LUSC_DMR.rds")

ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(M),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)

annotated$is.hyper <- (meanbeta_diff >= 0.2)
# 16100
annotated$is.hypo <- (meanbeta_diff <= -0.2)
# 25459
myAnnotation <- subset(annotated, is.sig == TRUE & abs(meanbeta_diff)>=0.2)
class(myAnnotation) <- "annot"
DMRs_test <- dmrcate_redefine(myAnnotation, lambda=250, C=2, min.cpgs = 2)
sum(DMRs_test$results$no.hyper)
# 3167
sum(DMRs_test$results$no.hypo)
# 4956

########### Region statistic part ###########
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
LUAD_sig_anno <- ann450k[match(as.character(ID_LUAD),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
LUSC_sig_anno <- ann450k[match(as.character(ID_LUSC),ann450k$Name),
                         c(1:4,12:19,24:ncol(ann450k))]
Overlap_sig_anno <- ann450k[match(intersect(as.character(ID_LUAD), as.character(ID_LUSC)), ann450k$Name),
                         c(1:4,12:19,24:ncol(ann450k))]
LUAD_hyper_anno <- ann450k[match(as.character(LUAD_hyper),ann450k$Name),
                         c(1:4,12:19,24:ncol(ann450k))]
LUSC_hyper_anno <- ann450k[match(as.character(LUSC_hyper),ann450k$Name),
                         c(1:4,12:19,24:ncol(ann450k))]
Overlap_hyper_anno <- ann450k[match(intersect(as.character(LUAD_hyper), as.character(LUSC_hyper)), ann450k$Name),
                              c(1:4,12:19,24:ncol(ann450k))]
LUAD_hypo_anno <- ann450k[match(as.character(LUAD_hypo),ann450k$Name),
                           c(1:4,12:19,24:ncol(ann450k))]
LUSC_hypo_anno <- ann450k[match(as.character(LUSC_hypo),ann450k$Name),
                           c(1:4,12:19,24:ncol(ann450k))]
Overlap_hypo_anno <- ann450k[match(intersect(as.character(LUAD_hypo), as.character(LUSC_hypo)), ann450k$Name),
                              c(1:4,12:19,24:ncol(ann450k))]

region_stat <- data.frame(Region = names(table(ann450k$Relation_to_Island)), stringsAsFactors = F)
region_stat$LUAD_sig <- table(LUAD_sig_anno$Relation_to_Island)[region_stat$Region]
region_stat$LUSC_sig <- table(LUSC_sig_anno$Relation_to_Island)[region_stat$Region]
region_stat$Overlap_sig <- table(Overlap_sig_anno$Relation_to_Island)[region_stat$Region]
region_stat$LUAD_hyper <- table(LUAD_hyper_anno$Relation_to_Island)[region_stat$Region]
region_stat$LUSC_hyper <- table(LUSC_hyper_anno$Relation_to_Island)[region_stat$Region]
region_stat$Overlap_hyper <- table(Overlap_hyper_anno$Relation_to_Island)[region_stat$Region]
region_stat$LUAD_hypo <- table(LUAD_hypo_anno$Relation_to_Island)[region_stat$Region]
region_stat$LUSC_hypo <- table(LUSC_hypo_anno$Relation_to_Island)[region_stat$Region]
region_stat$Overlap_hypo <- table(Overlap_hypo_anno$Relation_to_Island)[region_stat$Region]

region_stat_percent <- region_stat
for(i in 2:ncol(region_stat)){
  region_stat_percent[,i] <- round(region_stat_percent[,i]/sum(region_stat_percent[,i])*100, 2)
}
total <- region_stat[1,]
total$Region <- "total"
total[,2:ncol(total)] <- apply(region_stat[,2:ncol(total)], 2, sum)
region_stat <- rbind(region_stat, total)
region_stat_percent <- rbind(region_stat_percent, total)
write.table(region_stat, "../Results/updated/region_stat.txt", sep = "\t", quote = F, row.names = F)
write.table(region_stat_percent, "../Results/updated/region_stat_percent.txt", sep = "\t", quote = F, row.names = F)


########### DMRs statistic part ###########
LUAD_DMRs <- readRDS("../Results/updated/2013method_LUAD_DMR.rds")
LUSC_DMRs <- readRDS("../Results/updated/2013method_LUSC_DMR.rds")
length(intersect(LUAD_DMRs$input$ID, LUSC_DMRs$input$ID))
# 15069
length(intersect(LUAD_DMRs$input$ID[which(LUAD_DMRs$input$is.hyper)], LUSC_DMRs$input$ID[which(LUSC_DMRs$input$is.hyper)]))
# 4954
length(intersect(LUAD_DMRs$input$ID[which(LUAD_DMRs$input$is.hypo)], LUSC_DMRs$input$ID[which(LUSC_DMRs$input$is.hypo)]))
# 223
length(intersect(LUAD_DMRs$results$coord, LUSC_DMRs$results$coord))
# 9

results.ranges_LUAD <- extractRanges(LUAD_DMRs, genome = "hg19")
gene_LUAD <- unique(unlist(strsplit(results.ranges_LUAD$overlapping.promoters, split = ",")))
for(i in 1:length(gene_LUAD)){
  gene_LUAD[i] <- unlist(strsplit(gene_LUAD[i], split = "-"))[1]
}
gene_LUAD <- unique(gsub(" ", "", gene_LUAD))
gene_LUAD <- gene_LUAD[-is.na(gene_LUAD)]
length(gene_LUAD)

results.ranges_LUSC <- extractRanges(LUSC_DMRs, genome = "hg19")
gene_LUSC <- unique(unlist(strsplit(results.ranges_LUSC$overlapping.promoters, split = ",")))
for(i in 1:length(gene_LUSC)){
  gene_LUSC[i] <- unlist(strsplit(gene_LUSC[i], split = "-"))[1]
}
gene_LUSC <- unique(gsub(" ", "", gene_LUSC))
gene_LUSC <- gene_LUSC[-is.na(gene_LUSC)]
length(gene_LUSC)

########### get DMRs bed ###########
starts <- as.integer(unlist(strsplit(unlist(strsplit(LUAD_DMRs$results$coord, split = ":"))[seq(2,(2*nrow(LUAD_DMRs$results)),2)], 
                                     split = "-"))[seq(1,(nrow(LUAD_DMRs$results)*2-1),2)])
ends <- as.integer(unlist(strsplit(unlist(strsplit(LUAD_DMRs$results$coord, split = ":"))[seq(2,(2*nrow(LUAD_DMRs$results)),2)], 
                                   split = "-"))[seq(2,(nrow(LUAD_DMRs$results)*2),2)])
chr <- unlist(strsplit(LUAD_DMRs$results$coord, split = ":"))[seq(1,(2*nrow(LUAD_DMRs$results)-1),2)]

sum(ends-starts+1)
bed_LUAD <- data.frame(chr = chr, start = starts, end = (ends+1))
write.table(bed_LUAD, "../Results/updated/bed_DMR_2013method_LUAD.txt", sep = "\t", quote = F, row.names = F, col.names = F)

starts <- as.integer(unlist(strsplit(unlist(strsplit(LUSC_DMRs$results$coord, split = ":"))[seq(2,(2*nrow(LUSC_DMRs$results)),2)], 
                                     split = "-"))[seq(1,(nrow(LUSC_DMRs$results)*2-1),2)])
ends <- as.integer(unlist(strsplit(unlist(strsplit(LUSC_DMRs$results$coord, split = ":"))[seq(2,(2*nrow(LUSC_DMRs$results)),2)], 
                                   split = "-"))[seq(2,(nrow(LUSC_DMRs$results)*2),2)])
chr <- unlist(strsplit(LUSC_DMRs$results$coord, split = ":"))[seq(1,(2*nrow(LUSC_DMRs$results)-1),2)]

sum(ends-starts+1)
bed_LUSC <- data.frame(chr = chr, start = starts, end = (ends+1))
write.table(bed_LUSC, "../Results/updated/bed_DMR_2013method_LUSC.txt", sep = "\t", quote = F, row.names = F, col.names = F)


