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

Beta_o <- Beta
M_o <- M
targets_o <- targets

for(j in 1:100){
  cat(paste(j, "\n", sep = ""))
  picked <- sample(c(1:nrow(targets_o)), round(nrow(targets_o)*2/3, 0))
  Beta <- Beta_o[,picked]
  M <- M_o[,picked]
  targets <- targets_o[picked,]
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
  m <- match(rownames(M), rownames(tt))
  tt <- tt[m, ]
  stat <- tt$t
  annotated <- data.frame(ID = rownames(tt), stat = stat,
                          indfdr = tt$adj.P.Val, is.sig = tt$adj.P.Val < 0.001)
  annotated <- subset(annotated, is.sig == TRUE)
  N_index <- which(targets$NC == "Control")
  T_index <- which(targets$NC == "Case")
  meanbeta_N <- apply(Beta[annotated$ID,N_index], 1, mean)
  meanbeta_T <- apply(Beta[annotated$ID,T_index], 1, mean)
  SDbeta_T <- apply(Beta[,T_index], 1, sd)
  meanbeta_diff <- meanbeta_T - meanbeta_N
  percent_NoverT <- NA
  for(i in 1:length(meanbeta_T)){
    percent_NoverT[i] <- length(which(Beta[i,N_index] >= (meanbeta_T[i] + 2*SDbeta_T[i])))/length(N_index)
  }
  annotated$meanbeta_diff <- meanbeta_diff
  annotated$is.hypo <- (percent_NoverT >= 0.6 & meanbeta_diff <= -0.2)
  # 472
  myAnnotation <- subset(annotated, is.hypo == TRUE)
  if(j == 1){
    candidates <- data.frame(ID = myAnnotation$ID, round = j, stringsAsFactors = F)
  }else{
    temp <- data.frame(ID = myAnnotation$ID, round = j, stringsAsFactors = F)
    candidates <- rbind(candidates, temp)
  }
}

saveRDS(candidates, "../Results/updated/cv_2013method_LUAD.rds")

#######################################################

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

Beta_o <- Beta
M_o <- M
targets_o <- targets

for(j in 1:100){
  cat(paste(j, "\n", sep = ""))
  picked <- sample(c(1:nrow(targets_o)), round(nrow(targets_o)*2/3, 0))
  Beta <- Beta_o[,picked]
  M <- M_o[,picked]
  targets <- targets_o[picked,]
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
  m <- match(rownames(M), rownames(tt))
  tt <- tt[m, ]
  stat <- tt$t
  annotated <- data.frame(ID = rownames(tt), stat = stat,
                          indfdr = tt$adj.P.Val, is.sig = tt$adj.P.Val < 0.001)
  annotated <- subset(annotated, is.sig == TRUE)
  N_index <- which(targets$NC == "Control")
  T_index <- which(targets$NC == "Case")
  meanbeta_N <- apply(Beta[annotated$ID,N_index], 1, mean)
  meanbeta_T <- apply(Beta[annotated$ID,T_index], 1, mean)
  SDbeta_T <- apply(Beta[,T_index], 1, sd)
  meanbeta_diff <- meanbeta_T - meanbeta_N
  percent_NoverT <- NA
  for(i in 1:length(meanbeta_T)){
    percent_NoverT[i] <- length(which(Beta[i,N_index] >= (meanbeta_T[i] + 2*SDbeta_T[i])))/length(N_index)
  }
  annotated$meanbeta_diff <- meanbeta_diff
  annotated$is.hypo <- (percent_NoverT >= 0.6 & meanbeta_diff <= -0.2)
  # 472
  myAnnotation <- subset(annotated, is.hypo == TRUE)
  if(j == 1){
    candidates <- data.frame(ID = myAnnotation$ID, round = j, stringsAsFactors = F)
  }else{
    temp <- data.frame(ID = myAnnotation$ID, round = j, stringsAsFactors = F)
    candidates <- rbind(candidates, temp)
  }
}

saveRDS(candidates, "../Results/updated/cv_2013method_LUSC.rds")



LUAD_DMRs <- readRDS("../Results/updated/2013method_LUAD_DMR.rds")
candidates <- readRDS("../Results/updated/cv_2013method_LUAD.rds")
freq <- table(as.character(candidates$ID))
cv <- rownames(freq)[order(c(freq), decreasing = T)][1:548]
table(LUAD_DMRs$input$ID[which(LUAD_DMRs$input$is.hypo)] %in% cv)
length(which(freq[which(names(freq) %in% LUAD_DMRs$input$ID[which(LUAD_DMRs$input$is.hypo)])]>50))


candidates <- readRDS("../Results/updated/cv_2013method_LUSC.rds")
LUSC_DMRs <- readRDS("../Results/updated/2013method_LUSC_DMR.rds")
freq <- table(as.character(candidates$ID))
#cv <- rownames(freq)[order(c(freq), decreasing = T)][1:2058]
cv <- rownames(freq)[order(c(freq), decreasing = T)][1:2386]
table(LUSC_DMRs$input$ID[which(LUSC_DMRs$input$is.hypo)] %in% cv)
length(which(freq[which(names(freq) %in% LUSC_DMRs$input$ID[which(LUSC_DMRs$input$is.hypo)])]>50))

##################################
LUAD <- readRDS("../Results/updated/M_LUAD_top1000.rds")
LUSC <- readRDS("../Results/updated/M_LUSC_top1000.rds")
table(rownames(LUAD) %in% LUAD_DMRs$input$ID)
table(rownames(LUSC) %in% LUSC_DMRs$input$ID)
