library(glmnet)
library(ROCR)

# read Data
LUAD <- readRDS("Results/M_LUAD_top1000.rds")
targets_LUAD <- readRDS("Data/Clinical_LUAD.rds")
targets_LUAD$NC <- factor(targets_LUAD$NC, ordered = T, levels = rev(levels(factor(targets_LUAD$NC))))
targets_LUAD$age_qtile <- "over_65"
targets_LUAD$age_qtile[which(targets_LUAD$age < 65)] <- "under_65"
targets_LUAD$stage <- gsub("a", "", targets_LUAD$stage)
targets_LUAD$stage <- gsub("b", "", targets_LUAD$stage)
targets_LUAD <- subset(targets_LUAD, stage == "ii" | stage == "i" | NC == "Control")
targets_LUAD <- subset(targets_LUAD, Sample %in% colnames(LUAD))
LUAD <- LUAD[, targets_LUAD$Sample]

## set and save training testing dataset
LUAD <- t(LUAD)
training_LUAD <- sample(c(1:nrow(LUAD)), round(nrow(LUAD)*2/3, 0))
LUAD_training <- LUAD[training_LUAD,]
LUAD_testing <- LUAD[-training_LUAD,]
targets_LUAD_training <- targets_LUAD[training_LUAD,]
targets_LUAD_testing <- targets_LUAD[-training_LUAD,]
saveRDS(LUAD_training, "Data/LUAD_training.rds")
saveRDS(LUAD_testing, "Data/LUAD_testing.rds")
saveRDS(targets_LUAD_training, "Data/targets_LUAD_training.rds")
saveRDS(targets_LUAD_testing, "Data/targets_LUAD_testing.rds")

# if you already have grouped data files
LUAD_training <- readRDS("Data/LUAD_training.rds")
LUAD_testing <- readRDS("Data/LUAD_testing.rds")
targets_LUAD_training <- readRDS("Data/targets_LUAD_training.rds")
targets_LUAD_testing <- readRDS("Data/targets_LUAD_testing.rds")

## LASSO
for (i in 1:500){
  model.lasso <- cv.glmnet(
    LUAD_training, targets_LUAD_training$NC,
    family="binomial",
    alpha=1,
    nfolds=10
  )
  # plot(model.lasso)
  idx <- with(model.lasso, which.min(abs(glmnet.fit$lambda - lambda.1se)))
  coefs.lasso <- with(model.lasso$glmnet.fit, c(
    `(Intercept)`=unname(a0[idx]),
    beta[,idx]
  ))
  coefs.lasso <- coefs.lasso[abs(coefs.lasso) > .Machine$double.eps]
  if(i == 1){
    candidates <- names(coefs.lasso)
  }else{
    candidates <- c(candidates, names(coefs.lasso))
  }
  if (i %% 50 == 0){
    cat(paste(i, "predictions accomplished!\n", sep = " "))
  } 
}

# pick CpGs that was picked over 450 times in 500 iterations
predictors <- unique(candidates)
predictors <- predictors[which(predictors != "(Intercept)")]
cat(paste("there are", length(predictors), "predictors\n"))
LUAD_freq <- table(candidates)
predictors <- names(LUAD_freq)[which(LUAD_freq >= 450)[-1]]
saveRDS(LUAD_freq, "Results/lasso_LUAD_DMP.rds")
cpg <- names(LUAD_freq)[which(LUAD_freq >= 450)[-1]]
cat(paste(length(cpg), "predictors are picked\n"))

# use training data for model fitting and testing data for validation
glm_mat <- LUAD_training[,cpg]
glm_data <- as.data.frame(glm_mat)
colnames(glm_data) <- cpg
glm_data$group <- factor(targets_LUAD_training$NC)
f <- paste("group ~ ", cpg[1], sep = "")
if(length(cpg) != 1){
  for (j in 2:length(cpg)){
    f <- paste(f, " + ", cpg[j], sep = "")
  }
}
glm.model <- glm(f, family = "binomial", data = glm_data)
glm_mat <- LUAD_testing[,cpg]
glm_data <- as.data.frame(glm_mat)
colnames(glm_data) <- cpg
prediction <- predict(glm.model, glm_data, type="response")
perf <- performance(prediction(prediction, targets_LUAD_testing$NC), "tpr", "fpr")
# ROC plot
plot(perf)
# AUC
print(performance(prediction(prediction, targets_LUAD_testing$NC), "auc")@y.values)

