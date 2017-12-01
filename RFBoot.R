library(varSelRF)

# Read the data
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

# set and save training testing dataset
LUAD <- t(LUAD)
LUSC <- t(LUSC)
training_LUAD <- sample(c(1:nrow(LUAD)), round(nrow(LUAD)*2/3, 0))
LUAD_training <- LUAD[training_LUAD,]
LUAD_testing <- LUAD[-training_LUAD,]
targets_LUAD_training <- targets_LUAD[training_LUAD,]
targets_LUAD_testing <- targets_LUAD[-training_LUAD,]

saveRDS(LUAD_training, "Data/LUAD_training.rds")
saveRDS(LUAD_testing, "Data/LUAD_testing.rds")
saveRDS(targets_LUAD_training, "Data/targets_LUAD_training.rds")
saveRDS(targets_LUAD_testing, "Data/targets_LUAD_testing.rds")


# only training data are used
# doing RF as reference for bootstrap RF
# ntree is the number of trees built
# c.sd can be set to 0 or 1, if = 0, then in each tree, model with minumun OOB will be chosen; 
## if = 1, model with OOB with in min+-1sd and less variables will be chosen
# mtryFactor decides how may genes will be used at each split
# vars.drop.frac decides the fraction of variables dropped at each time
# By default, the least number of variables to be kept in each tree is 2.
# Example: 1000 variables, 5000 trees, vars.drop.frac = 0.3, c.sd = 1. 
## n st tree: 1000 > 667 > 467 > 327 > ...... > 4 > 3 > 2
## if OOB is: 12% > 11% > 10% > 10% > ...... > 10% > 8% > 7%,  sd = *% -> keep 2 variables
## if OOB is: 12% > 11% > 10% > 10% > ...... > 10% > 8% > 12%, sd = 1% -> keep 3 variables
## if OOB is: 12% > 11% > 10% > 10% > ...... > 10% > 8% > 12%, sd = 4% -> keep 2 variables
# Example: 1000 variables, 5000 trees, vars.drop.frac = 0.5, c.sd = 0.
## n st tree: 1000 > 500 > 250 > 125 > ...... > 8 > 4 > 2
## if OOB is: 12% > 11% > 10% > 10% > ...... > 10% > 8% > 7%,  sd = *% -> keep 2 variables
## if OOB is: 12% > 11% > 10% > 10% > ...... > 10% > 8% > 12%, sd = 1% -> keep 3 variables
## if OOB is: 12% > 11% > 10% > 10% > ...... > 10% > 8% > 12%, sd = 4% -> keep 3 variables

# varSelRF() is bagging
LUAD_RF <- varSelRF(LUAD_training, factor(targets_LUAD_training$NC, ordered = F), c.sd = 1, mtryFactor = 1, ntree = 5000,
                    ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2,
                    whole.range = FALSE, recompute.var.imp = FALSE, verbose = FALSE,
                    returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE)

# varSelRFBoot() is bootstrap. This should be run using Cluster, or it would be really slow
# it not then just comment the makeForkCluster and delete the cluster option in varSelRFBoot
core <- detectCores()
forkCL <- makeForkCluster(round(0.85*core, 0))
LUAD_RFBoot <- varSelRFBoot(LUAD_training, factor(targets_LUAD_training$NC, ordered = F), c.sd = 1,
                            mtryFactor = 1, ntree = 5000, ntreeIterat = 2000,
                            vars.drop.frac = 0.3, bootnumber = 1000,
                            whole.range = FALSE,
                            recompute.var.imp = FALSE, usingCluster = TRUE,
                            srf = LUAD_RF, TheCluster = forkCL)

saveRDS(LUAD_RFBoot, "Data/LUAD_RFBoot.rds")


# extract results
LUAD_RFBoot <- readRDS("Data/LUAD_RFBoot.rds")
# the variables chosen by RFBoot
cpg_LUAD <- LUAD_RFBoot$all.data.vars
# the prediction compared to observation
LUAD_rf <- LUAD_RFBoot$all.data.randomForest
table(LUAD_rf$predicted, LUAD_rf$y)

