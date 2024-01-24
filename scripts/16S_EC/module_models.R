#!/usr/bin/env Rscript

### Module models

# accept command line arguments and save them in a list called args
args = commandArgs(trailingOnly=TRUE)

# print task number
print(paste0('Hello! I am task number: ', args[1]))

# libraries
library(dplyr)
library(party)
library(tidyverse)
library(Boruta)
library(moreparty)
library(permimp)

myVars <- c('ace', 'SOM', 'activeC', 'resp', 'agg_stab', 'water_cap', 'ph', 'p', 'k', 'mg', 'fe', 'mn', 'zn')

################################################################################
#################### This chunk run once on local computer ####################
# # import ECs that will be assigned to network modules
# ECs <- read.csv("machine_learning/EC_RA_corr/ml_EC_RA_corr.csv", check.names = FALSE)
# ECs2 <- ECs[,4:2439]
# ECs2_t <- t(ECs2)
# ECs2_t <- data.frame(ECs2_t)
# 
# # make column names the sample names
# sample_names <- ECs$SampleID
# colnames(ECs2_t) <- sample_names
# 
# # make ECs a column for joining
# ECs2_t <- rownames_to_column(ECs2_t, var = "Predictor")
# 
# # import modules, join modules and ECs
# modules <- read.csv("scripts/16S_EC/res_node_table_0.8.csv")
# colnames(modules)[1] <- "Predictor"
# ECs_mod <- ECs2_t %>%
#   left_join(modules, by = join_by(Predictor))
# ECs_mod <- ECs_mod[,c(2:537,544)] # reduce to just enzymes and module number
# ECs_mod <- ECs_mod %>%
#   filter(module != "NA")
# ECs_mod <- ECs_mod %>%
#   group_by(module) %>%
#   summarise_all(., .funs = sum)
# 
# # format for joining with indicators
# ECs_mod <- column_to_rownames(ECs_mod, var = "module")
# ECs_mod_t <- t(ECs_mod)
# ECs_mod_t <- data.frame(ECs_mod_t)
# ECs_mod_t <- rownames_to_column(ECs_mod_t, var = "SampleID")
# 
# # add back indicator info using SampleID
# library(phyloseq)
# ECs_join <- sample_data(ECs[,c('SampleID','ace', 'SOM', 'activeC', 'resp', 'agg_stab', 'water_cap', 'ph', 'p', 'k', 'mg', 'fe', 'mn', 'zn','ClimateZ','clay','DNA')])
# ECs_mod_final <- ECs_mod_t %>%
#   left_join(ECs_join, by = "SampleID")
# 
# # convert M to MOD in module names for ease of column selection
# colnames(ECs_mod_final) <- gsub("^M","MOD", colnames(ECs_mod_final))
# # fix the SOM column
# colnames(ECs_mod_final)[which(names(ECs_mod_final) == "SOMOD")] <- "SOM"
#
#saveRDS(ECs_mod_final, 'machine_learning/16S_EC/ECs_mod_final.RDS')
################################################################################

ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ECs_mod_final.RDS")

for (myVar in myVars) {
  data <- ml_EC_16S %>% select(grep("MOD", colnames(ml_EC_16S)))
  df.myVar <- as.vector(ml_EC_16S[,myVar])
  df.climate <- as.vector(ml_EC_16S[,'ClimateZ'])
  df.clay <- as.vector(ml_EC_16S[,'clay'])
  df.dna <- as.vector(ml_EC_16S[,'DNA'])
  data <- cbind(df.myVar, df.climate, df.clay, df.dna, data)
  names(data)[1] <- myVar
  names(data)[2] <- 'ClimateZ'
  names(data)[3] <- 'clay'
  names(data)[4] <- 'DNA'
  
  myFormula <- as.formula(paste0(myVar, ' ~ .'))
  Boruta.res <- Boruta(myFormula, data=data)
  myFormula <- getNonRejectedFormula(Boruta.res)
  keep_X <- names(Boruta.res$finalDecision[Boruta.res$finalDecision != "Rejected"])
  keep_X
  
  print(paste0('Hello! I am task number ', args[1], " and I made it through Boruta"))
  
  final <- data %>%
    select(all_of(keep_X), myVar)
  
  if ('ClimateZ' %in% keep_X) {
    final$ClimateZ <- factor(final$ClimateZ)
  }
  if ('clay' %in% keep_X) {
    final$clay <- as.numeric(final$clay)
  }
  if ('DNA' %in% keep_X) {
    final$DNA <- as.numeric(final$DNA)
  }
  
  print(paste0("Starting run for indicator ", myVar))
  
  # split into train and test (4/5 proportion)
  final$id <- 1:nrow(final)
  train <- final %>% dplyr::sample_frac(0.80)
  test <- dplyr::anti_join(final, train, by = 'id')
  
  # get rid of id columns
  train <- train %>% select(-id)
  test <- test %>% select(-id)
  
  # cforest on training data
  my_cforest_control <- cforest_control(teststat = "quad",
                                        testtype = "Univ", mincriterion = 0, ntree = 500, 
                                        mtry = floor(length(keep_X)/3),
                                        replace = FALSE)
  cf.train <- cforest(myFormula, data = train,
                      controls = my_cforest_control)
  
  # predict the response
  cf.pred <- predict(cf.train, newdata = test, OOB = TRUE, type = "response")
  
  # observed vs predicted formatting
  colnames(cf.pred)[1] <- "pred"
  cf.pred <- data.frame(cf.pred)
  cf.pred <- rownames_to_column(cf.pred, var = "id")
  obs <- data.frame(test[,myVar])
  colnames(obs)[1] <- "obs"
  cf.pvso <- cbind(cf.pred, obs)
  more_stats <- moreparty::PerfsRegression(cf.pvso$pred, cf.pvso$obs)
  more_stats <- t(data.frame(more_stats))
  
  # observed vs predicted linear model
  res.lm <- lm(obs ~ pred, data = cf.pvso)
  
  # save R^2 and p-values to a file
  r2_val <- summary(res.lm)$adj.r.squared
  p_val <- summary(res.lm)$coef[2,4]
  write.table(cbind(args[1], r2_val, p_val, myVar, more_stats),
              file = paste0("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_modules/", myVar, "_modules_corr_stats", ".csv", sep = ""),
              col.names = FALSE, append = TRUE, sep = ",", row.names = FALSE)
  
  # variable importances
  imp <- permimp::permimp(cf.train, nperm=1, OOB=TRUE, scaled=FALSE,
                          conditional=FALSE, asParty=FALSE,
                          thresholdDiagnostics = FALSE, progressBar = TRUE)
  
  imp <- as.data.frame(imp$values)
  names(imp) <- 'vimp'
  imp <- imp %>% arrange(desc(vimp))
  imp <- tibble::rownames_to_column(imp, "Module")
  imp$Run <- args[1]
  
  # write variable importances to a file
  # writing each to their own file because scinet was messing up 
  write.table(imp,
              file = paste0("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_modules/", myVar, "_modules_varimp", args[1],".csv", sep = ""),
              col.names = FALSE, append = TRUE, sep = ",", row.names = FALSE)
  
  for (Module in seq(1,min(10,length(imp$Module)),1)) {
    print(paste0("Partial dependence for predictor ", Module, ": ", imp$Module[Module]))
    
    pd <- GetPartialData(cf.train, xnames=imp$Module[Module], 
                         quantiles=FALSE, grid.resolution = NULL,
                         parallel=TRUE)
    pd$run <- args[1]
    
    write.table(pd, 
                file = paste0("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_modules/", myVar, "_modules_pd", "_Predictor_", imp$Module[Module], ".csv", sep = ""),
                col.names=FALSE,
                append = TRUE, sep = ",", row.names = FALSE)
  }
}