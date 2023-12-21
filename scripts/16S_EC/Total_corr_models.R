#!/usr/bin/env Rscript

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

ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_Total_corr.RDS")

for (myVar in myVars) {
  data <- ml_EC_16S %>% select(grep("EC:", colnames(ml_EC_16S)))
  df.myVar <- as.vector(ml_EC_16S[,myVar])
  df.climate <- as.vector(ml_EC_16S[,'ClimateZ'])
  df.clay <- as.vector(ml_EC_16S[,'clay'])
  df.dna <- as.vector(ml_EC_16S[,'DNA'])
  data <- cbind(df.myVar, df.climate, df.clay, df.dna, data)
  names(data)[1] <- myVar
  names(data)[2] <- 'ClimateZ'
  names(data)[3] <- 'clay'
  names(data)[4] <- 'DNA'
  
  # format so : and . are replaced by _ (for varimp)
  names(data) <- gsub(":","_", names(data))
  names(data) <- gsub("\\.","_", names(data))
  
  myFormula <- as.formula(paste0(myVar, ' ~ .'))
  Boruta.res <- Boruta(myFormula, data=data)
  myFormula <- getNonRejectedFormula(Boruta.res)
  keep_X <- names(Boruta.res$finalDecision[Boruta.res$finalDecision != "Rejected"])
  keep_X
  
  final <- data %>%
    select(all_of(keep_X), myVar)
  
  if ('ClimateZ' %in% keep_X) {
    final$ClimateZ <- factor(final$ClimateZ)
  }
  if ('clay' %in% keep_X) {
    final$clay <- as.numeric(final$clay)
  }
  if ('DNA' %in% keep_X) {
    final$clay <- as.numeric(final$DNA)
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
              file = paste0("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_total_corr/", myVar, "_EC_total_corr_stats", ".csv", sep = ""),
              col.names = TRUE, append = TRUE, sep = ",", row.names = FALSE)
  
  # variable importances
  imp <- permimp::permimp(cf.train, nperm=1, OOB=TRUE, scaled=FALSE,
                          conditional=FALSE, asParty=FALSE,
                          thresholdDiagnostics = FALSE, progressBar = TRUE)
  
  imp <- as.data.frame(imp$values)
  names(imp) <- 'vimp'
  imp <- imp %>% arrange(desc(vimp))
  imp <- tibble::rownames_to_column(imp, "EC")
  imp$Run <- args[1]
  
  #write variable importances to a file
  write.table(imp,
              file = paste0("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_total_corr/", myVar, "_total_corr_varimp", ".csv", sep = ""),
              col.names = FALSE, append = TRUE, sep = ",", row.names = FALSE)
  
  for (EC in seq(1,min(10,length(imp$EC)),1)) {
    print(paste0("Partial dependence for predictor ", EC, ": ", imp$EC[EC]))
    
    pd <- GetPartialData(cf.train, xnames=imp$EC[EC], 
                         quantiles=FALSE, grid.resolution = NULL,
                         parallel=TRUE)
    pd$run <- args[1]
    
    write.table(pd, 
                file = paste0("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_total_corr/", myVar, "_total_corr_pd", "_Predictor_", imp$EC[EC], ".csv", sep = ""),
                col.names=FALSE,
                append = TRUE, sep = ",", row.names = FALSE)
  }
  
}
