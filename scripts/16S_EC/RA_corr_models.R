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

ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_RA_corr.RDS")

for (myVar in myVars) {
    data <- ml_EC_16S %>%
      select(grep("EC:", colnames(ml_EC_16S)), ClimateZ, clay, myVar)
    
    # format so : and . are replaced by _ (for varimp)
    names(data) <- gsub(":","_", names(data))
    names(data) <- gsub("\\.","_", names(data))
    
    myFormula <- as.formula(paste0(myVar, ' ~ .'))
    Boruta.res <- Boruta(myFormula, data=data)
    myFormula <- getConfirmedFormula(Boruta.res)
    keep_X <- names(Boruta.res$finalDecision[Boruta.res$finalDecision == "Confirmed"])
    keep_X
    
    final <- data %>%
      select(all_of(keep_X), myVar)
    
    if ('ClimateZ' %in% keep_X) {
      final$ClimateZ <- factor(final$ClimateZ)
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
      moreparty::PerfsRegression(cf.pvso$pred, cf.pvso$obs)
      saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_RA_corr/", myVar, "_EC_RA_corr_obs_vs_pred", args[1], ".RDS", sep = ""))
      
      # observed vs predicted linear model
      res.lm <- lm(obs ~ pred, data = cf.pvso)
      
      # save R^2 and p-values to a file
      r2_val <- summary(res.lm)$adj.r.squared
      p_val <- summary(res.lm)$coef[2,4]
      write.table(cbind(args[1], r2_val, p_val, myVar), 
                  file = paste0("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_RA_corr/", myVar, "_EC_RA_corr_stats", args[1], ".csv", sep = ""), 
                  col.names = FALSE, append = TRUE, sep = ",", row.names = FALSE)
      
      # variable importances
      # threshold = 1 identifies important features independent of the other predictors
      imp <- permimp::permimp(cf.train, nperm=1, OOB=TRUE, scaled=FALSE,
                             conditional=FALSE, threshold=1, asParty=FALSE,
                             thresholdDiagnostics = FALSE, progressBar = TRUE)

      imp <- as.data.frame(imp$values)
      names(imp) <- 'vimp'
      imp <- imp %>% arrange(desc(vimp))
      imp <- tibble::rownames_to_column(imp, "EC")
      imp$Run <- args[1]
      
      # write variable importances to a file
      write.table(imp,
                  file = paste0("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_RA_corr/", myVar, "_varimp", args[1], ".csv", sep = ""),
                  col.names = FALSE, append = TRUE, sep = ",", row.names = FALSE)
      
      # partial dependence plots
      pd <- GetPartialData(cf.train, xnames=imp$EC[2], quantiles=FALSE, grid.resolution = 10)
      ggForestEffects(pd, vline=0, xlabel="", ylabel="", main="")

      pd <- edarf::partial_dependence(cf.train, imp$EC[2], interaction=FALSE)
      p <- edarf::plot_pd(pd)
      ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/", myVar, "_model_results_RA_corr/", myVar, "_pd", args[1], ".png", sep = ""),
           unit = "in", width = 8, height = 8, dpi = 300, device = "png")
      
}
