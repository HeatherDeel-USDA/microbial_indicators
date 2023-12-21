### Final models

# libraries
library(dplyr)
library(party)
library(tidyverse)
library(Boruta)
library(moreparty)
library(permimp)

myVars <- c('ace', 'SOM', 'activeC', 'resp', 'agg_stab', 'water_cap', 'ph', 'p', 'k', 'mg', 'fe', 'mn', 'zn')

ml_EC_16S <- readRDS("machine_learning/EC_RA_corr/ml_EC_RA_corr.RDS")

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
  
  ########### before rerunning final models, take code from figures script for compiling the final enzymes and write description of making the final_varimps file
  #######
  # final_varimps.csv was made by manually compiling X (write better description after redone)
  varimp_df <- read.csv("machine_learning/16S_EC/final_varimps.csv", header = FALSE)
  
  colnames(varimp_df) <- c("Predictor","VIP","Run","Category","Class","Sub_class",
                           "Name","Short","Book","Missing","VIP_scaled","all_gibbs",
                           "indicator")
  
  myVars <- c('ace', 'SOM', 'activeC', 'resp', 'agg_stab', 'water_cap', 'ph', 'p', 'k', 'mg', 'fe', 'mn', 'zn')
  
  for (myVar in myVars) {
    varimp_df_filt <- filter(varimp_df, indicator == myVar)
    varimp_df_list <- varimp_df_filt %>% 
      group_by(Predictor) %>% 
      summarize(VIP_scaled_avg = mean(VIP_scaled))
    write.csv(varimp_df_list, file = paste0("machine_learning/16S_FINAL/", myVar,"_EC_GIBBs_final.csv"),
              row.names = FALSE)
  }
  ########
  
  # final enzymes
  final_list <- read.csv(paste0("machine_learning/16S_FINAL/",myVar,"_EC_GIBBs_final.csv"))
  keep_X <- final_list$Predictor
  myFormula <- as.formula(paste0(myVar, "~ ."))
  
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
  
  # not splitting into train/test so we can get pdps on all data
  # will report importances of the final model
  final$id <- 1:nrow(final)
  
  # get rid of id column
  final <- final %>% select(-id)

  # cforest on training data
  my_cforest_control <- cforest_control(teststat = "quad",
                                        testtype = "Univ", mincriterion = 0, ntree = 500, 
                                        mtry = floor(length(keep_X)/3),
                                        replace = FALSE)
  cf.final <- cforest(myFormula, data = final,
                      controls = my_cforest_control)
  
  # variable importances
  imp <- permimp::permimp(cf.final, nperm=1, OOB=TRUE, scaled=FALSE,
                          conditional=FALSE, asParty=FALSE,
                          thresholdDiagnostics = FALSE, progressBar = TRUE)
  
  imp <- as.data.frame(imp$values)
  names(imp) <- 'vimp'
  imp <- imp %>% arrange(desc(vimp))
  imp <- tibble::rownames_to_column(imp, "EC")
  
  # write variable importances to a file
  # when running these, check imp and make sure they're writing correctly
  write.table(imp,
              file = paste0("machine_learning/16S_FINAL/", myVar, "_model_results/", myVar, "_final_varimp", ".csv", sep = ""),
              col.names = FALSE, append = TRUE, sep = ",", row.names = FALSE)
  for (EC in seq(1,min(10,length(imp$EC)),1)) {
    print(paste0("Partial dependence for predictor ", EC, ": ", imp$EC[EC]))
    
    pd <- GetPartialData(cf.final, xnames=imp$EC[EC], 
                         quantiles=FALSE, grid.resolution = NULL,
                         parallel=TRUE)
    
    write.table(pd, 
                file = paste0("machine_learning/16S_FINAL/", myVar, "_model_results/", myVar, "_final_pd", "_Predictor_", imp$EC[EC], ".csv", sep = ""),
                col.names=FALSE,
                append = TRUE, sep = ",", row.names = FALSE)
  }
}