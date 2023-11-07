#!/usr/bin/env Rscript

# accept command line arguments and save them in a list called args
args = commandArgs(trailingOnly=TRUE)

# print task number
print(paste0('Hello! I am task number: ', args[1]))

### things to do:
# make sure the varimp output works (before doing anything else)
# need to spit varimps out to a file and save them
# spit R2 and p-values to a file, but make it so the task number is attached to it
# graphing of predicted vs observed of R2 values for everything, save to a file such that the task number is attached to it
# put together scripts for job array in scinet
# add cores = 12 (or something) in cforest controls

### libraries
library(dplyr)
library(party)

ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_16S.RDS")
ml_EC_16S <- readRDS("machine_learning/16S_EC/ml_EC_16S.RDS")

ml_EC_16S_ACE <- ml_EC_16S[,c(2:2445,2515,2549:2555,2522)]

# filter NAs
ml_EC_16S_ACE$clay <- as.numeric(ml_EC_16S_ACE$clay)
ml_EC_16S_ACE$ace <- as.numeric(ml_EC_16S_ACE$ace)
ml_EC_16S_ACE <- ml_EC_16S_ACE %>%
  filter(!is.na(ace))

# format so : and . are replaced by _ (for varimp)
names(ml_EC_16S_ACE) <- gsub(":","_", names(ml_EC_16S_ACE))
names(ml_EC_16S_ACE) <- gsub("\\.","_", names(ml_EC_16S_ACE))

# split into train and test (4/5 proportion)
ml_EC_16S_ACE$id <- 1:nrow(ml_EC_16S_ACE)
train <- ml_EC_16S_ACE %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_EC_16S_ACE, train, by = 'id')

# get rid of id columns
train <- train[,c(1:2453)]
test <- test[,c(1:2453)]

p = nrow(train)/3

# cforest on training data
cf.ace <- cforest(ace ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
cf.pred <- predict(cf.ace, newdata = test)

# variable importances
ace.imp <- party::varimp(object = cf.ace, conditional = TRUE)

