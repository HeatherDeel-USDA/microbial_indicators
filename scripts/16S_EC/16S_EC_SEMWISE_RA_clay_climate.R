#!/usr/bin/env Rscript

# accept command line arguments and save them in a list called args
args = commandArgs(trailingOnly=TRUE)

# print task number
print(paste0('Hello! I am task number: ', args[1]))

# libraries
library(dplyr)
library(party)
library(tidyverse)

ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/NAME OF RA FILE")

# change column numbers as necessary
# include all ECs, clay, ClimateZ, and SH_rating
ml_EC_16S_SEMWISE <- ml_EC_16S[,c(2:2445,2515,2547,2557)]

# filter NAs
ml_EC_16S_SEMWISE$clay <- as.numeric(ml_EC_16S_SEMWISE$clay)
ml_EC_16S_SEMWISE$SH_rating <- as.numeric(ml_EC_16S_SEMWISE$SH_rating)
ml_EC_16S_SEMWISE <- ml_EC_16S_SEMWISE %>%
  filter(!is.na(SH_rating))

# format so : and . are replSH_ratingd by _ (for varimp)
names(ml_EC_16S_SEMWISE) <- gsub(":","_", names(ml_EC_16S_SEMWISE))
names(ml_EC_16S_SEMWISE) <- gsub("\\.","_", names(ml_EC_16S_SEMWISE))

# split into train and test (4/5 proportion)
ml_EC_16S_SEMWISE$id <- 1:nrow(ml_EC_16S_SEMWISE)
train <- ml_EC_16S_SEMWISE %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_EC_16S_SEMWISE, train, by = 'id')

# get rid of id columns
train <- train[,c(1:2447)]
test <- test[,c(1:2447)]
train$ClimateZ <- as.factor(train$ClimateZ)
test$ClimateZ <- as.factor(test$ClimateZ)

p = nrow(train)/3

# cforest on training data
cf.SH_rating <- cforest(SH_rating ~ ., data = train,
                        controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.SH_rating, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/SEMWISE_model_results_RA_clay_climate/SEMWISE_EC_RA_clay_climate_cf", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.SH_rating, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "SH_rating.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.SH_rating <- data.frame(test[,2447])
colnames(test.SH_rating)[1] <- "SH_rating.obs"
test.SH_rating <- rownames_to_column(test.SH_rating, var = "id")
cf.pvso <- merge(cf.pred, test.SH_rating, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/SEMWISE_model_results_RA_clay_climate/SEMWISE_EC_RA_clay_climate_cf_obs_vs_pred", args[1], ".RDS", sep = ""))

SEMWISE_lm <- lm(SH_rating.obs ~ SH_rating.pred, data = cf.pvso)
p1 <- ggplot(SEMWISE_lm$model, aes(x = SH_rating.obs, y = SH_rating.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(SEMWISE_lm)$adj.r.squared, 2),
                     " P =",signif(summary(SEMWISE_lm)$coef[2,4], 2)),
       x = "Observed SEMWISE Rating", y = "Predicted SEMWISE Rating") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/SEMWISE_model_results_RA_clay_climate/SEMWISE_EC_RA_clay_climate_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(SEMWISE_lm)$adj.r.squared
p_val <- summary(SEMWISE_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/SEMWISE_model_results_RA_clay_climate/SEMWISE_EC_RA_clay_climate_stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
SH_rating.imp <- party::varimp(object = cf.SH_rating, conditional = TRUE)
write.csv(SH_rating.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/SEMWISE_model_results_RA_clay_climate/SEMWISE_EC_RA_clay_climate_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)
