#!/usr/bin/env Rscript

# accept command line arguments and save them in a list called args
args = commandArgs(trailingOnly=TRUE)

# print task number
print(paste0('Hello! I am task number: ', args[1]))

# libraries
library(dplyr)
library(party)
library(tidyverse)

ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_16S.RDS")

ml_EC_16S_DIST <- ml_EC_16S[,c(2:2445,2515,2549:2555,2546)]

# filter NAs
ml_EC_16S_DIST$clay <- as.numeric(ml_EC_16S_DIST$clay)
ml_EC_16S_DIST$dist <- as.numeric(ml_EC_16S_DIST$dist)
ml_EC_16S_DIST <- ml_EC_16S_DIST %>%
  filter(!is.na(dist))

# format so : and . are repldistd by _ (for varimp)
names(ml_EC_16S_DIST) <- gsub(":","_", names(ml_EC_16S_DIST))
names(ml_EC_16S_DIST) <- gsub("\\.","_", names(ml_EC_16S_DIST))

# split into train and test (4/5 proportion)
ml_EC_16S_DIST$id <- 1:nrow(ml_EC_16S_DIST)
train <- ml_EC_16S_DIST %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_EC_16S_DIST, train, by = 'id')

# get rid of id columns
train <- train[,c(1:2453)]
test <- test[,c(1:2453)]

p = nrow(train)/3

# cforest on training data
cf.dist <- cforest(dist ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.dist, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/DIST_model_results_clay_climate/cf.dist", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.dist, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "dist.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.dist <- data.frame(test[,2453])
colnames(test.dist)[1] <- "dist.obs"
test.dist <- rownames_to_column(test.dist, var = "id")
cf.pvso <- merge(cf.pred, test.dist, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/DIST_model_results_clay_climate/cf_obs_vs_pred", args[1], ".RDS", sep = ""))

DIST_lm <- lm(dist.obs ~ dist.pred, data = cf.pvso)
p1 <- ggplot(DIST_lm$model, aes(x = dist.obs, y = dist.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(DIST_lm)$adj.r.squared, 2),
                     " P =",signif(summary(DIST_lm)$coef[2,4], 2)),
       x = "Observed Disturbance Index", y = "Predicted Disturbance Index") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/DIST_model_results_clay_climate/DIST_EC_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(DIST_lm)$adj.r.squared
p_val <- summary(DIST_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/DIST_model_results_clay_climate/dist.ec.stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
dist.imp <- party::varimp(object = cf.dist, conditional = TRUE)
write.csv(dist.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/DIST_model_results_clay_climate/DIST_EC_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)
