#############################
#!/usr/bin/env Rscript

# accept command line arguments and save them in a list called args
args = commandArgs(trailingOnly=TRUE)

# print task number
print(paste0('Hello! I am task number: ', args[1]))

# libraries
library(dplyr)
library(party)
library(tidyverse)

ml_TAX_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ml_TAX_16S.RDS")

ml_TAX_16S_WATERCAP <- ml_TAX_16S[,c(70,102,71,152:7614)]

# filter NAs
ml_TAX_16S_WATERCAP$clay <- as.numeric(ml_TAX_16S_WATERCAP$clay)
ml_TAX_16S_WATERCAP$water_cap <- as.numeric(ml_TAX_16S_WATERCAP$water_cap)
ml_TAX_16S_WATERCAP <- ml_TAX_16S_WATERCAP %>%
  filter(!is.na(water_cap))

# format so : and . are replwater_capd by _ (for varimp)
names(ml_TAX_16S_WATERCAP) <- gsub(":","_", names(ml_TAX_16S_WATERCAP))
names(ml_TAX_16S_WATERCAP) <- gsub("\\.","_", names(ml_TAX_16S_WATERCAP))

# split into train and test (4/5 proportion)
ml_TAX_16S_WATERCAP$id <- 1:nrow(ml_TAX_16S_WATERCAP)
train <- ml_TAX_16S_WATERCAP %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_TAX_16S_WATERCAP, train, by = 'id')

# get rid of id columns
train <- train[,c(1:7466)]
test <- test[,c(1:7466)]
train$ClimateZ <- as.factor(train$ClimateZ)
test$ClimateZ <- as.factor(test$ClimateZ)

p = nrow(train)/3

# cforest on training data
cf.water_cap <- cforest(water_cap ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.water_cap, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/WATERCAP_model_results_clay_climate/cf.water_cap", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.water_cap, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "water_cap.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.water_cap <- data.frame(test[,3])
colnames(test.water_cap)[1] <- "water_cap.obs"
test.water_cap <- rownames_to_column(test.water_cap, var = "id")
cf.pvso <- merge(cf.pred, test.water_cap, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/WATERCAP_model_results_clay_climate/cf_obs_vs_pred", args[1], ".RDS", sep = ""))

WATERCAP_lm <- lm(water_cap.obs ~ water_cap.pred, data = cf.pvso)
p1 <- ggplot(WATERCAP_lm$model, aes(x = water_cap.obs, y = water_cap.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(WATERCAP_lm)$adj.r.squared, 2),
                     " P =",signif(summary(WATERCAP_lm)$coef[2,4], 2)),
       x = "Observed Water Capacity", y = "Predicted Water Capacity") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/WATERCAP_model_results_clay_climate/WATERCAP_TAX_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(WATERCAP_lm)$adj.r.squared
p_val <- summary(WATERCAP_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/WATERCAP_model_results_clay_climate/water_cap.ec.stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
water_cap.imp <- party::varimp(object = cf.water_cap, conditional = TRUE)
write.csv(water_cap.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/WATERCAP_model_results_clay_climate/WATERCAP_TAX_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)

