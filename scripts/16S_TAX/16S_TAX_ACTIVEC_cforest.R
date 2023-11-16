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

ml_TAX_16S_ACTIVEC <- ml_TAX_16S[,c(70,102,81,152:7614)]

# filter NAs
ml_TAX_16S_ACTIVEC$clay <- as.numeric(ml_TAX_16S_ACTIVEC$clay)
ml_TAX_16S_ACTIVEC$activeC <- as.numeric(ml_TAX_16S_ACTIVEC$activeC)
ml_TAX_16S_ACTIVEC <- ml_TAX_16S_ACTIVEC %>%
  filter(!is.na(activeC))

# format so : and . are replactiveCd by _ (for varimp)
names(ml_TAX_16S_ACTIVEC) <- gsub(":","_", names(ml_TAX_16S_ACTIVEC))
names(ml_TAX_16S_ACTIVEC) <- gsub("\\.","_", names(ml_TAX_16S_ACTIVEC))

# split into train and test (4/5 proportion)
ml_TAX_16S_ACTIVEC$id <- 1:nrow(ml_TAX_16S_ACTIVEC)
train <- ml_TAX_16S_ACTIVEC %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_TAX_16S_ACTIVEC, train, by = 'id')

# get rid of id columns
train <- train[,c(1:7466)]
test <- test[,c(1:7466)]
train$ClimateZ <- as.factor(train$ClimateZ)
test$ClimateZ <- as.factor(test$ClimateZ)

p = nrow(train)/3

# cforest on training data
cf.activeC <- cforest(activeC ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.activeC, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACTIVEC_model_results_clay_climate/cf.activeC", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.activeC, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "activeC.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.activeC <- data.frame(test[,3])
colnames(test.activeC)[1] <- "activeC.obs"
test.activeC <- rownames_to_column(test.activeC, var = "id")
cf.pvso <- merge(cf.pred, test.activeC, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACTIVEC_model_results_clay_climate/cf_obs_vs_pred", args[1], ".RDS", sep = ""))

ACTIVEC_lm <- lm(activeC.obs ~ activeC.pred, data = cf.pvso)
p1 <- ggplot(ACTIVEC_lm$model, aes(x = activeC.obs, y = activeC.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(ACTIVEC_lm)$adj.r.squared, 2),
                     " P =",signif(summary(ACTIVEC_lm)$coef[2,4], 2)),
       x = "Observed Active C", y = "Predicted Active C") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACTIVEC_model_results_clay_climate/ACTIVEC_TAX_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(ACTIVEC_lm)$adj.r.squared
p_val <- summary(ACTIVEC_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACTIVEC_model_results_clay_climate/activeC.ec.stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
activeC.imp <- party::varimp(object = cf.activeC, conditional = TRUE)
write.csv(activeC.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACTIVEC_model_results_clay_climate/ACTIVEC_TAX_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)

