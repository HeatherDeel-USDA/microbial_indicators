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

ml_TAX_16S_SHMI <- ml_TAX_16S[,c(70,102,149,152:7614)]

# filter NAs
ml_TAX_16S_SHMI$clay <- as.numeric(ml_TAX_16S_SHMI$clay)
ml_TAX_16S_SHMI$SHMI2_rating <- as.numeric(ml_TAX_16S_SHMI$SHMI2_rating)
ml_TAX_16S_SHMI <- ml_TAX_16S_SHMI %>%
  filter(!is.na(SHMI2_rating))

# format so : and . are replSHMI2_ratingd by _ (for varimp)
names(ml_TAX_16S_SHMI) <- gsub(":","_", names(ml_TAX_16S_SHMI))
names(ml_TAX_16S_SHMI) <- gsub("\\.","_", names(ml_TAX_16S_SHMI))

# split into train and test (4/5 proportion)
ml_TAX_16S_SHMI$id <- 1:nrow(ml_TAX_16S_SHMI)
train <- ml_TAX_16S_SHMI %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_TAX_16S_SHMI, train, by = 'id')

# get rid of id columns
train <- train[,c(1:7466)]
test <- test[,c(1:7466)]
train$ClimateZ <- as.factor(train$ClimateZ)
test$ClimateZ <- as.factor(test$ClimateZ)

p = nrow(train)/3

# cforest on training data
cf.SHMI2_rating <- cforest(SHMI2_rating ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.SHMI2_rating, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SHMI_model_results_clay_climate/cf.SHMI2_rating", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.SHMI2_rating, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "SHMI2_rating.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.SHMI2_rating <- data.frame(test[,3])
colnames(test.SHMI2_rating)[1] <- "SHMI2_rating.obs"
test.SHMI2_rating <- rownames_to_column(test.SHMI2_rating, var = "id")
cf.pvso <- merge(cf.pred, test.SHMI2_rating, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SHMI_model_results_clay_climate/cf_obs_vs_pred", args[1], ".RDS", sep = ""))

SHMI_lm <- lm(SHMI2_rating.obs ~ SHMI2_rating.pred, data = cf.pvso)
p1 <- ggplot(SHMI_lm$model, aes(x = SHMI2_rating.obs, y = SHMI2_rating.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(SHMI_lm)$adj.r.squared, 2),
                     " P =",signif(summary(SHMI_lm)$coef[2,4], 2)),
       x = "Observed SHMI Rating", y = "Predicted SHMI Rating") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SHMI_model_results_clay_climate/SHMI_TAX_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(SHMI_lm)$adj.r.squared
p_val <- summary(SHMI_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SHMI_model_results_clay_climate/SHMI2_rating.ec.stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
SHMI2_rating.imp <- party::varimp(object = cf.SHMI2_rating, conditional = TRUE)
write.csv(SHMI2_rating.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SHMI_model_results_clay_climate/SHMI_TAX_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)

