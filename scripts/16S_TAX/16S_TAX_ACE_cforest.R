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

ml_TAX_16S_ACE <- ml_TAX_16S[,c(70,102,77,152:7614)]

# filter NAs
ml_TAX_16S_ACE$clay <- as.numeric(ml_TAX_16S_ACE$clay)
ml_TAX_16S_ACE$ace <- as.numeric(ml_TAX_16S_ACE$ace)
ml_TAX_16S_ACE <- ml_TAX_16S_ACE %>%
  filter(!is.na(ace))

# format so : and . are replaced by _ (for varimp)
names(ml_TAX_16S_ACE) <- gsub(":","_", names(ml_TAX_16S_ACE))
names(ml_TAX_16S_ACE) <- gsub("\\.","_", names(ml_TAX_16S_ACE))

# split into train and test (4/5 proportion)
ml_TAX_16S_ACE$id <- 1:nrow(ml_TAX_16S_ACE)
train <- ml_TAX_16S_ACE %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_TAX_16S_ACE, train, by = 'id')

# get rid of id columns
train <- train[,c(1:7466)]
test <- test[,c(1:7466)]
train$ClimateZ <- as.factor(train$ClimateZ)
test$ClimateZ <- as.factor(test$ClimateZ)

p = nrow(train)/3

# cforest on training data
cf.ace <- cforest(ace ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.ace, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/cf.ace", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.ace, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "ace.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.ace <- data.frame(test[,3])
colnames(test.ace)[1] <- "ace.obs"
test.ace <- rownames_to_column(test.ace, var = "id")
cf.pvso <- merge(cf.pred, test.ace, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/cf_obs_vs_pred", args[1], ".RDS", sep = ""))

ACE_lm <- lm(ace.obs ~ ace.pred, data = cf.pvso)
p1 <- ggplot(ACE_lm$model, aes(x = ace.obs, y = ace.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(ACE_lm)$adj.r.squared, 2),
                     " P =",signif(summary(ACE_lm)$coef[2,4], 2)),
       x = "Observed ACE", y = "Predicted ACE") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ACE_TAX_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(ACE_lm)$adj.r.squared
p_val <- summary(ACE_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ace.ec.stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
ace.imp <- party::varimp(object = cf.ace, conditional = TRUE)
write.csv(ace.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ACE_TAX_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)

