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
#ml_TAX_16S <- readRDS("machine_learning/16S_TAX/ml_TAX_16S.RDS")

ml_TAX_16S_AGGSTAB <- ml_TAX_16S[,c(70,104:110,73,152:7614)]

# filter NAs
ml_TAX_16S_AGGSTAB$clay <- as.numeric(ml_TAX_16S_AGGSTAB$clay)
ml_TAX_16S_AGGSTAB$agg_stab <- as.numeric(ml_TAX_16S_AGGSTAB$agg_stab)
ml_TAX_16S_AGGSTAB <- ml_TAX_16S_AGGSTAB %>%
  filter(!is.na(agg_stab))

# format so : and . are replagg_stabd by _ (for varimp)
names(ml_TAX_16S_AGGSTAB) <- gsub(":","_", names(ml_TAX_16S_AGGSTAB))
names(ml_TAX_16S_AGGSTAB) <- gsub("\\.","_", names(ml_TAX_16S_AGGSTAB))

# split into train and test (4/5 proportion)
ml_TAX_16S_AGGSTAB$id <- 1:nrow(ml_TAX_16S_AGGSTAB)
train <- ml_TAX_16S_AGGSTAB %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_TAX_16S_AGGSTAB, train, by = 'id')

# get rid of id columns
train <- train[,c(1:7472)]
test <- test[,c(1:7472)]

p = nrow(train)/3

# cforest on training data
cf.agg_stab <- cforest(agg_stab ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.agg_stab, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/AGGSTAB_model_results_clay_climate/cf.agg_stab", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.agg_stab, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "agg_stab.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.agg_stab <- data.frame(test[,9])
colnames(test.agg_stab)[1] <- "agg_stab.obs"
test.agg_stab <- rownames_to_column(test.agg_stab, var = "id")
cf.pvso <- merge(cf.pred, test.agg_stab, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/AGGSTAB_model_results_clay_climate/cf_obs_vs_pred", args[1], ".RDS", sep = ""))

AGGSTAB_lm <- lm(agg_stab.obs ~ agg_stab.pred, data = cf.pvso)
p1 <- ggplot(AGGSTAB_lm$model, aes(x = agg_stab.obs, y = agg_stab.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(AGGSTAB_lm)$adj.r.squared, 2),
                     " P =",signif(summary(AGGSTAB_lm)$coef[2,4], 2)),
       x = "Observed Aggregate Stability", y = "Predicted Aggregate Stability") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/AGGSTAB_model_results_clay_climate/AGGSTAB_TAX_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(AGGSTAB_lm)$adj.r.squared
p_val <- summary(AGGSTAB_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/AGGSTAB_model_results_clay_climate/agg_stab.ec.stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
agg_stab.imp <- party::varimp(object = cf.agg_stab, conditional = TRUE)
write.csv(agg_stab.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/AGGSTAB_model_results_clay_climate/AGGSTAB_TAX_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)

