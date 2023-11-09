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

ml_TAX_16S_SOM <- ml_TAX_16S[,c(70,104:110,75,152:7614)]

# filter NAs
ml_TAX_16S_SOM$clay <- as.numeric(ml_TAX_16S_SOM$clay)
ml_TAX_16S_SOM$SOM <- as.numeric(ml_TAX_16S_SOM$SOM)
ml_TAX_16S_SOM <- ml_TAX_16S_SOM %>%
  filter(!is.na(SOM))

# format so : and . are replSOMd by _ (for varimp)
names(ml_TAX_16S_SOM) <- gsub(":","_", names(ml_TAX_16S_SOM))
names(ml_TAX_16S_SOM) <- gsub("\\.","_", names(ml_TAX_16S_SOM))

# split into train and test (4/5 proportion)
ml_TAX_16S_SOM$id <- 1:nrow(ml_TAX_16S_SOM)
train <- ml_TAX_16S_SOM %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_TAX_16S_SOM, train, by = 'id')

# get rid of id columns
train <- train[,c(1:7472)]
test <- test[,c(1:7472)]

p = nrow(train)/3

# cforest on training data
cf.SOM <- cforest(SOM ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.SOM, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SOM_model_results_clay_climate/cf.SOM", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.SOM, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "SOM.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.SOM <- data.frame(test[,9])
colnames(test.SOM)[1] <- "SOM.obs"
test.SOM <- rownames_to_column(test.SOM, var = "id")
cf.pvso <- merge(cf.pred, test.SOM, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SOM_model_results_clay_climate/cf_obs_vs_pred", args[1], ".RDS", sep = ""))

SOM_lm <- lm(SOM.obs ~ SOM.pred, data = cf.pvso)
p1 <- ggplot(SOM_lm$model, aes(x = SOM.obs, y = SOM.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(SOM_lm)$adj.r.squared, 2),
                     " P =",signif(summary(SOM_lm)$coef[2,4], 2)),
       x = "Observed SOM", y = "Predicted SOM") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SOM_model_results_clay_climate/SOM_TAX_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(SOM_lm)$adj.r.squared
p_val <- summary(SOM_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SOM_model_results_clay_climate/SOM.ec.stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
SOM.imp <- party::varimp(object = cf.SOM, conditional = TRUE)
write.csv(SOM.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/SOM_model_results_clay_climate/SOM_TAX_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)

