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

ml_TAX_16S_RESP <- ml_TAX_16S[,c(70,104:110,79,152:7614)]

# filter NAs
ml_TAX_16S_RESP$clay <- as.numeric(ml_TAX_16S_RESP$clay)
ml_TAX_16S_RESP$resp <- as.numeric(ml_TAX_16S_RESP$resp)
ml_TAX_16S_RESP <- ml_TAX_16S_RESP %>%
  filter(!is.na(resp))

# format so : and . are replrespd by _ (for varimp)
names(ml_TAX_16S_RESP) <- gsub(":","_", names(ml_TAX_16S_RESP))
names(ml_TAX_16S_RESP) <- gsub("\\.","_", names(ml_TAX_16S_RESP))

# split into train and test (4/5 proportion)
ml_TAX_16S_RESP$id <- 1:nrow(ml_TAX_16S_RESP)
train <- ml_TAX_16S_RESP %>% dplyr::sample_frac(0.80)
test <- dplyr::anti_join(ml_TAX_16S_RESP, train, by = 'id')

# get rid of id columns
train <- train[,c(1:7472)]
test <- test[,c(1:7472)]

p = nrow(train)/3

# cforest on training data
cf.resp <- cforest(resp ~ ., data = train,
                  controls = cforest_unbiased(mtry = p/3, ntree = 500))
saveRDS(cf.resp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/RESP_model_results_clay_climate/cf.resp", args[1], ".RDS", sep = ""))
cf.pred <- predict(cf.resp, newdata = test, OOB = TRUE, type = "response")

# observed vs predicted
colnames(cf.pred)[1] <- "resp.pred"
cf.pred <- data.frame(cf.pred)
cf.pred <- rownames_to_column(cf.pred, var = "id")
test.resp <- data.frame(test[,9])
colnames(test.resp)[1] <- "resp.obs"
test.resp <- rownames_to_column(test.resp, var = "id")
cf.pvso <- merge(cf.pred, test.resp, by = "id")

saveRDS(cf.pvso, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/RESP_model_results_clay_climate/cf_obs_vs_pred", args[1], ".RDS", sep = ""))

RESP_lm <- lm(resp.obs ~ resp.pred, data = cf.pvso)
p1 <- ggplot(RESP_lm$model, aes(x = resp.obs, y = resp.pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(RESP_lm)$adj.r.squared, 2),
                     " P =",signif(summary(RESP_lm)$coef[2,4], 2)),
       x = "Observed Respiration", y = "Predicted Respiration") +
  theme_bw()
p1

ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/RESP_model_results_clay_climate/RESP_TAX_pred_vs_obs", args[1], ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")

# save R^2 and p-values to a file
task_num <- args[1]
r2_val <- summary(RESP_lm)$adj.r.squared
p_val <- summary(RESP_lm)$coef[2,4]
print(paste0('Hello! I am task number: ', args[1]))
print(paste0('Hello! My R^2 value is: ', r2_val))
print(paste0('Hello! My p-value is: ', p_val))

write.table(cbind(task_num,r2_val,p_val), 
            file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/RESP_model_results_clay_climate/resp.ec.stats.csv", 
            col.names = c("task_number","r2_val","p_val"),append = TRUE, sep = ",", row.names = FALSE)

# variable importances
resp.imp <- party::varimp(object = cf.resp, conditional = TRUE)
write.csv(resp.imp, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/RESP_model_results_clay_climate/RESP_TAX_var_importance", args[1], ".csv", sep = ""), row.names = TRUE)

