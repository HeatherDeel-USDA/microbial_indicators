#################################################
# 16S TAX ACE RF Model
# Using cforest within partykit for unbiased variable importances
# Not tuning due to increased computational time, using typical default values
#################################################

### libraries
library(tidymodels)
library(workflows)
library(tune)
library(bonsai)
library(partykit)

### read in data
ml_TAX_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ml_TAX_16S.RDS")

### Predict ACE with clay and climate as predictors
# use regular ACE value
ml_TAX_16S_ACE <- ml_TAX_16S[,c(70,104:110,77,152:7614)]

# filter NAs
ml_TAX_16S_ACE$clay <- as.numeric(ml_TAX_16S_ACE$clay)
ml_TAX_16S_ACE$ace <- as.numeric(ml_TAX_16S_ACE$ace)
ml_TAX_16S_ACE <- ml_TAX_16S_ACE %>%
  filter(!is.na(ace))

set.seed(1)

for (i in 1:25) {
  soil_split <- initial_split(ml_TAX_16S_ACE, prop = 4/5)
  soil_split
  
  # extract the train and test sets
  soil_train <- training(soil_split)
  soil_test <- testing(soil_split)
  
  # cross validation
  soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)
  
  # define a recipe
  soil_recipe <- recipe(ace ~ ., data = ml_TAX_16S_ACE)
  soil_recipe
  
  # specify the model
  rf_model <- rand_forest() %>%
    set_args(mtry = 2490, min_n = 5, trees = 500) %>%
    set_engine("partykit") %>%
    set_mode("regression") %>%
    translate()
  rf_model
  
  # set the workflow
  rf_workflow <- workflow() %>%
    add_recipe(soil_recipe) %>%
    add_model(rf_model)
  
  # fit the model
  rf_fit <- rf_workflow %>%
    last_fit(soil_split)
  rf_fit
  
  # save the fit
  saveRDS(rf_fit, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ACE_TAX_fit", i, ".RDS", sep = ""))
  
  # see how well the model performs
  test_performance <- rf_fit %>% collect_metrics()
  test_performance
  
  # generate predictions from the test set
  test_predictions <- rf_fit %>% collect_predictions()
  test_predictions
  
  # graph a regression of predicted vs observed SH_rating values
  ACE_TAX_lm <- lm(ace ~ .pred, data = test_predictions)
  p1 <- ggplot(ACE_TAX_lm$model, aes(x = ace, y = .pred)) +
    geom_point() +
    stat_smooth(method = "lm", se = TRUE, level = 0.95) +
    labs(title = paste("Adj R2 =",signif(summary(ACE_TAX_lm)$adj.r.squared, 2),
                       " P =",signif(summary(ACE_TAX_lm)$coef[2,4], 2)),
         x = "Observed ACE", y = "Predicted ACE") +
    theme_bw()
  p1
  
  # save R^2 and p-values to files
  write.table(summary(ACE_TAX_lm)$adj.r.squared, file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ACE_TAX_r2_values.txt", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  write.table(summary(ACE_TAX_lm)$coef[2,4], file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ACE_TAX_p_values.txt", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # save plot
  ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ACE_TAX_pred_vs_obs", i, ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")
  
  # fitting the final model
  # uses all data that can be tested on a new data set
  final_model <- fit(rf_workflow, ml_TAX_16S_ACE)
  
  # save the final model
  saveRDS(final_model, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ACE_TAX_final_model", i, ".RDS", sep = ""))
  
  # variable importance
  ranger_obj <- pull_workflow_fit(final_model)$fit
  ranger_obj
  var_importance <- as.data.frame(ranger_obj$variable.importance)
  write.csv(var_importance, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_TAX/ACE_model_results_clay_climate/ACE_TAX_var_importance", i, ".csv", sep = ""), row.names = TRUE)
}
