#################################################
# 16S EC RICH RF Model
# Using cforest within partykit for unbiased variable importances
# Not tuning due to increased computational time, using typical default values
#################################################

### libraries
library(tidymodels)
library(workflows)
library(tune)
library(bonsai)
library(partykit)

# ### Predict RICH without clay and climate as predictors
# # read in data and subset to correct column
# ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_16S.RDS")
# 
# # use regular RICH value
# # column ranges are: EC numbers, prediction
# ml_EC_16S_RICH <- ml_EC_16S[,c(2:2445,2545)]
# 
# # filter NAs
# ml_EC_16S_RICH$rich <- as.numeric(ml_EC_16S_RICH$rich)
# ml_EC_16S_RICH <- ml_EC_16S_RICH %>% 
#   filter(!is.na(rich))

set.seed(1)

# for (i in 1:25) {
#   soil_split <- initial_split(ml_EC_16S_RICH, prop = 4/5)
#   soil_split
#   
#   # extract the train and test sets
#   soil_train <- training(soil_split)
#   soil_test <- testing(soil_split)
#   
#   # cross validation
#   soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)
#   
#   # define a recipe
#   soil_recipe <- recipe(rich ~ ., data = ml_EC_16S_RICH)
#   soil_recipe
#   
#   # specify the model
#   rf_model <- rand_forest() %>%
#     set_args(mtry = 815, min_n = 5, trees = 500) %>%
#     set_engine("partykit") %>%
#     set_mode("regression") %>%
#     translate()
#   rf_model
#   
#   # set the workflow
#   rf_workflow <- workflow() %>%
#     add_recipe(soil_recipe) %>%
#     add_model(rf_model)
#   
#   # fit the model
#   rf_fit <- rf_workflow %>%
#     last_fit(soil_split)
#   rf_fit
#   
#   # save the fit
#   saveRDS(rf_fit, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results/RICH_EC_fit", i, ".RDS", sep = ""))
#   
#   # see how well the model performs
#   test_performance <- rf_fit %>% collect_metrics()
#   test_performance
#   
#   # generate predictions from the test set
#   test_predictions <- rf_fit %>% collect_predictions()
#   test_predictions
#   
#   # graph a regression of predicted vs observed SH_rating values
#   RICH_EC_lm <- lm(rich ~ .pred, data = test_predictions)
#   p1 <- ggplot(RICH_EC_lm$model, aes(x = rich, y = .pred)) +
#     geom_point() +
#     stat_smooth(method = "lm", se = TRUE, level = 0.95) +
#     labs(title = paste("Adj R2 =",signif(summary(RICH_EC_lm)$adj.r.squared, 2),
#                        " P =",signif(summary(RICH_EC_lm)$coef[2,4], 2)),
#          x = "Observed Richness", y = "Predicted Richness") +
#     theme_bw()
#   p1
#   
#   # save R^2 and p-values to files
#   write.table(summary(RICH_EC_lm)$adj.r.squared, file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results/RICH_EC_r2_values.txt", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
#   
#   write.table(summary(RICH_EC_lm)$coef[2,4], file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results/RICH_EC_p_values.txt", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
#   
#   # save plot
#   ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results/RICH_EC_pred_vs_obs", i, ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")
#   
#   # fitting the final model
#   # uses all data that can be tested on a new data set
#   final_model <- fit(rf_workflow, ml_EC_16S_RICH)
#   
#   # save the final model
#   saveRDS(final_model, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results/RICH_EC_final_model", i, ".RDS", sep = ""))
#   
#   # variable importance
#   ranger_obj <- pull_workflow_fit(final_model)$fit
#   ranger_obj
#   var_importance <- as.data.frame(ranger_obj$variable.importance)
#   write.csv(var_importance, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results/RICH_EC_var_importance", i, ".csv", sep = ""), row.names = TRUE)
# }

### Predict RICH with clay and climate as predictors
# use regular RICH value
ml_EC_16S_RICH <- ml_EC_16S[,c(2:2445,2515,2549:2555,2545)]

# filter NAs
ml_EC_16S_RICH$clay <- as.numeric(ml_EC_16S_RICH$clay)
ml_EC_16S_RICH$rich <- as.numeric(ml_EC_16S_RICH$rich)
ml_EC_16S_RICH <- ml_EC_16S_RICH %>%
  filter(!is.na(rich))

for (i in 1:25) {
  soil_split <- initial_split(ml_EC_16S_RICH, prop = 4/5)
  soil_split
  
  # extract the train and test sets
  soil_train <- training(soil_split)
  soil_test <- testing(soil_split)
  
  # cross validation
  soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)
  
  # define a recipe
  soil_recipe <- recipe(rich ~ ., data = ml_EC_16S_RICH)
  soil_recipe
  
  # specify the model
  rf_model <- rand_forest() %>%
    set_args(mtry = 815, min_n = 5, trees = 500) %>%
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
  saveRDS(rf_fit, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results_clay_climate/RICH_EC_fit", i, ".RDS", sep = ""))
  
  # see how well the model performs
  test_performance <- rf_fit %>% collect_metrics()
  test_performance
  
  # generate predictions from the test set
  test_predictions <- rf_fit %>% collect_predictions()
  test_predictions
  
  # graph a regression of predicted vs observed SH_rating values
  RICH_EC_lm <- lm(rich ~ .pred, data = test_predictions)
  p1 <- ggplot(RICH_EC_lm$model, aes(x = rich, y = .pred)) +
    geom_point() +
    stat_smooth(method = "lm", se = TRUE, level = 0.95) +
    labs(title = paste("Adj R2 =",signif(summary(RICH_EC_lm)$adj.r.squared, 2),
                       " P =",signif(summary(RICH_EC_lm)$coef[2,4], 2)),
         x = "Observed Richness", y = "Predicted Richness") +
    theme_bw()
  p1
  
  # save R^2 and p-values to files
  write.table(summary(RICH_EC_lm)$adj.r.squared, file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results_clay_climate/RICH_EC_r2_values.txt", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  write.table(summary(RICH_EC_lm)$coef[2,4], file = "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results_clay_climate/RICH_EC_p_values.txt", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # save plot
  ggsave(paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results_clay_climate/RICH_EC_pred_vs_obs", i, ".pdf", sep = ""), unit = "in", width = 6, height = 6, dpi = 300, device = "pdf")
  
  # fitting the final model
  # uses all data that can be tested on a new data set
  final_model <- fit(rf_workflow, ml_EC_16S_RICH)
  
  # save the final model
  saveRDS(final_model, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results_clay_climate/RICH_EC_final_model", i, ".RDS", sep = ""))
  
  # variable importance
  ranger_obj <- pull_workflow_fit(final_model)$fit
  ranger_obj
  var_importance <- as.data.frame(ranger_obj$variable.importance)
  write.csv(var_importance, paste("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/RICH_model_results_clay_climate/RICH_EC_var_importance", i, ".csv", sep = ""), row.names = TRUE)
}
