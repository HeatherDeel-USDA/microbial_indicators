#################################################
# 16S EC tuning for ACE RF Model
# Using cforest methods as described in Strobl et al. 2007
#################################################

### submitted as bash script in Scinet

### libraries
library(tidymodels)
library(workflows)
library(tune)
library(bonsai)

### Predict ACE without clay and climate as predictors
# read in data and subset to correct column
ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_16S.RDS")

# use regular ACE value
# column ranges are: EC numbers, prediction
ml_EC_16S_ACE <- ml_EC_16S[,c(2:2445,2522)]

# filter NAs
ml_EC_16S_ACE$ace <- as.numeric(ml_EC_16S_ACE$ace)
ml_EC_16S_ACE <- ml_EC_16S_ACE %>% 
  filter(!is.na(ace))

soil_split <- initial_split(ml_EC_16S_ACE, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(ace ~ ., data = ml_EC_16S_ACE)
soil_recipe

# specify the model, tune
rf_model <- rand_forest() %>% 
  set_args(mtry = 611, min_n = 5, trees = 3) %>%
  set_engine("partykit", importance = "permutation", replace = FALSE) %>%
  set_mode("regression") %>% 
  translate()
rf_model

# set the workflow
rf_workflow <- workflow() %>%
  add_recipe(soil_recipe) %>%
  add_model(rf_model)

# tune the parameters
rf_grid <- expand.grid(mtry = c(611,1223,1834),
                       trees = c(100,250,500),
                       min_n = c(3,5,7))

# extract results
rf_tune_results <- rf_workflow %>%
  tune_grid(resamples = soil_cv, grid = rf_grid, metrics = metric_set(mae, rmse))

# save tune results
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ACE_model_results/ml_EC_16S_ACE_tune_results_cforest.RDS")

### Predict ACE with clay and climate as predictors
# use regular ACE value
ml_EC_16S_ACE <- ml_EC_16S[,c(2:2445,2515,2549:2555,2522)]

# filter NAs
ml_EC_16S_ACE$clay <- as.numeric(ml_EC_16S_ACE$clay)
ml_EC_16S_ACE$ace <- as.numeric(ml_EC_16S_ACE$ace)
ml_EC_16S_ACE <- ml_EC_16S_ACE %>% 
  filter(!is.na(ace))

soil_split <- initial_split(ml_EC_16S_ACE, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(ace ~ ., data = ml_EC_16S_ACE)
soil_recipe

# specify the model, tune
rf_model <- rand_forest() %>% 
  set_args(trees = tune()) %>%
  set_engine("partykit", importance = c("permutation","impurity"), replace = TRUE) %>%
  set_mode("regression") %>% 
  translate()
rf_model

# set the workflow
rf_workflow <- workflow() %>%
  add_recipe(soil_recipe) %>%
  add_model(rf_model)

# tune the parameters
# rf_grid <- expand.grid(mtry = c(613,1226,1839),
#                        trees = c(1,2,5),
#                        min_n = c(3,5,7))

rf_grid <- expand.grid(trees = c(1,2))

# extract results
rf_tune_results <- rf_workflow %>%
  tune_grid(resamples = soil_cv, grid = rf_grid, metrics = metric_set(mae, rmse))

# save tune results
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ACE_model_results/ml_EC_16S_ACE_tune_results_cforest_clay_climate.RDS")


#### erase later, just using for testing right now
test_performance <- rf_fit %>% collect_metrics()
test_performance

# generate predictions from the test set
test_predictions <- rf_fit %>% collect_predictions()
test_predictions

# graph a regression of predicted vs observed SH_rating values
ACE_EC_lm <- lm(ace ~ .pred, data = test_predictions)
p1 <- ggplot(ACE_EC_lm$model, aes(x = ace, y = .pred)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(title = paste("Adj R2 =",signif(summary(ACE_EC_lm)$adj.r.squared, 2),
                     " P =",signif(summary(ACE_EC_lm)$coef[2,4], 2)),
       x = "Observed ACE Rating", y = "Predicted ACE Rating") +
  theme_bw()
p1

final_model <- fit(rf_workflow, ml_EC_16S_ACE)

ranger_obj <- pull_workflow_fit(final_model)$fit
ranger_obj
var_importance <- as.data.frame(ranger_obj$variable.importance)

### questions
# which importance metric is it using? is it both or is it defaulting? can't find it
# is it actually sampling without replacement using replace = FALSE? need to find a way to test it

### ran the above code 3 times, twice with replace = FALSE and once with replace = TRUE, and the errors were the same every time
### I don't think it's reading the replace parameters, and I'm not sure it's giving both

### need to mess with this, something not right
