#################################################
# 16S EC tuning for SHMI RF Model
# Using cforest methods as described in Strobl et al. 2007
#################################################

### submitted as bash script in Scinet

### libraries
library(tidymodels)
library(workflows)
library(tune)
library(bonsai)

### Predict SHMI without clay and climate as predictors
# read in data and subset to correct column
ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_16S.RDS")

# use regular SHMI value
# column ranges are: EC numbers, prediction
ml_EC_16S_SHMI <- ml_EC_16S[,c(2:2445,2594)]

# filter NAs
ml_EC_16S_SHMI$SHMI2_rating <- as.numeric(ml_EC_16S_SHMI$SHMI2_rating)
ml_EC_16S_SHMI <- ml_EC_16S_SHMI %>% 
  filter(!is.na(SHMI2_rating))

soil_split <- initial_split(ml_EC_16S_SHMI, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(SHMI2_rating ~ ., data = ml_EC_16S_SHMI)
soil_recipe

# specify the model, tune
rf_model <- rand_forest() %>% 
  set_args(mtry = tune(), min_n = tune(), trees = tune()) %>%
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
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/SHMI_model_results/ml_EC_16S_SHMI_tune_results_cforest.RDS")

### Predict SHMI with clay and climate as predictors
# use regular SHMI value
ml_EC_16S_SHMI <- ml_EC_16S[,c(2:2445,2515,2549:2555,2594)]

# filter NAs
ml_EC_16S_SHMI$clay <- as.numeric(ml_EC_16S_SHMI$clay)
ml_EC_16S_SHMI$SHMI2_rating <- as.numeric(ml_EC_16S_SHMI$SHMI2_rating)
ml_EC_16S_SHMI <- ml_EC_16S_SHMI %>% 
  filter(!is.na(SHMI2_rating))

soil_split <- initial_split(ml_EC_16S_SHMI, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(SHMI2_rating ~ ., data = ml_EC_16S_SHMI)
soil_recipe

# specify the model, tune
rf_model <- rand_forest() %>% 
  set_args(mtry = tune(), min_n = tune(), trees = tune()) %>%
  set_engine("partykit", importance = "permutation", replace = FALSE) %>%
  set_mode("regression") %>% 
  translate()
rf_model

# set the workflow
rf_workflow <- workflow() %>%
  add_recipe(soil_recipe) %>%
  add_model(rf_model)

# tune the parameters
rf_grid <- expand.grid(mtry = c(613,1226,1839),
                       trees = c(100,250,500),
                       min_n = c(3,5,7))

# extract results
rf_tune_results <- rf_workflow %>%
  tune_grid(resamples = soil_cv, grid = rf_grid, metrics = metric_set(mae, rmse))

# save tune results
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/SHMI_model_results/ml_EC_16S_SHMI_tune_results_cforest_clay_climate.RDS")

