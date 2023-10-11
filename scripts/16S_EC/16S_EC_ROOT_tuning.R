#################################################
# 16S EC tuning for ROOT RF Model
#################################################

### submitted as bash script in Scinet

### libraries
library(tidymodels)
library(workflows)
library(tune)
library(ranger)

### Predict ROOT without clay and climate as predictors
# read in data and subset to correct column
ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_16S.RDS")

# use regular ROOT value
# column ranges are: EC numbers, prediction
ml_EC_16S_ROOT <- ml_EC_16S[,c(2:2445,2544)]

# filter NAs
ml_EC_16S_ROOT$root <- as.numeric(ml_EC_16S_ROOT$root)
ml_EC_16S_ROOT <- ml_EC_16S_ROOT %>% 
  filter(!is.na(root))

soil_split <- initial_split(ml_EC_16S_ROOT, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(root ~ ., data = ml_EC_16S_ROOT)
soil_recipe

# specify the model, tune
rf_model <- rand_forest() %>% 
  set_args(mtry = tune(), trees =tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("regression")

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
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ROOT_model_results/ml_EC_16S_ROOT_tune_results.RDS")

### Predict ROOT with clay and climate as predictors
# use regular ROOT value
ml_EC_16S_ROOT <- ml_EC_16S[,c(2:2445,2515,2549:2555,2544)]

# filter NAs
ml_EC_16S_ROOT$clay <- as.numeric(ml_EC_16S_ROOT$clay)
ml_EC_16S_ROOT$root <- as.numeric(ml_EC_16S_ROOT$root)
ml_EC_16S_ROOT <- ml_EC_16S_ROOT %>% 
  filter(!is.na(root))

soil_split <- initial_split(ml_EC_16S_ROOT, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(root ~ ., data = ml_EC_16S_ROOT)
soil_recipe

# specify the model, tune
rf_model <- rand_forest() %>% 
  set_args(mtry = tune(), trees =tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("regression")

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
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ROOT_model_results/ml_EC_16S_ROOT_tune_results_clay_climate.RDS")

