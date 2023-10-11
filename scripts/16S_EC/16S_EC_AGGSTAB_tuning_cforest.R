#################################################
# 16S EC tuning for AGGSTAB RF Model
# Using cforest methods as described in Strobl et al. 2007
#################################################

### submitted as bash script in Scinet

### libraries
library(tidymodels)
library(workflows)
library(tune)
library(bonsai)

### Predict AGGSTAB without clay and climate as predictors
# read in data and subset to correct column
ml_EC_16S <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/ml_EC_16S.RDS")

# use regular AGGSTAB value
# column ranges are: EC numbers, prediction
ml_EC_16S_AGGSTAB <- ml_EC_16S[,c(2:2445,2518)]

# filter NAs
ml_EC_16S_AGGSTAB$agg_stab <- as.numeric(ml_EC_16S_AGGSTAB$agg_stab)
ml_EC_16S_AGGSTAB <- ml_EC_16S_AGGSTAB %>% 
  filter(!is.na(agg_stab))

soil_split <- initial_split(ml_EC_16S_AGGSTAB, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(agg_stab ~ ., data = ml_EC_16S_AGGSTAB)
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
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/AGGSTAB_model_results/ml_EC_16S_AGGSTAB_tune_results_cforest.RDS")

### Predict AGGSTAB with clay and climate as predictors
# use regular AGGSTAB value
ml_EC_16S_AGGSTAB <- ml_EC_16S[,c(2:2445,2515,2549:2555,2518)]

# filter NAs
ml_EC_16S_AGGSTAB$clay <- as.numeric(ml_EC_16S_AGGSTAB$clay)
ml_EC_16S_AGGSTAB$agg_stab <- as.numeric(ml_EC_16S_AGGSTAB$agg_stab)
ml_EC_16S_AGGSTAB <- ml_EC_16S_AGGSTAB %>% 
  filter(!is.na(agg_stab))

soil_split <- initial_split(ml_EC_16S_AGGSTAB, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(agg_stab ~ ., data = ml_EC_16S_AGGSTAB)
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
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/16S_EC/AGGSTAB_model_results/ml_EC_16S_AGGSTAB_tune_results_cforest_clay_climate.RDS")

