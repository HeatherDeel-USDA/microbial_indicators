#!/bin/bash
#SBATCH --job-name="18S_EC"
#SBATCH -p debug
#SBATCH -N 1
#SBATCH --threads-per-core=1
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-user=heather.deel@usda.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "stdout.%j.%N"

### activate R module
module load r/4.3.0

### change directory so it can find R packages
cd /project/soil_micro_lab/micro_indicators/R_packages

### start R
Rserve(args="--no-save")

### libraries
library(tidymodels)
library(workflows)
library(tune)

### Predict Overall CASH rating
# read in data and subset to correct column
SHAI18S_ml <- readRDS("/project/soil_micro_lab/micro_indicators/machine_learning/18S_EC/SHAI18S_ml_EC.RDS")

SHAI18S_ml_CASH <- SHAI18S_ml[,c(97,154:1282)]

soil_split <- initial_split(SHAI18S_ml_CASH, prop = 4/5)
soil_split

# extract the train and test sets
soil_train <- training(soil_split)
soil_test <- testing(soil_split)

# cross validation
soil_cv <- vfold_cv(soil_train, v = 5, repeats = 10, strata = NULL)

# define the recipe
soil_recipe <- recipe(Overall ~ ., data = SHAI18S_ml_CASH)
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
rf_grid <- expand.grid(mtry = c(282,565,847),
                       trees = c(100,250,500),
                       min_n = c(3,5,7))

# extract results
rf_tune_results <- rf_workflow %>%
  tune_grid(resamples = soil_cv, grid = rf_grid, metrics = metric_set(mae, rmse))

# save tune results
saveRDS(rf_tune_results, "/project/soil_micro_lab/micro_indicators/machine_learning/18S_EC/CASH_model_results/SHAI18S_ml_CASH_tune_results.RDS")



