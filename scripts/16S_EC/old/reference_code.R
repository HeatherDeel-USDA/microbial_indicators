### code I wanted to keep for calculating percentages and cumulative sums across a list of data frames
#calculate percent importances
for (i in seq_along(ldf)) {
  ldf[[i]]$importance <- as.numeric(ldf[[i]]$importance)
  sum <- sum(ldf[[i]]$importance)
  ldf[[i]]$importance_percent <- ldf[[i]]$importance / sum * 100
}

ldf <- lapply(ldf, function(df){
  df[order(df$importance_percent),] #orders importance percentages
})

for (i in seq_along(ldf)) {
  ldf[[i]]$importance_percent_cumsum <- cumsum(ldf[[i]]$importance_percent) # cumulated sum of percentages
}

#ldf <- lapply(ldf, function(x) filter(x, importance_percent_cumsum > 50)) # filters cumulated sums < 50
# comment this line to include all features shared between all 25 models

### Make varimp list of GIBBs + important features for each indicator
####### can probably erase this after I run the final models - should be there and not here
####### probably won't do the final models, so can likely get rid of this later
```{r}
varimp_df <- read.csv("machine_learning/16S_EC/final_varimps.csv", header = FALSE)

colnames(varimp_df) <- c("Predictor","VIP","Run","Category","Class","Sub_class",
                         "Name","Short","Book","Missing","VIP_scaled","all_gibbs",
                         "indicator")

# old code for writing predictors to a file that could be used in "final" models for partial dependence plot generation
myVars <- c('ace', 'SOM', 'activeC', 'resp', 'agg_stab', 'water_cap', 'ph', 'p', 'k', 'mg', 'fe', 'mn', 'zn')

for (myVar in myVars) {
  varimp_df_filt <- filter(varimp_df, indicator == myVar)
  varimp_df_list <- varimp_df_filt %>% 
    group_by(Predictor) %>% 
    summarize(VIP_scaled_avg = mean(VIP_scaled))
  write.csv(varimp_df_list, file = paste0("machine_learning/16S_FINAL/", myVar,"_EC_GIBBs_final.csv"),
  row.names = FALSE)
}

# code for partial dependence plots
### Partial dependence plots of final models
### To do:
# rerun these with new model results
# haven't yet rerun final models because we may be reducing colinear features
# I also need to unscale clay when it is included
# may or may not do the final models, so not sure if we need this
# show response of important features to each indicator
# the nutrients aren't scaled, so process those in a separate loop
myVars1 <- c('ace', 'SOM', 'activeC', 'resp', 'agg_stab', 'water_cap')

myVar <- 'agg_stab'

# read in data.pred so I can use the attributes to unscale the indicators
data.pred <- readRDS("metadata/data.pred.SHMI.RDS")

for (myVar in myVars1) {
  # read in and format files
  filenames <- list.files(paste0("machine_learning/16S_FINAL/",myVar,"_model_results/"), pattern = paste0(myVar,"_final_pd_Predictor*"), full.names = TRUE)
  predictors <- lapply(filenames, read.csv, header = FALSE)
  colnames <- c("predictor","predictor_value","indicator_response")
  predictors <- lapply(predictors, setNames, colnames)
  
  # create df with scaling attributes
  scale_attr <- data.frame(attributes(data.pred[[myVar]]))
  scaled_scale <- scale_attr$scaled.scale[1]
  scaled_center <- scale_attr$scaled.center[1]
  
  # unscale
  for (i in seq_along(predictors)) {
    predictors[[i]]$unscaled <- (predictors[[i]]$indicator_response * scaled_scale) + scaled_center
  }
  
  # create functions for making pdps and bar charts
  make_graphs <- function(x) {
    ggplot(x, aes(x = predictor_value, y = unscaled, group = factor(predictor))) +
      geom_line() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      facet_wrap(~predictor)
  }
  make_bargraphs <- function(x) {
    ggplot(x, aes(x = predictor_value, y = unscaled, group = factor(predictor))) +
      geom_bar(stat = 'identity') +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      facet_wrap(~predictor)
  }
  
  # make pdps and climate bar charts
  for (i in seq_along(predictors)) {
    if (predictors[[i]]$predictor[1] == "ClimateZ") { 
      p1 = make_bargraphs(predictors[[i]])
    }
    if (predictors[[i]]$predictor[1] != "ClimateZ") { 
      predictors_EC <- lapply(predictors, function(x) filter(x, predictor != "ClimateZ"))
      predictors_EC <- predictors_EC[sapply(predictors_EC, function(x) dim(x)[1]) > 0]
      
      p2 = lapply(predictors_EC, make_graphs)
    }
  }
  
  # merge graphs and save
  p3 <- ggarrange(p1,p2[[1]],p2[[2]],p2[[3]],p2[[4]], 
                  p2[[5]],p2[[6]],p2[[7]],p2[[8]],p2[[9]],
                  nrow = 2, ncol = 5)
  p3
  annotate_figure(p3, left = text_grob(myVar, rot = 90),
                  bottom = text_grob("Feature Value"))
  #ggsave(paste0("figures/pdp/", myVar, "_pdp.png"), unit = "in", width = 10, height = 4, dpi = 300, device = "png")
}

myVars2 <- c('p', 'k', 'mg', 'fe', 'mn', 'zn','ph')

# read in ml_EC_RA_corr so I can use the attributes to unscale the nutrients and ph
ml_EC_16S <- readRDS("machine_learning/EC_RA_corr/ml_EC_RA_corr.RDS")

# create df with scaling attributes for each indicator
scale_attr <- data.frame(attributes(ml_EC_16S[[myVars2]]))
scaled_scale <- scale_attr$scaled.scale[1]
scaled_center <- scale_attr$scaled.center[1]

for (myVar in myVars1) {
  filenames <- list.files(paste0("machine_learning/16S_FINAL/",myVar,"_model_results/"), pattern = paste0(myVar,"_final_pd_Predictor*"), full.names = TRUE)
  predictors <- lapply(filenames, read.csv, header = FALSE)
  colnames <- c("predictor","predictor_value","indicator_response")
  predictors <- lapply(predictors, setNames, colnames)
  for (i in seq_along(predictors)) {
    predictors[[i]]$unscaled <- predictors[[i]]$indicator_response * scaled_scale + scaled_center
  }
  make_graphs <- function(x) {
    ggplot(x, aes(x = predictor_value, y = unscaled, group = factor(predictor))) +
      geom_line() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      facet_wrap(~predictor)
  }
  make_bargraphs <- function(x) {
    ggplot(x, aes(x = predictor_value, y = unscaled, group = factor(predictor))) +
      geom_bar(stat = 'identity') +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      facet_wrap(~predictor)
  }
  for (i in seq_along(predictors)) {
    if (predictors[[i]]$predictor[1] == "ClimateZ") { 
      p1 = make_bargraphs(predictors[[i]])
    }
    if (predictors[[i]]$predictor[1] != "ClimateZ") { 
      predictors_EC <- lapply(predictors, function(x) filter(x, predictor != "ClimateZ"))
      predictors_EC <- predictors_EC[sapply(predictors_EC, function(x) dim(x)[1]) > 0]
      
      p2 = lapply(predictors_EC, make_graphs)
    }
  }
  p3 <- ggarrange(p1,p2[[1]],p2[[2]],p2[[3]],p2[[4]], 
                  p2[[5]],p2[[6]],p2[[7]],p2[[8]],p2[[9]],
                  nrow = 2, ncol = 5)
  p3
  annotate_figure(p3, left = text_grob(myVar, rot = 90),
                  bottom = text_grob("Feature Value"))
  #ggsave(paste0("figures/pdp/", myVar, "_pdp.png"), unit = "in", width = 10, height = 4, dpi = 300, device = "png")
}


# might need to standardize the y axis limits at some point but this is good enough for now



