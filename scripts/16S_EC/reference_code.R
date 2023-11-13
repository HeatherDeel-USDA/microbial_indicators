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