### This file contains a bunch of code for things that I ended up not needing but didn't want to delete

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

### code for graphing relative abundance of modules
# create new TAX file, need full EC names
EC_full <- read_csv("ec_files_expasy/enzyme_R.csv")
EC_full <- column_to_rownames(EC_full, var = "Predictor")
EC_full <- as.matrix(EC_full)
TAX <- tax_table(EC_full)

# create sample file with just sample name, SHMI, SEMWISE, for merging with mod files later
EC_sam_SHMI_SEM <- EC_sam
EC_sam_SHMI_SEM <- rownames_to_column(EC_sam_SHMI_SEM, var = "Sample")
EC_sam_SHMI_SEM <- EC_sam_SHMI_SEM %>% 
  select(Sample, SH_rating, SHMI2_rating, SHMI2_bin, SH_bin)

# start with ECs2_t, filter out extra columns
ECs2_t2 <- ECs2_t[,1:537]
ECs2_t2 <- column_to_rownames(ECs2_t2, var = "ECs")

# make column names the sample names
sample_names <- ECs$SampleID
colnames(ECs2_t2) <- sample_names

# make ECs a column for joining
ECs2_t2 <- rownames_to_column(ECs2_t2, var = "Predictor")

# join modules and ECs
ECs_mod <- ECs2_t2 %>% 
  left_join(modules, by = join_by(Predictor))

ECs_mod <- column_to_rownames(ECs_mod, var = "Predictor")
ECs_mod <- ECs_mod[,c(1:536,543)]
ECs_mod <- ECs_mod %>% 
  filter(module != "NA")

# replace the M in module names so I can convert the column to numeric
ECs_mod$module <- gsub("M","", ECs_mod$module)
ECs_mod$module <- as.numeric(ECs_mod$module)

# create module 1 file, phyloseq, and df
ECs_mod8 <- ECs_mod %>% 
  filter(module == 8)
ECs_mod8 <- ECs_mod8[,1:536] # get rid of module col from otu table
OTU8 = otu_table(ECs_mod8, taxa_are_rows = TRUE)
p_mod8 <- phyloseq(OTU8, TAX, SAM)
p_mod8_melt <- psmelt(p_mod8)
p_mod8_persamp <- p_mod8_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod8_persamp$module <- "8"
p_mod8_persamp <- left_join(p_mod8_persamp, EC_sam_SHMI_SEM, by = "Sample")

##### quick graph of ECs_mod1, want to see if they're monotonic
# maybe do a correlation heatmap of everything to verify
library(ggcorrplot)
p.mat <- cor_pmat(ECs_mod8) 
ggcorrplot(ECs_mod8)
ggsave("figures/corr_heatmap_test_mod8.pdf", unit = "in", width = 30, height = 30, dpi = 300, device = "pdf")

# create module 126 file, phyloseq, and df
ECs_mod126 <- ECs_mod %>% 
  filter(module == 126)
ECs_mod126 <- ECs_mod126[,1:536] # get rid of module col from otu table
OTU126 = otu_table(ECs_mod126, taxa_are_rows = TRUE)
p_mod126 <- phyloseq(OTU126, TAX, SAM)
p_mod126_melt <- psmelt(p_mod126)
p_mod126_persamp <- p_mod126_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod126_persamp$module <- "126"
p_mod126_persamp <- left_join(p_mod126_persamp, EC_sam_SHMI_SEM, by = "Sample")

# create module 206 file, phyloseq, and df
ECs_mod206 <- ECs_mod %>% 
  filter(module == 206)
ECs_mod206 <- ECs_mod206[,1:536] # get rid of module col from otu table
OTU206 = otu_table(ECs_mod206, taxa_are_rows = TRUE)
p_mod206 <- phyloseq(OTU206, TAX, SAM)
p_mod206_melt <- psmelt(p_mod206)
p_mod206_persamp <- p_mod206_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod206_persamp$module <- "206"
p_mod206_persamp <- left_join(p_mod206_persamp, EC_sam_SHMI_SEM, by = "Sample")

# create module 23 file, phyloseq, and df
ECs_mod23 <- ECs_mod %>% 
  filter(module == 23)
ECs_mod23 <- ECs_mod23[,1:536] # get rid of module col from otu table
OTU23 = otu_table(ECs_mod23, taxa_are_rows = TRUE)
p_mod23 <- phyloseq(OTU23, TAX, SAM)
p_mod23_melt <- psmelt(p_mod23)
p_mod23_persamp <- p_mod23_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod23_persamp$module <- "23"
p_mod23_persamp <- left_join(p_mod23_persamp, EC_sam_SHMI_SEM, by = "Sample")

# create module 240 file, phyloseq, and df
ECs_mod240 <- ECs_mod %>% 
  filter(module == 240)
ECs_mod240 <- ECs_mod240[,1:536] # get rid of module col from otu table
OTU240 = otu_table(ECs_mod240, taxa_are_rows = TRUE)
p_mod240 <- phyloseq(OTU240, TAX, SAM)
p_mod240_melt <- psmelt(p_mod240)
p_mod240_persamp <- p_mod240_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod240_persamp$module <- "240"
p_mod240_persamp <- left_join(p_mod240_persamp, EC_sam_SHMI_SEM, by = "Sample")

# create module 272 file, phyloseq, and df
ECs_mod272 <- ECs_mod %>% 
  filter(module == 272)
ECs_mod272 <- ECs_mod272[,1:536] # get rid of module col from otu table
OTU272 = otu_table(ECs_mod272, taxa_are_rows = TRUE)
p_mod272 <- phyloseq(OTU272, TAX, SAM)
p_mod272_melt <- psmelt(p_mod272)
p_mod272_persamp <- p_mod272_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod272_persamp$module <- "272"
p_mod272_persamp <- left_join(p_mod272_persamp, EC_sam_SHMI_SEM, by = "Sample")

# create module 37 file, phyloseq, and df
ECs_mod37 <- ECs_mod %>% 
  filter(module == 37)
ECs_mod37 <- ECs_mod37[,1:536] # get rid of module col from otu table
OTU37 = otu_table(ECs_mod37, taxa_are_rows = TRUE)
p_mod37 <- phyloseq(OTU37, TAX, SAM)
p_mod37_melt <- psmelt(p_mod37)
p_mod37_persamp <- p_mod37_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod37_persamp$module <- "37"
p_mod37_persamp <- left_join(p_mod37_persamp, EC_sam_SHMI_SEM, by = "Sample")

# create module 5 file, phyloseq, and df
ECs_mod5 <- ECs_mod %>% 
  filter(module == 5)
ECs_mod5 <- ECs_mod5[,1:536] # get rid of module col from otu table
OTU5 = otu_table(ECs_mod5, taxa_are_rows = TRUE)
p_mod5 <- phyloseq(OTU5, TAX, SAM)
p_mod5_melt <- psmelt(p_mod5)
p_mod5_persamp <- p_mod5_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod5_persamp$module <- "5"
p_mod5_persamp <- left_join(p_mod5_persamp, EC_sam_SHMI_SEM, by = "Sample")

# create module 72 file, phyloseq, and df
ECs_mod72 <- ECs_mod %>% 
  filter(module == 72)
ECs_mod72 <- ECs_mod72[,1:536] # get rid of module col from otu table
OTU72 = otu_table(ECs_mod72, taxa_are_rows = TRUE)
p_mod72 <- phyloseq(OTU72, TAX, SAM)
p_mod72_melt <- psmelt(p_mod72)
p_mod72_persamp <- p_mod72_melt %>% 
  group_by(Sample) %>% 
  summarize(Abundance = sum(Abundance))
p_mod72_persamp$module <- "72"
p_mod72_persamp <- left_join(p_mod72_persamp, EC_sam_SHMI_SEM, by = "Sample")

# merge all module dfs
ace_EC_mods <- full_join(p_mod1_persamp, p_mod126_persamp)
ace_EC_mods <- full_join(ace_EC_mods, p_mod206_persamp)
ace_EC_mods <- full_join(ace_EC_mods, p_mod23_persamp)
ace_EC_mods <- full_join(ace_EC_mods, p_mod240_persamp)
ace_EC_mods <- full_join(ace_EC_mods, p_mod272_persamp)
ace_EC_mods <- full_join(ace_EC_mods, p_mod37_persamp)
ace_EC_mods <- full_join(ace_EC_mods, p_mod5_persamp)
ace_EC_mods <- full_join(ace_EC_mods, p_mod72_persamp)

# line graph, also playing with boxplot and SHMI2_bin - unsure which is best, but either way with ACE I'm not seeing much of a difference
# may need to put module 1 in it's own graph or something? Not seeing the trends per module like I wanted to
p_mod_ace <- ggplot(ace_EC_mods, aes(x = SHMI2_bin, y = Abundance, color = module)) +
  geom_boxplot()
p_mod_ace

### Relative abundance of ECs by subclasses and subsubclasses
ECs <- read.csv("machine_learning/EC_RA_corr/ml_EC_RA_corr.csv", check.names = FALSE)

ECs2 <- ECs[,4:2439]
ECs2_t <- t(ECs2)
ECs2_t <- data.frame(ECs2_t)
ECs2_t <- rownames_to_column(ECs2_t, var = "ECs")

# collapse EC list to subclasses
pattern <- "EC:(.*)\\.(.*)\\."
ECs2_t$ECs_subclass <- regmatches(ECs2_t$ECs, regexec(pattern, ECs2_t$ECs))
ECs2_t$ECs_subclass2 <- sapply(ECs2_t$ECs_subclass, "[[", 2)
ECs2_sc <- ECs2_t[,c(2:537,539)]
ECs2_sc <- ECs2_sc %>% 
  group_by(ECs_subclass2) %>% 
  summarise_all(., .funs = sum)

# add EC: back as prefix
ECs2_sc$ECs_subclass2 <- paste0("EC:", ECs2_sc$ECs_subclass2)

# format "otu" table for importing into phyloseq
ECs2_sc <- column_to_rownames(ECs2_sc, var = "ECs_subclass2")

# attach sample names to column names
sample_names <- ECs$SampleID
colnames(ECs2_sc) <- sample_names

OTU = otu_table(ECs2_sc, taxa_are_rows = TRUE)

# format sample data for importing into phyloseq
EC_sam <- ECs[,c(2:3,2440:2588)]
EC_sam <- column_to_rownames(EC_sam, var = "SampleID")
EC_sam$SHMI2_bin <- factor(EC_sam$SHMI2_bin, levels = c("very low","low","medium","high","very high"))
EC_sam$SH_bin <- factor(EC_sam$SH_bin, levels = c("very low","low","medium","high","very high"))

SAM <- sample_data(EC_sam)

# use the EC description as "taxa"
EC_desc <- read_xlsx("ec_files_expasy/enzclass_R.xlsx")
EC_desc <- column_to_rownames(EC_desc, var = "EC")
EC_desc <- as.matrix(EC_desc)
TAX <- tax_table(EC_desc)

# create phyloseq object
p_subclass <- phyloseq(OTU, TAX, SAM)

# make SHMI bar chart
p_df_subclass <- psmelt(p_subclass)
ggplot(p_df_subclass, aes(x = SHMI2_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "SHMI Bin") +
  scale_fill_discrete(name = "Subclass")
#ggsave("figures/EC_abund_SHMI_subclass.pdf", unit = "in", width = 6, height = 4, dpi = 300, device = "pdf")

# make SEMWISE bar chart
ggplot(p_df_subclass, aes(x = SH_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "SEMWISE Bin") +
  scale_fill_discrete(name = "Subclass")
#ggsave("figures/EC_abund_SEMWISE_subclass.pdf", unit = "in", width = 6, height = 4, dpi = 300, device = "pdf")

# graph each subclass on its own
# in p_df, create a column with just subclass and filter by that
p_df_subclass$subclass = substr(p_df_subclass$OTU, start = 4, stop = 4)
unique(p_df_subclass$subclass)

# create df for each subclass
# 1 = oxidoreductases
# 2 = transferases
# 3 = hydrolases
# 4 = lyases
# 5 = isomerases
# 6 = ligases

p_df_subclass_1 <- subset(p_df_subclass, subclass == "1")
p_df_subclass_2 <- subset(p_df_subclass, subclass == "2")
p_df_subclass_3 <- subset(p_df_subclass, subclass == "3")
p_df_subclass_4 <- subset(p_df_subclass, subclass == "4")
p_df_subclass_5 <- subset(p_df_subclass, subclass == "5")
p_df_subclass_6 <- subset(p_df_subclass, subclass == "6")

# SHMI bar charts for all levels, combine and save
p1 <- ggplot(p_df_subclass_1, aes(x = SHMI2_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Oxidoreductases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p2 <- ggplot(p_df_subclass_2, aes(x = SHMI2_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Transferases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p3 <- ggplot(p_df_subclass_3, aes(x = SHMI2_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Hydrolases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p4 <- ggplot(p_df_subclass_4, aes(x = SHMI2_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Lyases") +
  scale_fill_discrete(name = "Subclass") + 
  ylim(0,250)
p5 <- ggplot(p_df_subclass_5, aes(x = SHMI2_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Isomerases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p6 <- ggplot(p_df_subclass_6, aes(x = SHMI2_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Ligases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p_SHMI_all <- ggarrange(p1,p2,p3,p4,p5,p6,
                        nrow = 2, ncol = 3)
p_SHMI_all <- annotate_figure(p_SHMI_all, 
                              bottom = text_grob("SHMI Bin"),
                              left = text_grob("Abundance", rot = 90))
p_SHMI_all
#ggsave("figures/EC_abund_SHMI_subclass_split.pdf", unit = "in", width = 12, height = 8, dpi = 300, device = "pdf")

# SEMWISE bar charts for all levels, combine and save
p1 <- ggplot(p_df_subclass_1, aes(x = SH_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Oxidoreductases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p2 <- ggplot(p_df_subclass_2, aes(x = SH_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Transferases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p3 <- ggplot(p_df_subclass_3, aes(x = SH_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Hydrolases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p4 <- ggplot(p_df_subclass_4, aes(x = SH_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Lyases") +
  scale_fill_discrete(name = "Subclass") + 
  ylim(0,250)
p5 <- ggplot(p_df_subclass_5, aes(x = SH_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Isomerases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p6 <- ggplot(p_df_subclass_6, aes(x = SH_bin, y = Abundance, fill = Family)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Ligases") +
  scale_fill_discrete(name = "Subclass") +
  ylim(0,250)
p_SEMWISE_all <- ggarrange(p1,p2,p3,p4,p5,p6,
                           nrow = 2, ncol = 3)
p_SEMWISE_all <- annotate_figure(p_SEMWISE_all, 
                                 bottom = text_grob("SEMWISE Bin"),
                                 left = text_grob("Abundance", rot = 90))
p_SEMWISE_all
#ggsave("figures/EC_abund_SEMWISE_subclass_split.pdf", unit = "in", width = 12, height = 8, dpi = 300, device = "pdf")

# subsubclass level
# collapse EC list to subsubclasses
pattern <- "EC:(.*)\\.(.*)\\.(.*)\\."
ECs2_t$ECs_subsubclass <- regmatches(ECs2_t$ECs, regexec(pattern, ECs2_t$ECs))
ECs2_t$ECs_subsubclass2 <- sapply(ECs2_t$ECs_subsubclass, "[[", 1)
ECs2_ssc <- ECs2_t[,c(2:537,541)]
ECs2_ssc$ECs_subsubclass2 = substr(ECs2_ssc$ECs_subsubclass2, 1, nchar(ECs2_ssc$ECs_subsubclass2)-1)
ECs2_ssc <- ECs2_ssc %>% 
  group_by(ECs_subsubclass2) %>% 
  summarise_all(., .funs = sum)

# format "otu" table for importing into phyloseq
ECs2_ssc <- column_to_rownames(ECs2_ssc, var = "ECs_subsubclass2")

# attach sample names to column names
sample_names <- ECs$SampleID
colnames(ECs2_ssc) <- sample_names

OTU = otu_table(ECs2_ssc, taxa_are_rows = TRUE)

# format sample data for importing into phyloseq
EC_sam <- ECs[,c(2:3,2440:2588)]
EC_sam <- column_to_rownames(EC_sam, var = "SampleID")
EC_sam$SHMI2_bin <- factor(EC_sam$SHMI2_bin, levels = c("very low","low","medium","high","very high"))
EC_sam$SH_bin <- factor(EC_sam$SH_bin, levels = c("very low","low","medium","high","very high"))

SAM <- sample_data(EC_sam)

# use the EC description as "taxa"
EC_desc <- read_xlsx("ec_files_expasy/enzclass_R.xlsx")
EC_desc <- column_to_rownames(EC_desc, var = "EC")
EC_desc <- as.matrix(EC_desc)
TAX <- tax_table(EC_desc)

# create phyloseq object
p_subsubclass <- phyloseq(OTU, TAX, SAM)

# make SHMI bar chart
p_df_subsubclass <- psmelt(p_subsubclass)
ggplot(p_df_subsubclass, aes(x = SHMI2_bin, y = Abundance, fill = Genus)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "SHMI Bin") +
  scale_fill_discrete(name = "Sub-subclass")
#ggsave("figures/EC_abund_SHMI_subsubclass.pdf", unit = "in", width = 6, height = 4, dpi = 300, device = "pdf")

# make SEMWISE bar chart
ggplot(p_df_subsubclass, aes(x = SH_bin, y = Abundance, fill = Genus)) +
  geom_bar(stat = "summary", fun = "mean") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "SEMWISE Bin") +
  scale_fill_discrete(name = "Sub-subclass")
#ggsave("figures/EC_abund_SEMWISE_subsubclass.pdf", unit = "in", width = 6, height = 4, dpi = 300, device = "pdf")


