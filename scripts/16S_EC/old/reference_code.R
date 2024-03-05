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

### Feature importance distributions of models run using 0.8 modules
# Using indicators where models with all enzymes performed better
# # import gibbs
# gibbs <- read.csv('picrust2_files/GIBBs.EC.numbers.csv')
# gibbs <- gibbs %>%
#   filter(Missing != 'x')
# 
# # import modules
# modules <- read.csv("scripts/16S_EC/res_node_table_0.8.csv")
# modules <- modules[,2:12]
# colnames(modules)[1] <- "EC"
# 
# # assign gibbs ECs modules
# gibbs_mods <- gibbs %>% 
#   left_join(modules, by = "EC") %>% 
#   drop_na(module)
# gibbs_mods$module <- gsub("M","MOD", gibbs_mods$module)
# gibbs_mods$gibbs <- "G"
# gibbs_mods <- gibbs_mods[,c(14,19)]
# colnames(gibbs_mods)[1] <- "Predictor"
# 
# myVars <- c('ace', 'activeC', 'ph', 'p', 'k','fe')
# 
# for (myVar in myVars) {
#     filenames <- list.files(paste0("machine_learning/16S_EC/",myVar,"_model_results_modules/"), pattern = paste0(myVar,"_modules_varimp*"), full.names = TRUE)
#   varimps <- lapply(filenames, read.csv, header = FALSE)
#   varimps_full <- varimps %>% 
#     reduce(full_join)
#   
#   # changing myVar format for cleaner plot titles
#   if (myVar == 'ace'){
#     myVar <- "ACE"
#   }
#   if (myVar == 'activeC'){
#     myVar <- "ActiveC"
#   }
#   if (myVar == 'ph'){
#     myVar <- "pH"
#   }
#   if (myVar == 'p'){
#     myVar <- "P"
#   }
#   if (myVar == 'k'){
#     myVar <- "K"
#   }
#   if (myVar == 'fe'){
#     myVar <- "Fe"
#   }
# 
#   names(varimps_full) <- c('Predictor', 'VIP', 'Run')
#   
#   # left join with gibbs mods
#   varimps_full_gibbs <- varimps_full %>% 
#     left_join(gibbs_mods, by = "Predictor") 
#   
#   varimps_full_gibbs$gibbs <- varimps_full_gibbs$gibbs %>% 
#    replace_na("NG")
#   
#   # filter out poor indicators
#   temp <- varimps_full_gibbs %>% group_by(Predictor) %>%
#     summarise(Var = median(VIP)) %>%
#     top_n(10, Var) 
#   
#   gDF <- varimps_full_gibbs %>%
#     filter(Predictor %in% temp$Predictor) %>% 
#     mutate(Predictor=fct_reorder(Predictor, VIP, .fun=median, .desc=TRUE))
#   
#   colors <- paletteer_d(`"ggsci::category20_d3"`)
#   
#   means <- aggregate(VIP ~ Predictor, gDF, median)
#   means$Category <- 'Other'
#   means$VIP <- round(means$VIP, 3)
#   
#   # formatting for gibbs label
#   if (myVar == 'ACE'){
#     gDF$y <- 0.5
#   }
#   
#   if (myVar == 'ActiveC'){
#     gDF$y <- 0.35
#   }
#   
#   if (myVar == 'pH'){
#     gDF$y <- 0.35
#   }
#   
#   if (myVar == 'P'){
#     gDF$y <- 0.35
#   }
#   
#   if (myVar == 'K'){
#     gDF$y <- 0.62
#   }
#   
#   if (myVar == 'Fe'){
#     gDF$y <- 0.52
#   }
#   
#   p <- ggplot(gDF, aes(x=Predictor, y=VIP)) +
#     #scale_fill_manual(values=colors, na.value="white") + # maybe use this to fill GIBBs enzymes?
#     geom_boxplot(alpha = 1, fatten = 1) + 
#     theme_bw() +
#     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), plot.title = element_text(size = 12)) + 
#     labs(x="", y="", title = myVar) +
#     geom_text(inherit.aes=F, data=gDF, aes(x=Predictor, y=y, label=gibbs), size = 4)
#   p
#   assign(paste0(myVar,"_mods"),p)
#   #ggsave(paste0("figures/varimps/", myVar, "_module_varimp.png"), device='png', dpi=600, height=6, width=8)
# }
# 
# p_mods_merged <- ggarrange(ACE_mods, ActiveC_mods, Fe_mods, K_mods, P_mods, pH_mods,
#                             nrow = 2, ncol = 3)
# p_mods_merged
# 
# annotate_figure(p_mods_merged, left = text_grob("Variable Importance", rot = 90, vjust = 1, size = 14),
#                 bottom = text_grob("Model Predictors", size = 14))
# 
# ggsave("figures/mods_varimps.png", device='png', dpi=600, height=8, width=12)

### Feature importance distributions of GIBBs models
# separating by positively and negatively correlated features
# gibbs <- read.csv('picrust2_files/GIBBs.EC.numbers.csv')
# gibbs <- gibbs %>%
#   filter(Missing != 'x')
# 
# myVars <- c('SOM', 'resp', 'agg_stab', 'water_cap', 'mg', 'mn', 'zn', 'ace', 'activeC', 'ph', 'p', 'k','fe')
# 
# for (myVar in myVars) {
#      filenames <- list.files(paste0("machine_learning/16S_EC/",myVar,"_model_results_RA_corr_GIBBs/"), pattern = paste0(myVar,"_RA_corr_GIBBs_varimp*"), full.names = TRUE)
#   varimps <- lapply(filenames, read.csv, header = FALSE)
#   varimps_full <- varimps %>% 
#     reduce(full_join)
#   
#   names(varimps_full) <- c('Predictor', 'VIP', 'Run')
#   varimps_full$Predictor <- gsub("EC_","EC:", varimps_full$Predictor)
#   varimps_full$Predictor <- gsub("_","\\.", varimps_full$Predictor)
# 
#  gDF <- varimps_full %>%
#     left_join(gibbs, by=c("Predictor"="EC")) %>% 
#     mutate(Predictor=fct_reorder(Predictor, VIP, .fun=median, .desc=TRUE))
#  
#  # calculate positive and negative correlations between enzymes and indicators
#  ml_EC_RA_corr <- readRDS('machine_learning/EC_RA_corr/ml_EC_RA_corr.RDS')
#  
#  ml_indic <- ml_EC_RA_corr %>% 
#    select(gDF$Predictor, paste0(myVar)) %>% 
#    subset(select=-c(ClimateZ, clay, DNA))
#  
#  ml_indic <- as.matrix(ml_indic)
#  ec.mat <- rcorr(ml_indic)
#  ec.r <- ec.mat$r
#  ec.indic <- ec.r %>% 
#    subset(select = paste0(myVar))
#  ec.indic <- data.frame(ec.indic)
# 
#  # get positive and negative correlation values
#  ec.indic.pos <- ec.indic %>% 
#    filter(ec.indic[,1] > 0)
#  ec.indic.neg <- ec.indic %>% 
#    filter(ec.indic[,1] < 0)
#  ec.indic.pos <- rownames_to_column(ec.indic.pos, var = "Predictor")
#  ec.indic.neg <- rownames_to_column(ec.indic.neg, var = "Predictor")
#  ec.indic.pos <- ec.indic.pos[ec.indic.pos$Predictor != myVar,]
#  ec.indic.neg <- ec.indic.neg[ec.indic.neg$Predictor != myVar,]
# 
#   # filter out poor indicators in pos and neg correlations
#   # add this if I want to include climate, clay, and DNA:
#   # | Predictor %in% c('clay','ClimateZ','DNA')
#  gDF.pos <- gDF %>% 
#    filter(Predictor %in% ec.indic.pos$Predictor) 
#  gDF.neg <- gDF %>% 
#    filter(Predictor %in% ec.indic.neg$Predictor) 
#  
#   temp.pos <- gDF.pos %>% group_by(Predictor) %>%
#     summarise(Var = median(VIP)) %>%
#     top_n(10, Var) 
#   temp.neg <- gDF.neg %>% group_by(Predictor) %>%
#     summarise(Var = median(VIP)) %>%
#     top_n(10, Var) 
#   
#   gDF.pos <- gDF.pos %>%
#     filter(Predictor %in% temp.pos$Predictor)
#   gDF.neg <- gDF.neg %>%
#     filter(Predictor %in% temp.neg$Predictor)
#   
#   colors <- c("Biocontrol" = "#1F77B4FF",
#               "C cycling" = "#FF7F0EFF",
#               "C/N cycling" = "#2CA02CFF",
#               "N cycling" = "#D62728FF",
#               "P cycling" = "#9467BDFF",
#               "S cycling" = "#8C564BFF",
#               "Stress" = "#E377C2FF",
#               "Siderophore" = "#7F7F7FFF",
#               "NA" = "white")
#   
#   # changing myVar format for cleaner plot titles
#   if (myVar == 'resp'){
#     myVar <- "Resp"
#   }
#   if (myVar == 'agg_stab'){
#     myVar <- "AggStab"
#   }
#   if (myVar == 'water_cap'){
#     myVar <- "WaterCap"
#   }
#   if (myVar == 'mg'){
#     myVar <- "Mg"
#   }
#   if (myVar == 'mn'){
#     myVar <- "Mn"
#   }
#   if (myVar == 'zn'){
#     myVar <- "Zn"
#   }
#   if (myVar == 'ace'){
#     myVar <- "ACE"
#   }
#   if (myVar == 'activeC'){
#     myVar <- "ActiveC"
#   }
#   if (myVar == 'ph'){
#     myVar <- "pH"
#   }
#   if (myVar == 'p'){
#     myVar <- "P"
#   }
#   if (myVar == 'k'){
#     myVar <- "K"
#   }
#   if (myVar == 'fe'){
#     myVar <- "Fe"
#   }
#   
#   # positive correlation graph
#   means <- aggregate(VIP ~ Predictor, gDF.pos, median)
#   means$Category <- 'Other'
#   means$VIP <- round(means$VIP, 3)
#   p.pos <- ggplot(gDF.pos, aes(x=Predictor, y=VIP, fill=Category)) +
#     scale_fill_manual(values=colors, na.value="white") +
#     geom_boxplot(alpha = 1, fatten = 1) + 
#     theme_bw() +
#     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
#           legend.position = "none", plot.title = element_text(size = 12)) +
#     labs(x="", y="", title = myVar)
#   p.pos
#   assign(paste0(myVar,"_gibbs_pos"),p.pos)
#   ggsave(paste0("figures/varimps/", myVar, "_RA_corr_GIBBs_varimp.pos.png"), device='png', dpi=600, height=8, width=8)
#   
#   # negative correlation graph
#   means <- aggregate(VIP ~ Predictor, gDF.neg, median)
#   means$Category <- 'Other'
#   means$VIP <- round(means$VIP, 3)
#   p.neg <- ggplot(gDF.neg, aes(x=Predictor, y=VIP, fill=Category)) +
#     scale_fill_manual(values=colors, na.value="white") +
#     geom_boxplot(alpha = 1, fatten = 1) + 
#     theme_bw() +
#     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
#           legend.position = "none", plot.title = element_text(size = 12)) +
#     labs(x="", y="", title = myVar)
#   p.neg
#   assign(paste0(myVar,"_gibbs_neg"),p.neg)
#   ggsave(paste0("figures/varimps/", myVar, "_RA_corr_GIBBs_varimp.neg.png"), device='png', dpi=600, height=8, width=8)
# }
# 
# # create legend of all GIBBs colors to add to merged plot
# p_legend <- ggplot(gibbs, aes(x = EC, fill = Category)) +
#   scale_fill_manual(values = colors) +
#   geom_boxplot()
# p_legend
# p_legend2 <- get_legend(p_legend)
# p_legend3 <- as_ggplot(p_legend2)
# p_legend3
# 
# # positive correlations graph
# p_gibbs_merged_pos <- ggarrange(ACE_gibbs_pos, ActiveC_gibbs_pos, AggStab_gibbs_pos, 
#                                 Fe_gibbs_pos, K_gibbs_pos, Mg_gibbs_pos, Mn_gibbs_pos, 
#                                 P_gibbs_pos, pH_gibbs_pos, Resp_gibbs_pos, 
#                                 SOM_gibbs_pos, WaterCap_gibbs_pos, Zn_gibbs_pos,
#                                 p_legend3, nrow = 2, ncol = 7)
# p_gibbs_merged_pos
# 
# annotate_figure(p_gibbs_merged_pos, left = text_grob("Variable Importance", rot = 90, vjust = 1, size = 14),
#                 bottom = text_grob("Positively Correlated Model Predictors", size = 14),
#                 top = text_grob("GIBBs Enzymes", size = 14))
# 
# ggsave("figures/gibbs_varimps_pos.png", device='png', dpi=600, height=6, width=12)
# 
# # negative correlations graph
# p_gibbs_merged_neg <- ggarrange(ACE_gibbs_neg, ActiveC_gibbs_neg, AggStab_gibbs_neg, 
#                                 Fe_gibbs_neg, K_gibbs_neg, Mg_gibbs_neg, Mn_gibbs_neg, 
#                                 P_gibbs_neg, pH_gibbs_neg, Resp_gibbs_neg, 
#                                 SOM_gibbs_neg, WaterCap_gibbs_neg, Zn_gibbs_neg,
#                                 p_legend3, nrow = 2, ncol = 7)
# p_gibbs_merged_neg
# 
# annotate_figure(p_gibbs_merged_neg, left = text_grob("Variable Importance", rot = 90, vjust = 1, size = 14),
#                 bottom = text_grob("Negatively Correlated Model Predictors", size = 14),
#                 top = text_grob("GIBBs Enzymes", size = 14))
# 
# ggsave("figures/gibbs_varimps_neg.png", device='png', dpi=600, height=6, width=12)


# create list of ECs counts across all data - which ECs are highly functionally redundant?
# hops2018_count <- table(unlist(lapply(all_merged$Hops2018_ECs, unique)))
# hops2018_count <- data.frame(hops2018_count)
# hopsARS_count <- table(unlist(lapply(all_merged$hopsARS_ECs, unique)))
# hopsARS_count <- data.frame(hopsARS_count)
# NRCS1_count <- table(unlist(lapply(all_merged$NRCS1_ECs, unique)))
# NRCS1_count <- data.frame(NRCS1_count)
# NRCS2_count <- table(unlist(lapply(all_merged$NRCS2_ECs, unique)))
# NRCS2_count <- data.frame(NRCS2_count)
# rangeland_count <- table(unlist(lapply(all_merged$rangeland_ECs, unique)))
# rangeland_count <- data.frame(rangeland_count)
# 
# EC_counts <- hops2018_count %>% 
#   full_join(hopsARS_count, by = "Var1") %>% 
#   full_join(NRCS1_count, by = "Var1") %>% 
#   full_join(NRCS2_count, by = "Var1") %>% 
#   full_join(rangeland_count, by = "Var1")
# EC_counts <- EC_counts %>% 
#   mutate(Freq_total = select(., 2:6) %>% rowSums(na.rm = TRUE))
# EC_counts <- EC_counts[,c(1,7)]
# colnames(EC_counts)[1] <- "EC"
# saveRDS(EC_counts, "picrust2_files/EC_freq.RDS")

### Functional redundancy in taxonomic models - old code
### import important ECs
top_ECs <- read.csv("machine_learning/16S_EC/top_ECs_HLD.csv", header = FALSE)
colnames(top_ECs) <- c("EC","mean_VIP","indicator")
top_ECs$EC <- gsub("EC:", "EC.", top_ECs$EC) 
top_ECs_unique <- top_ECs %>% # create list without duplicates
  distinct(EC)

### import important taxa
#need to import varimps first
myVars <- c('SOM', 'resp', 'agg_stab', 'water_cap', 'mg', 'mn', 'zn', 'ace', 'activeC', 'ph', 'p', 'k','fe')

for (myVar in myVars) {
  filenames <- list.files(paste0("machine_learning/16S_TAX/",myVar,"_model_results_RA_TAX/"), pattern = paste0(myVar,"_RA_TAX_varimp*"), full.names = TRUE)
  varimps <- lapply(filenames, read.csv, header = FALSE)
  varimps_full <- varimps %>%
    reduce(full_join)
  
  names(varimps_full) <- c('Predictor', 'VIP', 'Run')
  
  gDF <- varimps_full %>%
    mutate(Predictor=fct_reorder(Predictor, VIP, .fun=median, .desc=TRUE))
  
  gDF <- gDF %>%
    filter(gDF$Predictor != "ClimateZ" & gDF$Predictor != "clay" & gDF$Predictor != "DNA")
  
  temp <- gDF %>% group_by(Predictor) %>%
    summarise(Var = median(VIP)) %>%
    top_n(10, Var)
  
  gDF <- gDF %>%
    filter(Predictor %in% temp$Predictor)
  
  # write important variables to a csv file
  gDF.grouped <- gDF %>%
    group_by(Predictor) %>%
    summarise(mean_VIP = mean(VIP)) %>%
    mutate(indicator = myVar)
  
  write.table(gDF.grouped, "machine_learning/16S_TAX/top_TAX_HLD.csv", sep = ",",
              append = TRUE, row.names = FALSE, col.names = FALSE)
}

# create list without duplicates
top_TAX <- read.csv("machine_learning/16S_TAX/top_TAX_HLD.csv", header = FALSE)
colnames(top_TAX) <- c("species","mean_VIP","indicator")
top_TAX_unique <- top_TAX %>%
  distinct(species)

### import picrust2 copy number data for each data set
# hops 2018
hops2018 <- read_tsv("picrust2_files/Hops.2018/EC_predicted_and_nsti.tsv")
hops2018$sequence <- gsub("OTU", "", hops2018$sequence)  # get rid of "OTU" in taxon column
colnames(hops2018)[1] <- "OTU_ID"

# associate OTU IDs with species names
hops2018_otuid <- read.delim("16S_dada2_output/Hops.2018.taxonomy.otuid.csv",
                             sep = ",")

hops2018_otuid$OTU_ID <- as.character(hops2018_otuid$OTU_ID)
hops2018_sp <- hops2018 %>%
  left_join(hops2018_otuid, by = "OTU_ID")
hops2018_sp <- hops2018_sp %>%
  select(starts_with("EC.") | starts_with("species"))

# get rid of unclassified species
hops2018_sp <- hops2018_sp[!hops2018_sp$species == "",]

# split taxonomy to get just species
hops2018_sp2 <- separate_wider_delim(hops2018_sp, cols = species, delim = ";",
                                     names = c("kingdom", "phylum", "class",
                                               "order", "family", "genus",
                                               "species"))

# filter to only species and replace spaces with _ and . with _
# also format some specific species so they are identical to that in top_TAX_unique list (for filtering later)
hops2018_sp3 <- hops2018_sp2 %>%
  select(starts_with("EC.") | starts_with("species"))
hops2018_sp3$species <- gsub(" ", "_", hops2018_sp3$species)
hops2018_sp3$species <- gsub("\\.", "_", hops2018_sp3$species)
hops2018_sp3$species <- gsub("Bacillus_sp__X1\\(2014\\)", "Bacillus_sp__X1_2014_", hops2018_sp3$species)
hops2018_sp3$species <- gsub("Bacillus_sp__HBCD-sjtu", "Bacillus_sp__HBCD_sjtu", hops2018_sp3$species)
hops2018_sp3$species <- gsub("Burkholderiales_bacterium_GJ-E10", "Burkholderiales_bacterium_GJ_E10", hops2018_sp3$species)
hops2018_sp3$species <- gsub("Rhodoplanes_sp__Z2-YC6860", "Rhodoplanes_sp__Z2_YC6860", hops2018_sp3$species)

# filter by important ECs and important species
hops2018_sp4 <- hops2018_sp3 %>%
  group_by(species) %>%
  summarise(across(everything(), mean))
hops2018_sp4 <- column_to_rownames(hops2018_sp4, var = "species")
hops2018_sp5 <- hops2018_sp4[,names(hops2018_sp4) %in% top_ECs_unique$EC]
hops2018_sp5 <- rownames_to_column(hops2018_sp5, "species")
hops2018_sp6 <- hops2018_sp5[hops2018_sp5$species %in% top_TAX_unique$species,]

# convert to long format for merging with other data sets
hops2018_sp6_long <- hops2018_sp6 %>%
  pivot_longer(names_to = "EC",
               values_to = "ECs_per_species",
               cols = starts_with("EC."))

## hops ARS
hopsARS <- read.delim("picrust2_files/Hops.ARS/EC_predicted_and_nsti.tsv.gz", sep = "\t")
hopsARS$sequence <- gsub("OTU", "", hopsARS$sequence)
colnames(hopsARS)[1] <- "OTU_ID"
hopsARS_otuid <- read.delim("16S_dada2_output/Hops.ARS.taxonomy.otuid.csv",
                            sep = ",")

hopsARS_otuid$OTU_ID <- as.character(hopsARS_otuid$OTU_ID)
hopsARS_sp <- hopsARS %>%
  left_join(hopsARS_otuid, by = "OTU_ID")
hopsARS_sp <- hopsARS_sp %>%
  select(starts_with("EC.") | starts_with("species"))
hopsARS_sp <- hopsARS_sp[!hopsARS_sp$species == "",]
hopsARS_sp2 <- separate_wider_delim(hopsARS_sp, cols = species, delim = ";",
                                    names = c("kingdom", "phylum", "class",
                                              "order", "family", "genus",
                                              "species"))
hopsARS_sp3 <- hopsARS_sp2 %>%
  select(starts_with("EC.") | starts_with("species"))
hopsARS_sp3$species <- gsub(" ", "_", hopsARS_sp3$species)
hopsARS_sp3$species <- gsub("\\.", "_", hopsARS_sp3$species)
hopsARS_sp3$species <- gsub("Bacillus_sp__X1\\(2014\\)", "Bacillus_sp__X1_2014_", hopsARS_sp3$species)
hopsARS_sp3$species <- gsub("Bacillus_sp__HBCD-sjtu", "Bacillus_sp__HBCD_sjtu", hopsARS_sp3$species)
hopsARS_sp3$species <- gsub("Burkholderiales_bacterium_GJ-E10", "Burkholderiales_bacterium_GJ_E10", hopsARS_sp3$species)
hopsARS_sp3$species <- gsub("Rhodoplanes_sp__Z2-YC6860", "Rhodoplanes_sp__Z2_YC6860", hopsARS_sp3$species)
hopsARS_sp4 <- hopsARS_sp3 %>%
  group_by(species) %>%
  summarise(across(everything(), mean))
hopsARS_sp4 <- column_to_rownames(hopsARS_sp4, var = "species")
hopsARS_sp5 <- hopsARS_sp4[,names(hopsARS_sp4) %in% top_ECs_unique$EC]
hopsARS_sp5 <- rownames_to_column(hopsARS_sp5, "species")
hopsARS_sp6 <- hopsARS_sp5[hopsARS_sp5$species %in% top_TAX_unique$species,]
hopsARS_sp6_long <- hopsARS_sp6 %>%
  pivot_longer(names_to = "EC",
               values_to = "ECs_per_species",
               cols = starts_with("EC."))

## NRCS
NRCS <- read.delim("picrust2_files/NRCS/EC_predicted_and_nsti.tsv.gz", sep = "\t")
NRCS$sequence <- gsub("OTU", "", NRCS$sequence)
colnames(NRCS)[1] <- "OTU_ID"
NRCS_otuid <- read.delim("16S_dada2_output/NRCS.taxonomy.otuid.csv",
                         sep = ",")

NRCS_otuid$OTU_ID <- as.character(NRCS_otuid$OTU_ID)
NRCS_sp <- NRCS %>%
  left_join(NRCS_otuid, by = "OTU_ID")
NRCS_sp <- NRCS_sp %>%
  select(starts_with("EC.") | starts_with("species"))
NRCS_sp <- NRCS_sp[!NRCS_sp$species == "",]
NRCS_sp2 <- separate_wider_delim(NRCS_sp, cols = species, delim = ";",
                                 names = c("kingdom", "phylum", "class",
                                           "order", "family", "genus",
                                           "species"))
NRCS_sp3 <- NRCS_sp2 %>%
  select(starts_with("EC.") | starts_with("species"))
NRCS_sp3$species <- gsub(" ", "_", NRCS_sp3$species)
NRCS_sp3$species <- gsub("\\.", "_", NRCS_sp3$species)
NRCS_sp3$species <- gsub("Bacillus_sp__X1\\(2014\\)", "Bacillus_sp__X1_2014_", NRCS_sp3$species)
NRCS_sp3$species <- gsub("Bacillus_sp__HBCD-sjtu", "Bacillus_sp__HBCD_sjtu", NRCS_sp3$species)
NRCS_sp3$species <- gsub("Burkholderiales_bacterium_GJ-E10", "Burkholderiales_bacterium_GJ_E10", NRCS_sp3$species)
NRCS_sp3$species <- gsub("Rhodoplanes_sp__Z2-YC6860", "Rhodoplanes_sp__Z2_YC6860", NRCS_sp3$species)
NRCS_sp4 <- NRCS_sp3 %>%
  group_by(species) %>%
  summarise(across(everything(), mean))
NRCS_sp4 <- column_to_rownames(NRCS_sp4, var = "species")
NRCS_sp5 <- NRCS_sp4[,names(NRCS_sp4) %in% top_ECs_unique$EC]
NRCS_sp5 <- rownames_to_column(NRCS_sp5, "species")
NRCS_sp6 <- NRCS_sp5[NRCS_sp5$species %in% top_TAX_unique$species,]
NRCS_sp6_long <- NRCS_sp6 %>%
  pivot_longer(names_to = "EC",
               values_to = "ECs_per_species",
               cols = starts_with("EC."))

## rangeland
rangeland <- read.delim("picrust2_files/Rangeland/EC_predicted_and_nsti.tsv.gz", sep = "\t")
rangeland$sequence <- gsub("OTU", "", rangeland$sequence)  # get rid of "OTU" in taxon column
colnames(rangeland)[1] <- "OTU_ID"
rangeland_otuid <- read.delim("16S_dada2_output/Rangeland.taxonomy.otuid.csv",
                              sep = ",")

rangeland_otuid$OTU_ID <- as.character(rangeland_otuid$OTU_ID)
rangeland_sp <- rangeland %>%
  left_join(rangeland_otuid, by = "OTU_ID")
rangeland_sp <- rangeland_sp %>%
  select(starts_with("EC.") | starts_with("species"))
rangeland_sp <- rangeland_sp[!rangeland_sp$species == "",]
rangeland_sp2 <- separate_wider_delim(rangeland_sp, cols = species, delim = ";",
                                      names = c("kingdom", "phylum", "class",
                                                "order", "family", "genus",
                                                "species"))
rangeland_sp3 <- rangeland_sp2 %>%
  select(starts_with("EC.") | starts_with("species"))
rangeland_sp3$species <- gsub(" ", "_", rangeland_sp3$species)
rangeland_sp3$species <- gsub("\\.", "_", rangeland_sp3$species)
rangeland_sp3$species <- gsub("Bacillus_sp__X1\\(2014\\)", "Bacillus_sp__X1_2014_", rangeland_sp3$species)
rangeland_sp3$species <- gsub("Bacillus_sp__HBCD-sjtu", "Bacillus_sp__HBCD_sjtu", rangeland_sp3$species)
rangeland_sp3$species <- gsub("Burkholderiales_bacterium_GJ-E10", "Burkholderiales_bacterium_GJ_E10", rangeland_sp3$species)
rangeland_sp3$species <- gsub("Rhodoplanes_sp__Z2-YC6860", "Rhodoplanes_sp__Z2_YC6860", rangeland_sp3$species)
rangeland_sp4 <- rangeland_sp3 %>%
  group_by(species) %>%
  summarise(across(everything(), mean))
rangeland_sp4 <- column_to_rownames(rangeland_sp4, var = "species")
rangeland_sp5 <- rangeland_sp4[,names(rangeland_sp4) %in% top_ECs_unique$EC]
rangeland_sp5 <- rownames_to_column(rangeland_sp5, "species")
rangeland_sp6 <- rangeland_sp5[rangeland_sp5$species %in% top_TAX_unique$species,]
rangeland_sp6_long <- rangeland_sp6 %>%
  pivot_longer(names_to = "EC",
               values_to = "ECs_per_species",
               cols = starts_with("EC."))

### merge, get rid of duplicates, and save
all_merged <- rbind(hops2018_sp6_long, hopsARS_sp6_long, NRCS_sp6_long, rangeland_sp6_long)
all_merged_unique <- all_merged %>%
  distinct()

# check that all_merged_unique includes all of the top species
length(unique(all_merged_unique$species))
length(top_TAX_unique$species)
# both include 93 species

saveRDS(all_merged_unique, "picrust2_files/Imp_EC_copies_per_imp_species_HLD.RDS")
all_merged_unique <- readRDS("picrust2_files/Imp_EC_copies_per_imp_species_HLD.RDS")

# shorten length of Pseudarthrobacter_sp__NIBRBAC000502772 so the plots aren't so squished
# write in figure description somewhere that this was changed
all_merged_unique$species <- gsub("Pseudarthrobacter_sp__NIBRBAC000502772", "Pseudarthrobacter_sp", all_merged_unique$species)
top_TAX$species <- gsub("Pseudarthrobacter_sp__NIBRBAC000502772", "Pseudarthrobacter_sp", top_TAX$species)


### graphs
myVars <- c("ace", "activeC", "agg_stab", "fe", "k", "mg", "mn", "p", "ph", "resp", "SOM", "water_cap", "zn")
for (myVar in myVars) {
  top_ECs_indic <- top_ECs %>% # don't need to use unique list since it's per indicator (and I can keep mean_VIP)
    filter(indicator == myVar)
  top_TAX_indic <- top_TAX %>% 
    filter(indicator == myVar)
  
  # filter all_merged down to just top taxa for each indicator
  all_merged_imp <- all_merged_unique %>% 
    filter(species %in% top_TAX_indic$species) 
  
  # filter by ECs in top ECs list
  all_merged_imp2 <- all_merged_imp %>% 
    filter(EC %in% top_ECs_indic$EC)
  
  # calculate number of species each EC occurs in
  all_merged_imp3 <- all_merged_imp2 %>% 
    group_by(species, EC) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    group_by(EC) %>% 
    summarise(species_freq = n())
  
  # sort ECs based on importance
  all_merged_imp3 <- all_merged_imp3 %>% 
    left_join(top_ECs_indic, by = "EC")
  all_merged_imp3 <- all_merged_imp3[order(all_merged_imp3$mean_VIP, decreasing = TRUE),]
  all_merged_imp3$EC <- reorder(all_merged_imp3$EC, all_merged_imp3$mean_VIP)
  
  # changing myVar format for cleaner plot titles
  if (myVar == 'resp'){
    myVar <- "Resp"
  }
  if (myVar == 'agg_stab'){
    myVar <- "AggStab"
  }
  if (myVar == 'water_cap'){
    myVar <- "WaterCap"
  }
  if (myVar == 'mg'){
    myVar <- "Mg"
  }
  if (myVar == 'mn'){
    myVar <- "Mn"
  }
  if (myVar == 'zn'){
    myVar <- "Zn"
  }
  if (myVar == 'ace'){
    myVar <- "ACE"
  }
  if (myVar == 'activeC'){
    myVar <- "ActiveC"
  }
  if (myVar == 'ph'){
    myVar <- "pH"
  }
  if (myVar == 'p'){
    myVar <- "P"
  }
  if (myVar == 'k'){
    myVar <- "K"
  }
  if (myVar == 'fe'){
    myVar <- "Fe"
  }
  
  # make bar chart of number of species for each EC
  # ECs increase in importance from left to right
  # just counts the important species now. Do we want to measure the total number of enzymes in each species? Unsure
  p <- ggplot(all_merged_imp3, aes(x = EC, y = species_freq)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.6)) +
    labs(x = "", y = "", title = myVar)
  p
  ggsave(paste0("../figures/species_per_enzyme/species_freq_per_ez_",myVar,".png"), device='png', dpi=600, height=4, width=6)
  assign(paste0("bar_",myVar),p)
  
  # reorder ECs by importance
  all_merged_imp4 <- all_merged_imp2 %>% 
    left_join(top_ECs_indic, by = "EC")
  all_merged_imp4 <- all_merged_imp4[order(all_merged_imp4$mean_VIP),]
  all_merged_imp4$EC <- reorder(all_merged_imp4$EC, all_merged_imp4$mean_VIP)
  
  # reorder species by importance
  all_merged_imp5 <- all_merged_imp4 %>% 
    left_join(top_TAX_indic, by = "species")
  colnames(all_merged_imp5)[4] <- "mean_VIP_EC"
  colnames(all_merged_imp5)[6] <- "mean_VIP_TAX"
  all_merged_imp5$species <- reorder(all_merged_imp5$species, all_merged_imp5$mean_VIP_TAX)
  
  # create heatmaps
  all_merged_imp5[all_merged_imp5$ECs_per_species == 0,"ECs_per_species"] <- NA
  all_merged_imp5 <- all_merged_imp5[all_merged_imp5$indicator.x == tolower(myVar),]
  p2 <- ggplot(all_merged_imp5, aes(x = EC, y = species, fill = ECs_per_species)) +
    geom_tile(color='grey20') +
    geom_text(aes(label=round(ECs_per_species,3))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          legend.text = element_text(size = 8)) +
    scale_fill_binned(low = "white", high = "red", limits = c(0,10),
                      breaks = c(1,2,3,4,5,6,7,8,9,10)) +
    labs(x = "", y = "", title = myVar, fill = "Enzyme Copies")
  ggsave(paste0("../figures/species_per_enzyme/species_freq_per_ez_heatmap_",myVar,".png"), device='png', dpi=600, height=4, width=6)
  assign(paste0("heatmap_",myVar),p2)
}

# create merged plot of species per EC
p_bar_all <- ggarrange(bar_ACE, bar_ActiveC, bar_AggStab,bar_Fe, bar_K, bar_Mg, 
                       bar_Mn,bar_P, bar_pH, bar_Resp, bar_SOM,bar_WaterCap, 
                       bar_Zn, ncol = 4, nrow = 4)
p_bar_all

annotate_figure(p_bar_all, left = text_grob("Species Frequency", rot = 90, vjust = 1, size = 14),
                bottom = text_grob("Top Enzymes", size = 14))
ggsave("../figures/species_per_enzyme/species_freq_per_ez.png", device='png', dpi=600, height=12, width=10)

# create merged plot of heatmaps
p_heatmap_all <- ggarrange(heatmap_ACE, heatmap_ActiveC, heatmap_AggStab,heatmap_Fe, heatmap_K, heatmap_Mg, 
                           heatmap_Mn,heatmap_P, heatmap_pH, heatmap_Resp, heatmap_SOM,heatmap_WaterCap, 
                           heatmap_Zn, ncol = 4, nrow = 4, common.legend = TRUE)
p_heatmap_all

annotate_figure(p_heatmap_all, left = text_grob("Top Species", rot = 90, vjust = 1, size = 14),
                bottom = text_grob("Top Enzymes", size = 14))
ggsave("../figures/species_per_enzyme/ez_copies_per_species_heatmap.png", device='png', dpi=600, height=14, width=14)





