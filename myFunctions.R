## ---------------------------
##
## Script name: Custom R functions
##
## Purpose of script: Load reusable common reusable scripts
##
## Author: Dan Manter
##
## Date Created: March 6 2022
##
## ## Email: daniel.manter@usda.gov
##
## ---------------------------
##
## Notes: List of functions
##   1) set_panel_size
##   2) tax.clean
##   3) emu_to_phyloseq 
##   4) gibbs_to_phyloseq 
##   5) log2FC_mean
##
## ---------------------------

# function to set the panel size of a multi-panel ggplot object
set_panel_size <- function (p=NULL, g=ggplotGrob(p), 
                            file=NULL, margin = unit(1,"mm"), 
                            width=unit(4, "cm"), height=unit(4, "cm")) {
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  if(getRversion() < "3.3.0") {
    g$widths <- grid:::unit.list(g$widths)
    g$heights <- grid:::unit.list(g$heights)
    g$widths[panel_index_w] <-  rep(list(width),  nw)
    g$heights[panel_index_h] <- rep(list(height), nh)
  } else {
    g$widths[panel_index_w] <-  rep(width,  nw)
    g$heights[panel_index_h] <- rep(height, nh)
  }
  if (!is.null(file)) {
    ggsave(file, g,
           width = convertWidth(sum(g$widths) + margin, unitTo = "in", valueOnly = TRUE),
           height = convertHeight(sum(g$heights) + margin, unitTo = "in", valueOnly = TRUE))
  }
  invisible(g)  
}

# this function changes blanks in a taxonomy table to the lowest common ancestor
tax.clean <- function (tax) {
  for (i in 1:nrow(tax)) {
    if (tax[i,7] == "") {
      if (tax[i,6] == "") {
        if (tax[i,5] == "") {
          if (tax[i,4] == "") {
            if (tax[i,3] == "") {
              if (tax[i,2] == "") {
                if (tax[i,1] == "") {
                  tax[i,7] <- "<NA>"
                } else {
                  tax[i,7] <- paste0("s_of_", tax[i,1])
                }
              } else {
                tax[i,7] <- paste0("s_of_", tax[i,2])
              }
            } else {
              tax[i,7] <- paste0("s_of_", tax[i,3])
            }
          } else {
            tax[i,7] <- paste0("s_of_", tax[i,4])
          }
        } else{
          tax[i,7] <- paste0("s_of_", tax[i,5])
        }
      } else {
        tax[i,7] <- paste0("s_of_", tax[i,6])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,6] == "") {
      if (tax[i,5] == "") {
        if (tax[i,4] == "") {
          if (tax[i,3] == "") {
            if (tax[i,2] == "") {
              if (tax[i,1] == "") {
                tax[i,6] <- "<NA>"
              } else {
                tax[i,6] <- paste0("g_of_", tax[i,1])
              }
            } else {
              tax[i,6] <- paste0("g_of_", tax[i,2])
            }
          } else {
            tax[i,6] <- paste0("g_of_", tax[i,3])
          }
        } else {
          tax[i,6] <- paste0("g_of_", tax[i,4])
        }
      } else {
        tax[i,6] <- paste0("g_of_", tax[i,5])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,5] == "") {
      if (tax[i,4] == "") {
        if (tax[i,3] == "") {
          if (tax[i,2] == "") {
            if (tax[i,1] == "") {
              tax[i,5] <- "<NA>"
            } else {
              tax[i,5] <- paste0("f_of_", tax[i,1])
            }
          } else {
            tax[i,5] <- paste0("f_of_", tax[i,2])
          }
        } else {
          tax[i,5] <- paste0("f_of_", tax[i,3])
        }
      } else {
        tax[i,5] <- paste0("f_of_", tax[i,4])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,4] == "") {
      if (tax[i,3] == "") {
        if (tax[i,2] == "") {
          if (tax[i,1] == "") {
            tax[i,4] <- "<NA>"
          } else {
            tax[i,4] <- paste0("o_of_", tax[i,1])
          }
        } else {
          tax[i,4] <- paste0("o_of_", tax[i,2])
        }
      } else {
        tax[i,4] <- paste0("o_of_", tax[i,3])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,3] == "") {
      if (tax[i,2] == "") {
        if (tax[i,1] == "") {
          tax[i,3] <- "<NA>"
        } else {
          tax[i,3] <- paste0("c_of_", tax[i,1])
        }
      } else {
        tax[i,3] <- paste0("c_of_", tax[i,2])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,2] == "") {
      if (tax[i,1] == "") {
        tax[i,2] <- "<NA>"
      } else {
        tax[i,2] <- paste0("p_of_", tax[i,1])
      }
    }
  }
  return(tax)  
}


# custom function to import emu output files into phyloseq
emu_to_phyloseq <- function (RA_file=NULL, meta_file=NULL, sheet=NULL, 
                             range=NULL, sample_names=NULL, run_name=NULL) {
  df <- read.csv(RA_file, header=T,  
                 row.names='tax_id')
  tax <- df[,2:8]
  tax <- tax.clean(tax)
  TAX <- tax_table(as.matrix(tax))
  
  otu <- df[,9:ncol(df)]
  otu[is.na(otu)] <- 0
  
  meta <- data.frame(read_excel(meta_file, sheet=sheet, range=range))
  common <- Reduce(intersect, list(meta$barcode_key, names(otu)))
  meta <- meta[meta$barcode_key %in% common,]
  SAMP <- sample_data(meta)
  sample_names(SAMP) <- paste(meta$barcode_key, run_name, sep='_')
  
  otu <- otu[,common]
  names(otu) <- paste(names(otu), run_name, sep='_')
  OTU <- otu_table(otu, taxa_are_rows=T)
  
  phy_obj <- phyloseq(SAMP, OTU, TAX )
  
  return(phy_obj)
}

# custom function to import GIBBs output files into phyloseq
gibbs_to_phyloseq <- function (RA_file=NULL, meta_file=NULL, sheet=NULL, 
                             range=NULL, sample_names=NULL, run_name=NULL) {
  otu <- read.csv(RA_file, header=T,  
                  row.names='X')
  
  tax <- read.csv('GIBBs.EC.numbers.csv', row.names='EC')
  TAX <- tax_table(as.matrix(tax))
  
  meta <- data.frame(read_excel(meta_file, sheet=sheet, range=range))
  common <- Reduce(intersect, list(meta$barcode_key, names(otu)))
  meta <- meta[meta$barcode_key %in% common,]
  SAMP <- sample_data(meta)
  sample_names(SAMP) <- paste(meta$barcode_key, run_name, sep='_')
  
  otu <- otu[,common]
  names(otu) <- paste(names(otu), run_name, sep='_')
  OTU <- otu_table(otu, taxa_are_rows=T)
  
  phy_obj <- phyloseq(SAMP, OTU, TAX )
  
  return(phy_obj)
}

# custom function for metacoder to calculate log2FC of the means
log2FC_mean <- function(abund_1, abund_2) {
  log_ratio <- log2(median(abund_1)) - log2(median(abund_2))
  if (is.nan(log_ratio)) {
    log_ratio <- 0
  }
  list(log2_FC = log_ratio,
       median_diff = median(abund_1) - median(abund_2),
       mean_diff = mean(abund_1) - mean(abund_2),
       wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
}

# custom ggplot theme
theme_dan <- function (base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # Titles
      plot.title = element_text(size = rel(1), face = "bold", margin = margin(0,0,5,0), hjust = 0),
      plot.subtitle = element_text(size = rel(0.85), face = "bold", margin = margin(0,0,5,0), hjust = 0),
      plot.caption = element_text(size = rel(0.70), hjust=0),
      
      # panel
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      
      # Axes
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.line = element_line(color = "black", arrow = arrow(length = unit(0.3, "lines"), type = "closed")),
      
      # Legend
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.70), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(1.5, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      
      # Facets
      strip.background = element_rect(fill = "#17252D", color = "#17252D"),
      strip.text = element_text(size = rel(0.85), face = "bold", color = "white", margin = margin(5,0,5,0))
    )
}