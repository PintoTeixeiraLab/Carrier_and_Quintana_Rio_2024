
# Load necessary libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)

# Set working directory
setwd("Your_directory")

# Read in the Seurat object
T4T5_all_stages <- readRDS("T4T5_all_stages.rds")

# Modify metadata to change Time labels
T4T5_all_stages@meta.data <- T4T5_all_stages@meta.data %>%
  mutate(time = as.factor(str_remove(MergeOrigin, "_T45"))) %>%  # Create 'time' column by removing '_T45' from 'MergeOrigin'
  # mutate(time = fct_recode(as.factor(time), "10hapf" = 'P15',"30hapf" = 'P30')) %>% # Optional: Rename time labels (e.g., from P to hours)
  mutate(time = forcats::fct_relevel(time, "Adult", after = Inf)) %>% # Reorder factor levels to put 'Adult' at the end
  mutate(subtype = Idents(T4T5_all_stages)) # Create 'subtype' column from Seurat object identities

# Define list of genes of interest
my_genes = c("side","CG42313","CG34113","CG14372","CG34371","CG34114","CG12950","CG12484","beat-Ia","beat-Ib","beat-Ic","beat-IIa","beat-IIb","beat-IIIa","beat-IIIb","beat-IIIc","beat-IV","beat-Va","beat-Vb","beat-Vc","beat-VI","beat-VII")

# Check if the genes in 'my_genes' are present in the dataset
T4T5_all_stages_rows <- T4T5_all_stages@assays$RNA@data %>% rownames
my_genes %in% T4T5_all_stages_rows 

# Logical vector indicating which genes are present in the dataset
genes_in_data <- my_genes %in% rownames(T4T5_all_stages) 

# Create a new column for combined time and subtype information
T4T5_all_stages$newcol <- str_c(T4T5_all_stages$time, T4T5_all_stages$subtype, sep = ".")

# Calculate average expression for each group defined by the new column
avgexp.obj = AverageExpression(T4T5_all_stages, return.seurat = T, group.by = 'newcol')

# Subset the Seurat object to include only the specified genes
avgexp.obj.sub <- subset(avgexp.obj, features = my_genes[genes_in_data]) 

# Transform the average expression data for plotting
avgexp <- avgexp.obj.sub@assays$RNA@data %>% as.data.frame() %>% 
  tibble::rownames_to_column("Feature") %>%
  pivot_longer(-Feature, names_to = "Identity", values_to = "Expression") %>%
  separate(Identity, into = c("time", "subtype"), sep = "\\.") %>%
  mutate(time = forcats::fct_relevel(time, "Adult", after = Inf)) %>%
  # Rename genes that have CG prefixes
  mutate(Feature = fct_recode(as.factor(Feature), "side-II" = 'CG42313', "side-III" = 'CG34113', "side-IV" = 'CG14372', "side-V" = 'CG34371', "side-VI" = 'CG34114', "side-VII" = 'CG12950', "side-VIII" = 'CG12484')) %>%
  # Specify the order of genes in the plot
  mutate(Feature = fct_relevel(Feature, c("side", "side-II", "side-III", "side-IV", "side-V", "side-VI", "side-VII", "side-VIII", "beat-Ia", "beat-Ib", "beat-Ic", "beat-IIa", "beat-IIb", "beat-IIIa", "beat-IIIb", "beat-IIIc", "beat-IV", "beat-Va", "beat-Vb", "beat-Vc", "beat-VI", "beat-VII")))

# Define subtypes to plot
vector_subtypes_to_plot = c("T4c", "T5c", "T4d", "T5d")

# Create heatmap plot
Heatmap <- avgexp %>%
  dplyr::filter(subtype %in% vector_subtypes_to_plot) %>%
  mutate(subtype = factor(subtype, levels = vector_subtypes_to_plot)) %>%
  ggplot(aes(x = Feature, y = time, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) +
  facet_wrap(~subtype, ncol = 1, strip.position = "right") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  coord_cartesian(expand = F) +
  theme(legend.position = "bottom", legend.justification = "right", legend.box = "horizontal") +
  theme(text = element_text(size = 11)) +
  theme(text = element_text(family = "Helvetica"))
