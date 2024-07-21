# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(tidyverse)

# Set working directory
setwd("~/Google Drive/My Drive/Documents Gdrive/T4T5 Sequencing/Zipursky FULL dataset_V2")

# Load dataset from Yoo et al., 2023
dataset_V2 <- read_rds("data_V1.1a.rds")

# Define gene list
my_genes <- c("side", "side-II", "side-III", "side-IV", "side-V", "side-VI", "side-VII", "side-VIII",
              "beat-Ia", "beat-Ib", "beat-Ic", "beat-IIa", "beat-IIb", "beat-IIIa", "beat-IIIb",
              "beat-IIIc", "beat-IV", "beat-Va", "beat-Vb", "beat-Vc", "beat-VI", "beat-VII")

# Create a new column combining time and subtype2
dataset_V2$newcol <- str_c(dataset_V2$time, dataset_V2$subtype2, sep = ".")

# Compute average expression and create Seurat object
avgexp.obj <- AverageExpression(dataset_V2, return.seurat = TRUE, group.by = 'newcol')

# Subset Seurat object using genes_in_data
avgexp.obj.sub <- subset(avgexp.obj, features = my_genes[genes_in_data])

# Prepare data for plotting
avgexp <- avgexp.obj.sub@assays$RNA@data %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Feature") %>%
  pivot_longer(-Feature, names_to = "Identity", values_to = "Expression") %>%
  separate(Identity, into = c("time", "subtype2"), sep = "h\\.")

# Convert columns to factors with appropriate levels and labels
avgexp$subtype2 <- factor(avgexp$subtype2, levels = c("LPC2", "LLPC2", "LLPC3"))
avgexp$Feature <- factor(avgexp$Feature, levels = my_genes)
avgexp$time <- factor(avgexp$time, levels = c("24", "36", "48", "60", "72", "84", "96"),
                      labels = c("P24", "P36", "P48", "P60", "P72", "P84", "P96"))

# Filter data for plotting
plot_data <- avgexp %>%
  filter(Feature %in% genes_present, subtype2 %in% c("LPC2", "LLPC2", "LLPC3"))

# Generate the heatmap plot
heatmap <- ggplot(plot_data, aes(x = Feature, y = time, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))) +
  facet_wrap(~subtype2, ncol = 1, strip.position = "right") +
  scale_y_discrete(name = "Time") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  coord_cartesian(expand = FALSE) +
  theme(legend.position = "none")

# Print the plot
print(heatmap)
