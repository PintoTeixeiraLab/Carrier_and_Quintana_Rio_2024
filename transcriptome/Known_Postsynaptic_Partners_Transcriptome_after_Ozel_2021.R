# Load necessary libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)

# Set working directory to the specified path
setwd("YOUR_DIRECTORY")

# Load the Seurat object containing the RNA data
full_rna = readRDS('Ozel_2021_Integrated.rds')

# Rename cell identities in the Seurat object
full_rna <- RenameIdents(full_rna,
                         '40' = "Lpi4-3", 
                         "120" = "Lpi3-4",
                         '137' = 'TmY14',
                         '44' = 'Y3-like', 
                         '66' = 'LLPC1',
                         '42' = 'TmY5a', 
                         '75' = 'TmY4',
                         '148' = 'LPLC1',
                         '150' = 'LPLC2')

# Update metadata with   time points names and reorder the factor levels
full_rna@meta.data <- full_rna@meta.data %>% 
  mutate(time = fct_recode(as.factor(MergeOrigin), 
                           "P15" = 'GSE142787_P15',
                           "P30" = 'GSE142787_P30',
                           'P40' = 'GSE142787_P40',
                           'P50' = 'GSE142787_P50',
                           'P70' = 'GSE142787_P70',
                           'Adult' = 'Adult_Annotated-003')) %>%
  mutate(time = forcats::fct_relevel(time, "Adult", after = Inf)) %>%  # Reorder factors to place "Adult" at the end
  mutate(celltype = Idents(full_rna))  # Add cell type information to metadata



# Determine whether to plot all cell types or a subset
full_celltype = FALSE

if (full_celltype) {
  my_celltypes = unique(full_rna@meta.data$celltype)
  subset_rna = full_rna
} else {
  my_celltypes = c("Lpi4-3", "Lpi3-4", "TmY14", "Y3-like", "LLPC1","LPLC1","LPLC2", "TmY5a", "TmY4")
  subset_rna = subset(full_rna, idents = my_celltypes)
}

# Define a list of genes of interest
my_genes = c("side", "CG42313", "CG34113", "CG14372", "CG34371", "CG34114", "CG12950", "CG12484", 
             "beat-Ia", "beat-Ib", "beat-Ic", "beat-IIa", "beat-IIb", "beat-IIIa", "beat-IIIb", 
             "beat-IIIc", "beat-IV", "beat-Va", "beat-Vb", "beat-Vc", "beat-VI", "beat-VII")

# Check if the genes are present in the dataset
subset_rna_rows <- subset_rna@assays$RNA@data %>% rownames
my_genes %in% subset_rna_rows

# Map over the genes to extract metadata and create a combined data frame
names(my_genes) <- my_genes
subset_rna_df <- my_genes %>% map(function(gene) {
  message(gene)
  subset_rna[gene, ]@meta.data %>% 
    dplyr::select(nCount_RNA, time, celltype) %>% 
    tibble::rownames_to_column("cell")
}) %>% bind_rows(.id = "gene")

# Rename and reorder gene factors
subset_rna_df = subset_rna_df %>% 
  mutate(gene = fct_recode(as.factor(gene), 
                           "side-II" = 'CG42313', 
                           "side-III" = 'CG34113', 
                           "side-IV" = 'CG14372', 
                           "side-V" = 'CG34371', 
                           "side-VI" = 'CG34114', 
                           "side-VII" = 'CG12950', 
                           "side-VIII" = 'CG12484')) %>%
  mutate(gene = fct_relevel(gene, 
                            c("side", "side-II", "side-III", "side-IV", "side-V", "side-VI", 
                              "side-VII", "side-VIII", "beat-Ia", "beat-Ib", "beat-Ic", 
                              "beat-IIa", "beat-IIb", "beat-IIIa", "beat-IIIb", "beat-IIIc", 
                              "beat-IV", "beat-Va", "beat-Vb", "beat-Vc", "beat-VI", "beat-VII")))

# Check which genes are present in the dataset
genes_in_data <- my_genes %in% rownames(subset_rna)

# Create a new column combining time and celltype for grouping
subset_rna$newcol <- str_c(subset_rna$time, subset_rna$celltype, sep = ".")

# Calculate average expression per group and subset the data
avgexp.obj = AverageExpression(subset_rna, return.seurat = TRUE, group.by = 'newcol')
avgexp.obj.sub <- subset(avgexp.obj, features = my_genes[genes_in_data])

# Reshape the average expression data for heatmap plotting
avgexp <- avgexp.obj.sub@assays$RNA@data %>% as.data.frame() %>% 
  tibble::rownames_to_column("Feature") %>%
  pivot_longer(-Feature, names_to = "Identity", values_to = "Expression") %>% 
  separate(Identity, into = c("time", "celltype"), sep = "\\.") %>%
  mutate(time = forcats::fct_relevel(time, "Adult", after = Inf)) %>%
  mutate(Feature = fct_recode(as.factor(Feature), 
                              "side-II" = 'CG42313', 
                              "side-III" = 'CG34113', 
                              "side-IV" = 'CG14372', 
                              "side-V" = 'CG34371', 
                              "side-VI" = 'CG34114', 
                              "side-VII" = 'CG12950', 
                              "side-VIII" = 'CG12484')) %>%
  mutate(Feature = fct_relevel(Feature, 
                               c("side", "side-II", "side-III", "side-IV", "side-V", "side-VI", 
                                 "side-VII", "side-VIII", "beat-Ia", "beat-Ib", "beat-Ic", 
                                 "beat-IIa", "beat-IIb", "beat-IIIa", "beat-IIIb", "beat-IIIc", 
                                 "beat-IV", "beat-Va", "beat-Vb", "beat-Vc", "beat-VI", "beat-VII")))

# Plot heatmap of gene expression across cell types and time points
Heatmap <- avgexp %>%
  mutate(celltype = factor(celltype, levels = my_celltypes)) %>%
  ggplot(aes(x = Feature, y = time, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) +
  facet_wrap(~celltype, ncol = 1, strip.position = "right") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  coord_cartesian(expand = FALSE)
