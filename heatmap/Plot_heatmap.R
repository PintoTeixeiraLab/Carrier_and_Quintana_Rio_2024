
# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Set the working directory
setwd("~/Library/CloudStorage/GoogleDrive-filipepts@gmail.com/My Drive/Documents Gdrive/Pinto-Teixeira Lab/Lop manuscript_2023/T4T5_LOP_Connectivity")

# Read data
data <- read_excel("T4T5_OUTPUT_NEURONS_SYNAPSES.xlsx", sheet = 1)

# Rename columns
colnames(data) <- c("cell_name", "Subtype", "Output", "Synapses")

# Identify rows with empty values in the "Output" column
missing_output_rows <- data[!complete.cases(data$Output), ]

# Print rows with empty values in the "Output" column
view(missing_output_rows)

# Remove rows with empty values in the "Output" column
data <- data[complete.cases(data$Output), ]

# Continue with the rest of the script
# Define the desired order of the Subtype levels
desired_order <- c("T4c", "T5c", "T4d", "T5d")

# Convert the "Subtype" column to a factor with the desired order
data$Subtype <- factor(data$Subtype, levels = desired_order)

# Filter the data to include only the specified subtypes
data <- data[data$Subtype %in% desired_order, ]

# Order the data based on the "Subtype" column
data <- data[order(data$Subtype), ]

# If you want to plot all outputs
plot_all_outputs <- FALSE

# Select outputs based on the plotting choice
if (plot_all_outputs) {
  selected_outputs <- unique(data$Output)
} else {
  selected_outputs <- c("VS", "LPi3-4", "LPi4-3", "LPi34-12", "LPi3a", "Tlp12", "Tlp13", "Y3", "Y11", "Y12", "LLPC2", "LPC2", "LLPC3", "LPLC2", "TmY14", "TmY5a","TmY4")
}

# Convert "Output" to a factor with the desired order
data$Output <- factor(data$Output, levels = selected_outputs)

# Filter the data based on the selected outputs
filtered_data <- data %>% filter(Output %in% selected_outputs)

# Sum the number of connections for each combination of cell_name and Output
summed_data <- filtered_data %>%
  group_by(cell_name, Output) %>%
  summarise(Sum_Synapses = sum(Synapses))

# Arrange the columns in the specified order
summed_data <- summed_data %>%
  select(cell_name, Output, Sum_Synapses)

# Ensure the order of cell_name corresponds to the order of Subtype
summed_data$cell_name <- factor(summed_data$cell_name, levels = unique(data$cell_name))

# Create the heatmap
heatmap_cellname_output <- ggplot(summed_data, aes(x = Output, y = cell_name, fill = Sum_Synapses)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = brewer.pal(n = 9, name = "Blues"),
    na.value = "white",
    guide = guide_colorbar(title = "Synapse Sum")
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    legend.position = "right",
    legend.justification = "right"
  ) +
  geom_hline(yintercept = seq(0.5, nlevels(summed_data$cell_name) - 0.5), color = "black", size = 0.2) +
  geom_vline(xintercept = seq(0.5, length(unique(summed_data$Output)) - 0.5), color = "black", size = 0.2) +
  theme(axis.line = element_line(size = 0.2))

# Print the heatmap
print(heatmap_cellname_output)


