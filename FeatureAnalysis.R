# Overview ----
# Title: SPIAT feature extraction analysis.
# Date: 20/05/2023
# Author: Cyril P
# Director: Carlos de Andrea
# Institution: University of Navarra
# Libraries ----
library(gridExtra)
library(BBmisc)
library(stats)
library(readxl)
library(car)
library(moments)
library(reshape2)
library(pbapply)
library(ggridges)
library(FSA)
library(ggsignif)
library(corrplot)
library(spatstat)
library(RColorBrewer)
library(ggExtra)
library(ggplotify)
library(grid)
library(plotly)
library(tidyverse)
library(SPIAT)
library(parallel)
library(parallelly)
library(pheatmap)
library(igraph)
library(heatmaply)
library(dplyr)
library(magrittr)
library(fields)
library(viridis)

# Load and format data ----
# Raw measurements exported from QuPath
raw_measurements <- read.table("Data/raw_measurements.tsv", header = TRUE, sep = "\t")

# Clinical data
responses <- read.csv("Data/Responses.csv", header = TRUE, sep = ",")

## Formatting
# Treatment response dataset
responses <- responses[, -c(4, 6:8)] # Eliminate unnecessary columns
colnames(responses) <- c("ID", "Age", "Sex", "Response") # Set name to columns
responses$ID <- sub("^0", "", responses$ID) # Eliminate 0 suffix to compare to eda_df
responses$Response <- factor(responses$Response)
responses$Sex <- factor(responses$Sex)

# Exploratory dataset
eda_df <- raw_measurements
colnames(eda_df) <- c("Image", "X", "Y", "DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")
rownames(eda_df) <- str_c("Cell", as.character(1: nrow(eda_df)), sep = "_")

# Correct NAs
summary(eda_df) # Summarise data
NAs <- is.na(eda_df) # Find NAs
images_with_NAs <- unique(eda_df$Image[rowSums(NAs) > 0]) # Find images NAs come from
eda_df <- eda_df[!eda_df$Image %in% images_with_NAs, ] # Remove images with NAs

# Match response groups to exploratory dataset
eda_df$ID <- gsub(".*?(\\d{3}|\\d{4})\\.tif", "\\1", eda_df$Image) # Extract patient ID pt 1/2
eda_df$ID <- as.numeric(eda_df$ID) # Extract patient ID pt 2/2
response_vector <- responses$Response[match(eda_df$ID, responses$ID)] # Match response to patient ID
eda_df$Response <- response_vector

# Drop NE group
responses$ID <- as.numeric(responses$ID)
eda_df <- eda_df[eda_df$Response != "Not evaluable/NE", ]
responses <- responses[responses$Response != "Not evaluable/NE", ]

eda_df <- eda_df[, -ncol(eda_df)]
new_order <- c("Image", "X", "Y", "ID", "DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")
eda_df <- eda_df[, new_order]

# Counts cells per TIFF and binds it to responses df
cell_counts <- table(eda_df$ID)
cell_counts<- as.data.frame(cell_counts)
colnames(cell_counts) <- c("ID", "Cell_Count")
cell_count_vector <- cell_counts$Cell_Count[match(responses$ID, cell_counts$ID)]
responses$total_cells <- cell_count_vector

# Final retouches
# Create convenient lists
SPIAT_tifs <- unique(eda_df$Image) # Contains all curated images
SPIAT_px <- unique(eda_df$ID) # Contains all patient IDs

responses <- responses %>%
  left_join(
    eda_df %>%
    group_by(ID) %>%
    summarise(n_images = n_distinct(Image)))

responses <- responses %>%
  left_join(
    eda_df %>%
      group_by(ID) %>%
      summarise(image_names = list(unique(Image))))

responses$SPIAT_indices <- lapply(responses$image_names, function(x) match(x, SPIAT_tifs))

responses$Response <- droplevels(responses$Response, "Not evaluable/NE")

levels(responses$Response) <- c("CR", "PR", "PD", "SD")

# Cleaning memory
suppressWarnings(rm(list = c("cell_counts", "NAs", "cell_count_vector", "images_with_NAs", "new_order", "response_vector")))
cat("\014")
# SPIAT package modules ----
# 1. Create a spatial object ----
# Select image subset
selected_image <- SPIAT_tifs[65]
selected_data <- subset(raw_measurements, Image %in% selected_image)
selected_data <- selected_data %>% # Remove "Image" name column
  select(-1)

# Data formating 
x_cord <- t(selected_data)[1 , ] 
y_cord <- t(selected_data)[2 , ]

intensity_matrix <- t(selected_data[ , -c(1:2)])
intensity_matrix <- as.data.frame(intensity_matrix)

colnames(intensity_matrix) <- str_c("Cell", as.character(1: ncol(intensity_matrix)), sep = "_")
rownames(intensity_matrix) <- c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

phenotype = phenotype = rep(NA, ncol(intensity_matrix))

# Create spatial experiment (S4 object)
general_format_image <- format_image_to_spe(
  format = "general",
  phenotypes = phenotype,
  intensity_matrix = intensity_matrix,
  coord_x = x_cord,
  coord_y = y_cord)

# Inspecting the spatial experiment object
dim(general_format_image)

# Phenotypes and cell properties
colData(general_format_image)[1:5, ]

# Cell coordinates
spatialCoords(general_format_image)[1:5, ]

# Marker intensities
assay(general_format_image)[, 1:5]

# Predict cell phenotypes
predicted_image <- predict_phenotypes(
  spe_object = general_format_image,
  thresholds = NULL,
  tumour_marker = "CK",
  baseline_markers = c("PD1", "CD8", "CD3", "TIM3", "LAG3"),
  nuclear_marker = "DAPI",
  reference_phenotypes = FALSE,
  plot_distribution = TRUE)

# 2. Quality control and visualization ----
# Boxplot of marker intensities
marker_names <- predicted_image@rowRanges@partitioning@NAMES
marker_names <- marker_names[marker_names != "DAPI"]
plots_list <- list()

for (marker in marker_names) {
  tryCatch({
    marker_plot <- marker_intensity_boxplot(predicted_image, marker)
    plots_list[[marker]] <- marker_plot
  }, error = function(e) {
    message(paste("Error for marker", marker, ":", e$message))
  })
}

combined_plots <- do.call(grid.arrange, c(plots_list, ncol = 2))
print(combined_plots)

# Scatter plots of marker levels
plots_list <- list()

for (marker in marker_names) {
  tryCatch({
    marker_plot <- plot_cell_marker_levels(predicted_image, marker)
    plots_list[[marker]] <- marker_plot
  }, error = function(e) {
    message(paste("Error for marker", marker, ":", e$message))
  })
}

combined_plots <- do.call(grid.arrange, c(plots_list, ncol = 2))
print(combined_plots)

# Heatmaps of marker levels
plots_list <- list()

for (marker in marker_names) {
  tryCatch({
    marker_plot <- plot_marker_level_heatmap(predicted_image, num_splits = 100, marker)
    plots_list[[marker]] <- marker_plot
  }, error = function(e) {
    message(paste("Error for marker", marker, ":", e$message))
  })
}

combined_plots <- do.call(grid.arrange, c(plots_list, ncol = 2))
print(combined_plots)

# Identifying incorrect phenotypes
unique(predicted_image$Phenotype)

g <- dimensionality_reduction_plot(
  predicted_image, 
  plot_type = "TSNE",
  feature_colname = "Phenotype",
  perplexity = 20)

plotly::ggplotly(g) 

# Specifying cell types
formatted_image <- define_celltypes(
  predicted_image, 
  categories = c("CK",
                 "PD1,CD8,CD3,TIM3,LAG3", 
                 "CD3", 
                 "CD8,CD3",
                 "CD8,CD3,TIM3",
                 "PD1,CD8,CD3", 
                 "CD8,CD3,LAG3"), 
  category_colname = "Phenotype", 
  names = c("Tumour", "PanExhaust", "T-cell", "Active T-cell", "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3"),
  new_colname = "Cell.Type")

cleaned_image <- select_celltypes(
  formatted_image,
  "Undefined",
  feature_colname = "Cell.Type",
  keep = FALSE)

# Dimentionality reduction to identify misclassified cells
cleaned_g <- dimensionality_reduction_plot(
  predicted_image, 
  plot_type = "TSNE",
  feature_colname = "Cell.Type")

plotly::ggplotly(cleaned_g) 

# Comparison between the two TSNE
g_plotly <- plotly::ggplotly(g)
cleaned_g_plotly <- plotly::ggplotly(cleaned_g)

combined_plot <- plotly::subplot(g_plotly, cleaned_g_plotly, nrows = 1) %>%
  layout(title = list(text = "Raw vs Formatted t-SNE"))

combined_plot

# Categorical dot plot
p <- plot_cell_categories(
  predicted_image, 
  c("None", "CK", "CD3,CK", "CD3,TIM3,LAG3,CK", "TIM3,CK", "PD1,CD3,TIM3,CK", "PD1,CD3,CK", "LAG3,CK", "CD3,LAG3,CK", "CD3,TIM3","CD3","CD8,CK", "PD1,CD8,CD3,TIM3", "CD8,CD3,TIM3", "PD1,CD3,TIM3", "CD3,TIM3,CK", 
    "TIM3", "PD1,CD3", "CD8", "CD8,CD3", "PD1,CK", "CD8,TIM3,CK", "CD8,CD3,TIM3,CK"),
  color_vector <- c("grey", "yellow", "cyan", "red", "darkred", "purple", "darkgreen", "blue", "orange", "brown", "pink", "magenta", "limegreen", "gold", "tan", "steelblue", "darkgrey", "darkblue", "darkorange", "darkmagenta", "green4", "skyblue", "darkcyan"), 
  "Phenotype",
  layered = TRUE)
ggplotly(p) 

# 3D surface plot
# Warning: Uses rownames(SummarizedExperiment::assay(spe_object)) to compare markers so only works with original markers and not phenotypes
marker_names <- c("PD1",  "LAG3", "CK")
marker_surface_plot_stack(predicted_image, num_splits=50, marker_names)

# 3. Basic analyses ----
# Cell percentages
p_cells <- calculate_cell_proportions(predicted_image, 
                                      reference_celltypes = NULL, 
                                      feature_colname = "Phenotype",
                                      celltypes_to_exclude = "Others",
                                      plot.image = TRUE)

p_cells

plot_cell_percentages(
  cell_proportions = p_cells, 
  cells_to_exclude = "Tumour", 
  cellprop_colname = "Proportion_name")

# Cell distances
distances <- calculate_pairwise_distances_between_celltypes(
  spe_object = formatted_image, 
  cell_types_of_interest = c("Tumour", "PanExhaust", "T-cell", "Active T-cell", "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3"),
  feature_colname = "Cell.Type")

plot_cell_distances_violin(distances)

summary_distances <- calculate_summary_distances_between_celltypes(distances)

summary_distances

plot_distance_heatmap(phenotype_distances_result = summary_distances, metric = "median")

# Minimum cell distances
min_dist <- calculate_minimum_distances_between_celltypes(
  spe_object = formatted_image, 
  cell_types_of_interest = c("Tumour", "PanExhaust", "T-cell", "Active T-cell", "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3"),
  feature_colname = "Cell.Type")

plot_cell_distances_violin(cell_to_cell_dist = min_dist)

min_summary_dist <- calculate_summary_distances_between_celltypes(min_dist)

min_summary_dist

plot_distance_heatmap(phenotype_distances_result = min_summary_dist, metric = "mean")

# 4. Quantifying cell co-localisation ----
# Cells in Neighbourhood (CIN)
average_percentage_of_cells_within_radius(spe_object = cleaned_image, 
                                          reference_celltype = "Tumour", 
                                          target_celltype = "ExhaustedLAG3", 
                                          radius=100, 
                                          feature_colname = "Cell.Type")

marker_names <- predicted_image@rowRanges@partitioning@NAMES
plots_list <- list()

for (marker in marker_names) {
  tryCatch({
    marker_plot <- plot_average_intensity(spe_object=predicted_image, reference_marker="CK", 
                                          target_marker=marker, radii=c(15, 25, 35, 50, 75, 100))
    
    plots_list[[marker]] <- marker_plot
  }, error = function(e) {
    message(paste("Error for marker", marker, ":", e$message))
  })
}

combined_plots <- do.call(grid.arrange, c(plots_list, ncol = 2))
print(combined_plots)

# Mixing score (MS) and normalised mixing score (NMS)
mixing_score_summary(spe_object = formatted_image, 
                     reference_celltype = "Tumour", 
                     target_celltype = "PanExhaust", 
                     radius = 100, 
                     feature_colname = "Cell.Type")

# Cross K function
marker_names <- c("PanExhaust", "T-cell", "Active T-cell", "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3")
plots_list <- list()

for (marker in marker_names) {
  tryCatch({
    df_cross <- calculate_cross_functions(formatted_image, 
                                          method = "Kcross", 
                                          cell_types_of_interest = c("Tumour", marker), 
                                          feature_colname ="Cell.Type",
                                          dist = 100)
    plots_list[[marker]] <- df_cross
  }, error = function(e) {
    message(paste("Error for marker", marker, ":", e$message))
  })
}

AUC_of_cross_function(df_cross) # Mixed
crossing_of_crossK(df_cross) # Ring

# 5. Spatial heterogeneity ----
# Localised enthropy
calculate_entropy(formatted_image, 
                  cell_types_of_interest = c("PanExhaust","ExhaustedLAG3"), 
                  feature_colname = "Cell.Type")

# Fishnet grid
grid <- grid_metrics(formatted_image, 
                     FUN = calculate_entropy, 
                     n_split = 20,
                     cell_types_of_interest = c("Tumour", "PanExhaust" , "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3"), 
                     feature_colname = "Cell.Type")

calculate_percentage_of_grids(grid, 
                              threshold = 0.75, 
                              above = TRUE)

calculate_spatial_autocorrelation(grid, metric = "globalmoran")

# Gradients (concentric circles)
gradient_positions <- c(15, 50, 100)

gradient_entropy <- 
  compute_gradient(formatted_image, 
                   radii = gradient_positions, 
                   FUN = calculate_entropy,  
                   cell_types_of_interest = c("Tumour", "PanExhaust" , "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3"),
                   feature_colname = "Cell.Type")

length(gradient_entropy)

gradient_pos <- seq(50, 500, 50)
gradient_results <- entropy_gradient_aggregated(formatted_image, 
                                                cell_types_of_interest = c("Tumour", "ExhaustedLAG3"),
                                                feature_colname = "Cell.Type", 
                                                radii = gradient_pos)

plot(1:10,gradient_results$gradient_df[1, 3:12])


# 6. Characterising tissue structures ----
# Determining whether there is a clear tumour margin
R_BC(formatted_image, cell_type_of_interest = "Tumour", "Cell.Type")

# Automatic identification of the tumour margin
formatted_border <- identify_bordering_cells(formatted_image, 
                                             reference_cell = "Tumour", 
                                             feature_colname = "Cell.Type")
attr(formatted_border, "n_of_clusters")

# Classification of cells relative to the margin
formatted_distance <- calculate_distance_to_margin(formatted_border)

names_of_immune_cells <- c("Tumour", "PanExhaust", "T-cell", "Active T-cell", "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3")

formatted_structure <- define_structure(
  formatted_distance, 
  cell_types_of_interest = names_of_immune_cells, 
  feature_colname = "Cell.Type", 
  n_margin_layers = 5)

categories <- unique(formatted_structure$Structure)

plot_cell_categories(formatted_structure, feature_colname = "Structure")

immune_proportions <- calculate_proportions_of_cells_in_structure(
  spe_object = formatted_structure, 
  cell_types_of_interest = names_of_immune_cells, feature_colname = "Cell.Type")

immune_proportions

immune_distances <- calculate_summary_distances_of_cells_to_borders(
  spe_object = formatted_structure, 
  cell_types_of_interest = names_of_immune_cells, feature_colname = "Cell.Type")

immune_distances

# 7. Cellular neighborhood (CIN) ----
# Cellular neighborhood
average_minimum_distance(formatted_image)

clusters <- identify_neighborhoods(
  formatted_image, 
  method = "hierarchical", 
  min_neighborhood_size = 50,
  cell_types_of_interest =  c("Tumour", "PanExhaust", "T-cell", "Active T-cell", "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3"), 
  radius = 10, 
  feature_colname = "Cell.Type")

neighorhoods_vis <- 
  composition_of_neighborhoods(clusters, feature_colname = "Cell.Type")
neighorhoods_vis <- 
  neighorhoods_vis[neighorhoods_vis$Total_number_of_cells >=5,]

plot_composition_heatmap(neighorhoods_vis, feature_colname="Cell.Type")

# Average nearest neighbour index (ANNI)
average_nearest_neighbor_index(clusters, 
                               reference_celltypes=c("Cluster_1"), 
                               feature_colname="Neighborhood", 
                               p_val = 0.05)

average_nearest_neighbor_index(
  formatted_image, 
  reference_celltypes = c("Tumour", "PanExhaust", "T-cell", "Active T-cell", "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3"), 
  feature_colname="Cell.Type", 
  p_val = 0.05)

# F1 | Predicted phenotypes -----
## PROCESS - Extraction and formatting
spiat_predicted_phenotypes <- list()

# Custom extraction function
extract_data_from_spatial_experiment <- function(spe_object) {
  spatial_coords <- spe_object@int_colData@listData$spatialCoords
  phenotype <- spe_object@colData@listData$Phenotype
  
  return(list(spatial_coords = spatial_coords, phenotype = phenotype))
}

# Extract all predictions using pblapply
spiat_predicted_phenotypes <- pblapply(SPIAT_tifs, function(selected_image) {
  
  # Filter data based on the current image
  selected_data <- subset(raw_measurements, Image %in% selected_image)
  selected_data <- selected_data %>% 
    select(-1) # Remove "Image" name column
  
  # Extract coordinates and intensity matrix
  x_cord <- t(selected_data)[1 , ] 
  y_cord <- t(selected_data)[2 , ]
  
  intensity_matrix <- t(selected_data[ , -c(1:2)])
  intensity_matrix <- as.data.frame(intensity_matrix)
  
  colnames(intensity_matrix) <- paste0("Cell_", as.character(1: ncol(intensity_matrix)))
  rownames(intensity_matrix) <- c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")
  
  # Extract phenotype for the current image
  phenotype = rep(NA, ncol(intensity_matrix))
  
  # Create spatial experiment (S4 object) for the current image
  general_format_image <- format_image_to_spe(
    format = "general",
    phenotypes = phenotype,
    intensity_matrix = intensity_matrix,
    coord_x = x_cord,
    coord_y = y_cord
  )
  
  # Predict cell phenotypes for the current image
  predicted_image <- predict_phenotypes(
    spe_object = general_format_image,
    thresholds = NULL,
    tumour_marker = "CK",
    baseline_markers = c("PD1", "CD8", "CD3", "TIM3", "LAG3"),
    nuclear_marker = "DAPI",
    reference_phenotypes = FALSE,
    plot_distribution = FALSE
  )
  
  # Extract the required data from the SpatialExperiment object
  extracted_data <- extract_data_from_spatial_experiment(predicted_image)
  
  # Return the data along with the image name
  return(list(image_name = selected_image, spatial_coords = extracted_data$spatial_coords, phenotype = extracted_data$phenotype))
})

# Name each entry with their image name
names(spiat_predicted_phenotypes) <- sapply(spiat_predicted_phenotypes, function(x) x$image_name)

## Obtain predicted phenotype counts
# Step 1: Extract phenotypes
phenotype_vector <- unlist(lapply(spiat_predicted_phenotypes, function(x) x$phenotype))

# Step 2: Get unique phenotypes of all images
unique_phenotypes <- unique(phenotype_vector)

# Step 3: Create df with image names and columns named based on unique phenos
spiat_pheno_counts <- data.frame(image_name = names(spiat_predicted_phenotypes),
                                 matrix(NA, ncol = length(unique_phenotypes), nrow = length(spiat_predicted_phenotypes)))

colnames(spiat_pheno_counts) <- c("Image", unique_phenotypes)

# Step 4: Iterate through images, count phenotypes, and updat df
for (i in seq_along(spiat_predicted_phenotypes)) {
  image <- spiat_predicted_phenotypes[[i]]
  image_name <- names(spiat_predicted_phenotypes[i])
  image_phenotype_counts <- table(image$phenotype)
  
  # Finding the row index for the current image in the dataframe
  row_index <- which(spiat_pheno_counts$Image == image_name)
  
  # For each unique phenotype in the count table, update the dataframe
  phenotypes_in_image <- names(image_phenotype_counts)
  
  for (phenotype in phenotypes_in_image) {
    # Finding the column index for the current phenotype
    col_index <- which(colnames(spiat_pheno_counts) == phenotype)
    
    # Updating the corresponding cell with the count from image_phenotype_counts
    spiat_pheno_counts[row_index, col_index] <- image_phenotype_counts[phenotype]
  }
}

## Re-ordering columns
# Add ID and response columns
spiat_pheno_counts$ID <- eda_df$ID[match(spiat_pheno_counts$Image, eda_df$Image)]
spiat_pheno_counts$Response <- responses$Response[match(spiat_pheno_counts$ID, responses$ID)]
current_cols <- setdiff(colnames(spiat_pheno_counts), c("Image","ID", "Response", "None", "CK"))

# Assign primary weights (prefixes)
get_primary_weight <- function(col) {
  if (startsWith(col, "CD3")) return(1)
  if (startsWith(col, "CD8")) return(2)
  if (startsWith(col, "PD1")) return(3)
  if (startsWith(col, "LAG3")) return(4)
  if (startsWith(col, "TIM3")) return(5)
  
  # Return a high weight for unmatched columns so they appear at the end
  return(100)
}

# Assign secondary weight (number of immune markers)
get_secondary_weight <- function(col) {
  
  # Count number of commas and add 1 to determine the number of markers
  return(1 + sum(nchar(gsub("[^,]", "", col))))
}

# Calculate weights for columns
primary_weights <- sapply(current_cols, get_primary_weight)
secondary_weights <- sapply(current_cols, get_secondary_weight)

# Combine primary and secondary weights to get the order
combined_weights <- primary_weights * 1000 + secondary_weights

# Sort columns based on combined weights
ordered_cols <- current_cols[order(combined_weights)]

# Add "Image", "None", and "CK" columns at the start
final_order <- c("Image","ID", "Response", "None", "CK", ordered_cols)

# Re-arrange columns (final order)
spiat_pheno_counts <- spiat_pheno_counts[, final_order]

## PROCESS - Convert to proportions
# Reshape data
spiat_long_format <- spiat_pheno_counts %>%
  select(-Image) %>%
  gather(key = 'spiat_Phenotype', value = 'spiat_Count', -ID, -Response)

# Group data by ID and summarize counts
spiat_total_cells_per_ID <- spiat_long_format %>%
  group_by(ID) %>%
  summarise(spiat_TotalCells = sum(spiat_Count, na.rm = TRUE))

# Mutate to long format
spiat_long_with_props <- spiat_long_format %>%
  group_by(ID, spiat_Phenotype) %>%
  summarise(spiat_PhenotypeCount = sum(spiat_Count, na.rm = TRUE)) %>%
  inner_join(spiat_total_cells_per_ID, by = "ID") %>%
  mutate(spiat_Prop = spiat_PhenotypeCount / spiat_TotalCells)

# Assign responses
spiat_long_with_props$spiat_Response <- responses$Response[match(spiat_long_with_props$ID, responses$ID)]

spiat_long_with_props$spiat_Response <- factor(
  spiat_long_with_props$spiat_Response,
  levels = c("PD", "SD", "PR", "CR")
)

# Aggregate by Phenotype and Response
spiat_aggregated_pheno <- spiat_long_with_props %>%
  group_by(spiat_Phenotype, spiat_Response) %>%
  summarise(spiat_PhenotypeSum = sum(spiat_PhenotypeCount)) %>%
  ungroup() %>%
  left_join(
    spiat_long_with_props %>%
      group_by(spiat_Response) %>%
      summarise(spiat_TotalSum = sum(spiat_PhenotypeCount)),
    by = "spiat_Response"
  )

# Calculate proportions for each group
spiat_aggregated_pheno <- spiat_aggregated_pheno %>%
  mutate(spiat_GroupProp = spiat_PhenotypeSum / spiat_TotalSum)

# Absent markers counted as 0 instead of ignored
spiat_pheno_counts[,4:67] <- apply(spiat_pheno_counts[,4:67], 2, function(x) ifelse(is.na(x), 0, x))

# Aggregation by summation (since phenotypes are counts, we test proportions instead of counts to lessen the effects of repeated paired measures)
spiat_pheno_temp <- spiat_pheno_counts %>%
  select(-Image) %>%  # Exclude the Image column
  group_by(ID, Response) %>%
  summarise(across(everything(), sum)) %>%
  ungroup()

spiat_pheno_temp$total_cells <- rowSums(spiat_pheno_temp[, -(1:2)])

spiat_pheno_temp[, -(1:2)] <- spiat_pheno_temp[, -(1:2)] / spiat_pheno_temp$total_cells

spiat_pheno_temp$total_cells <- NULL

spiat_phenoKW <- spiat_pheno_temp %>%
  pivot_longer(cols = -c(Response, ID), names_to = "Phenotype", values_to = "Value")

phenos <- unique(spiat_phenoKW$Phenotype)

# Storage
spiat_pheno_levene <- data.frame(
  Phenotype = character(0), 
  levene = numeric(0)
)

spiat_phenoKWfdr <- data.frame(
  Phenotype = character(0), 
  p.value = numeric(0)
)

spiat_pheno_dunn <- data.frame(
  Phenotype = character(0),
  Comparison = character(0),
  Z = numeric(0),
  P.unadj = numeric(0),
  P.adj = numeric(0)
)

# Levene's test
for (pheno in phenos) {
  lev_test <- leveneTest(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == pheno, ])$`Pr(>F)`[1]
  spiat_pheno_levene <- rbind(spiat_pheno_levene, data.frame(Phenotype = pheno, levene = lev_test))
}

# KW test
for (pheno in phenos) {
  p_val <- kruskal.test(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == pheno, ])$p.value
  spiat_phenoKWfdr <- rbind(spiat_phenoKWfdr, data.frame(Phenotype = pheno, p.value = p_val))
}

# Dunn's test for KW p < 0.05
for (pheno in spiat_phenoKWfdr$Phenotype[spiat_phenoKWfdr$p.value < 0.05]) {
  dunn_results <- dunnTest(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == pheno, ], method = "bh")
  dunn_df <- data.frame(
    Phenotype = rep(pheno, nrow(dunn_results$res)),
    Comparison = dunn_results$res$Comparison,
    Z = dunn_results$res$Z,
    P.unadj = dunn_results$res$`P.unadj`,
    P.adj = dunn_results$res$`P.adj`
  )
  spiat_pheno_dunn <- rbind(spiat_pheno_dunn, dunn_df)
}

## VISUALS - Relative phenotype proportion per patient (Barplot)
# Merge datasets
merged_data <- merge(spiat_long_with_props, responses[, c("ID", "Response")], by = "ID")

# Order IDs based on response
ordered_ids <- merged_data %>% 
  arrange(factor(Response, levels = c("CR", "PR", "SD", "PD"))) %>%
  pull(ID) %>%
  unique()

# Re-level factor levels of ID
merged_data$ID <- factor(merged_data$ID, levels = ordered_ids)

gg <- ggplot(merged_data, aes(x=factor(ID), y=spiat_Prop, fill=spiat_Phenotype)) +
  geom_bar(stat="identity", position="fill") + 
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)))

gg

## VISUALS - Relative phenotype proportion per response group (Barplot)
spiat_plot <- ggplot(spiat_aggregated_pheno, aes(x = spiat_Response, y = spiat_GroupProp, fill = spiat_Phenotype)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Proportion of Each Phenotype by Response Group (SPIAT)",
    x = "Response Group",
    y = "Proportion",
    fill = "Phenotype"
  )

plot(spiat_plot)

# Cleaning space from memory
suppressWarnings(
  rm(col_index, current_cols, i, image, image_name, image_phenotype_counts, ordered_cols, phenotype, phenotype_vector, phenotypes_in_image, row_index, 
     unique_phenotypes, primary_weights, secondary_weights, combined_weights, get_primary_weight, get_secondary_weight, final_order, spiat_pheno_temp, 
     marker, markers, p_val, dunn_df, dunn_results, lev_test, spiat_long_format, spiat_phenoKW, spiat_total_cells_per_ID, pheno, spiat_aggregated_pheno,
     spiat_pheno_levene, spiat_plot, gg, ordered_ids, phenos, spiat_long_with_props)
)

# Clean everything in the console
cat("\014")
## STATS - KW comparison & post-hoc print
print(
  spiat_phenoKWfdr %>% 
    filter(p.value <= 0.05) %>% 
    arrange(p.value)
)

print(
  spiat_pheno_dunn %>% 
    filter(P.adj <= 0.05) %>% 
    arrange(P.adj)
)

## DUNN'S TEST RESULTS
filtered_data <- spiat_pheno_dunn %>% 
  filter(P.adj <= 0.05) %>% 
  arrange(P.adj)

ggplot(filtered_data, aes(x = P.adj, y = reorder(Phenotype, P.adj), color = Comparison)) +
  geom_point(size = 4) +
  geom_text(aes(label = round(Z, 2)), hjust = -0.5, vjust = 0.5) +
  labs(
    title = "Dunn's Post-hoc Test Results",
    x = "Adjusted P-value",
    y = "Phenotypes",
    color = "Comparison Groups"
  ) +
  theme_minimal()
ggplotly(
  p = ggplot2::last_plot())

# F2 | Scatterplots -----
## PROCESS - Extract scatterplots
spiat_scattercount <- list()

# Extraction function
extract_data_from_spatial_experiment <- function(plots_list) {
  
  for (marker in names(plots_list)) {
    X <- plots_list[[marker]]$data$Cell.X.Position
    Y <- plots_list[[marker]]$data$Cell.Y.Position
    PD1 <- plots_list[[marker]]$data$PD1
    CD8 <- plots_list[[marker]]$data$CD8
    CD3 <- plots_list[[marker]]$data$CD3
    LAG3 <- plots_list[[marker]]$data$LAG3
    TIM3 <- plots_list[[marker]]$data$TIM3
    CK <- plots_list[[marker]]$data$CK
    Pheno <- plots_list[[marker]]$data$Phenotype
    L <- list(P = Pheno, X = X, Y = Y, PD1 = PD1, CD8 = CD8, CD3 = CD3, LAG3 = LAG3, TIM3 = TIM3, CK = CK)
  }
  
  return(list(P = L$P, X = L$X, Y = L$Y, PD1 = L$PD1, CD8 = L$CD8, CD3 = L$CD3, LAG3 = L$LAG3, TIM3 = L$TIM3, CK = L$CK))
}

# Loop through each image
spiat_scattercount <- pblapply(SPIAT_tifs, function(selected_image) {
  
  # Filter data
  selected_data <- subset(raw_measurements, Image %in% selected_image)
  selected_data <- selected_data %>% 
    select(-1) # Remove "Image" name column
  
  # Extract coordinates and intensity matrix
  x_cord <- t(selected_data)[1 , ] 
  y_cord <- t(selected_data)[2 , ]
  
  intensity_matrix <- t(selected_data[ , -c(1:2)])
  intensity_matrix <- as.data.frame(intensity_matrix)
  
  colnames(intensity_matrix) <- paste0("Cell_", as.character(1: ncol(intensity_matrix)))
  rownames(intensity_matrix) <- c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")
  
  # Extract phenotype
  phenotype = rep(NA, ncol(intensity_matrix))
  
  # Create spatial experiment (S4 object)
  general_format_image <- format_image_to_spe(
    format = "general",
    phenotypes = phenotype,
    intensity_matrix = intensity_matrix,
    coord_x = x_cord,
    coord_y = y_cord
  )
  
  # Predict cell phenotypes
  predicted_image <- predict_phenotypes(
    spe_object = general_format_image,
    thresholds = NULL,
    tumour_marker = "CK",
    baseline_markers = c("PD1", "CD8", "CD3", "TIM3", "LAG3"),
    nuclear_marker = "DAPI",
    reference_phenotypes = FALSE,
    plot_distribution = FALSE
  )
  
  marker_names <- predicted_image@rowRanges@partitioning@NAMES
  marker_names <- marker_names[marker_names != "DAPI"]
  plots_list <- list()
  
  for (marker in marker_names) {
    tryCatch({
      marker_plot <- plot_cell_marker_levels(predicted_image, marker)
      plots_list[[marker]] <- marker_plot
    }, error = function(e) {
      message(paste("Error for marker", marker, ":", e$message))
    })
  }
  
  # Extract the required data
  extracted_data <- extract_data_from_spatial_experiment(plots_list)

  # Return the data along with the image name
  return(list(Image = selected_image, P = extracted_data$P, X = extracted_data$X, Y = extracted_data$Y, PD1 = extracted_data$PD1, CD8 = extracted_data$CD8, CD3 = extracted_data$CD3, LAG3 = extracted_data$LAG3, TIM3 = extracted_data$TIM3, CK = extracted_data$CK))
  
})

# Name entries
names(spiat_scattercount) <- sapply(spiat_scattercount, function(x) x$Image)

# Cat matching response and ID
spiat_scattercount <- pblapply(spiat_scattercount, function(x) {
  
  # Find matching ID based on image name
  matching_ID <- unique(eda_df$ID[eda_df$Image == x$Image])
  
  if (length(matching_ID) == 1) {
    x$ID <- matching_ID
    
    # Find matching response based on ID
    matching_Response <- responses$Response[responses$ID == matching_ID]
    
    if (length(matching_Response) == 1) {
      x$Response <- matching_Response
      
    } else {
      warning(paste("No matching Response found for ID:", matching_ID))
    }
    
  } else {
    warning(paste("No matching ID found for Image:", x$Image))
  }
  
  return(x)
})
cat("\014")

nq <- 15 # Select n of quadrats (<90) ----
## PROCESS - Creates PPP, density maps and quadrats
# Generate standardized window of observation
allX_list <- lapply(spiat_scattercount, function(x) x$X)
allY_list <- lapply(spiat_scattercount, function(x) x$Y)
center_x <- median(unlist(allX_list))
center_y <- median(unlist(allY_list))
distances <- sqrt((unlist(allX_list) - center_x)^2 + (unlist(allY_list) - center_y)^2)
radius <- max(distances)
w <- disc(radius=radius, centre=c(center_x, center_y)) # Window w

immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Storage lists
pointpatterns <- vector("list", length = length(names(spiat_scattercount))) # PPP objects
density_maps <- vector("list", length = length(names(spiat_scattercount))) # Density maps
qtest_results <- vector("list", length = length(names(spiat_scattercount))) # MC simulations
quadrats <- vector("list", length = length(names(spiat_scattercount))) # Quadratcount

skipped_images <- data.frame(ImageName = character(),
                             ID = character(),
                             Marker = character(),
                             Response = character(),
                             stringsAsFactors = FALSE)

# Name them
names(pointpatterns) <- names(spiat_scattercount)
names(density_maps) <- names(spiat_scattercount)
names(qtest_results) <- names(spiat_scattercount)
names(quadrats) <- names(spiat_scattercount)

## Most important line, does all plots and tests
for (image_name in names(spiat_scattercount)) {
  
  # Subset image from list
  image_data <- spiat_scattercount[[image_name]]
  image_data_df <- as.data.frame(image_data)
  
  # Create lists for storage
  marker_ppps <- list()
  marker_density_maps <- list()
  q_count <- list()
  qtest_images <- list()
  
  # Loop through each marker
  for (marker in immune_markers) {
    
    q_test <- NULL
    
    # Subset the data to only include rows where "P" contains the marker
    subset_data <- image_data_df[grepl(marker, image_data_df$P, fixed = TRUE), ]
    
    # Create a point pattern object
    P <- ppp(subset_data$X, subset_data$Y, 
                         checkdup = FALSE, 
                         window = w)
    
    Q <- quadratcount(P, 
                      nx = nq, 
                      keepempty = T)
    
    # Estimate density
    D <- density.ppp(P,
                     sigma = 85,
                     edge = TRUE,
                     diggle = TRUE)
    
    # Debug
    if (nrow(subset_data) >= nq) {
      
      # Calculate quadrat test
      q_test <- quadrat.test.ppp(P,
                                 alternative = "clustered",
                                 method = "M",
                                 conditional = T,
                                 nsim = 999)
      
      
    } else {
      
      # Check if ID and Response are available for image_name
      id_available <- unlist(responses$ID[sapply(responses$image_names, function(x) any(x == image_name))])
      response_available <- unlist(responses$Response[sapply(responses$image_names, function(x) any(x == image_name))])
      
      if (length(id_available) > 0 && length(response_available) > 0) {
        
        # Store in the skipped_images list
        skipped_images <- rbind(skipped_images, data.frame(ImageName = image_name,
                                                           ID = id_available,
                                                           Marker = marker,
                                                           Response = response_available,
                                                           stringsAsFactors = FALSE))
      }
    } # Monte Carlo requires at least as many cells as quadrats. This else skips images with lower cells and stores them in skipped_image
    
    marker_ppps[[marker]] <- P # Store PPPs
    marker_density_maps[[marker]] <- D # Store density maps
    q_count[[marker]] <- Q # Store quadrat counts
    qtest_images[[marker]] <- q_test # Store quadrat tests
    
  }
  
  # Store images and results
  pointpatterns[[image_name]] <- marker_ppps
  density_maps[[image_name]] <- marker_density_maps
  quadrats[[image_name]] <- q_count
  qtest_results[[image_name]] <- qtest_images
  
  print(paste("Process complete for image:", image_name))
  
}

cat("\014")
print(skipped_images)
readline(prompt=paste(cat("\014"), "Press enter to continue the execution."))

## PROCESS - Aggregation by ID
transformed_quadrats <- list()

# Image's bins to proportions relative to total patient count
for (image_name in names(quadrats)) {
  
  image_data <- quadrats[[image_name]]
  
  # Storage
  transformed_image_data <- list()
  
  for (marker in names(image_data)) {
    
    # Extract counts for current marker
    marker_counts <- image_data[[marker]]
    
    # Calculate total count for current marker
    greped_ID <- which(unlist(sapply(responses$image_names, function(x) image_name %in% x)))
    total_count <- responses$total_cells[greped_ID]
    
    # Transform into relative proportions
    relative_proportions <- round((marker_counts / total_count)*100, 3)
    
    # Store transformed data
    transformed_image_data[[marker]] <- relative_proportions
  }
  
  # Store transformed data for current
  transformed_quadrats[[image_name]] <- transformed_image_data
}

# Aggregate values by average based on: marker, and bin index
grouped_quadrats <- list()

for (patient_id in unique(responses$ID)) {
  
  # Find image names corresponding to current ID
  image_names <- unlist(responses$image_names[responses$ID == patient_id])
  
  grouped_patient_data <- list()
  
  for (image_name in image_names) {
    
    # Extract data for current image
    image_data <- transformed_quadrats[[image_name]]

    for (marker in names(image_data)) {
      
      # Append counts for current marker and image
      grouped_patient_data[[marker]] <- append(grouped_patient_data[[marker]], list(image_data[[marker]]))
    }
  }
  
  # Now, aggregate each marker by the mean
  for (marker in names(grouped_patient_data)) {
    
    # Combine all the counts for this marker into a matrix
    marker_matrix <- do.call(rbind, grouped_patient_data[[marker]])
    
    # Replace NaN with 0
    marker_matrix[is.nan(marker_matrix)] <- 0
    
    # Calculate mean for each bin (column)
    mean_values <- apply(marker_matrix, 2, mean, na.rm = TRUE)
    
    # Create a new quadratcount object with the mean values but retaining the owin object and tess attributes
    new_quadrat <- grouped_patient_data[[marker]][[1]]
    attributes(new_quadrat)$quadratcount <- mean_values
    
    # Store
    grouped_patient_data[[marker]] <- new_quadrat
  }
  
  # Store grouped and aggregated data for the current patient
  grouped_quadrats[[as.character(patient_id)]] <- grouped_patient_data
}
cat("\014")

## PROCESS - Aggregation by response group
grouped_by_response <- list()

for (response_group in unique(responses$Response)) {
  
  grouped_response_data <- list()
  
  # Find IDs corresponding to current response group
  patient_ids <- unique(responses$ID[responses$Response == response_group])
  
  for (patient_id in patient_ids) {
    
    # Extract data for the current patient
    patient_data <- grouped_quadrats[[as.character(patient_id)]]
    
    for (marker in names(patient_data)) {
      
      # Append counts for current marker and patient
      grouped_response_data[[marker]] <- append(grouped_response_data[[marker]], list(patient_data[[marker]]))
    }
  }
  
  # Now, aggregate each marker by mean
  for (marker in names(grouped_response_data)) {
    
    # Combine all the counts for this marker
    marker_matrix <- do.call(rbind, lapply(grouped_response_data[[marker]], function(x) attributes(x)$quadratcount))
    
    # Calculate mean for each quadrant
    mean_val <- apply(marker_matrix, 2, mean, na.rm = TRUE)
    
    # Create a new quadratcount object with the mean values but retaining the owin object and tess attributes
    new_quadrat <- grouped_response_data[[marker]][[1]]
    attributes(new_quadrat)$quadratcount <- mean_val
    
    # Store the new quadrat object
    grouped_response_data[[marker]] <- new_quadrat
  }
  
  # Store the grouped and aggregated data for the current response group
  grouped_by_response[[response_group]] <- grouped_response_data
}
cat("\014")
selected_image <- sample(SPIAT_tifs, 1) #----
# Density map sampler per group ----
# Order of response groups and immune markers
response_order <- c("SD", "PD", "CR", "PR")
marker_order <- c("CD8", "PD1", "CK", "TIM3", "CD3", "LAG3")

group_ppp <- list()
group_densities <- list()

for (response_group in response_order) {
  
  # Filter responses data frame to only include patients in the current response group
  group_images <- unlist(responses$image_names[responses$Response == response_group])
  
  # Sample an image
  sampled_images <- sample(group_images, size = min(1, length(group_images)), replace = FALSE)

  for (image_name in sampled_images) {
    
    # Extract data
    gppp <- pointpatterns[[image_name]]
    dmap <- density_maps[[image_name]]
    
    # Store
    index_name <- paste(response_group, image_name, sep = "_")
    group_ppp[[index_name]] <- gppp
    group_densities[[index_name]] <- dmap
  }
}

# Store min and max density
marker_min_max = list()

for (image_name in names(group_densities)) {
  
  for (marker in names(group_densities[[image_name]])) {
    
    current_density_values = group_densities[[image_name]][[marker]]$v
    
    # Update min and max if current are more extreme
    marker_min_max[[marker]]$min = min(marker_min_max[[marker]]$min, min(current_density_values, na.rm = TRUE))
    marker_min_max[[marker]]$max = max(marker_min_max[[marker]]$max, max(current_density_values, na.rm = TRUE))
  }
}

# Plot
par(mfrow = c(6, length(group_ppp)), mar = c(1,1,1,1), bty = "n")

for (marker in marker_order) {
  for (response_group in response_order) {
    for (image_name in names(group_ppp)) {
      
      if (startsWith(image_name, response_group)) {
        
        f_ppp <- group_ppp[[image_name]][[marker]]
        f_densities <- group_densities[[image_name]][[marker]]
        
        # Use zlim from marker_min_max
        zlim_values = c(0, marker_min_max[[marker]]$max)
        brks <- seq(0, marker_min_max[[marker]]$max, length.out = 6)
        
        plot(f_densities, 
             col = colorRampPalette(c("white", "cyan", "limegreen", "yellow", "firebrick1")), 
             zlim = zlim_values,
             main = ifelse(marker == "CD8", unfactor(responses$Response[sapply(responses$ID, function(x) any(x == paste(gsub(".*_([0-9]+)\\.tif$", "\\1", image_name))))]), paste(gsub(".*_([0-9]+)\\.tif$", "\\1", image_name))))
        
        contour(f_densities,
                add = T,
                labels = NULL,
                drawlabels = FALSE,
                nlevels = 5)
        
        if (image_name == names(group_ppp)[1]) {
          mtext(marker, side = 2, line = 0, cex = 1)
        }
      }
    }
  }
}

# Quadrat index comparison per group ----
summary_table <- data.frame()

for (marker in names(grouped_quadrats[[1]])) {
  
  # Store data for KW
  all_data <- list()
  
  for (response_group in unique(responses$Response)) {
    
    # IDs of current response
    patient_ids <- unique(responses$ID[responses$Response == response_group])
    
    # Extract data
    marker_data <- lapply(patient_ids, function(id) grouped_quadrats[[as.character(id)]][[marker]])
    
    # Bind
    marker_matrix <- do.call(rbind, marker_data)
    
    # Store data for KW
    all_data[[as.character(response_group)]] <- marker_matrix
  }
  
  # Combine
  combined_matrix <- do.call(rbind, all_data)
  
  # Group vector for KW
  group_vector <- as.vector(rep(names(all_data), sapply(all_data, nrow)))
  
  # KW to each bin (column in combined_matrix)
  for (bin_index in 1:ncol(combined_matrix)) {
    
    # Data of current bin
    bin_data <- combined_matrix[, bin_index]
    
    # Perform Kruskal-Wallis test
    kw_result <- kruskal.test(bin_data ~ group_vector)$p.value
    
    # Temp df to store results
    temp_df <- data.frame(
      Marker = marker,
      Bin = bin_index,
      KW_p_value = kw_result
    )
    
    # Append to the summary table
    summary_table <- rbind(summary_table, temp_df)
  }
}

summary_table <- summary_table %>% 
  arrange(KW_p_value) %>%
  filter(KW_p_value <= 0.05)

# Dunn test
dunn_results <- data.frame()

for (i in 1:nrow(summary_table)) {
  
  # Filter
  if (!is.na(summary_table$KW_p_value[i]) && summary_table$KW_p_value[i] <= 0.05) {
    marker_bin <- summary_table[i, c("Marker", "Bin"), drop = FALSE]
    
    # Extract marker and bin index
    marker <- marker_bin$Marker
    bin_index <- marker_bin$Bin
    
    # Extract current bin
    bin_data <- combined_matrix[, bin_index]
    
    # Temp df for data and bin
    temp_data <- data.frame(Value = bin_data, Response = as.factor(group_vector))
    
    # Dunn's test
    dunn_result <- dunnTest(Value ~ Response, data = temp_data, method = "bh")
    
    # Temp df
    temp_df <- data.frame(
      Marker = rep(marker, nrow(dunn_result$res)),
      Bin = rep(bin_index, nrow(dunn_result$res)),
      Comparison = dunn_result$res$Comparison,
      Z = dunn_result$res$Z,
      P = dunn_result$res$`P.unadj`,
      Adj_P = dunn_result$res$`P.adj`
    )
    
    # Append to the Dunn's table
    dunn_results <- rbind(dunn_results, temp_df)
  }
}

# Sort by FDR adj
dunn_results <- dunn_results %>% 
  arrange(Adj_P) %>%
  filter(Adj_P <= 0.05)

## VISUALS - Phenotype intensity per group subsetted by marker
# Order of response groups and immune markers
response_order <- c("SD", "PD", "CR", "PR")
marker_order <- c("CD8", "PD1", "CK", "TIM3", "CD3", "LAG3")

# Min and max of each marker
marker_min_max <- marker_order %>%
  lapply(function(marker) {
    counts <- response_order %>%
      lapply(function(response_group) {
        Q <- grouped_by_response[[response_group]][[marker]]
        as.vector(round(attributes(Q)$quadratcount, 3))
      }) %>%
      unlist()
    c(min = min(counts), max = max(counts))
  }) %>%
  set_names(marker_order)

# Plot layout
par(mfrow = c(6, 4), mar = c(0, 0, 2, 2), oma = c(0, 0, 0, 1))  # Added room on the right for the colorbar

marker_order %<>%
  map(function(marker) {
    
    # Color palette
    marker_min <- marker_min_max[[marker]]['min']
    marker_max <- marker_min_max[[marker]]['max']
    cols <- colorRampPalette(c("white", "yellow", "gold", "orange", "firebrick1", "orchid1", "orchid4", "black"))(1000)
    
    response_order %<>%
      map(function(response_group) {
        
        # Retrieve quadrats
        Q <- grouped_by_response[[response_group]][[marker]] # Class quadratcount
        
        # Plot window
        plot(attr(Q, "tess")$window, main = NULL, border = NA)
        
        # Draw and color each tile
        counts <- as.vector(round(attributes(Q)$quadratcount, 3))
        
        tiles <- attr(Q, "tess")$tiles
        tiles %<>%
          map2(counts, function(tile, count) {
            col_idx <- findInterval(count, seq(marker_min, marker_max, length.out = length(cols)))
            plot(tile, add = TRUE, col = cols[col_idx], border = FALSE)
          })
        
        if (response_group == "SD") {
          mtext(marker, side = 2, line = -1, cex = 1)
          
        } else if (response_group == "PR") {
          par(xpd = NA)  # Allow plotting in outer margins
          
          # Normalized colorbar
          color_range <- seq(marker_min, marker_max, length.out = length(cols))
          colorbar <- as.matrix(color_range)
          image.plot(legend.only = TRUE, zlim = c(marker_min, marker_max),
                     col = cols, horizontal = FALSE,
                     axis.args = list(at = seq(marker_min, marker_max, by = (marker_max - marker_min) / 4),
                                      labels = round(seq(marker_min, marker_max, by = (marker_max - marker_min) / 4), 2)))
          
          par(xpd = FALSE)  # Restore clipping
        }
      })
  })

# Adjust column_positions to align text
column_positions <- seq(1/8, 1, by = 1/4)
mtext(response_order, side = 3, line = -2, outer = TRUE, at = column_positions, cex = 1)

cat("\014")
dunn_results

# Complete Spatial Randomness (CSR) #----
SpaFx <- data.frame(
  Image = character(0),
  PD1_CSR = numeric(0),
  CD8_CSR = numeric(0),
  CD3_CSR = numeric(0),
  TIM3_CSR = numeric(0),
  LAG3_CSR = numeric(0),
  CK_CSR = numeric(0)
)

pval_df <- data.frame(
  Image = character(0),
  PD1_pval = numeric(0),
  CD8_pval = numeric(0),
  CD3_pval = numeric(0),
  TIM3_pval = numeric(0),
  LAG3_pval = numeric(0),
  CK_pval = numeric(0)
)

# Compare all ppps against inhomogeneous Poisson distribution with MC simulations
for (i in seq_along(pointpatterns)) {
  image_name <- names(pointpatterns[i])
  pointpattern <- pointpatterns[[i]]
  csr_list <- list()
  pval_list <- list()
  
  for (marker in names(pointpattern)) {
    ppp_obj <- pointpattern[[marker]]
    
    if (ppp_obj$n > 0) {
      csr_test <- quadrat.test(ppp_obj, alternative = "clustered", method = "M", conditional = T, nsim = 999)
      csr_list[[marker]] <- csr_test$statistic
      pval_list[[marker]] <- csr_test$p.value 
    } else {
      csr_list[[marker]] <- NA
      pval_list[[marker]] <- NA 
    }
  }
  
  new_row <- data.frame(
    Image = image_name,
    PD1_CSR = csr_list$PD1,
    CD8_CSR = csr_list$CD8,
    CD3_CSR = csr_list$CD3,
    TIM3_CSR = csr_list$TIM3,
    LAG3_CSR = csr_list$LAG3,
    CK_CSR = csr_list$CK
  )
  SpaFx <- rbind(SpaFx, new_row)
  row.names(SpaFx) <- NULL 
  
  new_pval_row <- data.frame(
    Image = image_name,
    PD1_pval = pval_list$PD1,
    CD8_pval = pval_list$CD8,
    CD3_pval = pval_list$CD3,
    TIM3_pval = pval_list$TIM3,
    LAG3_pval = pval_list$LAG3,
    CK_pval = pval_list$CK
  )
  
  pval_df <- rbind(pval_df, new_pval_row)
  row.names(pval_df) <- NULL 
}

unique_eda_df <- eda_df[!duplicated(eda_df$Image), ]
SpaFx$ID <- unique_eda_df$ID[match(SpaFx$Image, unique_eda_df$Image)]
pval_df$ID <- unique_eda_df$ID[match(pval_df$Image, unique_eda_df$Image)]

SpaFx$R <- as.character(NA)

# Cat matching responses
for (i in 1:nrow(SpaFx)) {
  matched_id <- SpaFx$ID[i]
  matching_rows <- which(responses$ID == matched_id)
  
  
  if (length(matching_rows) >= 1) {
    SpaFx$R[i] <- as.character(responses$Response[matching_rows[1]])
  }
}

cat("\014")
summary(SpaFx)
summary(pval_df)

# KW test function
perform_kw_test <- function(df, markers, response_col = "R") {
  kw_p_values <- list()
  
  for(marker in markers) {
    df_for_kw <- df[!is.na(df[, marker]), c(marker, response_col)]
    kw_test_result <- kruskal.test(as.formula(paste0(marker, "~ ", response_col)), data = df_for_kw)
    kw_p_values[[marker]] <- kw_test_result$p.value
  }

  data.frame(Marker = names(kw_p_values), P_Value = as.numeric(kw_p_values))
}

# Dunn test function
perform_dunn_test <- function(df, markers, response_col = "R") {
  dunn_results <- list()
  
  for(marker in markers) {
    df_for_dunn <- df[!is.na(df[, marker]), c(marker, response_col)]
    dunn_test_result <- dunnTest(df_for_dunn[, 1], g = as.factor(df_for_dunn[, response_col]), method = "bh")
    dunn_results[[marker]] <- dunn_test_result
  }
  
  dunn_adj_p_values_df <- do.call(rbind, lapply(names(dunn_results), function(marker) {
    res_df <- dunn_results[[marker]]$res
    data.frame(Marker = marker, Comparison = res_df$Comparison, P_adj = res_df$P.adj, Z = res_df$Z)
  }))
  
  dunn_adj_p_values_df[order(dunn_adj_p_values_df$P_adj),]
}

markers <- colnames(SpaFx)[2:7]
response_col <- "R"

# KW test
kw_test_results <- perform_kw_test(SpaFx, markers, response_col)

# Dunn test
dunn_test_results <- perform_dunn_test(SpaFx, markers, response_col)

## VISUALS - CSR (Heatmap)
long_SpaFx <- SpaFx %>%
  gather(key = "Marker", value = "CSR", PD1_CSR:CK_CSR) %>%
  filter(!is.na(CSR))  # Remove NAs

long_SpaFx$Ordered_R <- factor(long_SpaFx$R, levels = c("CR", "PR", "SD", "PD"), ordered = TRUE)  # New Line

# Sort df by ordered treatment response and Image
long_SpaFx <- long_SpaFx %>%
  arrange(Ordered_R, Image)

long_SpaFx$Ordered_Image <- factor(long_SpaFx$Image, levels = unique(long_SpaFx$Image), ordered = TRUE)

# Compute min and max CSR
agg_data <- long_SpaFx %>%
  group_by(Marker) %>%
  summarise(min_CSR = min(CSR, na.rm = TRUE),
            max_CSR = max(CSR, na.rm = TRUE))

# Merge with long_SpaFx
long_SpaFx <- left_join(long_SpaFx, agg_data, by = "Marker")

# Normalize CSR between 0 and 1
long_SpaFx <- long_SpaFx %>%
  mutate(Norm_CSR = (CSR - min_CSR) / (max_CSR - min_CSR))

# Compute global mean of the CSR
global_median <- mean(long_SpaFx$CSR, na.rm = TRUE)

# Generate labels for y-axis
unique_images <- unique(long_SpaFx$Ordered_Image)
response_groups <- unique(long_SpaFx$Ordered_R)
response_labels <- rep("", length(unique_images))

for (r in response_groups) {
  indices <- which(long_SpaFx[!duplicated(long_SpaFx$Ordered_Image), "Ordered_R"] == r)
  response_labels[indices] <- r
}

# Wide-format df for heatmap
heatmap_data <- long_SpaFx %>% 
  select(Ordered_Image, Marker, CSR) %>% 
  spread(key = Marker, value = CSR)

rows <- response_labels
cols <- str_replace(colnames(heatmap_data)[-1], "_CSR", "")

mat <- as.matrix(heatmap_data[, -1])

row_annot <- data.frame(Response = factor(rows, levels = c("CR", "PR", "SD", "PD")))

rownames(row_annot) <- rownames(mat)

# Remove _csr suffix
colnames(mat) <- gsub("_CSR", "", colnames(mat))
dunn_test_results$Marker <- gsub("_CSR", "", dunn_test_results$Marker)

crosses <- dunn_test_results %>%
  filter(dunn_test_results$P_adj <= 0.05)

# Plot dunn test results
p <- ggplot(dunn_test_results, aes(x=Comparison, y=Marker, fill=P_adj)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(11, "RdYlBu"), name = "Adj.p") +
  geom_point(data = crosses, aes(x=Comparison, y=Marker), 
             shape=3, size=6, stroke=3, color="white")

p

# Plot CSR heatmap
heatmaply(mat,
          row_side_colors = row_annot$Response,
          row_side_palette = c("CR" = "chartreuse", "PR" = "tan1", "SD" = "tomato", "PD" = "firebrick4"),
          colors = turbo, 
          na.value = "black",
          main = "",
          scale = "column",
          titleX = FALSE,
          titleY = FALSE,
          subplot_widths = c(0.95, 0.05),  
          margins = c(20, 20, 20, 100),  
          hide_colorbar = F,
          labRow = NULL,  
          labCol = colnames(mat), 
          column_text_angle = 0,  
          showticklabels = c(TRUE, FALSE),
          grid_color = NA, 
          dendrogram = "none",
          width = 800,  
          height = 800,
          plot_method = "ggplot")

dev.set(which = dev.prev())


# Clean everything in the console
cat("\014")

kw_test_results %>%
  filter(P_Value <= 0.05) %>%
  arrange(P_Value)

dunn_test_results %>%
  filter(P_adj <= 0.05) %>%
  arrange(P_adj)

# Optimal number of quadrats ----
# User-provided data and information
allX_list <- lapply(spiat_scattercount, function(x) x$X)
allY_list <- lapply(spiat_scattercount, function(x) x$Y)

center_x <- mean(unlist(allX_list))
center_y <- mean(unlist(allY_list))

# Calculate distances from center to each point
distances <- sqrt((unlist(allX_list) - center_x)^2 + (unlist(allY_list) - center_y)^2)
radius <- max(distances)
w <- disc(radius=max_radius, centre=c(center_x, center_y))
markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Store overall RSS for each number of bins
overall_bin_metrics <- data.frame(BinNumber = numeric(), OverallRSS = numeric(), stringsAsFactors = FALSE)

# Function to calculate RSS for a given n_bins
calc_RSS <- function(n_bins) {
  overall_RSS <- 0
  for (image in names(spiat_scattercount)) {
    image_data <- spiat_scattercount[[image]]
    image_data_df <- as.data.frame(image_data)
    for (marker in markers) {
      sub_df <- image_data_df[grepl(marker, image_data_df$P, fixed = TRUE),]
      if (nrow(sub_df) == 0) {
        next
      }
      point_pattern <- ppp(sub_df$X, sub_df$Y, window = w, checkdup = F)
      bin_counts <- quadratcount(point_pattern, nx = n_bins, keepempty = T)
      observed_counts <- as.numeric(bin_counts)
      total_points = nrow(sub_df)
      expected_counts = rep(total_points / (n_bins^2), n_bins^2)
      RSS = sum((observed_counts - expected_counts)^2)
      overall_RSS <- overall_RSS + RSS
    }
    print(paste("Finished processing image:", image, "for", n_bins, "bins."))
  }
  print(paste("Finished processing for n_bins:", n_bins))
  return(data.frame(BinNumber = n_bins, OverallRSS = overall_RSS))
}

# Use parallel processing
n_cores <- detectCores(logical = T)  # Number of physical cores

# If you want to leave one core free for other tasks, use (n_cores - 1)
cl <- makeCluster(n_cores)

# Export variables to the cluster
clusterExport(cl, c('spiat_scattercount', 'markers', 'w'))

# Load the required packages on each node in the cluster
clusterEvalQ(cl, {
  library(spatstat)
})

# Apply function in parallel
result_list <- parLapply(cl, 5:21, calc_RSS)

# Combine results
overall_bin_metrics <- do.call(rbind, result_list)

# Stop the cluster
stopCluster(cl)

## VISUALS - Bin tests (Elbow plot)
# Plot the results
plot(overall_bin_metrics$BinNumber, overall_bin_metrics$OverallRSS, type="b", xlab="Number of Bins", ylab="Overall RSS")
cat("\014")
# Bandwidth testing ----
# Function to calculate optimal bandwidths for each image
calculate_bandwidths <- function(image_name, df) {
  # Subset data for the current image
  image_data <- df[df$Image == image_name,]
  center_x_image <- mean(image_data$X)
  center_y_image <- mean(image_data$Y)
  w <- disc(radius=max_radius, centre=c(center_x_image, center_y_image))
  
  # Create a point pattern
  ppp_image <- ppp(image_data$X, image_data$Y, w)
  
  # Calculate bandwidths
  bw_diggle_val <- bw.diggle(ppp_image)
  ppl_val <- bw.ppl(ppp_image)
  CvL_val <- bw.CvL(ppp_image)
  adapt_val <- bw.CvL.adaptive(ppp_image)
  abraham_val <- bw.abram.ppp(ppp_image)
  
  # Store results in a list
  result <- list(bw_diggle = bw_diggle_val, ppl = ppl_val, cvl = CvL_val, adpt = adapt_val, abraham = abraham_val)
  print(paste("Bandwidths processed for image:", image_name))
  
  return(result)
}

cl <- makeCluster(parallel::detectCores())
clusterExport(cl, list("calculate_bandwidths", "max_radius", "eda_df"))
clusterEvalQ(cl, library(spatstat))
bandwidth_results <- parallel::parLapply(cl, image_names, function(image_name) {
  calculate_bandwidths(image_name, eda_df)
})

## VISUALS - Bandwidth tests (Histograms)
# Plot distribution of optimal bandwidths
par(mfrow = c(5, 1))
hist(sapply(bandwidth_results, function(x) x$bw_diggle), main = NULL, xlab = 'Bandwidth', breaks = 20)
hist(sapply(bandwidth_results, function(x) x$ppl), main = NULL, xlab = 'Bandwidth', breaks = 20)
hist(sapply(bandwidth_results, function(x) x$cvl), main = NULL, xlab = 'Bandwidth', breaks = 20)
hist(sapply(bandwidth_results, function(x) x$adpt), main = NULL, xlab = 'Bandwidth', breaks = 20)
hist(unlist(sapply(bandwidth_results, function(x) x$abraham)), main = NULL, xlab = 'Bandwidth', breaks = 20)
cat("\014")
############################################ F3 | Heatmaps -----
# EXTRACT - Extract intensity heatmaps ----
spiat_scattercount <- list()

# Extraction function
extract_data_from_spatial_experiment <- function(plots_list) {
  
  for (marker in names(plots_list)) {
    X <- plots_list[[marker]]$data$Cell.X.Position
    Y <- plots_list[[marker]]$data$Cell.Y.Position
    PD1 <- plots_list[[marker]]$data$PD1
    CD8 <- plots_list[[marker]]$data$CD8
    CD3 <- plots_list[[marker]]$data$CD3
    LAG3 <- plots_list[[marker]]$data$LAG3
    TIM3 <- plots_list[[marker]]$data$TIM3
    CK <- plots_list[[marker]]$data$CK
    Pheno <- plots_list[[marker]]$data$Phenotype
    L <- list(P = Pheno, X = X, Y = Y, PD1 = PD1, CD8 = CD8, CD3 = CD3, LAG3 = LAG3, TIM3 = TIM3, CK = CK)
  }
  
  return(list(P = L$P, X = L$X, Y = L$Y, PD1 = L$PD1, CD8 = L$CD8, CD3 = L$CD3, LAG3 = L$LAG3, TIM3 = L$TIM3, CK = L$CK))
}

# Loop through each image
spiat_scattercount <- pblapply(SPIAT_tifs, function(selected_image) {
  
  # Filter data
  selected_data <- subset(raw_measurements, Image %in% selected_image)
  selected_data <- selected_data %>% 
    select(-1) # Remove "Image" name column
  
  # Extract coordinates and intensity matrix
  x_cord <- t(selected_data)[1 , ] 
  y_cord <- t(selected_data)[2 , ]
  
  intensity_matrix <- t(selected_data[ , -c(1:2)])
  intensity_matrix <- as.data.frame(intensity_matrix)
  
  colnames(intensity_matrix) <- paste0("Cell_", as.character(1: ncol(intensity_matrix)))
  rownames(intensity_matrix) <- c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")
  
  # Extract phenotype
  phenotype = rep(NA, ncol(intensity_matrix))
  
  # Create spatial experiment (S4 object)
  general_format_image <- format_image_to_spe(
    format = "general",
    phenotypes = phenotype,
    intensity_matrix = intensity_matrix,
    coord_x = x_cord,
    coord_y = y_cord
  )
  
  # Predict cell phenotypes
  predicted_image <- predict_phenotypes(
    spe_object = general_format_image,
    thresholds = NULL,
    tumour_marker = "CK",
    baseline_markers = c("PD1", "CD8", "CD3", "TIM3", "LAG3"),
    nuclear_marker = "DAPI",
    reference_phenotypes = FALSE,
    plot_distribution = FALSE
  )
  
  marker_names <- predicted_image@rowRanges@partitioning@NAMES
  marker_names <- marker_names[marker_names != "DAPI"]
  plots_list <- list()
  
  for (marker in marker_names) {
    tryCatch({
      marker_plot <- plot_cell_marker_levels(predicted_image, marker)
      plots_list[[marker]] <- marker_plot
    }, error = function(e) {
      message(paste("Error for marker", marker, ":", e$message))
    })
  }
  
  # Extract the required data
  extracted_data <- extract_data_from_spatial_experiment(plots_list)
  
  # Return the data along with the image name
  return(list(Image = selected_image, P = extracted_data$P, X = extracted_data$X, Y = extracted_data$Y, PD1 = extracted_data$PD1, CD8 = extracted_data$CD8, CD3 = extracted_data$CD3, LAG3 = extracted_data$LAG3, TIM3 = extracted_data$TIM3, CK = extracted_data$CK))
  
})

# Name entries
names(spiat_scattercount) <- sapply(spiat_scattercount, function(x) x$Image)

# Cat matching response and ID
spiat_scattercount <- pblapply(spiat_scattercount, function(x) {
  
  # Find matching ID based on image name
  matching_ID <- unique(eda_df$ID[eda_df$Image == x$Image])
  
  if (length(matching_ID) == 1) {
    x$ID <- matching_ID
    
    # Find matching response based on ID
    matching_Response <- responses$Response[responses$ID == matching_ID]
    
    if (length(matching_Response) == 1) {
      x$Response <- matching_Response
      
    } else {
      warning(paste("No matching Response found for ID:", matching_ID))
    }
    
  } else {
    warning(paste("No matching ID found for Image:", x$Image))
  }
  
  return(x)
})

# Generate standardized window of observation
allX_list <- lapply(spiat_scattercount, function(x) x$X)
allY_list <- lapply(spiat_scattercount, function(x) x$Y)
center_x <- median(unlist(allX_list))
center_y <- median(unlist(allY_list))
distances <- sqrt((unlist(allX_list) - center_x)^2 + (unlist(allY_list) - center_y)^2)
radius <- max(distances)
w <- disc(radius=radius, centre=c(center_x, center_y)) # Window w

immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# PROCESS - Creation of a PPP Hyperframe ----
first_img_data <- spiat_scattercount[[1]]
first_img_data_df <- as.data.frame(first_img_data)
first_ppp_list <- list()

for (phenotype in unique(first_img_data[["P"]])) {
  subset_data <- subset(first_img_data_df, P == phenotype)
  
  if(nrow(subset_data) > 0) {
    ppp <- ppp(x = subset_data$X, y = subset_data$Y,
               marks = data.frame(PD1 = subset_data$PD1, CD8 = subset_data$CD8,
                                  CD3 = subset_data$CD3, LAG3 = subset_data$LAG3,
                                  TIM3 = subset_data$TIM3, CK = subset_data$CK), 
               window = w, checkdup = FALSE)
    first_ppp_list[[phenotype]] <- ppp
  }
}

hf <- hyperframe(I = unique(first_img_data$Image), ppp = list(first_ppp_list), ID = unique(first_img_data$ID), Response = unique(first_img_data$Response))

for (i in 2:length(spiat_scattercount)) {
  
  img_data <- spiat_scattercount[[i]]
  img_data_df <- as.data.frame(img_data)
  
  ppp_list <- list()
  
  # Loop through phenotypes
  for (phenotype in unique(img_data[["P"]])) {
    subset_data <- subset(img_data_df, P == phenotype)
    
    if(nrow(subset_data) > 0) {
      ppp <- ppp(x = subset_data$X, y = subset_data$Y,
                 marks = data.frame(PD1 = subset_data$PD1, CD8 = subset_data$CD8,
                                    CD3 = subset_data$CD3, LAG3 = subset_data$LAG3,
                                    TIM3 = subset_data$TIM3, CK = subset_data$CK), 
                 window = w, checkdup = FALSE)
      ppp_list[[phenotype]] <- ppp
    }
  }
  
  # Append to hyperframe
  new_row <- hyperframe(I = unique(img_data$Image), ppp = list(ppp_list), ID = unique(img_data$ID), Response = unique(img_data$Response))
  hf <- rbind.hyperframe(hf, new_row)
}

plot(hf$ppp[[sample(1:90, 1)]]$CK, clipwin = w, use.marks = F, pch = 20) # Plots the ppp object of the 2oth image from the hyperframe.

# VISUALS - Intratumor variation (Heatmap) ----
p <- ggplot(eda_df[eda_df$ID == unique(eda_df$ID)[1],], aes(x = X, y = Y)) +
  stat_summary_2d(aes(z = LAG3, fill = after_stat(value)), fun = median, bins = 20) +
  facet_wrap(~ Image) +
  scale_fill_gradient(low="yellow", high="purple") +
  theme_minimal() +
  labs(fill = "Average PD1 Intensity")

ggplotly(p)

bin_n <- 90 # Select bin number for tiling ----
marker_of_interest <- "PD1" # Select marker of interest ----
# PROCESS - Continuous mark aggregation----
# Global binning based on entire data set
global_bin_width <- (max(eda_df$X) - min(eda_df$X)) / bin_n
global_min_x <- min(eda_df$X)
global_min_y <- min(eda_df$Y)

# Convert df to long format
eda_df_long <- eda_df %>%
  pivot_longer(cols = c(PD1, CD8, CD3, TIM3, LAG3, CK),
               names_to = "marker",
               values_to = "value")

# Store results for each patient
tiled_markerheatmap_xID <- list()

# Get unique IDs
unique_ids <- unique(eda_df$ID)

for(id in unique_ids) {
  
  # Creates bins for the current patient
  bin_data <- eda_df_long[eda_df_long$ID == id,] %>%  
    mutate(
      bin_x = floor((X - global_min_x) / global_bin_width),
      bin_y = floor((Y - global_min_y) / global_bin_width)
    )
  
  # Group by bin, Image, and marker
  resultado <- suppressMessages({
    bin_data %>%
      group_by(Image, marker, bin_x, bin_y) %>%
      summarise(avg_value = median(value, na.rm = TRUE))
  })
  
  # Aggregate across images by bin and marker for the current patient
  resultado <- suppressMessages({
    resultado %>%
      group_by(marker, bin_x, bin_y) %>%
      summarise(median_value_across_images = median(avg_value, na.rm = TRUE))
  })
  
  resultado <- resultado %>%
    mutate(
      X = global_min_x + (bin_x + 0.5) * global_bin_width,
      Y = global_min_y + (bin_y + 0.5) * global_bin_width,
    )
  
  # Store result
  tiled_markerheatmap_xID[[as.character(id)]] <- resultado
  print(paste("Processing patient", id, "data."))
}

tiled_markerheatmap_xResponse <- list()

for (response in unique(responses$Response)) {
  
  print(paste("Retrieving", response, "IDs."))
  
  # Get IDs associated with current response
  matching_ids <- as.character(responses$ID[responses$Response == response])
  
  print(paste("Matching IDs:", cat(matching_ids)))
  
  # Filter tiled_markerheatmap_xResponse by the matching IDs
  tiled_markerheatmap_xResponse[[response]] <- tiled_markerheatmap_xID[matching_ids]
}

# 1. Process the data
tiled_medmarkerint_xResponse <- list()

for (response in names(tiled_markerheatmap_xResponse)) {
  patient_list <- tiled_markerheatmap_xResponse[[response]]
  
  # Combine all data per response group
  combined_data <- bind_rows(patient_list) %>% 
    filter(marker == marker_of_interest)
  
  # Compute mean of each bin
  response_data <- combined_data %>%
    group_by(bin_x, bin_y, X, Y) %>%
    summarise(median_value = mean(median_value_across_images, na.rm = TRUE)) %>%
    mutate(marker = marker_of_interest) %>%
    ungroup()
  
  # Store in the results list
  tiled_medmarkerint_xResponse[[response]] <- response_data
}

# 2. Combine results for plotting
all_data <- bind_rows(lapply(names(tiled_medmarkerint_xResponse), function(response) {
  df <- tiled_medmarkerint_xResponse[[response]]
  df$Response <- response
  return(df)
}))

# VISUALS - Marker intensity per response group (Raster) ----
p_resultado <- ggplot(all_data, aes(x = X, y = Y, fill = log1p(median_value))) +
  geom_tile() +
  scico::scale_fill_scico(palette = "lajolla") +
  theme_minimal() +
  labs(title = paste("Median Intensity Of", marker_of_interest, "Across Response")) +
  facet_wrap(~ Response, ncol = 4)

ggplotly(p_resultado)

suppressWarnings({
  rm(bin_data, resultado, combined_data, response_data, tiled_markerheatmap_xID)
})

############################################ F4 | Distances -----
# EXTRACT - Cell distances ----
# PROCESS - Picking phenotypes ----
# Filtering phenotypes to chose
# Calculate mean proportion
mean_proportions <- merged_data %>%
  group_by(spiat_Phenotype) %>%
  summarise(mean_prop = mean(spiat_Prop, na.rm = TRUE),
            n = n_distinct(ID)) 

# Filter out phenotypes
all_patient_count <- length(unique(merged_data$ID))
filtered_phenotypes <- mean_proportions %>%
  filter(n == all_patient_count)

# Sort by mean proportion
sorted_phenotypes <- filtered_phenotypes %>%
  arrange(desc(mean_prop))

sorted_phenotypes <- sorted_phenotypes %>%
  mutate(mean_prop_percent = mean_prop * 100)

print(sorted_phenotypes)

# PROCESS - Main extraction ----
# Iterate over all images
distance_features <- list()

# Loop through each image
distance_features <- pblapply(SPIAT_tifs, function(selected_image) {
  
  # Filter data based on the current image
  selected_data <- subset(raw_measurements, Image %in% selected_image)
  selected_data <- selected_data %>% 
    select(-1) # Remove "Image" name column
  
  # Extract coordinates and intensity matrix
  x_cord <- t(selected_data)[1 , ] 
  y_cord <- t(selected_data)[2 , ]
  
  intensity_matrix <- t(selected_data[ , -c(1:2)])
  intensity_matrix <- as.data.frame(intensity_matrix)
  
  colnames(intensity_matrix) <- paste0("Cell_", as.character(1: ncol(intensity_matrix)))
  rownames(intensity_matrix) <- c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")
  
  # Extract phenotype
  phenotype = rep(NA, ncol(intensity_matrix))
  
  # Create spatial experiment (S4 object)
  general_format_image <- format_image_to_spe(
    format = "general",
    phenotypes = phenotype,
    intensity_matrix = intensity_matrix,
    coord_x = x_cord,
    coord_y = y_cord
  )
  
  # Predict cell phenotypes
  predicted_image <- predict_phenotypes(
    spe_object = general_format_image,
    thresholds = NULL,
    tumour_marker = "CK",
    baseline_markers = c("PD1", "CD8", "CD3", "TIM3", "LAG3"),
    nuclear_marker = "DAPI",
    reference_phenotypes = FALSE,
    plot_distribution = FALSE
  )
  
  marker_names <- predicted_image@rowRanges@partitioning@NAMES
  marker_names <- marker_names[marker_names != "DAPI"]

  # Cell distances
  distances <- calculate_pairwise_distances_between_celltypes(
    spe_object = predicted_image, 
    cell_types_of_interest = c("CD3", "CD8,CD3", "CD8,CD3,TIM3", "PD1,CD8,CD3", "CD3,TIM3", "PD1,CD3", "TIM3"),
    feature_colname = "Phenotype")
  summary_distances <- calculate_summary_distances_between_celltypes(distances)
  
  # Minimum cell distances
  min_dist <- calculate_minimum_distances_between_celltypes(
    spe_object = predicted_image, 
    cell_types_of_interest = c("CD3", "CD8", "CD8,CD3", "CD8,CD3,TIM3", "PD1,CD8,CD3", "CD3,TIM3", "PD1,CD3", "TIM3"),
    feature_colname = "Phenotype")
  min_summary_dist <- calculate_summary_distances_between_celltypes(min_dist)
  
  # Return the data along with the image name
  return(list(Image = selected_image, PDists = distances, PDistsSum = summary_distances, MinDists = min_dist, MindistsSum = min_summary_dist))
  
})

# Name list entries
names(distance_features) <- sapply(distance_features, function(x) x$Image)

## Matching response and ID
distance_features <- pblapply(distance_features, function(x) {
  
  # Find the matching ID from eda_df based on Image
  matching_ID <- unique(eda_df$ID[eda_df$Image == x$Image])
  
  if (length(matching_ID) == 1) {
    x$ID <- matching_ID
    
    # Find the matching Response from responses based on ID
    matching_Response <- responses$Response[responses$ID == matching_ID]
    
    if (length(matching_Response) == 1) {
      x$Response <- matching_Response
    } else {
      warning(paste("No matching Response found for ID:", matching_ID))
    }
  } else {
    warning(paste("No matching ID found for Image:", x$Image))
  }
  
  return(x)
})

# PROCESS - Main extraction (parallel) ----
# Iterate over all images
distance_features <- list()

# Set up the cluster
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)

# Export variables and libraries to the cluster
clusterExport(cl, c("SPIAT_tifs", "raw_measurements", "eda_df"))
clusterEvalQ(cl, {
  library(dplyr)
  library(pbapply)
  library(SPIAT)
})

# Run the parallelized loop
distance_features <- parLapply(cl, SPIAT_tifs, function(selected_image) {
  
  # Filter data based on the current image
  selected_data <- subset(raw_measurements, Image %in% selected_image)
  selected_data <- selected_data %>% 
    select(-1) # Remove "Image" name column
  
  # Extract coordinates and intensity matrix
  x_cord <- t(selected_data)[1 , ] 
  y_cord <- t(selected_data)[2 , ]
  
  intensity_matrix <- t(selected_data[ , -c(1:2)])
  intensity_matrix <- as.data.frame(intensity_matrix)
  
  colnames(intensity_matrix) <- paste0("Cell_", as.character(1: ncol(intensity_matrix)))
  rownames(intensity_matrix) <- c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")
  
  # Extract phenotype
  phenotype = rep(NA, ncol(intensity_matrix))
  
  # Create spatial experiment (S4 object)
  general_format_image <- format_image_to_spe(
    format = "general",
    phenotypes = phenotype,
    intensity_matrix = intensity_matrix,
    coord_x = x_cord,
    coord_y = y_cord
  )
  
  # Predict cell phenotypes
  predicted_image <- predict_phenotypes(
    spe_object = general_format_image,
    thresholds = NULL,
    tumour_marker = "CK",
    baseline_markers = c("PD1", "CD8", "CD3", "TIM3", "LAG3"),
    nuclear_marker = "DAPI",
    reference_phenotypes = FALSE,
    plot_distribution = FALSE
  )
  
  marker_names <- predicted_image@rowRanges@partitioning@NAMES
  marker_names <- marker_names[marker_names != "DAPI"]
  
  # Cell distances
  distances <- calculate_pairwise_distances_between_celltypes(
    spe_object = predicted_image, 
    cell_types_of_interest = c("CK", "CD3", "CD8", "CD8,CD3", "CD8,CD3,TIM3", "PD1,CD8,CD3", "CD3,TIM3", "PD1,CD3", "TIM3"),
    feature_colname = "Phenotype")
  summary_distances <- calculate_summary_distances_between_celltypes(distances)
  
  # Minimum cell distances
  min_dist <- calculate_minimum_distances_between_celltypes(
    spe_object = predicted_image, 
    cell_types_of_interest = c("CK", "CD3", "CD8", "CD8,CD3", "CD8,CD3,TIM3", "PD1,CD8,CD3", "CD3,TIM3", "PD1,CD3", "TIM3"),
    feature_colname = "Phenotype")
  min_summary_dist <- calculate_summary_distances_between_celltypes(min_dist)
  
  # Return the data along with the image name
  return(list(Image = selected_image, PDists = distances, PDistsSum = summary_distances, MinDists = min_dist, MindistsSum = min_summary_dist))
  
})

stopCluster(cl)


# PROCESS - Convert extractions to summary dataframes ----
format_summary <- function(df, id, response, image_name){
  df <- df %>%
    select(Pair, Mean, Min, Max, Median, Std.Dev) %>%
    mutate(Image = image_name,
           ID = id,
           Response = response)
  
  return(df)
}

pdsum_list <- list()
mdsum_list <- list()

for(i in 1:length(distance_features)){
  pdsum_list[[i]] <- format_summary(
    distance_features[[i]][['PDistsSum']],
    distance_features[[i]][['ID']],
    distance_features[[i]][['Response']],
    names(distance_features)[i]
  )
  
  mdsum_list[[i]] <- format_summary(
    distance_features[[i]][['MindistsSum']],
    distance_features[[i]][['ID']],
    distance_features[[i]][['Response']],
    names(distance_features)[i]
  )
}

# Combine all data frames into single data frames
pdsum_df <- bind_rows(pdsum_list)
mdsum_df <- bind_rows(mdsum_list)

# Rearrange columns
pdsum_df <- pdsum_df %>% select(Image, ID, Response, Pair, Min, Mean, Median, Max, Std.Dev)
mdsum_df <- mdsum_df %>% select(Image, ID, Response, Pair, Min, Mean, Median, Max, Std.Dev)

# PROCESS - Aggregate by patient ID ----
# Aggregate by ID
agg_pdsum_df <- pdsum_df %>%
  group_by(ID, Pair) %>%
  summarise(
    Min = min(Min, na.rm = TRUE),
    Max = max(Max, na.rm = TRUE),
    Mean = mean(Mean, na.rm = TRUE),
    Response = dplyr::first(as.character(Response))
  )

agg_mdsum_df <- mdsum_df %>%
  group_by(ID, Pair) %>%
  summarise(
    Min = min(Min, na.rm = TRUE),
    Max = max(Max, na.rm = TRUE),
    Mean = mean(Mean, na.rm = TRUE),
    Response = dplyr::first(as.character(Response))
  )

# STATS - KW ----
# Perform Kruskal-Wallis test on the Mean distance grouped by Response
# Splitting the data by 'Pair'
# KW for pairwise distances
list_of_dfs_by_pair <- split(agg_pdsum_df, agg_pdsum_df$Pair)

kruskal_results <- data.frame(Pair=character(), Kruskal_p_value=numeric())

# Kruskal-Wallis test for each pair
for (pair in names(list_of_dfs_by_pair)) {
  df_pair <- list_of_dfs_by_pair[[pair]]
  kruskal_test_result <- kruskal.test(Mean ~ Response, data = df_pair)
  kruskal_results <- rbind(kruskal_results, data.frame(Pair = pair, Kruskal_p_value = kruskal_test_result$p.value))
}

print(kruskal_results[order(kruskal_results$Kruskal_p_value), ])

# Dunn test for pairwise distances
# Filter out pairs with Kruskal p-value <= 0.05
significant_kruskal_results <- kruskal_results[kruskal_results$Kruskal_p_value <= 0.05, ]

dunn_test_results_df <- data.frame()

# Dunn Test for the filtered pairs
for (pair in significant_kruskal_results$Pair) {
  df_pair <- list_of_dfs_by_pair[[pair]]
  dunn_test_result <- dunnTest(Mean ~ Response, data = df_pair, method = "bh")
  
  dunn_df_temp <- data.frame(
    Pair = rep(pair, nrow(dunn_test_result$res)),
    Comparison = dunn_test_result$res$Comparison,
    Z = dunn_test_result$res$Z,
    P.unadj = dunn_test_result$res$`P.unadj`,
    P.adj = dunn_test_result$res$`P.adj`
  )
  
  # Append to main results data frame
  dunn_test_results_df <- rbind(dunn_test_results_df, dunn_df_temp)
}

# Filter and sort the results based on adjusted p-value
filtered_sorted_dunn_results <- dunn_test_results_df %>% 
  arrange(P.adj)

# Print the filtered and sorted results
print(filtered_sorted_dunn_results)

# KW for min distances
# Splitting the data by 'Pair'
list_of_dfs_by_pair <- split(agg_mdsum_df, agg_mdsum_df$Pair)
kruskal_results <- data.frame(Pair=character(), Kruskal_p_value=numeric())

# Kruskal-Wallis test for each pair
for (pair in names(list_of_dfs_by_pair)) {
  df_pair <- list_of_dfs_by_pair[[pair]]
  kruskal_test_result <- kruskal.test(Mean ~ Response, data = df_pair)
  kruskal_results <- rbind(kruskal_results, data.frame(Pair = pair, Kruskal_p_value = kruskal_test_result$p.value))
}

print(kruskal_results[order(kruskal_results$Kruskal_p_value), ])

# Dunn test for pairwise distances
# Filter out pairs with Kruskal p-value <= 0.05
significant_kruskal_results <- kruskal_results[kruskal_results$Kruskal_p_value <= 0.05, ]

dunn_test_results_df <- data.frame()

# Dunn Test for the filtered pairs
for (pair in significant_kruskal_results$Pair) {
  df_pair <- list_of_dfs_by_pair[[pair]]
  dunn_test_result <- dunnTest(Mean ~ Response, data = df_pair, method = "bh")
  
  dunn_df_temp <- data.frame(
    Pair = rep(pair, nrow(dunn_test_result$res)),
    Comparison = dunn_test_result$res$Comparison,
    Z = dunn_test_result$res$Z,
    P.unadj = dunn_test_result$res$`P.unadj`,
    P.adj = dunn_test_result$res$`P.adj`
  )
  
  # Append to main results data frame
  dunn_test_results_df <- rbind(dunn_test_results_df, dunn_df_temp)
}

# Filter and sort the results based on adjusted p-value
filtered_sorted_dunn_results <- dunn_test_results_df %>% 
  arrange(P.adj)

# Print the filtered and sorted results
print(filtered_sorted_dunn_results)


# Session info ----