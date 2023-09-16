# LIBRARIES ----
library(gridExtra)
library(BBmisc)
library(stats)
library(reshape2)
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

############################################ DATA READING AND FORMATTING ---------------------------------------------
## 1.1. READING IN DATA ----
raw_measurements <- read.table("Data/raw_measurements.tsv", header = TRUE, sep = "\t")
responses <- read.csv("Data/Responses.csv", header = TRUE, sep = ",")

## 1.2. FORMATTING FOR EDA ----
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

# Delete groups that are not evaluable
responses$ID <- as.numeric(responses$ID)
eda_df <- eda_df[eda_df$Response != "Not evaluable/NE", ]
responses <- responses[responses$Response != "Not evaluable/NE", ]

# Finish formatting
eda_df <- eda_df[, -ncol(eda_df)]
new_order <- c("Image", "X", "Y", "ID", "DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")
eda_df <- eda_df[, new_order]


# Counting the number of cells per TIFF and binding it to the response dataset
cell_counts <- table(eda_df$ID)
cell_counts<- as.data.frame(cell_counts)
colnames(cell_counts) <- c("ID", "Cell_Count")
cell_count_vector <- cell_counts$Cell_Count[match(responses$ID, cell_counts$ID)]
responses$total_cells <- cell_count_vector

# Final retouches
# Create a convenient list with all images
SPIAT_tifs <- unique(eda_df$Image)
SPIAT_px <- unique(eda_df$ID)

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

# Limpieza de la memoria
responses$Response <- droplevels(responses$Response, "Not evaluable/NE")
levels(responses$Response) <- c("CR", "PR", "PD", "SD")

rm(list = c("cell_counts", "NAs", "cell_count_vector", "images_with_NAs", "new_order", "response_vector"))

############################################ SPIAT GLIMPSE (SPATIAL IMAGE ANALYSIS OF TISSUES) ----
## 4.1 CREATE SPATIAL EXPERIMENT ----
# Select image subset
selected_image <- SPIAT_tifs[1]
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

#my_phenotype = final_df$Phenotype[final_df$Image == selected_image] 
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

## 4.2. QUALITY CONTROL AND VISUALISATION ----
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
plot_cell_categories(
  predicted_image, 
  c("Undefined", "Tumour", "PanExhaust", "T-cell", "Active T-cell", "ExhaustedTIM3", "ExhaustedPD1", "ExhaustedLAG3"),
  c("grey", "yellow","cyan","red", "darkred","purple", "purple", "purple"), 
  "Cell.Type",
  layered = TRUE)

# 3D surface plot
# Warning: Uses rownames(SummarizedExperiment::assay(spe_object)) to compare markers so only works with original markers and not phenotypes
marker_names <- c("PD1",  "LAG3", "CK")
marker_surface_plot_stack(predicted_image, num_splits=50, marker_names)

## 4.3. BASIC ANALYSES WITH SPIAT ----
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

## 4.4. QUANTIFYING CELL COLOCALISATION WITH SPIAT ----
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

## 4.5. SPATIAL HETEROGENEITY ----
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


## 4.6. CHARACTERISING TISSUE STRUCTURE ----
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

## 4.7. CELLULAR NEIGHBORHOOD ----
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

########################################### FEATURE 1 - PHENOTYPES -----
# EXTRACT - Predicted phenotypes ----
## PROCESS - Extraction and re-ordering ----
spiat_predicted_phenotypes <- list()

# Create a custom function for extracting the necessary details from the SpatialExperiment object
extract_data_from_spatial_experiment <- function(spe_object) {
  spatial_coords <- spe_object@int_colData@listData$spatialCoords
  phenotype <- spe_object@colData@listData$Phenotype
  
  return(list(spatial_coords = spatial_coords, phenotype = phenotype))
}

# Use pblapply (from pbapply package) to loop through each image in SPIAT_tifs
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

# Name each list entry with its corresponding image name
names(spiat_predicted_phenotypes) <- sapply(spiat_predicted_phenotypes, function(x) x$image_name)

## OBTAIN PREDICTED PHENOTYPE COUNTS
# Step 1: Extract phenotypes and store them in a vector
phenotype_vector <- unlist(lapply(spiat_predicted_phenotypes, function(x) x$phenotype))

# Step 2: Get unique phenotypes across all images
unique_phenotypes <- unique(phenotype_vector)

# Step 3: Create a dataframe with image names and columns named based on unique phenotypes
spiat_pheno_counts <- data.frame(image_name = names(spiat_predicted_phenotypes),
                                 matrix(NA, ncol = length(unique_phenotypes), nrow = length(spiat_predicted_phenotypes)))
colnames(spiat_pheno_counts) <- c("Image", unique_phenotypes)

# Step 4 (Continuation): Iterate through images, count phenotypes, and update the dataframe
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

# Step 5: Re-ordering columns
# Add the ID and response columns
spiat_pheno_counts$ID <- eda_df$ID[match(spiat_pheno_counts$Image, eda_df$Image)]
spiat_pheno_counts$Response <- responses$Response[match(spiat_pheno_counts$ID, responses$ID)]
current_cols <- setdiff(colnames(spiat_pheno_counts), c("Image","ID", "Response", "None", "CK"))

# Assign primary weights based on prefixes
get_primary_weight <- function(col) {
  if (startsWith(col, "CD3")) return(1)
  if (startsWith(col, "CD8")) return(2)
  if (startsWith(col, "PD1")) return(3)
  if (startsWith(col, "LAG3")) return(4)
  if (startsWith(col, "TIM3")) return(5)
  
  # Return a high weight for unmatched columns so they appear at the end
  return(100)
}

# Assign secondary weight based on the number of immune markers
get_secondary_weight <- function(col) {
  # Count number of commas and add 1 to determine the number of markers
  return(1 + sum(nchar(gsub("[^,]", "", col))))
}

# Calculate weights for columns
primary_weights <- sapply(current_cols, get_primary_weight)
secondary_weights <- sapply(current_cols, get_secondary_weight)

# Combine primary and secondary weights to get the order
# Use primary weight * a big number (like 1000) + secondary weight
combined_weights <- primary_weights * 1000 + secondary_weights

# Sort columns based on combined weights
ordered_cols <- current_cols[order(combined_weights)]

# Add the "Image", "None", and "CK" columns at the start
final_order <- c("Image","ID", "Response", "None", "CK", ordered_cols)

# Re-arrange columns of the dataframe based on final_order
spiat_pheno_counts <- spiat_pheno_counts[, final_order]

## PROCESS - Convert to proportions ----
# Step 1: Reshape the data
spiat_long_format <- spiat_pheno_counts %>%
  select(-Image) %>%
  gather(key = 'spiat_Phenotype', value = 'spiat_Count', -ID, -Response)

# Group the data by ID and summarize the counts
spiat_total_cells_per_ID <- spiat_long_format %>%
  group_by(ID) %>%
  summarise(spiat_TotalCells = sum(spiat_Count, na.rm = TRUE))

# Long Data with Proportions
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

# So that abscent markers are not ignored but counted as 0
spiat_pheno_counts[,4:67] <- apply(spiat_pheno_counts[,4:67], 2, function(x) ifelse(is.na(x), 0, x))

# Aggregation by summation (since phenotypes are count, we test proportions instead of counts to lessen the effects of repeated paired measures)
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

# Initialize data frames
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

# Step 1: Perform Levene's test
for (pheno in phenos) {
  lev_test <- leveneTest(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == pheno, ])$`Pr(>F)`[1]
  spiat_pheno_levene <- rbind(spiat_pheno_levene, data.frame(Phenotype = pheno, levene = lev_test))
}

# Step 2: Perform Kruskal-Wallis test for Levene's test (p >= 0.05)
for (pheno in phenos) {
  p_val <- kruskal.test(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == pheno, ])$p.value
  spiat_phenoKWfdr <- rbind(spiat_phenoKWfdr, data.frame(Phenotype = pheno, p.value = p_val))
}

# Adjust p-values
spiat_phenoKWfdr$adj.p <- p.adjust(spiat_phenoKWfdr$p.value, method = "fdr")

# Step 3: Perform Dunn's test for Kruskal-Wallis p < 0.05
for (pheno in spiat_phenoKWfdr$Phenotype[spiat_phenoKWfdr$p.value < 0.05]) {
  dunn_results <- dunnTest(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == pheno, ], method = "bonferroni")
  dunn_df <- data.frame(
    Phenotype = rep(pheno, nrow(dunn_results$res)),
    Comparison = dunn_results$res$Comparison,
    Z = dunn_results$res$Z,
    P.unadj = dunn_results$res$`P.unadj`,
    P.adj = dunn_results$res$`P.adj`
  )
  spiat_pheno_dunn <- rbind(spiat_pheno_dunn, dunn_df)
}

# Cleaning up
rm(col_index, current_cols, i, image, image_name, image_phenotype_counts, ordered_cols, phenotype, phenotype_vector, phenotypes_in_image, row_index, 
   unique_phenotypes, primary_weights, secondary_weights, combined_weights, get_primary_weight, get_secondary_weight, final_order, spiat_pheno_temp, 
   marker, markers, p_val, dunn_df, dunn_results, lev_test, spiat_long_format, spiat_phenoKW, spiat_total_cells_per_ID)

## RESULTS - Statistics (Aqui falta point plot con p valores) ----
print(
  spiat_phenoKWfdr %>% 
    filter(p.value < 0.05) %>% 
    arrange(p.value)
)

print(
  spiat_pheno_dunn %>% 
    filter(P.adj < 0.05) %>% 
    arrange(P.adj)
)

## VISUALS - Horizontal stacked barplots (PATIENTS) ----
# Generate the base plot
# Merge the datasets
merged_data <- merge(spiat_long_with_props, responses[, c("ID", "Response")], by = "ID")

# Order IDs based on Response categories
ordered_ids <- merged_data %>% 
  arrange(factor(Response, levels = c("CR", "PR", "SD", "PD"))) %>%
  pull(ID) %>%
  unique()

# Relevel factor levels of ID
merged_data$ID <- factor(merged_data$ID, levels = ordered_ids)

# Plot
gg <- ggplot(merged_data, aes(x=factor(ID), y=spiat_Prop, fill=spiat_Phenotype)) +
  geom_bar(stat="identity", position="fill") + 
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)))

gg

## VISUALS - Horizontal stacked barplots (RESPONSES) ----
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

########################################### FEATURE 2 - SCATTERPLOT (S/T) -----
# EXTRACT - Scatterplot data ----
# Iterate over all images
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

# Name list entries
names(spiat_scattercount) <- sapply(spiat_scattercount, function(x) x$Image)

### Matching response and ID
# Update each element in spiat_scattercount
spiat_scattercount <- pblapply(spiat_scattercount, function(x) {
  
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

nq <- 7 # Select n of quadrats (<90) ----
# PROCESSING - Creates PPP, density maps and quadrats (ALL IMAGES) ---- 
allY_list <- lapply(spiat_scattercount, function(x) x$Y)
center_y <- median(unlist(allY_list))

allX_list <- lapply(spiat_scattercount, function(x) x$X)
center_x <- median(unlist(allX_list))

# Calculate distances from center to each point
distances <- sqrt((unlist(allX_list) - center_x)^2 + (unlist(allY_list) - center_y)^2)
max_radius <- max(distances)

w <- disc(radius=max_radius, centre=c(center_x, center_y))

immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Create lists for the results
pointpatterns <- vector("list", length = length(names(spiat_scattercount)))
density_maps <- vector("list", length = length(names(spiat_scattercount)))
qtest_results <- vector("list", length = length(names(spiat_scattercount)))
quadrats <- vector("list", length = length(names(spiat_scattercount)))

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

# Most important line, does all plots and tests
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
    
    # Estimate the density
    D <- density.ppp(P,
                     sigma = 85,
                     edge = TRUE,
                     diggle = TRUE)
    
    # Check if marker is present and if the number of counts is above the threshold
    if (nrow(subset_data) >= nq) {
      
      # Calculate quadrat test
      q_test <- quadrat.test.ppp(P,
                                 alternative = "clustered",
                                 method = "M",
                                 conditional = T,
                                 nsim = 999)
      
      
    } else {
      
      # Check if ID and Response are available for the image_name
      id_available <- unlist(responses$ID[sapply(responses$image_names, function(x) any(x == image_name))])
      response_available <- unlist(responses$Response[sapply(responses$image_names, function(x) any(x == image_name))])
      
      if (length(id_available) > 0 && length(response_available) > 0) {
        # Store the image name, ID, marker, and response group in the skipped_images list
        skipped_images <- rbind(skipped_images, data.frame(ImageName = image_name,
                                                           ID = id_available,
                                                           Marker = marker,
                                                           Response = response_available,
                                                           stringsAsFactors = FALSE))
      }
    }
    
    # Store each image's marker and quadrat test
    marker_ppps[[marker]] <- P
    marker_density_maps[[marker]] <- D
    q_count[[marker]] <- Q
    qtest_images[[marker]] <- q_test
    
  }
  
  # Store each image and results
  pointpatterns[[image_name]] <- marker_ppps
  density_maps[[image_name]] <- marker_density_maps
  quadrats[[image_name]] <- q_count
  qtest_results[[image_name]] <- qtest_images
  
  print(paste("Process complete for image:", image_name))
  
}

print(skipped_images)

# PROCESSING - Aggregation by ID ----
# (IMAGES) Transform each bin to relative proportions
transformed_quadrats <- list()

# Image's bins to proportions relative to total patient count
for (image_name in names(quadrats)) {
  
  # Extract current image data
  image_data <- quadrats[[image_name]]
  
  # Create empty list for storage
  transformed_image_data <- list()
  
  # Loop through each marker in the image data
  for (marker in names(image_data)) {
    
    # Extract the counts for the current marker
    marker_counts <- image_data[[marker]]
    
    # Calculate the total count for the current marker in the current image
    greped_ID <- which(unlist(sapply(responses$image_names, function(x) image_name %in% x)))
    total_count <- responses$total_cells[greped_ID]
    
    # Transform the counts into relative proportions
    relative_proportions <- round((marker_counts / total_count)*100, 3)
    
    # Store the transformed data in the list
    transformed_image_data[[marker]] <- relative_proportions
  }
  
  # Store the transformed data for the current image in the main list
  transformed_quadrats[[image_name]] <- transformed_image_data
}

# (PATIENT ID) Aggregate values by average based on: marker, and bin index
grouped_quadrats <- list()
for (patient_id in unique(responses$ID)) {
  
  # Find the image names corresponding to the current patient ID
  image_names <- unlist(responses$image_names[responses$ID == patient_id])
  
  # Initialize an empty list to store the grouped data for the current patient
  grouped_patient_data <- list()
  
  # Loop through each image name for the current patient
  for (image_name in image_names) {
    
    # Extract the data for the current image
    image_data <- transformed_quadrats[[image_name]]
    
    # Loop through each marker in the image data
    for (marker in names(image_data)) {
      
      # If the marker list doesn't exist, initialize it
      if (is.null(grouped_patient_data[[marker]])) {
        grouped_patient_data[[marker]] <- list()
      }
      
      # Append the counts for the current marker and image to the list
      grouped_patient_data[[marker]] <- append(grouped_patient_data[[marker]], list(image_data[[marker]]))
    }
  }
  
  # Now, aggregate each marker by the mean
  for (marker in names(grouped_patient_data)) {
    
    # Combine all the counts for this marker into a matrix
    marker_matrix <- do.call(rbind, grouped_patient_data[[marker]])
    
    # Replace NaN with 0
    marker_matrix[is.nan(marker_matrix)] <- 0
    
    # Calculate the mean for each bin (column)
    mean_values <- apply(marker_matrix, 2, mean, na.rm = TRUE)
    
    # Create a new quadratcount object with the mean values but retaining the owin object and tess attributes
    new_quadrat <- grouped_patient_data[[marker]][[1]]
    attributes(new_quadrat)$quadratcount <- mean_values
    
    # Store the new quadrat object in the list
    grouped_patient_data[[marker]] <- new_quadrat
  }
  
  # Store the grouped and aggregated data for the current patient in the main list
  grouped_quadrats[[as.character(patient_id)]] <- grouped_patient_data
}

# PROCESSING - Aggregation by response group ----
grouped_by_response <- list()

for (response_group in unique(responses$Response)) {
  
  # Initialize an empty list to store the grouped data for the current response group
  grouped_response_data <- list()
  
  # Find the patient IDs corresponding to the current response group
  patient_ids <- unique(responses$ID[responses$Response == response_group])
  
  # Loop through each patient ID for the current response group
  for (patient_id in patient_ids) {
    
    # Extract the data for the current patient
    patient_data <- grouped_quadrats[[as.character(patient_id)]]
    
    # Loop through each marker in the patient data
    for (marker in names(patient_data)) {
      
      # Append the counts for the current marker and patient to the list
      grouped_response_data[[marker]] <- append(grouped_response_data[[marker]], list(patient_data[[marker]]))
    }
  }
  
  # Now, aggregate each marker by the median
  for (marker in names(grouped_response_data)) {
    
    # Combine all the counts for this marker into a matrix
    marker_matrix <- do.call(rbind, lapply(grouped_response_data[[marker]], function(x) attributes(x)$quadratcount))
    
    # Calculate the median for each quadrant
    mean_val <- apply(marker_matrix, 2, mean, na.rm = TRUE)
    
    # Create a new quadratcount object with the median values but retaining the owin object and tess attributes
    new_quadrat <- grouped_response_data[[marker]][[1]]
    attributes(new_quadrat)$quadratcount <- mean_val
    
    # Store the new quadrat object in the list
    grouped_response_data[[marker]] <- new_quadrat
  }
  
  # Store the grouped and aggregated data for the current response group in the main list
  grouped_by_response[[response_group]] <- grouped_response_data
}

selected_image <- SPIAT_tifs[1] #----
# VISUALS - Sampler density maps ----
# Define the order of response groups and immune markers
response_order <- c("SD", "PD", "CR", "PR")
marker_order <- c("CD8", "PD1", "CK", "TIM3", "CD3", "LAG3")

# Initialize an empty list to store the point patterns, densities, and contours
group_ppp <- list()
group_densities <- list()

# Loop through each unique response group
for (response_group in response_order) {
  
  # Filter the responses data frame to only include patients in the current response group
  group_images <- unlist(responses$image_names[responses$Response == response_group])
  
  # Sample an image from the current response group (if available)
  sampled_images <- sample(group_images, size = min(1, length(group_images)), replace = FALSE)
  
  # Loop through each sampled image name
  for (image_name in sampled_images) {
    
    # Extract the data for the current image
    gppp <- pointpatterns[[image_name]]
    dmap <- density_maps[[image_name]]
    
    # Store the point patterns, densities, and contours in the list, indexed by response group and image name
    index_name <- paste(response_group, image_name, sep = "_")
    group_ppp[[index_name]] <- gppp
    group_densities[[index_name]] <- dmap
  }
}

# Initialize a list to store min and max density for each marker
marker_min_max = list()

# Loop through each image to get the min and max density for each marker
for (image_name in names(group_densities)) {
  
  for (marker in names(group_densities[[image_name]])) {
    
    current_density_values = group_densities[[image_name]][[marker]]$v
    
    # Update the min and max if the current density values are more extreme
    marker_min_max[[marker]]$min = min(marker_min_max[[marker]]$min, min(current_density_values, na.rm = TRUE))
    marker_min_max[[marker]]$max = max(marker_min_max[[marker]]$max, max(current_density_values, na.rm = TRUE))
  }
}

# Your plotting code remains largely the same
par(mfrow = c(6, length(group_ppp)), mar = c(1,1,1,1), bty = "n")

for (marker in marker_order) {
  for (response_group in response_order) {
    for (image_name in names(group_ppp)) {
      
      if (startsWith(image_name, response_group)) {
        
        f_ppp <- group_ppp[[image_name]][[marker]]
        f_densities <- group_densities[[image_name]][[marker]]
        
        # Use zlim parameter from marker_min_max
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
        
        # Add titles and labels
        if (image_name == names(group_ppp)[1]) {
          mtext(marker, side = 2, line = 0, cex = 1)
        }
      }
    }
  }
}

# VISUALS - Quadrat plots ----
par(mfrow = c(2, 3))
for (marker in immune_markers) {
  
  # Retrieve the density maps for the selected image
  Q <- quadrats[[selected_image]][[marker]]
  
  plot(Q, main = print(paste(marker, "-", sum(Q))))
  
}
# VISUALS - Quadrat test plots ----
par(mfrow = c(2, 3))
for (marker in immune_markers) {
  plot(qtest_results[[selected_image]][[marker]], main = marker)
}

# VISUALS - Quadrat proportion plots ----
par(mfrow = c(2, 3))
for (marker in immune_markers) {
  
  # Retrieve the density maps for the selected image
  Q <- transformed_quadrats[[selected_image]][[marker]]
  
  plot(Q, main = print(paste(marker, "-", round(sum(Q), 2))))
  
}

# VISUALS - Results ----
# Define the order of response groups and immune markers
response_order <- c("SD", "PD", "CR", "PR")
marker_order <- c("LAG3", "CD3", "TIM3", "CK", "PD1","CD8")

# Calculate min and max for each marker
marker_min_max <- lapply(marker_order, function(marker) {
  counts <- unlist(lapply(response_order, function(response_group) {
    Q <- grouped_by_response[[response_group]][[marker]]
    as.vector(round(attributes(Q)$quadratcount, 3))
  }))
  c(min = min(counts), max = max(counts))
})

names(marker_min_max) <- marker_order

# Initialize plotting layout
par(mfrow = c(6, 4), mar = c(0, 2, 2, 2), oma = c(0, 0, 0, 3))  # Added room on the right for the colorbar

# Loop through each immune marker
for (marker in marker_order) {
  
  # Generate the color palette
  marker_min <- marker_min_max[[marker]]['min']
  marker_max <- marker_min_max[[marker]]['max']
  cols <- colorRampPalette(c("white", "cyan", "limegreen", "yellow", "firebrick2"))(100)
  
  # Loop through each response group
  for (response_group in response_order) {
    
    # Retrieve the density maps
    Q <- grouped_by_response[[response_group]][[marker]]
    
    # Plot the window
    plot(attr(Q, "tess")$window, main = NULL, box = FALSE)
    
    # Loop to draw and color each tile
    counts <- as.vector(round(attributes(Q)$quadratcount, 3))
    tiles <- attr(Q, "tess")$tiles
    for (i in 1:length(tiles)) {
      tile <- tiles[[i]]
      count <- counts[i]
      col_idx <- findInterval(count, seq(marker_min, marker_max, length.out = length(cols)))
      plot(tile, add = TRUE, col = cols[col_idx], border = "black")
    }
    
    if (response_group == "SD") {
      mtext(marker, side = 2, line = 0, cex = 1)
      
    } else if (response_group == "PR") {
      par(xpd = NA)  # Allow plotting in the outer margins
      legend("right", inset = c(-0.1, 0),  # Adjust inset to move legend
             legend = rev(round(seq(marker_min, marker_max, length.out = 5), 3)),
             fill = rev(colorRampPalette(c("white", "cyan", "limegreen", "yellow", "firebrick2"))(5)), 
             bty = "n", cex = 0.8, title = NULL, box.lty=0)
      par(xpd = FALSE)  # Restore clipping
    }
  }
}

# Adjust column_positions to better align text
column_positions <- seq(1/8, 1, by = 1/4)
mtext(response_order, side = 3, line = -2, outer = TRUE, at = column_positions, cex = 1)

# STATS - Kruskall Wallis ----
# Initialize an empty data frame to store the summary statistics
summary_table <- data.frame()

# Loop through each marker
for (marker in names(grouped_quadrats[[1]])) {
  
  # Initialize an empty list to store the data for Kruskal-Wallis and Levene's tests
  all_data <- list()
  
  # Loop through each unique response
  for (response_group in unique(responses$Response)) {
    
    # Get IDs belonging to the current response group
    patient_ids <- unique(responses$ID[responses$Response == response_group])
    
    # Extract the data for the current marker and response group
    marker_data <- lapply(patient_ids, function(id) grouped_quadrats[[as.character(id)]][[marker]])
    
    # Assuming each element in marker_data is a vector, bind them into a matrix
    marker_matrix <- do.call(rbind, marker_data)
    
    # Store the data for later use in Kruskal-Wallis test
    all_data[[as.character(response_group)]] <- marker_matrix
  }
  
  # Combine all the matrices into one large matrix
  combined_matrix <- do.call(rbind, all_data)
  
  # Create a group vector for Kruskal-Wallis test
  group_vector <- as.vector(rep(names(all_data), sapply(all_data, nrow)))
  
  # Perform Kruskal-Wallis and Levene's tests for each bin (column in combined_matrix)
  for (bin_index in 1:ncol(combined_matrix)) {
    
    # Data for the current bin
    bin_data <- combined_matrix[, bin_index]
    
    # Perform Levene's test
    levene_result <- leveneTest(y = bin_data, group = group_vector)$`Pr(>F)`[1]
    
    # Perform Kruskal-Wallis test
    kw_result <- kruskal.test(bin_data ~ group_vector)$p.value
    
    # Create a temporary data frame to store the results
    temp_df <- data.frame(
      Marker = marker,
      Bin = bin_index,
      Levene_p_value = levene_result,
      KW_p_value = kw_result
    )
    
    # Append to the summary table
    summary_table <- rbind(summary_table, temp_df)
  }
}

# Perform FDR correction for Kruskal-Wallis p-values
summary_table <- summary_table %>% 
  mutate(FDR_adj_p_value = p.adjust(KW_p_value, method = "fdr"))

# Sort the summary table by FDR-adjusted p-values
summary_table <- summary_table %>% 
  arrange(FDR_adj_p_value)

# Print the summary table
print(summary_table)
# STATS - Dunn test ----
# Initialize an empty data frame to store the Dunn's test results
for (i in 1:nrow(summary_table)) {
  if (!is.na(summary_table$KW_p_value[i]) && summary_table$KW_p_value[i] <= 0.05) {
    marker_bin <- summary_table[i, c("Marker", "Bin"), drop = FALSE]
    
    # Extract the marker and bin index
    marker <- marker_bin$Marker
    bin_index <- marker_bin$Bin
    
    # Extract the data for the current bin
    bin_data <- combined_matrix[, bin_index]
    
    # Temp df for data and bin
    temp_data <- data.frame(Value = bin_data, Response = as.factor(group_vector))
    
    # Perform Dunn's test
    dunn_result <- dunnTest(Value ~ Response, data = temp_data, method = "bh")
    
    # Temp df to store results
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
    
  } else {
    if (is.na(summary_table$KW_p_value[i])) {
      print(paste("NA encountered in:", i))
    } else if (summary_table$KW_p_value[i] >= 0.05) {
      print(paste("P-value >= 0.05 encountered in:", i))
    }
  }
}

# Sort by FDR adj
dunn_results <- dunn_results %>% 
  arrange(Adj_P)

print(dunn_results)

# STATS - Optimal number of bins ----
# Function to calculate RSS for a given number of bins
calculate_RSS <- function(df, n_bins) {
  # Initialize counts for each bin
  bin_counts <- table(df$bin_x, df$bin_y, df$Freq)
  
  # Convert table to data frame
  bin_counts_df <- as.data.frame(as.table(bin_counts))
  if(ncol(bin_counts_df) == length(c("bin_x", "bin_y", "Freq"))){
    names(bin_counts_df) <- c("bin_x", "bin_y", "Freq")
  } else {
    print("Eidolon help me please")
  }
  names(bin_counts_df) <- c("bin_x", "bin_y", "Freq")
  
  # Make sure bin_counts includes zeros for empty bins
  all_possible_bins <- expand.grid(bin_x = unique(df$bin_x), bin_y = unique(df$bin_y))
  bin_counts_df <- merge(all_possible_bins, bin_counts_df, all.x = TRUE, by = c("bin_x", "bin_y"))
  bin_counts_df[is.na(bin_counts_df)] <- 0
  
  # Calculate the observed proportions
  total_points <- sum(bin_counts_df$Freq)
  observed_proportions <- bin_counts_df$Freq / total_points
  if(length(observed_proportions) == length(expected_proportions)) {
    residuals <- observed_proportions - expected_proportions
  } else {
    print("This wrong too gah")
  }
    
  
  # Calculate the expected counts and proportions
  base_count <- floor(total_points / n_bins)
  extra_counts = total_points %% n_bins
  
  # Distribute the "extra counts" randomly
  set.seed(123) # for reproducibility
  extra_indices <- sample(1:n_bins, extra_counts)
  
  expected_counts = rep(base_count, n_bins)
  expected_counts[extra_indices] = expected_counts[extra_indices] + 1
  
  expected_proportions = expected_counts / total_points
  
  # Calculate residuals for each bin
  residuals <- observed_proportions - expected_proportions
  
  # Calculate RSS
  RSS <- sum(residuals^2)
  
  return(RSS)
}

# Calculate the center and largest radius of all data
all_data_center_x <- mean(eda_df$X)
all_data_center_y <- mean(eda_df$Y)
largest_radius <- max(sqrt((eda_df$X - all_data_center_x)^2 + (eda_df$Y - all_data_center_y)^2))

# Make the observation window 5% larger
window_radius <- largest_radius * 1.05

# Define the radius of the influence region (this should be the same as your optimal bandwidth)
radius <- 40  # Replace with your optimal bandwidth if different

markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Initialize the bin_metrics dataframe
bin_metrics <- data.frame(I = character(), M = character(), stringsAsFactors = FALSE)

# Loop through each number of bins from 15 to 100
for (n_bins in 15:100) {
  
  # Temporary dataframe to store the current bin_metrics
  temp_bin_metrics <- data.frame(I = character(), M = character(), RSS = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each image and marker
  for (image in names(spiat_scattercount)) {
    for (marker in markers) {
      
      # Subset data for current image and marker
      image_data <- spiat_scattercount[[image]]
      image_data_df <- as.data.frame(image_data)
      sub_df <- image_data_df[grepl(marker, image_data_df$P, fixed = TRUE), ]
      
      # # Filter data points that fall within the circular window
      # sub_df <- sub_df %>%
      #   filter(sqrt((X - all_data_center_x)^2 + (Y - all_data_center_y)^2) <= window_radius)
      
      # Bin the data within the circular window
      sub_df <- sub_df %>%
        mutate(
          bin_x = floor((X - all_data_center_x) / radius),
          bin_y = floor((Y - all_data_center_y) / radius)
        )
      
      # Calculate RSS
      RSS <- calculate_RSS(sub_df, n_bins)
      
      # Append to the temporary dataframe
      temp_bin_metrics <- rbind(temp_bin_metrics, data.frame(I = image, M = marker, RSS = RSS))
      print(marker)
    }
    print(image)
  }
  
  # Rename the RSS column to the current n_bins
  names(temp_bin_metrics)[names(temp_bin_metrics) == "RSS"] <- as.character(n_bins)
  
  # Merge with the existing bin_metrics dataframe
  if (nrow(bin_metrics) == 0) {
    bin_metrics <- temp_bin_metrics
  } else {
    bin_metrics <- merge(bin_metrics, temp_bin_metrics, by = c("I", "M"), all = TRUE)
  }
  
  print(paste(n_bins, "bins test completed."))
}

# STATS - Bandwidth testing ----
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

# VISUALS - Bandwidth tests ----
# Plot the distribution of optimal bandwidths
par(mfrow = c(5, 1))
hist(sapply(bandwidth_results, function(x) x$bw_diggle), main = NULL, xlab = 'Bandwidth', breaks = 20)
hist(sapply(bandwidth_results, function(x) x$ppl), main = NULL, xlab = 'Bandwidth', breaks = 20)
hist(sapply(bandwidth_results, function(x) x$cvl), main = NULL, xlab = 'Bandwidth', breaks = 20)
hist(sapply(bandwidth_results, function(x) x$adpt), main = NULL, xlab = 'Bandwidth', breaks = 20)
hist(unlist(sapply(bandwidth_results, function(x) x$abraham)), main = NULL, xlab = 'Bandwidth', breaks = 20)

########################################### FEATURE 3 - MARKER HEATMAP (S/T) -----
# EXTRACT - MARKER INTENSITY HEATMAP (SPATSTAT continuous?) ----
# VISUALS - 1 patient variation heatmap ----
p <- ggplot(eda_df[eda_df$ID == unique(eda_df$ID)[1],], aes(x = X, y = Y)) +
  stat_summary_2d(aes(z = LAG3, fill = after_stat(value)), fun = median, bins = 20) +
  facet_wrap(~ Image) +
  scale_fill_gradient(low="yellow", high="purple") +
  theme_minimal() +
  labs(fill = "Average PD1 Intensity")

ggplotly(p)

bin_n <- 100 # Select bin number for tiling ----
marker_of_interest <- "PD1" # Select marker of interest ----
## Extraction of the data ----
#! PROBLEM - Handling intensity on borders of tumor !#
# Global binning based on entire data set
global_bin_width <- (max(eda_df$X) - min(eda_df$X)) / bin_n
global_min_x <- min(eda_df$X)
global_min_y <- min(eda_df$Y)

# Convert the dataframe to long format
eda_df_long <- eda_df %>%
  pivot_longer(cols = c(PD1, CD8, CD3, TIM3, LAG3, CK),
               names_to = "marker",
               values_to = "value")

# Create an empty list to store the results for each patient
tiled_markerheatmap_xID <- list()

# Get unique patient IDs
unique_ids <- unique(eda_df$ID)

# Iterate over each unique patient ID
for(id in unique_ids) {
  
  # Create bins for the current patient based on global bins
  bin_data <- eda_df_long[eda_df_long$ID == id,] %>%  ########### ODIO ESTA LINEA!!!!!!
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
  
  # Store the result in the list
  tiled_markerheatmap_xID[[as.character(id)]] <- resultado
  print(paste("Processing patient", id, "data."))
}

# Initialize an empty list to store the results for each response
tiled_markerheatmap_xResponse <- list()

# Iterate through the unique response groups
for (response in unique(responses$Response)) {
  
  print(paste("Retrieving", response, "IDs."))
  
  # Get the IDs associated with the current response group
  matching_ids <- as.character(responses$ID[responses$Response == response])
  
  print(paste("Matching IDs:", cat(matching_ids)))
  
  # Filter tiled_markerheatmap_xResponse by the matching IDs but maintain the individual structure
  tiled_markerheatmap_xResponse[[response]] <- tiled_markerheatmap_xID[matching_ids]
}

# 1. Process the data
tiled_medmarkerint_xResponse <- list()

for (response in names(tiled_markerheatmap_xResponse)) {
  patient_list <- tiled_markerheatmap_xResponse[[response]]
  
  # Combine all patient data for current response
  combined_data <- bind_rows(patient_list) %>% 
    filter(marker == marker_of_interest)
  
  # Compute median for each bin location across all patients of current response
  response_data <- combined_data %>%
    group_by(bin_x, bin_y, X, Y) %>%
    summarise(median_value = median(median_value_across_images, na.rm = TRUE)) %>%
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

## VISUALS - Results ----
# Plot
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

## STATISTICS ----
########################################### FEATURE 4 (S/T) -----
# EXTRACT - Categorical dotplot locations ----
# Initialize empty list
list_of_dataframes <- list()

# Iterate through each image
for (i in seq_along(spiat_predicted_phenotypes)) {
  
  # Extract data
  image_name <- spiat_predicted_phenotypes[[i]]$image_name
  coords <- spiat_predicted_phenotypes[[i]]$spatial_coords
  phenotype <- spiat_predicted_phenotypes[[i]]$phenotype
  
  # Convert data to dataframe
  df <- data.frame(
    Image = image_name,
    X = coords[, "Cell.X.Position"],
    Y = coords[, "Cell.Y.Position"],
    Phenotype = phenotype
  )
  
  # Append
  list_of_dataframes[[i]] <- df
  
}

# Combine
spiat_pheno_xy <- do.call(rbind, list_of_dataframes)

# Add the ID and response columns
spiat_pheno_xy$ID <- eda_df$ID[match(spiat_pheno_xy$Image, eda_df$Image)]
spiat_pheno_xy$Response <- responses$Response[match(spiat_pheno_xy$ID, responses$ID)]

# Rearrange
spiat_pheno_xy <- spiat_pheno_xy[, c("Image", "ID", "Phenotype", "Response", "X", "Y")]

rm(coords, df, i, image_name, phenotype, list_of_dataframes)

########################################### FEATURE 5 (S/T) -----
# EXTRACT - Cell distances ----
# Aqui tengo que hacer un histograma antes de cuales son las celulas de todos estos hp que mas influyen tal vez PCA?
########################################### FEATURE 6 (S/T) -----
## Custom tiling ----
plot_sample <- function(spiat_scattercount, image_name, marker, bins) {
  
  # Get the data for the specific image
  image_data_list <- spiat_scattercount[[image_name]]
  
  # Convert the list to a data frame
  original <- as.data.frame(image_data_list)
  
  # Global binning based on entire data set
  global_bin_width <- (max(original$X) - min(original$X)) / bins
  
  # Subset the data to only include rows where "P" contains the specified marker
  image_data <- original[grepl(marker, original$P, fixed = TRUE), ]
  
  # Create the ggplot
  gg <- ggplot() +
    geom_bin2d(data = image_data, 
               aes(x = X, y = Y, fill = after_stat(count)),  # Update as per the warning
               binwidth = c(50, 50), 
               alpha = 0.9) +
    geom_point(data = image_data, 
               aes(x = X, y = Y, alpha = 0.99)) +
    scale_fill_distiller(palette = "RdYlBu") +
    coord_fixed() +
    theme(panel.background = element_blank())
  
  print(gg)
}

plot_sample(spiat_scattercount, SPIAT_tifs[1], "LAG3", 10)

## Custom plot symbol con ggplot ----
plot_sample <- function(spiat_scattercount, image_name, marker, r) {
  
  image_data <- as.data.frame(spiat_scattercount[[image_name]])
  image_data <- image_data[grepl(marker, image_data$P, fixed = TRUE), ]
  
  # Calculate IQR for X and Y
  iqr_x <- IQR(image_data$X)
  iqr_y <- IQR(image_data$Y)
  
  # Use a fraction of the smaller IQR as radius
  radius <- min(iqr_x, iqr_y) * r
  
  bin_counts <- data.frame(X = numeric(), Y = numeric(), Count = numeric())
  
  for (i in 1:nrow(image_data)) {
    x_center <- image_data$X[i]
    y_center <- image_data$Y[i]
    
    within_radius <- (image_data$X - x_center)^2 + (image_data$Y - y_center)^2 <= radius^2
    count <- sum(within_radius)
    
    # Debugging print
    if(count == 1) {
      print(paste("Only one point at X:", x_center, ", Y:", y_center))
    }
    
    bin_counts <- rbind(bin_counts, data.frame(X = x_center, Y = y_center, Count = count))
  }
  
  # Sort by Count to ensure smaller points are plotted first
  bin_counts <- bin_counts[order(bin_counts$Count), ]
  
  # Create base R plot
  plot(bin_counts$X, bin_counts$Y, type = "n", xlab = "X", ylab = "Y", main = marker)  # Empty plot for setting axes
  
  # Add points with color and size based on Count
  color_range <- colorRampPalette(c("blue", "red"))(max(bin_counts$Count))
  point_colors <- color_range[bin_counts$Count]
  point_sizes <- sqrt(bin_counts$Count)  # Adjust as needed
  
  points(bin_counts$X, bin_counts$Y, col = point_colors, cex = point_sizes, pch = 19, alpha = 0.9)
  
  return(NULL)
}

gg_custom_scatter <- plot_sample(spiat_scattercount, SPIAT_tifs[1], "LAG3", 0.2)

## Tesselation plots ---- 
# Create a subset of the data
subset_data <- spiat_scattercount[[SPIAT_tifs[[1]]]]

# Create ppp object
Z <- ppp(subset_data$X, subset_data$Y, window = w, marks = subset_data$LAG3)

# Extract the marks ("LAG3" values)
marks_values <- marks(Z)

# Calculate quartiles on the marks
b <- quantile(marks_values, probs = (0:4)/4)

# Cut the data based on quartiles
Zcut <- cut(Z, breaks = b, labels = 1:4)

# Create tessellation
V <- tess(image = Zcut, window = w)

# Plot tessellation
plot(V)