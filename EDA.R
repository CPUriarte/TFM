# Overview ----
# Title: Exploratory data analysis of advanced urothelial carcinoma samples segmented with DeepCell Mesmer and exported with QuPath.
# Date: 20/05/2023
# Author: Cyril P
# Director: Carlos de Andrea
# Institution: University of Navarra
# Libraries ----
library(tidyverse)
library(gridExtra)
library(plotly)
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
library(igraph)
library(ggdist) # Applies genlog scale to ggplot values and colors
library(yarrr) # Pirate plots
# Load and format data ----
cat("\014")

# Raw measurements exported from QuPath
raw_measurements <- read.table("Data/raw_measurements.tsv", header = TRUE, sep = "\t")

# Raw associated clinical data
responses <- read.csv("Data/Responses.csv", header = TRUE, sep = ",")

# Formatting ----
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
responses$n_cells <- cell_count_vector

# Dropping non-evaluable patients
responses$Response <- droplevels(responses$Response, "Not evaluable/NE")
levels(responses$Response) <- c("CR", "PR", "PD", "SD")
SPIAT_tifs <- unique(eda_df$Image) # Unique list of SPIAT curated images

# Clear intermediary variables
rm(list = c("cell_counts", "NAs", "cell_count_vector", "images_with_NAs", "new_order", "response_vector"))

# Clear console for results
cat("\014")
summary(responses)
summary(eda_df[,c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")])

# S1 | Number of cells ----
cat("\014")

# Number of cells per image (QQPLOT)
eda_df %>%
  group_by(Image) %>%
  summarise(cell_count = n()) %>%
  {qqPlot(.$cell_count, dist = "norm", xlab = "Quantiles", ylab = "Number of Cells", main="Number of Cells per Image")}

# Number of cells per ID (QQPLOT)
qqPlot(responses$n_cells, dist = "norm", xlab = "Quantiles", ylab = "Number of Cells", main = "Number of Cells Per Patient")

# S2 | Age ----
age_stats <- list()
age_stats$residuals <- residuals(aov(Age ~ Response, data=responses))
shapiro.test(age_stats$residuals)

age_stats$levene <- leveneTest(Age ~ Response, data=responses)

age_stats$levene$`Pr(>F)`[1]

age_stats$aov <- aov(Age ~ Response, data=responses)

responses %>%
  group_by(ID) %>%
  ggplot(aes(x = Age)) +
  geom_histogram(aes(y = ..density..), fill = "steelblue", color = "white", bins = 30) +
  geom_density(aes(y = ..density..), color = "red") +  
  labs(title = "Distribution of Age",
       x = "Age",
       y = "Count")

cat("\014")
summary(age_stats$aov)

# S3 | Sex ----
# Fisher exact test
age_fishertest <- fisher.test(table(responses$Sex, responses$Response))

# Hypergeometric test - LaTeX: P(X=k) = \frac{\binom{N}{n} \binom{M}{k} \binom{N-M}{n-k}}{\binom{N}{n}}
cat("\014")
print("Probability of having 10 males in a group given female to male ratio is 1:3: \n")
choose(31, 10) / choose(46, 10)

print("Probability of having 10 males in a group given female to male ratio is 1:3 and all other distributions among groups:")
age_fishertest$p.value

# S4 | Marker intensity ----
eda_by_ID <- eda_df %>%
  group_by(ID) %>%
  summarise(
    DAPI = median(DAPI, na.rm = TRUE),
    PD1 = median(PD1, na.rm = TRUE),
    CD8 = median(CD8, na.rm = TRUE),
    CD3 = median(CD3, na.rm = TRUE),
    TIM3 = median(TIM3, na.rm = TRUE),
    LAG3 = median(LAG3, na.rm = TRUE),
    CK = median(CK, na.rm = TRUE)
  )

eda_by_ID <- eda_by_ID %>%
  mutate(Response = responses$Response[match(ID, responses$ID)]) %>%
  select(ID, Response, everything())

# STATISTICS - Marker intensity comparison across response
## KRUSKAL WALLIS TEST
markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

test_for_marker <- function(marker) {
  
  # Kruskal-Wallis test
  kw_p_val <- kruskal.test(get(marker) ~ Response, data = eda_by_ID)$p.value
  
  # Return combined data frame
  return(data.frame(marker = marker, 
                    kw_p_value = kw_p_val))
}

# Run tests and bind results
mIntensity_KWresults <- bind_rows(lapply(markers, test_for_marker))

## POST-HOC ANALYSIS WITH DUNN'S TEST
# CD3 Dunn test
CD3_dun <- dunnTest(CD3 ~ Response, data = eda_by_ID, method = "bh")

CD3_mIntensity_dunn <- data.frame(
  Comparison = CD3_dun$res$Comparison,
  Z = CD3_dun$res$Z,
  P.unadj = CD3_dun$res$`P.unadj`,
  P.adj = CD3_dun$res$`P.adj`
)

## CD8 Dunn test
CD8_dun <- dunnTest(CD8 ~ Response, data = eda_by_ID, method = "bh")

CD8_mIntensity_dunn <- data.frame(
  Comparison = CD8_dun$res$Comparison,
  Z = CD8_dun$res$Z,
  P.unadj = CD8_dun$res$`P.unadj`,
  P.adj = CD8_dun$res$`P.adj`
)

CD8_mIntensity_dunn <- CD8_mIntensity_dunn %>%
  arrange(P.adj)

# Violin plots
cd3_plot <- ggplot(eda_by_ID, aes(x=Response, y=CD3)) +
  geom_violin() +
  stat_summary(fun.y=median, geom="point", size=3, color="red") +
  geom_signif(comparisons = list(c("PD", "PR")), annotations = "**", y_position = c(0.75), size = 1.15) +
  geom_signif(comparisons = list(c("PR", "SD")), annotations = "*", y_position = c(1.05), size = 1.15) +
  theme_classic()

cd8_plot <- ggplot(eda_by_ID, aes(x=Response, y=CD8)) +
  geom_violin() +
  stat_summary(fun.y=median, geom="point", size=3, color="red") +
  geom_signif(comparisons = list(c("PD", "PR")), annotations = "*", y_position = c(0.4), size = 1.15) +
  theme_classic()

combined_plot <- cowplot::plot_grid(cd3_plot, cd8_plot, ncol = 2)
print(combined_plot)

cat("\014")
CD3_mIntensity_dunn %>%
  filter(P.adj <= 0.05) %>% 
  arrange(P.adj)

CD8_mIntensity_dunn %>%
  filter(P.adj <= 0.05) %>% 
  arrange(P.adj)

## Z-score plot
# Filter significant results from Dunn's test for CD3 and CD8
sig_CD3_dunn <- CD3_mIntensity_dunn %>%
  filter(P.adj <= 0.05) %>% 
  arrange(P.adj)

sig_CD8_dunn <- CD8_mIntensity_dunn %>%
  filter(P.adj <= 0.05) %>% 
  arrange(P.adj)

# Combine significant results into a single dataframe
sig_dunn_combined <- bind_rows(
  mutate(sig_CD3_dunn, marker = "CD3"),
  mutate(sig_CD8_dunn, marker = "CD8")
)

# Plotting Z-scores from significant results
zscore_plot <- ggplot(sig_dunn_combined, aes(x = Comparison, y = Z, color = marker, shape = marker)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_classic()

print(zscore_plot)

# S5 | Custom phenotype prediction ----
cat("\014")
## Thresholding - V1 ----
# Threshold dataset
img_threshs <- eda_df %>% 
  dplyr::select(Image, ID) %>% 
  distinct() %>%
  mutate(
    PD1 = NA_real_,
    LAG3 = NA_real_,
    TIM3 = NA_real_,
    CD3 = NA_real_,
    CD8 = NA_real_,
    CK = NA_real_
  )

rownames(img_threshs) <- NULL

# Custom optimal cutoff function using valleys of density
optimal_cutoff_valley <- function(intensity_values) {
  
  intensity_density <- stats::density(intensity_values, na.rm=TRUE)
  
  # Peaks
  peak_ycords <- pracma::findpeaks(intensity_density$y)[,1]
  peak_xcords <- intensity_density$x[match(peak_ycords, intensity_density$y)]
  
  # Valleys
  valley_ycords <- pracma::findpeaks(-intensity_density$y)[,1] * -1
  valley_xcords <- intensity_density$x[match(valley_ycords, intensity_density$y)]
  
  # Sort peaks to find the highest ones
  peaks_df <- data.frame(peak_xcords, peak_ycords) %>%
    arrange(desc(peak_ycords))
  
  # Get highest peak
  highest_peak <- peaks_df$peak_ycords[1]
  
  # Check for second peak as high as 10% of the highest peak
  if (nrow(peaks_df) >= 2 && peaks_df$peak_ycords[2] >= 0.10 * highest_peak) {
    
    # Find the valley between these two peaks
    first_peak_xcord <- peaks_df$peak_xcords[1]
    second_peak_xcord <- peaks_df$peak_xcords[2]
    valley_between <- valley_xcords[which(valley_xcords > min(first_peak_xcord, second_peak_xcord) & valley_xcords < max(first_peak_xcord, second_peak_xcord))]
    if(length(valley_between) > 0) {
      return(valley_between[1])
    }
  }
  
  # Determine the cutoff based on the first valley after the highest peak
  ycord_max_density <- max(intensity_density$y)
  xcord_max_density <- intensity_density$x[match(ycord_max_density, intensity_density$y)]
  valley_df <- data.frame(valley_xcords, valley_ycords) %>%
    filter(valley_xcords >= xcord_max_density, 
           valley_ycords <= 0.25 * ycord_max_density)
  
  return(valley_df$valley_xcords[1])
}

# Summarize and compute the optimal cutoff based on valleys
img_threshs <- eda_df %>%
  group_by(Image) %>%
  summarise(across(c('PD1', 'CD8', 'CD3', 'TIM3', 'LAG3', 'CK'), 
                   ~optimal_cutoff_valley(.x),
                   .names = "cutoff_{.col}")) %>%
  left_join(eda_df %>% distinct(Image, ID), by = "Image")

original_cols <- colnames(eda_df)
merged_df <- merge(eda_df, img_threshs, by="Image", all.x = TRUE)
names(merged_df)[names(merged_df) == "ID.x"] <- "ID"
merged_df$ID.y <- NULL

cat("\014")
## Thresholding - V2 ----
img_threshs <- eda_df %>% 
  select(Image, ID) %>% 
  distinct() %>%
  mutate(
    PD1 = NA_real_,
    LAG3 = NA_real_,
    TIM3 = NA_real_,
    CD3 = NA_real_,
    CD8 = NA_real_,
    CK = NA_real_
  )

rownames(img_threshs) <- NULL

optimal_cutoff_valley <- function(intensity_values) {
  
  intensity_density <- stats::density(intensity_values, na.rm=TRUE)
  
  # Compute differences
  nabla_f <- diff(intensity_density$y)
  nabla_2f <- diff(nabla_f)
  
  # Compute the w_v(x) metric
  min_val <- min(nabla_2f)
  max_val <- max(nabla_2f)
  w_v <- (nabla_2f - min_val) / (max_val - min_val)
  
  # Identify valleys using peaks in w_v
  valleys <- pracma::findpeaks(w_v)

  # Sort peaks to find the highest ones
  peaks_df <- data.frame(x = intensity_density$x[-c(1, length(intensity_density$x))], 
                         w_v) %>%
    arrange(desc(w_v))
  
  # Get highest peak (corresponds to most prominent valley in original)
  highest_valley <- peaks_df$w_v[1]
  
  if (length(peaks_df)<2) {
    intensity_density <- stats::density(intensity_values, na.rm=TRUE)
    
    # Peaks
    peak_ycords <- pracma::findpeaks(intensity_density$y)[,1]
    peak_xcords <- intensity_density$x[match(peak_ycords, intensity_density$y)]
    
    # Valleys
    valley_ycords <- pracma::findpeaks(-intensity_density$y)[,1] * -1
    valley_xcords <- intensity_density$x[match(valley_ycords, intensity_density$y)]
    
    # Sort peaks to find the highest ones
    peaks_df <- data.frame(peak_xcords, peak_ycords) %>%
      arrange(desc(peak_ycords))
    
    # Get highest peak
    highest_peak <- peaks_df$peak_ycords[1]
    
    # Check for second peak as high as 10% of the highest peak
    if (nrow(peaks_df) >= 2 && peaks_df$peak_ycords[2] >= 0.10 * highest_peak) {
      
      # Find the valley between these two peaks
      first_peak_xcord <- peaks_df$peak_xcords[1]
      second_peak_xcord <- peaks_df$peak_xcords[2]
      valley_between <- valley_xcords[which(valley_xcords > min(first_peak_xcord, second_peak_xcord) & valley_xcords < max(first_peak_xcord, second_peak_xcord))]
      if(length(valley_between) > 0) {
        return(valley_between[1])
      }
    }
    
    # Determine cutoff based on first valley after highest peak
    ycord_max_density <- max(intensity_density$y)
    xcord_max_density <- intensity_density$x[match(ycord_max_density, intensity_density$y)]
    valley_df <- data.frame(valley_xcords, valley_ycords) %>%
      filter(valley_xcords >= xcord_max_density, 
             valley_ycords <= 0.25 * ycord_max_density)
    
    return(valley_df$valley_xcords[1])
  }else{
    
  return(peaks_df$x[1])
  }
}

img_threshs <- eda_df %>%
  group_by(Image) %>%
  summarise(across(c('PD1', 'CD8', 'CD3', 'TIM3', 'LAG3', 'CK'), 
                   ~optimal_cutoff_valley(.x),
                   .names = "cutoff_{.col}")) %>%
  left_join(eda_df %>% distinct(Image, ID), by = "Image")

original_cols <- colnames(eda_df)
merged_df <- merge(eda_df, img_threshs, by="Image", all.x = TRUE)
names(merged_df)[names(merged_df) == "ID.x"] <- "ID"
merged_df$ID.y <- NULL

## Phenotype assignment ----
## PROCESSING - Applying cutoff (Step 1/3)
immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

for(marker in immune_markers) {
  cutoff_col <- paste0("cutoff_", marker)
  merged_df[, marker] <- ifelse(merged_df[, marker] >= merged_df[, cutoff_col], 1, 0)
}

## PROCESSING - Handling superposition cases (Step 2/3)
# Add a column for duplicate identification
# CK adjustments
result <- pbapply(merged_df, 1, function(row) {
  
  # If CK is 1 and both CD3 and CD8 are 0, CK is left as 1.
  if(row["CK"] == 1 & (row["CD3"] == 0 & row["CD8"] == 0)) {
    row["CK"] <- 1
  }
  
  # If CK is 1 but either CD3 or CD8 are 1 too:
  if(row["CK"] == 1 & sum(as.numeric(row[immune_markers])) != 0 & (row["CD3"] == 1 | row["CD8"] == 1)) {
    
    # Duplicate the row
    duplicate_row <- row
    
    # For the duplicate cell, make CK 1 and all immune markers 0
    duplicate_row["CK"] <- 1
    duplicate_row[setdiff(immune_markers, "CK")] <- 0
    
    # For the original cell, make CK 0.
    row["CK"] <- 0
    
    # Return both original and duplicate rows
    return(list(original = row, duplicate = duplicate_row))
  } else {
    return(list(original = row))
  }
  
}, cl = 4)

# Convert the list back to dataframe
original_rows <- do.call(rbind, lapply(result, `[[`, "original"))
duplicate_rows <- do.call(rbind, lapply(result, function(x) if ("duplicate" %in% names(x)) x$duplicate else NULL))

final_df <- rbind.data.frame(original_rows, duplicate_rows)
final_df <- final_df[, c("Image", "X", "Y", "ID", "DAPI", immune_markers)]
final_df[,-1] <- lapply(final_df[,-1], as.numeric)

## PROCESSING - Assigning phenotypes (Step 3/3)
# Immune_markers vector
immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Assigning phenotypes
cat("\014")
phenotype_strings <- pbapply(final_df[, immune_markers], 1, function(row) {
  paste(names(row)[as.logical(row)], collapse = ",")
}, cl = 4)

# Handle 'Tumor' and 'Undefined' cases
phenotype_strings[phenotype_strings == "CK"] <- "Tumor" # Only name "Tumor" if "CK" is the sole marker
phenotype_strings[phenotype_strings == ""] <- "Undefined"

# Add phenotype column to the dataframe
final_df$Phenotype <- phenotype_strings

# Join final_df and responses
final_df_joined <- final_df %>%
  left_join(responses, by = "ID")

# Remove intermediate variables
rm(cutoff_col, duplicate_rows, merged_df, original_rows, original_cols, result, phenotype_strings)
cat("\014")

## Sample image with predicted phenotypes ----
IMG <- sample(SPIAT_tifs, size = 1)

p <- final_df_joined %>%
  filter(Image == IMG) %>%
  ggplot(aes(x = X, y = Y, color = Phenotype)) +
  geom_point(alpha = 0.7) +
  labs(title = paste(IMG, "|", unique(final_df_joined$Response[final_df_joined$Image==IMG])), x  = "X coordinate", y = "Y coordinate") +
  theme_minimal() +
  scale_color_discrete(name = "Phenotypes") +
  coord_cartesian(xlim = c(0, 1300), ylim = c(0, 1300))

ggplotly(p)

## VISUAL - Cutoff quality (IMAGES)
long_eda_df <- eda_df[eda_df$Image == IMG, ] %>%
  pivot_longer(cols = colnames(eda_df[, 6:11]), names_to = "marker", values_to = "intensity")

long_eda_df$cutoff <- sapply(1:nrow(long_eda_df), function(i) {
  img_threshs[img_threshs$Image == IMG, paste0("cutoff_", long_eda_df$marker[i])]
})

long_eda_df$cutoff <- as.numeric(long_eda_df$cutoff)

ggplot(long_eda_df, aes(x = log(intensity+1))) +
  geom_density(fill = "cyan") +
  geom_vline(aes(xintercept = cutoff), color = "firebrick3", linetype = "dashed") +
  labs(title = paste("Density distribution of markers for", IMG)) +
  facet_wrap(~marker, scales = "free") +
  theme_minimal()

## Resulting phenotype proportions ----
# Median Proportion of Phenotypes per Response Group (PATIENTS)
# Step 1: Get unique phenotypes across all images
cp_unique_phenos <- unique(final_df$Phenotype)

# Step 2: Create a dataframe with image names and columns named based on unique phenotypes
cp_pheno_counts <- data.frame(Image = unique(final_df$Image),
                                 matrix(NA, ncol = length(cp_unique_phenos), nrow = length(SPIAT_tifs)))

colnames(cp_pheno_counts) <- c("Image", cp_unique_phenos)

cp_pheno_counts$ID <- eda_df$ID[match(cp_pheno_counts$Image, eda_df$Image)]
cp_pheno_counts$Response <- responses$Response[match(cp_pheno_counts$ID, responses$ID)]
cp_current_cols <- setdiff(colnames(cp_pheno_counts), c("Image", "ID", "Response", "Undefined", "Tumor"))

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
primary_weights <- sapply(cp_current_cols, get_primary_weight)
secondary_weights <- sapply(cp_current_cols, get_secondary_weight)

# Combine primary and secondary weights to get the order
# Use primary weight * a big number (like 1000) + secondary weight
combined_weights <- primary_weights * 1000 + secondary_weights

# Sort columns based on combined weights
cp_ordered_cols <- cp_current_cols[order(combined_weights)]

# Add the "Image", "None", and "CK" columns at the start
cp_final_order <- c("Image","ID", "Response", "Undefined", "Tumor", cp_ordered_cols)

# Re-arrange columns of the dataframe based on final_order
cp_pheno_counts <- cp_pheno_counts[, cp_final_order]

# Count the occurrences of each phenotype for each image
phenotype_counts <- aggregate(. ~ Image + Phenotype, final_df, function(x) sum(!is.na(x)))

# Step 2 and 3: Match and replace the NA values in cp_pheno_counts
for(i in 1:nrow(cp_pheno_counts)) {
  for(j in 2:ncol(cp_pheno_counts)) {
    # Fetch the image and phenotype we're trying to match
    current_image <- cp_pheno_counts$Image[i]
    current_pheno <- colnames(cp_pheno_counts)[j]
    
    # Find the matching row in phenotype_counts
    matching_row <- which(phenotype_counts$Image == current_image & 
                            phenotype_counts$Phenotype == current_pheno)
    
    # If a match is found, replace NA with the count
    if(length(matching_row) > 0) {
      cp_pheno_counts[i, j] <- phenotype_counts[matching_row, 3]
    }
  }
}

# FINAL RESHAPE
cp_long_format <- cp_pheno_counts %>%
  select(-Image) %>%
  gather(key = 'cp_Phenotype', value = 'cp_Count', -ID, -Response)

# Group data by ID summarize counts
cp_total_cells_per_ID <- cp_long_format %>%
  group_by(ID) %>%
  summarise(cp_TotalCells = sum(cp_Count, na.rm = TRUE))

# Long Data with Proportions
cp_long_with_props <- cp_long_format %>%
  group_by(ID, cp_Phenotype) %>%
  summarise(cp_PhenotypeCount = sum(cp_Count, na.rm = TRUE)) %>%
  inner_join(cp_total_cells_per_ID, by = "ID") %>%
  mutate(cp_Prop = cp_PhenotypeCount / cp_TotalCells)

# Assign responses
cp_long_with_props$cp_Response <- responses$Response[match(cp_long_with_props$ID, responses$ID)]

cp_long_with_props$cp_Response <- factor(
  cp_long_with_props$cp_Response,
  levels = c("PD", "SD", "PR", "CR")
)

# Aggregate by Phenotype and Response
cp_aggregated_pheno <- cp_long_with_props %>%
  group_by(cp_Phenotype, cp_Response) %>%
  summarise(cp_PhenotypeSum = sum(cp_PhenotypeCount)) %>%
  ungroup() %>%
  left_join(
    cp_long_with_props %>%
      group_by(cp_Response) %>%
      summarise(cp_TotalSum = sum(cp_PhenotypeCount)),
    by = "cp_Response"
  )

# Calculate proportions for each group
cp_aggregated_pheno <- cp_aggregated_pheno %>%
  mutate(cp_GroupProp = cp_PhenotypeSum / cp_TotalSum)

# Generate Plot
cp_plot <- ggplot(cp_aggregated_pheno, aes(x = cp_Response, y = cp_GroupProp, fill = cp_Phenotype)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Proportion of Each Phenotype by Response Group",
    x = "Response Group",
    y = "Proportion",
    fill = "Phenotype"
  )

# Create interactive plot
ggplotly(cp_plot)
cat("\014")
# Session info ----
# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 22621)
# 
# Matrix products: default
# 
# Locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# Attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# Other attached packages:
# [1] ggdist_3.3.0           igraph_1.5.1           spatstat_3.0-6         spatstat.linnet_3.1-1  spatstat.model_3.2-4   rpart_4.1.19          
# [7] spatstat.explore_3.2-1 nlme_3.1-162           spatstat.random_3.1-5  spatstat.geom_3.2-1    spatstat.data_3.0-1    corrplot_0.92         
# [13] ggsignif_0.6.4         FSA_0.9.4              ggridges_0.5.4         pbapply_1.7-0          reshape2_1.4.4         moments_0.14.1        
# [19] car_3.1-2              carData_3.0-5          readxl_1.4.2           BBmisc_1.13            plotly_4.10.2          gridExtra_2.3         
# [25] lubridate_1.9.2        forcats_1.0.0          stringr_1.5.0          dplyr_1.1.2            purrr_1.0.1            readr_2.1.4           
# [31] tidyr_1.3.0            tibble_3.2.1           ggplot2_3.4.3          tidyverse_2.0.0       
# 
# Loaded via a namespace (and not attached):
# [1] httr_1.4.6            jsonlite_1.8.4        viridisLite_0.4.2     splines_4.2.2         distributional_0.3.2  cellranger_1.1.0     
# [7] yaml_2.3.7            pillar_1.9.0          backports_1.4.1       lattice_0.20-45       glue_1.6.2            digest_0.6.33        
# [13] polyclip_1.10-4       checkmate_2.2.0       colorspace_2.1-0      htmltools_0.5.6       Matrix_1.5-4.1        plyr_1.8.8           
# [19] spatstat.sparse_3.0-1 pkgconfig_2.0.3       scales_1.2.1          tensor_1.5            spatstat.utils_3.0-3  tzdb_0.4.0           
# [25] pracma_2.4.2          timechange_0.2.0      mgcv_1.8-42           farver_2.1.1          generics_0.1.3        ellipsis_0.3.2       
# [31] withr_2.5.0           lazyeval_0.2.2        cli_3.6.1             crayon_1.5.2          magrittr_2.0.3        deldir_1.0-9         
# [37] fansi_1.0.4           dunn.test_1.3.5       tools_4.2.2           data.table_1.14.8     hms_1.1.3             lifecycle_1.0.3      
# [43] munsell_0.5.0         compiler_4.2.2        rlang_1.1.1           grid_4.2.2            rstudioapi_0.14       htmlwidgets_1.6.2    
# [49] goftest_1.2-3         crosstalk_1.2.0       labeling_0.4.2        gtable_0.3.4          abind_1.4-5           R6_2.5.1             
# [55] fastmap_1.1.1         utf8_1.2.3            stringi_1.7.12        parallel_4.2.2        Rcpp_1.0.10           vctrs_0.6.3          
# [61] tidyselect_1.2.0 