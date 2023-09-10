############################################ LIBRARIES ----
library(SPIAT)
library(tidyverse)
library(gridExtra)
library(plotly)
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
############################################ LOADING DATA ----
# 1. DATA READING AND FORMATTING
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
responses$n_cells <- cell_count_vector

# Dropping non-evaluable patients
responses$Response <- droplevels(responses$Response, "Not evaluable/NE")
levels(responses$Response) <- c("CR", "PR", "PD", "SD")

# Clear intermediary variables
rm(list = c("cell_counts", "NAs", "cell_count_vector", "images_with_NAs", "new_order", "response_vector"))

SPIAT_tifs <- unique(eda_df$Image)

# Number of images per patient and TIF index
tif_index <- eda_df %>%
  group_by(ID) %>%
  summarise(
    Count = n_distinct(Image),
    Images = list(unique(Image))
  )

tif_index$Matched_Indices <- sapply(tif_index$Images, function(images_list) which(SPIAT_tifs %in% images_list))

############################################ EXPLORATORY DATA ANALYSIS (EDA) ----
# Retrieve SPIAT image vector indexes based on ID ----
ID <- 901 # ----
tif_index_retriever <- tif_index[tif_index$ID == ID, "Matched_Indices"][[1]]
print(paste("SPIAT_tifs indexes for patient", ID, "are:", tif_index_retriever))
# 1. IMAGES PER PATIENT ----
summary(responses)
summary(eda_df[,c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")])

ggplot(tif_index, aes(x = as.factor(Count))) +
  geom_bar(aes(fill = responses$Response[match(ID, responses$ID)])) +
  labs(
    title = "Distribution of Patients based on Number of Images and Treatment Response",
    x = "Number of Images",
    y = "Count of Patients",
    fill = "Response to Treatment"
  )

## 1.1. NUMBER OF CELLS
# 2. NUMBER OF CELLS
## 1.1.1. VISUALS - QQplot of number of cells distribution(IMG & ID) ----
# Number of cells per image (DENSITY)
eda_df %>%
  group_by(Image) %>%
  tally() %>%
  ggplot(aes(x = n)) +
  geom_histogram(aes(y = ..density..), fill = "steelblue", color = "white", bins = 30) +
  geom_density(aes(y = ..density..), color = "red") +
  labs(title = "Distribution of Number of Cells per Image",
       y = "Density",
       x = "Number of Cells")

result_dagostino <- eda_df %>% # Normality test (AGOSTINO PEARSON)
  group_by(Image) %>%
  summarise(cell_count = n()) %>%
  {agostino.test(.$cell_count)}

result_dagostino

# Number of cells per image (QQplot)
eda_df %>%
  group_by(Image) %>%
  summarise(cell_count = n()) %>%
  {qqPlot(.$cell_count, dist = "norm", xlab = "Quantiles", ylab = "Number of Cells", main="Number of Cells per Image")}

# Number of cells per patient (DENSITY)
shapiro.test(responses$n_cells) # Normality test
responses %>%
  group_by(ID) %>%
  ggplot(aes(x = n_cells)) +
  geom_histogram(aes(y = ..density..), fill = "steelblue", color = "white", bins = 30) +
  geom_density(aes(y = ..density..), color = "red") +
  labs(title = "Distribution of Number of Cells per Patient",
       y = "Density",
       x = "Number of Cells")

# Number of cells (QQplot)
qqPlot(responses$n_cells, dist = "norm", xlab = "Quantiles", ylab = "Number of Cells", main = "Number of Cells Per Patient")

## 1.1.2. VISUALS - Barplot of age distribution (ID) ----
# Overall age distribution
responses %>%
  group_by(ID) %>%
  ggplot(aes(x = Age)) +
  geom_histogram(aes(y = ..density..), fill = "steelblue", color = "white", bins = 30) +
  geom_density(aes(y = ..density..), color = "red") +  
  labs(title = "Distribution of Age",
       x = "Age",
       y = "Count")

# 2. AGE
# 2. Age distribution ----
age_stats <- list()
age_stats$residuals <- residuals(aov(Age ~ Response, data=responses))
shapiro.test(age_stats$residuals)

age_stats$levene <- leveneTest(Age ~ Response, data=responses)

age_stats$levene$`Pr(>F)`[1]

age_stats$aov <- aov(Age ~ Response, data=responses)

summary(age_stats$aov)

## 2.3. SEX (Explicado en papers) ----  # RARO
binom.test(sum(responses$Sex == "Male"), nrow(responses), p = 0.3)
ggplot(responses, aes(x = Sex)) +
  geom_bar(fill = "#0072B2", color = "black") + 
  labs(title = "Distribution of Sex", x = "Sex", y = "Count") +
  theme_minimal()


# 3. SEX
# 3. Sex distribution ----
# LaTeX: P(X=k) = \frac{\binom{N}{n} \binom{M}{k} \binom{N-M}{n-k}}{\binom{N}{n}}
choose(31, 10) / choose(46, 10)

# Fisher exact test
age_fishertest <- fisher.test(table(responses$Sex, responses$Response))
age_fishertest$p.value

# 4. Marker intensity ----
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

## 4.1. STATISTICS - Normality test ----
markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Shapiro-Wilk test for normality
shapiro_test_by_group <- function(marker) {
  eda_by_ID %>%
    group_by(Response) %>%
    summarise(p_value = shapiro.test(get(marker))$p.value) %>%
    mutate(marker = marker)
}

eda_mIntensity_shapiro_byID <- bind_rows(lapply(markers, shapiro_test_by_group))

## 4.2. STATISTICS - Homoscedasticity and KW test ----
# Levene's test for homogeneity
test_for_marker <- function(marker) {
  
  # Levene's test
  levene_p_val <- leveneTest(get(marker) ~ Response, data = eda_by_ID)$`Pr(>F)`[1]
  
  # Kruskal-Wallis test
  kw_p_val <- kruskal.test(get(marker) ~ Response, data = eda_by_ID)$p.value
  
  # Return combined data frame
  return(data.frame(marker = marker, 
                    levene_p_value = levene_p_val,
                    kw_p_value = kw_p_val))
}

# Running tests and binding results
eda_mIntensity_levene_KWfdr <- bind_rows(lapply(markers, test_for_marker))

# Applying Benjamini-Hochberg adjustment
eda_mIntensity_levene_KWfdr$adj_kw_p_value <- p.adjust(eda_mIntensity_levene_KWfdr$kw_p_value, method = "BH")

## 4.3. STATISTICS - Dunn test ----
dunn_results <- dunnTest(CD3 ~ Response, data = eda_by_ID, method = "bh")

eda_mIntensity_dunn <- data.frame(
  Comparison = dunn_results$res$Comparison,
  Z = dunn_results$res$Z,
  P.unadj = dunn_results$res$`P.unadj`,
  P.adj = dunn_results$res$`P.adj`
)

# Visualize
p <- ggplot(eda_by_ID, aes(x = Response, y = CD3)) +
  geom_violin() + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(title = "Mean CD3 Intensity by Response Group", y = NULL, x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

# Significance bars with star notations
p + geom_signif(comparisons = list(c("PD", "PR"), c("PR", "SD")), 
                y_position = c(max(eda_by_ID$CD3) * 1.08, max(eda_by_ID$CD3) * 1.14), 
                map_signif_level = TRUE)

# Sorting
eda_mIntensity_shapiro_byID <- eda_mIntensity_shapiro_byID %>%
  arrange(p_value)
print(eda_mIntensity_shapiro_byID)

eda_mIntensity_levene_KWfdr <- eda_mIntensity_levene_KWfdr %>%
  arrange(adj_kw_p_value)
print(eda_mIntensity_levene_KWfdr)

eda_mIntensity_dunn <- eda_mIntensity_dunn %>%
  arrange(P.adj)
print(eda_mIntensity_dunn)

## 4.4. VISUALS - Results ----
# Correlation plot by ID
combined_data <- merge(eda_df, responses, by = "ID")

# Subset
subset_data <- combined_data[, c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")]

# Calculate
correlation_matrix <- cor(subset_data, method = "spearman", use = "complete.obs")

# Correlation matrix
corrplot(correlation_matrix, method = "color", addrect = NULL)

# Marker intensities distribution accross images (BOXPLOT)
ggplot(eda_df %>% pivot_longer(cols = PD1:CK, names_to = "Marker", values_to = "Intensity"), 
       aes(x = log1p(Intensity), fill = Marker)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ Marker, scales = "free_x") +
  labs(title = "Density Plot of Log-transformed Marker Intensities", 
       y = "Density", 
       x = "Log-transformed Intensity (log(1+Intensity))") +
  theme_minimal()

# Distribution of marker intensity accross patients per response groups
# Pivoting eda_by_ID (consists of aggregated patients by median)
long_eda_df <- eda_by_ID %>% 
  pivot_longer(cols = PD1:CK, names_to = "Marker", values_to = "Value")

long_eda_df <- left_join(long_eda_df, responses[,1:3], by = "ID")

# Log-scaled
shift_constant <- 0.01

ggplot(long_eda_df, aes(x = Value + shift_constant, y = Response, fill = Response)) +
  geom_density_ridges2(alpha = 0.5) +
  facet_wrap(~ Marker, scales = "free_x") +
  scale_x_log10() +
  labs(title = "Log-scaled Density of Marker Intensities by Response") +
  theme_minimal()

############################################ PHENOTYPING ATTEMPT ----
# 1. Determining cutoff V1.0 (default) ----
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

# Define the optimal cutoff function using valleys of density
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

# Use dplyr to group, summarize and compute the optimal cutoff based on valleys
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

## 1.1. Determining cutoff V3.0 (not working currently) ----
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
  
  # Compute the differences
  nabla_f <- diff(intensity_density$y)
  nabla_2f <- diff(nabla_f)
  
  # Compute the w_v(x) metric
  min_val <- min(nabla_2f)
  max_val <- max(nabla_2f)
  w_v <- (nabla_2f - min_val) / (max_val - min_val)
  
  # Identify valleys in the original histogram using peaks in w_v
  valleys <- pracma::findpeaks(w_v)
  
  # Assuming you still want to work with the highest peaks logic:
  # Sort peaks to find the highest ones
  peaks_df <- data.frame(x = intensity_density$x[-c(1, length(intensity_density$x))], 
                         w_v) %>%
    arrange(desc(w_v))
  
  # Get highest peak (now corresponds to the most prominent valley in original)
  highest_valley <- peaks_df$w_v[1]
  
  # You can adjust this logic based on how you want to use this metric for thresholding
  # For this example, we're returning the x corresponding to the highest value in w_v
  return(peaks_df$x[1])
}

# Use dplyr to group, summarize and compute the optimal cutoff based on valleys
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

# 2. PROCESSING - Applying cutoff ----
immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

for(marker in immune_markers) {
  cutoff_col <- paste0("cutoff_", marker)
  merged_df[, marker] <- ifelse(merged_df[, marker] >= merged_df[, cutoff_col], 1, 0)
}

## 2.1. PROCESSING - Handling superposition ---- 
# Add a column for duplicate identification
# CK adjustment
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

## 2.2. PROCESSING - Assigning phenotypes ----
# Immune_markers vector
immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Paste values
phenotype_strings <- pbapply(final_df[, immune_markers], 1, function(row) {
  paste(names(row)[as.logical(row)], collapse = ",")
}, cl = 4)

# Handle 'Tumor' and 'Undefined' cases
phenotype_strings[phenotype_strings == "CK"] <- "Tumor" # Only name "Tumor" if "CK" is the sole marker
phenotype_strings[phenotype_strings == ""] <- "Undefined"

# Add phenotype column to the dataframe
final_df$Phenotype <- phenotype_strings

random <- function() {
  return(sample(1:90, 1))
}

# Remove intermediate variables to clear memory
rm(cutoff_col, duplicate_rows, merged_df, original_rows, original_cols, result, phenotype_strings)

# 3. VISUALS - Phenotype dotplot (IMAGES) ----
selected_image <- SPIAT_tifs[random()] # 1 to 90 possible images ----

# Join final_df and responses
final_df_joined <- final_df %>%
  left_join(responses, by = "ID")

# Plot
p <- final_df_joined %>%
  filter(Image == selected_image) %>%
  ggplot(aes(x = X, y = Y, color = Phenotype)) +
  geom_point(alpha = 0.7) +
  labs(
    title = paste(
      selected_image,
      "(", first(unique(final_df_joined$Response[final_df_joined$Image == selected_image])), ")",
      sep = " "
    ),
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_minimal() +
  scale_color_discrete(name = "Phenotypes") +
  coord_cartesian(xlim = c(0, 1300), ylim = c(0, 1300))

ggplotly(p)

## 3.1. VISUALS - Marker distribution and tresholding (IMAGES) ----

long_eda_df <- eda_df[eda_df$Image == selected_image, ] %>%
  pivot_longer(cols = colnames(eda_df[, 6:11]), names_to = "marker", values_to = "intensity")

long_eda_df$cutoff <- sapply(1:nrow(long_eda_df), function(i) {
  img_threshs[img_threshs$Image == selected_image, paste0("cutoff_", long_eda_df$marker[i])]
})

long_eda_df$cutoff <- as.numeric(long_eda_df$cutoff)

ggplot(long_eda_df, aes(x = intensity)) +
  geom_density(fill = "aquamarine") +
  geom_vline(aes(xintercept = cutoff), color = "red", linetype = "dashed") +
  labs(title = paste("Density distribution of markers for", selected_image),
       subtitle = paste("Cutoff values are displayed as red dashed lines"),
       x = "Intensity",
       y = "Density") +
  facet_wrap(~marker, scales = "free") +
  theme_minimal()

# Marker intensity dotplot
eda_df %>%
  filter(Image == selected_image) %>%
  ggplot(aes(x=X, y=Y, color=CK)) +
  geom_point(alpha=0.6) +
  scale_color_gradientn(n.breaks = 3, colors = terrain.colors(93)) +
  coord_cartesian(xlim=c(0, 1300), ylim=c(0, 1300)) +
  theme_minimal()


## 3.2. VISUALS - Median phenotype proportion per response group (RESPONSE) ----
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
# Step 1: Reshape the data
cp_long_format <- cp_pheno_counts %>%
  select(-Image) %>%
  gather(key = 'cp_Phenotype', value = 'cp_Count', -ID, -Response)

# Group the data by ID and then summarize the counts
cp_total_cells_per_ID <- cp_long_format %>%
  group_by(ID) %>%
  summarise(cp_TotalCells = sum(cp_Count, na.rm = TRUE))

# Create Long Data with Proportions
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

# SESSION INFO ----
# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 22621)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] spatstat_3.0-6              spatstat.linnet_3.1-1       spatstat.model_3.2-4        rpart_4.1.19                spatstat.explore_3.2-1     
#  [6] nlme_3.1-162                spatstat.random_3.1-5       spatstat.geom_3.2-1         spatstat.data_3.0-1         corrplot_0.92              
# [11] ggsignif_0.6.4              FSA_0.9.4                   ggridges_0.5.4              pbapply_1.7-0               moments_0.14.1             
# [16] car_3.1-2                   carData_3.0-5               readxl_1.4.2                reshape2_1.4.4              BBmisc_1.13                
# [21] plotly_4.10.2               gridExtra_2.3               lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0              
# [26] dplyr_1.1.2                 purrr_1.0.1                 readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1               
# [31] ggplot2_3.4.3               tidyverse_2.0.0             SPIAT_1.0.4                 SpatialExperiment_1.8.1     SingleCellExperiment_1.20.1
# [36] SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0             
# [41] S4Vectors_0.36.2            BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_0.63.0         
# 
# loaded via a namespace (and not attached):
#  [1] colorspace_2.1-0          rjson_0.2.21              deldir_1.0-9              ellipsis_0.3.2            scuttle_1.8.4            
#  [6] XVector_0.38.0            rstudioapi_0.14           farver_2.1.1              fansi_1.0.4               splines_4.2.2            
# [11] codetools_0.2-19          R.methodsS3_1.8.2         sparseMatrixStats_1.10.0  dunn.test_1.3.5           polyclip_1.10-4          
# [16] jsonlite_1.8.4            R.oo_1.25.0               HDF5Array_1.26.0          spatstat.sparse_3.0-1     compiler_4.2.2           
# [21] httr_1.4.6                dqrng_0.3.0               backports_1.4.1           Matrix_1.5-4.1            fastmap_1.1.1            
# [26] lazyeval_0.2.2            limma_3.54.2              cli_3.6.1                 htmltools_0.5.6           tools_4.2.2              
# [31] gtable_0.3.4              glue_1.6.2                GenomeInfoDbData_1.2.9    Rcpp_1.0.10               cellranger_1.1.0         
# [36] vctrs_0.6.3               rhdf5filters_1.10.1       crosstalk_1.2.0           DelayedMatrixStats_1.20.0 beachmat_2.14.2          
# [41] timechange_0.2.0          lifecycle_1.0.3           goftest_1.2-3             edgeR_3.40.2              zlibbioc_1.44.0          
# [46] scales_1.2.1              hms_1.1.3                 spatstat.utils_3.0-3      parallel_4.2.2            rhdf5_2.42.1             
# [51] yaml_2.3.7                stringi_1.7.12            checkmate_2.2.0           BiocParallel_1.32.6       rlang_1.1.1              
# [56] pkgconfig_2.0.3           bitops_1.0-7              pracma_2.4.2              lattice_0.20-45           tensor_1.5               
# [61] Rhdf5lib_1.20.0           labeling_0.4.2            htmlwidgets_1.6.2         tidyselect_1.2.0          plyr_1.8.8               
# [66] magrittr_2.0.3            R6_2.5.1                  magick_2.7.4              generics_0.1.3            DelayedArray_0.23.2      
# [71] mgcv_1.8-42               pillar_1.9.0              withr_2.5.0               abind_1.4-5               RCurl_1.98-1.12          
# [76] crayon_1.5.2              DropletUtils_1.18.1       utf8_1.2.3                tzdb_0.4.0                locfit_1.5-9.7           
# [81] grid_4.2.2                data.table_1.14.8         digest_0.6.33             R.utils_2.12.2            munsell_0.5.0            
# [86] viridisLite_0.4.2        
