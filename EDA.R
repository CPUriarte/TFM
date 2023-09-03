# LIBRARIES ----
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
responses$n_cells <- cell_count_vector

# Limpieza de la memoria
responses$Response <- droplevels(responses$Response, "Not evaluable/NE")
levels(responses$Response) <- c("CR", "PR", "PD", "SD")

rm(list = c("cell_counts", "NAs", "cell_count_vector", "images_with_NAs", "new_order", "response_vector"))

SPIAT_tifs <- unique(eda_df$Image)

############################################ EXPLORATORY DATA ANALYSIS (EDA) -----------------------------------------
# TIF RETRIEVER (retrieve tifs with patient ID) ---- 
# Number of images per patient and TIF index
tif_index <- eda_df %>%
  group_by(ID) %>%
  summarise(
    Count = n_distinct(Image),
    Images = list(unique(Image))
  )

tif_index$Matched_Indices <- sapply(tif_index$Images, function(images_list) which(SPIAT_tifs %in% images_list))

# To retrieve indices for a specific ID:
tif_index_retriever <- tif_index[tif_index$ID == 3501, "Matched_Indices"][[1]]
print(tif_index_retriever)

# 2. EDA ----
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

## 2.1. NUMBER OF CELLS (IMG & ID) ----
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

## 2.2. AGE (ID) ----
# Overall age distribution
responses %>%
  group_by(ID) %>%
  ggplot(aes(x = Age)) +
  geom_histogram(aes(y = ..density..), fill = "steelblue", color = "white", bins = 30) +
  geom_density(aes(y = ..density..), color = "red") +  
  labs(title = "Distribution of Age",
       x = "Age",
       y = "Count")

### 2.2.1. Age distribution statistics ----
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

### 2.3.1. Sex distribution statistics ----
# LaTeX: P(X=k) = \frac{\binom{N}{n} \binom{M}{k} \binom{N-M}{n-k}}{\binom{N}{n}}
choose(31, 10) / choose(46, 10)

# Fisher exact test
age_fishertest <- fisher.test(table(responses$Sex, responses$Response))
age_fishertest$p.value

## 2.4 MARKER INTENSITY AVERAGED BY ID ----
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

### 2.4.1. Normality test ----
markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Shapiro-Wilk test for normality
shapiro_test_by_group <- function(marker) {
  eda_by_ID %>%
    group_by(Response) %>%
    summarise(p_value = shapiro.test(get(marker))$p.value) %>%
    mutate(marker = marker)
}

eda_mIntensity_shapiro_byID <- bind_rows(lapply(markers, shapiro_test_by_group))

### 2.4.2. Homoscedasticity and KW test ----
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

### 2.4.3. Dunn test ----
dunn_results <- dunnTest(CD3 ~ Response, data = eda_by_ID, method = "bh")

eda_mIntensity_dunn <- data.frame(
  Comparison = dunn_results$res$Comparison,
  Z = dunn_results$res$Z,
  P.unadj = dunn_results$res$`P.unadj`,
  P.adj = dunn_results$res$`P.adj`
)

# Visualize
p <- ggplot(eda_by_ID, aes(x = Response, y = CD3)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(title = "CD3 Intensity by Response Group", y = NULL, x = NULL) +
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

### 2.4.4. Graphics and plots ----
# Correlation plot by ID
combined_data <- merge(eda_df, responses, by = "ID")

# Subset
subset_data <- combined_data[, c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK", "X", "Y", "Age")]

# Calculate
correlation_matrix <- cor(subset_data, method = "spearman", use = "pairwise.complete.obs")

# Correlation matrix
corrplot(correlation_matrix, method = "number", type = "lower", order = "FPC")

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
## 3.1. DETERMINE CUTOFF ----
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

## 3.1. DETERMINE CUTOFF V3 ----
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

## 3.2 APPLY CUTOFF ----
immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

for(marker in immune_markers) {
  cutoff_col <- paste0("cutoff_", marker)
  merged_df[, marker] <- ifelse(merged_df[, marker] >= merged_df[, cutoff_col], 1, 0)
}

### 3.2.3. Version 3.0 ---- 
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

## 3.3. ASSIGN PHENOTYPES ----
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

# Remove intermediate variables to clear memory
rm(cutoff_col, duplicate_rows, merged_df, original_rows, original_cols, result, phenotype_strings)

## 3.4. SPATIAL PHENOTYPES PLOT ----
# Phenotype distribution dotplot
# Imagenes interesantes (CR: 53, 89 | Other: 65)
tif <- SPIAT_tifs[53]

p <- final_df %>%
  filter(Image == tif) %>%
  ggplot(aes(x=X, y=Y, color=Phenotype)) + 
  geom_point(alpha=0.7) + 
  labs(title=paste("Dotplot for", tif), x = "X coordinate", y = "Y coordinate") + 
  theme_minimal() + 
  scale_color_discrete(name="Phenotypes") +
  coord_cartesian(xlim=c(0, 1300), ylim=c(0, 1300))

ggplotly(p)

## 3.5. THRESHOLDING QUALITY CONTROL ----

long_eda_df <- eda_df[eda_df$Image == tif, ] %>%
  pivot_longer(cols = colnames(eda_df[, 6:11]), names_to = "marker", values_to = "intensity")

long_eda_df$cutoff <- sapply(1:nrow(long_eda_df), function(i) {
  img_threshs[img_threshs$Image == tif, paste0("cutoff_", long_eda_df$marker[i])]
})

long_eda_df$cutoff <- as.numeric(long_eda_df$cutoff)

ggplot(long_eda_df, aes(x = intensity)) +
  geom_density(fill = "aquamarine") +
  geom_vline(aes(xintercept = cutoff), color = "red", linetype = "dashed") +
  labs(title = paste("Density distribution of markers for", tif),
       subtitle = paste("Cutoff values are displayed as red dashed lines"),
       x = "Intensity",
       y = "Density") +
  facet_wrap(~marker, scales = "free") +
  theme_minimal()

# Marker intensity dotplot
eda_df %>%
  filter(Image == tif) %>%
  ggplot(aes(x=X, y=Y, color=CK)) +
  geom_point(alpha=0.6) +
  scale_color_gradientn(n.breaks = 3, colors = terrain.colors(93)) +
  coord_cartesian(xlim=c(0, 1300), ylim=c(0, 1300)) +
  theme_minimal()


## 3.6. PHENOTYPE GROUP ID AVGED PROPORTIONS BARPLOT ----
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

########################################### SPIAT FEATURE EXTRACTION -----
# EXTRACT - PREDICTED PHENOTYPES ----
# Initialize a list to store processed data for each image
spiat_predicted_phenotypes <- list()

# Create a custom function for extracting the necessary details from the SpatialExperiment object
extract_data_from_spatial_experiment <- function(spe_object) {
  spatial_coords <- spe_object@int_colData@listData$spatialCoords
  phenotype <- spe_object@colData@listData$Phenotype
  
  return(list(spatial_coords = spatial_coords, phenotype = phenotype))
}

# Use pblapply (from pbapply package) to loop through each image in SPIAT_tifs with a progress bar
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

# If you want to name each list entry with its corresponding image name (optional)
names(spiat_predicted_phenotypes) <- sapply(spiat_predicted_phenotypes, function(x) x$image_name)

## Phenotype counts ---- 
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

# Cleaning up
rm(col_index, current_cols, i, image, image_name, image_phenotype_counts,
   ordered_cols, phenotype, phenotype_vector, phenotypes_in_image, row_index, 
   unique_phenotypes, primary_weights, secondary_weights, combined_weights, 
   get_primary_weight, get_secondary_weight, final_order)

## Phenotype image sampled proportion barplot ----
# Sample
sampled_images <- spiat_pheno_counts %>%
  sample_n(4) %>%
  pull(Image)

# Prepare
plot_data <- spiat_pheno_counts %>%
  filter(Image %in% sampled_images) %>%
  select(-ID) %>%
  gather(Phenotype, Count, -Image, -Response) %>%
  filter(!is.na(Count)) %>%
  group_by(Image) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  mutate(Label = paste(Response, "(", Image, ")", sep=""))  # New label

# Plot
spiat_hbarplot_proposample <- plot_data %>%
  ggplot(aes(x = Label, y = Proportion, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Proportion of Phenotypes for Randomly Sampled Images", y = "Proportion", x = "Response (Image Name)") +
  scale_fill_viridis_d()

ggplotly(spiat_hbarplot_proposample, tooltip = c("x", "y", "fill"))

## Predicted phenotype results visualization ----
# Averaged Proportion of Phenotypes per Response Group (IMAGES)
avg_plot_data <- spiat_pheno_counts %>%
  select(-Image, -ID) %>%
  gather(Phenotype, Count, -Response) %>%
  filter(!is.na(Count)) %>%
  group_by(Response, Phenotype) %>%
  summarize(TotalCount = sum(Count, na.rm = TRUE)) %>%
  mutate(Proportion = TotalCount / sum(TotalCount)) %>%
  ungroup()

spiat_hbarplot_proporesponse <- avg_plot_data %>%
  ggplot(aes(x = Response, y = Proportion, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Averaged Proportion of Phenotypes from Images per Response Group", y = "Proportion", x = "Response") +
  scale_fill_viridis_d()

plot(spiat_hbarplot_proporesponse)

ggplotly(spiat_hbarplot_proporesponse, tooltip = c("x", "y", "fill"))

# Averaged Proportion of Phenotypes per Response Group (PATIENTS)
# Step 1: Reshape the data
spiat_long_format <- spiat_pheno_counts %>%
  select(-Image) %>%
  gather(key = 'spiat_Phenotype', value = 'spiat_Count', -ID, -Response)

# Group the data by ID and then summarize the counts
spiat_total_cells_per_ID <- spiat_long_format %>%
  group_by(ID) %>%
  summarise(spiat_TotalCells = sum(spiat_Count, na.rm = TRUE))

# Create Long Data with Proportions
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

# Generate Plot
spiat_plot <- ggplot(spiat_aggregated_pheno, aes(x = spiat_Response, y = spiat_GroupProp, fill = spiat_Phenotype)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Proportion of Each Phenotype by Response Group",
    x = "Response Group",
    y = "Proportion",
    fill = "Phenotype"
  )

# Create interactive plot
ggplotly(spiat_plot)

## Statistics - Predicted phenotypes (Summation aggregation)----
spiat_pheno_counts[,4:67] <- apply(spiat_pheno_counts[,4:67], 2, function(x) ifelse(is.na(x), 0, x))

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

markers <- unique(spiat_phenoKW$Phenotype)

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

# Step 1: Perform Levene's test first for all markers and store results
for (marker in markers) {
  lev_test <- leveneTest(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == marker, ])$`Pr(>F)`[1]
  spiat_pheno_levene <- rbind(spiat_pheno_levene, data.frame(Phenotype = marker, levene = lev_test))
}

# Print or save Levene's test results if needed
print(spiat_pheno_levene)

# Step 2: Perform Kruskal-Wallis test for those markers that pass Levene's test (p >= 0.05)
for (marker in markers) {
  if (spiat_pheno_levene[spiat_pheno_levene$Phenotype == marker, "levene"] >= 0.05) {
    p_val <- kruskal.test(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == marker, ])$p.value
    spiat_phenoKWfdr <- rbind(spiat_phenoKWfdr, data.frame(Phenotype = marker, p.value = p_val))
  }
}

# Adjust p-values
spiat_phenoKWfdr$adj.p <- p.adjust(spiat_phenoKWfdr$p.value, method = "fdr")

# Step 3: Perform Dunn's test for those with Kruskal-Wallis p < 0.05
for (marker in spiat_phenoKWfdr$Phenotype[spiat_phenoKWfdr$p.value < 0.05]) {
  dunn_results <- dunnTest(Value ~ Response, data = spiat_phenoKW[spiat_phenoKW$Phenotype == marker, ], method = "bonferroni")
  dunn_df <- data.frame(
    Phenotype = rep(marker, nrow(dunn_results$res)),
    Comparison = dunn_results$res$Comparison,
    Z = dunn_results$res$Z,
    P.unadj = dunn_results$res$`P.unadj`,
    P.adj = dunn_results$res$`P.adj`
  )
  spiat_pheno_dunn <- rbind(spiat_pheno_dunn, dunn_df)
}

rm(spiat_pheno_temp, temp_df, marker, markers, p_val, dunn_df, dunn_results, lev_test)

# EXTRACT - PHENOTYPE LOCATIONS ----
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

# EXTRACT - MARKER CUTOFF P/N COUNTS ----
# This is the thing to iterate over all images
# Initialize a list to store processed data for each image
spiat_PNcounts <- list()

# Create a custom function for extracting the necessary details
extract_data_from_spatial_experiment <- function(plots_list) {
  
  for (marker in names(plots_list)) {
    P_count <- sum(plots_list[[marker]]$data$intensity == "P")
    N_count <- sum(plots_list[[marker]]$data$intensity == "N")
    Marker_name <- marker
    L <- list(M = Marker_name, P = P_count, N = N_count)
  }
  
  return(list(M = L$M, P = L$P, N = L$N))
}

# Use pblapply (from pbapply package) to loop through each image in SPIAT_tifs with a progress bar
spiat_PNcounts <- pblapply(SPIAT_tifs, function(selected_image) {
  
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
  
  marker_names <- predicted_image@rowRanges@partitioning@NAMES[predicted_image@rowRanges@partitioning@NAMES != "DAPI"]
  plots_list <- list()
  
  for (marker in marker_names) {
    tryCatch({
      marker_plot <- marker_intensity_boxplot(predicted_image, marker)
      plots_list[[marker]] <- marker_plot
    }, error = function(e) {
      message(paste("Error for marker", marker, ":", e$message))
    })
  }
  
  # Extract the required data from the SpatialExperiment object
  extracted_data <- extract_data_from_spatial_experiment(plots_list)
  
  # Return the data along with the image name
  return(list(Image = selected_image, M = extracted_data$M, P = extracted_data$P, N = extracted_data$N))
})

# If you want to name each list entry with its corresponding image name (optional)
names(spiat_predicted_phenotypes) <- sapply(spiat_predicted_phenotypes, function(x) x$Image)

# EXTRACT - MARKER HEATMAP ----
p <- ggplot(eda_df[eda_df$ID == unique(eda_df$ID)[1],], aes(x = X, y = Y)) +
  stat_summary_2d(aes(z = LAG3, fill = after_stat(value)), fun = median, bins = 200) +
  facet_wrap(~ Image) +
  scale_fill_gradient(low="yellow", high="purple") +
  theme_minimal() +
  labs(fill = "Average PD1 Intensity")

ggplotly(p)
 #! PROBLEM - Handling intensity on borders of tumor !#
## Extraction/Aggregation for 1 patient ----
bin_n <- 200
bin_width <- (max(eda_df[eda_df$ID == unique(eda_df$ID)[1],]$X) - min(eda_df[eda_df$ID == unique(eda_df$ID)[1],]$X)) / bin_n

# Convert the dataframe to long format
eda_df_long <- eda_df %>%
  pivot_longer(cols = c(PD1, CD8, CD3, TIM3, LAG3, CK),
               names_to = "marker",
               values_to = "value")

# Step 1: Create bins
bin_example <- eda_df_long[eda_df_long$ID == unique(eda_df_long$ID)[1],] %>%
  mutate(
    bin_x = floor(X / bin_width),
    bin_y = floor(Y / bin_width)
  )

# Step 2: Group by bin, Image, and marker
resultado <- bin_example %>%
  group_by(Image, marker, bin_x, bin_y) %>%
  summarise(avg_value = median(value, na.rm = TRUE))

# Step 3: Aggregate across images by bin and marker
resultado <- resultado %>%
  group_by(marker, bin_x, bin_y) %>%
  summarise(median_value_across_images = median(avg_value, na.rm = TRUE))

resultado <- resultado %>%
  mutate(
    X = (bin_x + 0.5) * bin_width,
    Y = (bin_y + 0.5) * bin_width,
  )

# Plot using ggplot2
p_resultado <- ggplot(resultado, aes(x = X, y = Y, fill = median_value_across_images)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "purple", name = NULL) +
  theme_minimal() +
  labs(title = paste("Aggregated Median Intensity Across Images Of Patient", unique(eda_df_long$ID)[1])) +
  facet_wrap(~ marker, ncol = 3)

ggplotly(p_resultado)

rm(bin_width, resultado, bin_example, eda_df_long)

## Extraction/Aggregation for all patients ----
# Define the bin number
bin_n <- 4

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
  
  # Store the result in the list
  tiled_markerheatmap_xID[[as.character(id)]] <- resultado
}

# Initialize an empty list to store the results for each response
tiled_markerheatmap_xResponse <- list()

# Iterate through the unique response groups
for (response in unique(responses$Response)) {
  # Get the IDs associated with the current response group
  matching_ids <- as.character(responses$ID[responses$Response == response])
  
  # Filter tiled_markerheatmap_xID by the matching IDs but maintain the individual structure
  tiled_markerheatmap_xResponse[[response]] <- tiled_markerheatmap_xID[matching_ids]
}

## Pick Marker ----
marker_of_interest <- "CK"
## Visualize aggregated median intensity per response groups ----
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

# Plot
p_resultado <- ggplot(all_data, aes(x = X, y = Y, fill = median_value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", name = NULL) +
  theme_minimal() +
  labs(title = paste("Median Intensity Of", marker_of_interest, "Across Response")) +
  facet_wrap(~ Response, ncol = 2)

ggplotly(p_resultado)

suppressWarnings({
  rm(bin_data, resultado, combined_data, response_data, tiled_markerheatmap_xID)
})

## Statistics - Marker heatmaps
# ----
# Reshape the data for ggplot2
long_df <- eda_df[eda_df$Image == SPIAT_tifs[1], ] %>% 
  select(-Image, -ID, -X, -Y) %>% 
  gather(key = "Marker", value = "Value", -TIM3)

# Create the plot
ggplot(long_df, aes(x = TIM3, y = Value)) +
  geom_point(aes(color = Marker)) +
  geom_smooth(method = "lm", se = TRUE, aes(color = "black")) +  # se=TRUE adds the confidence interval
  facet_wrap(~ Marker) +
  theme_minimal()