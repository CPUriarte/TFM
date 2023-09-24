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
library(igraph)
############################################ LOADING DATA ----
# Raw measurements exported from QuPath
raw_measurements <- read.table("Data/raw_measurements.tsv", header = TRUE, sep = "\t")

# Raw associated clinical data
responses <- read.csv("Data/Responses.csv", header = TRUE, sep = ",")

# PROCESSING ----
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

############################################ EDA ----
# Image retriever by ID ----
ID <- 303 # Other interesting cases(1101, 502, 303) ----
tif_index_retriever <- tif_index[tif_index$ID == ID, "Matched_Indices"][[1]]
print(paste("SPIAT_tifs indexes for patient", ID, "are:", tif_index_retriever))
############################################  SECTION 1 / Images and cells ----
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
# VISUALS - Cells per img and ID (QQ-plots) ----
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

############################################  SECTION 2 / Age ----
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
# VISUALS - Age per ID (Barplot) ----
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
############################################  SECTION 3 / Sex ----
# LaTeX: P(X=k) = \frac{\binom{N}{n} \binom{M}{k} \binom{N-M}{n-k}}{\binom{N}{n}}
choose(31, 10) / choose(46, 10)

# Fisher exact test
age_fishertest <- fisher.test(table(responses$Sex, responses$Response))
age_fishertest$p.value

############################################  SECTION 4 / Marker intensity ----
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

# STATISTICS -  Shapiro, Levene, KW, Dunn ----
markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Shapiro-Wilk
shapiro_test_by_group <- function(marker) {
  eda_by_ID %>%
    group_by(Response) %>%
    summarise(p_value = shapiro.test(get(marker))$p.value) %>%
    mutate(marker = marker)
}
eda_mIntensity_shapiro_byID <- bind_rows(lapply(markers, shapiro_test_by_group))

# Levene and KW test
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

# Run tests and bind results
eda_mIntensity_levene_KWfdr <- bind_rows(lapply(markers, test_for_marker))

# CD3 Dunn test
dunn_results <- dunnTest(CD3 ~ Response, data = eda_by_ID, method = "bh")

eda_mIntensity_dunn <- data.frame(
  Comparison = dunn_results$res$Comparison,
  Z = dunn_results$res$Z,
  P.unadj = dunn_results$res$`P.unadj`,
  P.adj = dunn_results$res$`P.adj`
)

ggplot(eda_mIntensity_dunn, aes(x = reorder(Comparison, -Z), y = Z)) +
  geom_bar(stat = "identity", aes(fill = P.adj <= 0.05)) + 
  geom_text(aes(label = sprintf("%.2f", P.adj)), vjust = -0.5) +
  scale_fill_manual(values = c("red", "gray"), 
                    name = "Significance", 
                    labels = c("p <= 0.05", "p > 0.05")) +
  labs(title = "Dunn's Test Results for CD3",
       x = "Group Comparison",
       y = "Z-Score") +
  theme_minimal()

eda_mIntensity_levene_KWfdr <- eda_mIntensity_levene_KWfdr %>%
  arrange(kw_p_value)
print(eda_mIntensity_levene_KWfdr)

eda_mIntensity_dunn <- eda_mIntensity_dunn %>%
  arrange(P.adj)
print(eda_mIntensity_dunn)

# CD8 Dunn test
dunn_results <- dunnTest(CD8 ~ Response, data = eda_by_ID, method = "bh")

eda_mIntensity_dunn <- data.frame(
  Comparison = dunn_results$res$Comparison,
  Z = dunn_results$res$Z,
  P.unadj = dunn_results$res$`P.unadj`,
  P.adj = dunn_results$res$`P.adj`
)

eda_mIntensity_dunn <- eda_mIntensity_dunn %>%
  arrange(P.adj)

ggplot(eda_mIntensity_dunn, aes(x = reorder(Comparison, -Z), y = Z)) +
  geom_bar(stat = "identity", aes(fill = P.adj <= 0.05)) + 
  geom_text(aes(label = sprintf("%.2f", P.adj)), vjust = -0.5) +
  scale_fill_manual(values = c("red", "gray"), 
                    name = "Significance", 
                    labels = c("p <= 0.05", "p > 0.05")) +
  labs(title = "Dunn's Test Results for CD8",
       x = "Group Comparison",
       y = "Z-Score") +
  theme_minimal()

print(eda_mIntensity_dunn)

# VISUALS - Results ----
# Correlation plot by ID
combined_data <- merge(eda_df, responses, by = "ID")

# Subset markers
subset_data <- combined_data[, c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")]

# Compute with spearman's
correlation_matrix <- cor(subset_data, method = "spearman", use = "complete.obs")

# Correlation matrix
corrplot(correlation_matrix, type = "lower", method = "color", addCoef.col = 'black', tl.pos = 'd')
corrplot(correlation_matrix, type = "upper", method = "ellipse", diag = F, tl.pos = 'n', cl.pos = 'n', add=T)

# Network analysis
cor_matrix <- matrix(correlation_matrix, nrow=6, byrow=TRUE)

rownames(cor_matrix) <- colnames(cor_matrix) <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

# Create list on arbitrary treshold
edges <- which(abs(cor_matrix) > 0.5 & upper.tri(cor_matrix), arr.ind = TRUE)
edge_list <- data.frame(from = rownames(cor_matrix)[edges[, 1]], 
                        to = rownames(cor_matrix)[edges[, 2]], 
                        weight = abs(cor_matrix[edges]))

g <- graph_from_data_frame(edge_list, directed = FALSE)

# Plot graph with scaled edge width
plot(g, vertex.label=V(g)$name,
     vertex.size=30,
     edge.width=(E(g)$weight - min(E(g)$weight)) * 50,
     vertex.color="firebrick3", edge.color='orange', vertex.label.color = "ivory")

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
# (Version 1) - Determining cutoff ----
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

# (Version 2) - Determining cutoff ----
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
    
    # Determine the cutoff based on the first valley after the highest peak
    ycord_max_density <- max(intensity_density$y)
    xcord_max_density <- intensity_density$x[match(ycord_max_density, intensity_density$y)]
    valley_df <- data.frame(valley_xcords, valley_ycords) %>%
      filter(valley_xcords >= xcord_max_density, 
             valley_ycords <= 0.25 * ycord_max_density)
    
    return(valley_df$valley_xcords[1])
  }else{
  
  # You can adjust this logic based on how you want to use this metric for thresholding
  return(peaks_df$x[1])
  }
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

# PROCESSING - Applying cutoff ----
immune_markers <- c("PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")

for(marker in immune_markers) {
  cutoff_col <- paste0("cutoff_", marker)
  merged_df[, marker] <- ifelse(merged_df[, marker] >= merged_df[, cutoff_col], 1, 0)
}

# PROCESSING - Handling superposition ---- 
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

# PROCESSING - Assigning phenotypes ----
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

# VISUALS - Predicted phenotypes (Dotplot) ----
selected_image <- SPIAT_tifs[random()] # 1 to 90 possible images ----
selected_image <- SPIAT_tifs[65] # Select a specific image
# Join final_df and responses
final_df_joined <- final_df %>%
  left_join(responses, by = "ID")

p <- final_df_joined %>%
  filter(Image == selected_image) %>%
  ggplot(aes(x = X, y = Y, color = Phenotype)) +
  geom_point(alpha = 0.7) +
  labs(
    title = paste(
      selected_image,
      "(", first(unique(final_df_joined$Response[final_df_joined$Image == 1])), ")",
      sep = " "
    ),
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_minimal() +
  scale_color_discrete(name = "Phenotypes") +
  coord_cartesian(xlim = c(0, 1300), ylim = c(0, 1300))

ggplotly(p)

# VISUALS - Cutoff quality (IMAGES) ----

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


# VISUALS - Phenotype proportion per response group (Stacked Barplot) ----
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
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
# [1] igraph_1.5.1                spatstat_3.0-6             
# [3] spatstat.linnet_3.1-1       spatstat.model_3.2-4       
# [5] rpart_4.1.19                spatstat.explore_3.2-1     
# [7] nlme_3.1-162                spatstat.random_3.1-5      
# [9] spatstat.geom_3.2-1         spatstat.data_3.0-1        
# [11] corrplot_0.92               ggsignif_0.6.4             
# [13] FSA_0.9.4                   ggridges_0.5.4             
# [15] pbapply_1.7-0               moments_0.14.1             
# [17] car_3.1-2                   carData_3.0-5              
# [19] readxl_1.4.2                reshape2_1.4.4             
# [21] BBmisc_1.13                 plotly_4.10.2              
# [23] gridExtra_2.3               lubridate_1.9.2            
# [25] forcats_1.0.0               stringr_1.5.0              
# [27] dplyr_1.1.2                 purrr_1.0.1                
# [29] readr_2.1.4                 tidyr_1.3.0                
# [31] tibble_3.2.1                ggplot2_3.4.3              
# [33] tidyverse_2.0.0             SPIAT_1.0.4                
# [35] SpatialExperiment_1.8.1     SingleCellExperiment_1.20.1
# [37] SummarizedExperiment_1.28.0 Biobase_2.58.0             
# [39] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
# [41] IRanges_2.32.0              S4Vectors_0.36.2           
# [43] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
# [45] matrixStats_0.63.0         
# 
# loaded via a namespace (and not attached):
# [1] colorspace_2.1-0          rjson_0.2.21             
# [3] deldir_1.0-9              ellipsis_0.3.2           
# [5] scuttle_1.8.4             XVector_0.38.0           
# [7] rstudioapi_0.14           farver_2.1.1             
# [9] scico_1.5.0               fansi_1.0.4              
# [11] splines_4.2.2             codetools_0.2-19         
# [13] R.methodsS3_1.8.2         sparseMatrixStats_1.10.0 
# [15] dunn.test_1.3.5           polyclip_1.10-4          
# [17] jsonlite_1.8.4            R.oo_1.25.0              
# [19] HDF5Array_1.26.0          spatstat.sparse_3.0-1    
# [21] compiler_4.2.2            httr_1.4.6               
# [23] dqrng_0.3.0               backports_1.4.1          
# [25] Matrix_1.5-4.1            fastmap_1.1.1            
# [27] lazyeval_0.2.2            limma_3.54.2             
# [29] cli_3.6.1                 htmltools_0.5.6          
# [31] tools_4.2.2               gtable_0.3.4             
# [33] glue_1.6.2                GenomeInfoDbData_1.2.9   
# [35] Rcpp_1.0.10               cellranger_1.1.0         
# [37] vctrs_0.6.3               rhdf5filters_1.10.1      
# [39] crosstalk_1.2.0           DelayedMatrixStats_1.20.0
# [41] beachmat_2.14.2           timechange_0.2.0         
# [43] lifecycle_1.0.3           goftest_1.2-3            
# [45] edgeR_3.40.2              zlibbioc_1.44.0          
# [47] scales_1.2.1              hms_1.1.3                
# [49] spatstat.utils_3.0-3      parallel_4.2.2           
# [51] rhdf5_2.42.1              yaml_2.3.7               
# [53] stringi_1.7.12            checkmate_2.2.0          
# [55] BiocParallel_1.32.6       rlang_1.1.1              
# [57] pkgconfig_2.0.3           bitops_1.0-7             
# [59] pracma_2.4.2              lattice_0.20-45          
# [61] tensor_1.5                Rhdf5lib_1.20.0          
# [63] labeling_0.4.2            htmlwidgets_1.6.2        
# [65] tidyselect_1.2.0          plyr_1.8.8               
# [67] magrittr_2.0.3            R6_2.5.1                 
# [69] magick_2.7.4              generics_0.1.3           
# [71] DelayedArray_0.23.2       mgcv_1.8-42              
# [73] pillar_1.9.0              withr_2.5.0              
# [75] abind_1.4-5               RCurl_1.98-1.12          
# [77] crayon_1.5.2              DropletUtils_1.18.1      
# [79] utf8_1.2.3                tzdb_0.4.0               
# [81] locfit_1.5-9.7            grid_4.2.2               
# [83] data.table_1.14.8         digest_0.6.33            
# [85] R.utils_2.12.2            munsell_0.5.0            
# [87] viridisLite_0.4.2