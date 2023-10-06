# Overview ----
# Title: Data pre-processing.
# Date: 20/05/2023
# Author: Cyril P
# Director: Carlos de Andrea
# Institution: University of Navarra
# Load raw data ----
# Raw measurements exported from QuPath
raw_measurements <- read.table("Data/raw_measurements.tsv", header = TRUE, sep = "\t")

# Raw associated clinical data
responses <- read.csv("Data/Responses.csv", header = TRUE, sep = ",")

# Format datasets ----
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
responses$Total_cells <- cell_count_vector

# Dropping non-evaluable patients
responses$Response <- droplevels(responses$Response, "Not evaluable/NE")
levels(responses$Response) <- c("CR", "PR", "PD", "SD")
SPIAT_tifs <- unique(eda_df$Image) # Unique list of SPIAT curated images

IntensityMatrix <- eda_df
ClinicalData <- responses

# Clear intermediary variables
rm(list = c("cell_counts", "NAs", "cell_count_vector", "images_with_NAs", "new_order", "response_vector", "eda_df"))

# Clear console for results
cat("\014")
summary(ClinicalData)
summary(IntensityMatrix[,c("DAPI", "PD1", "CD8", "CD3", "TIM3", "LAG3", "CK")])
