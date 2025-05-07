library(tidyverse)

attributes_url <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
attributes <- read_tsv(attributes_url)

phenotypes_url <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
phenotypes <- read_tsv(phenotypes_url)


brain_samples <- attributes |>
  ggplot(brain_samples, aes(x = SMTSD, y = SMRIN)) +
  geom_boxplot() +
  labs(title = "RNA Integrity by Brain Region", x = "Tissue", y = "RIN Score")

heart_samples <- attributes |>
  filter(SMTSD %in% c("Heart - Atrial Appendage", "Heart - Left Ventricle"))
ggplot(heart_samples, aes(x = SMTSD, y = SMRIN)) +
  geom_boxplot() +
  labs(title = "RNA Integrity by Heart Region", x = "Tissue", y = "RIN Score")


library(ggplot2)
library(readxl)

# --- read gct data
read_gct <- function(file) {
  lines <- readLines(file)
  
  # First line is version
  version <- lines[1]
  
  # Second line: number of rows and columns
  dims <- strsplit(lines[2], "\t")[[1]]
  n_rows <- as.integer(dims[1])
  n_cols <- as.integer(dims[2])
  
  # Read the actual data
  df <- read.delim(file, skip = 2, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract expression matrix
  expr <- as.matrix(df[, -(1:2)])  # Remove Name and Description
  rownames(expr) <- df$Name
  
  return(list(
    version = version,
    data = df,
    expression = expr
  ))
}

gtex_data <- read_gct("~/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_reads.gct")
gtex_attr <- read_xlsx("C:/Users/bdanc/Downloads/GTEx_Analysis_v10_Annotations_SampleAttributesDD.xlsx")
url <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
destfile <- "GTEx_SampleAttributesDS.txt"
download.file(url, destfile, method = "auto")
library(readr)
gtex_metadata <- read_tsv("GTEx_SampleAttributesDS.txt")

brain_IDs <- gsub("-", ".", gtex_metadata$SAMPID[gtex_metadata$SMTS=="Brain"])
heart_IDs <- gsub("-", ".", gtex_metadata$SAMPID[gtex_metadata$SMTS=="Heart"])
gtex_expr_bh <- gtex_data$expression[,colnames(gtex_data$expression) %in% c(brain_IDs,heart_IDs)]
tissue_df <- rbind(data.frame(tissue="Brain",sampleID = brain_IDs), data.frame(tissue="Heart",sampleID=heart_IDs) )
tissue_df <- tissue_df[match(colnames(gtex_expr_bh), tissue_df$sampleID), ]
gtex_expr_bh_filtered <- gtex_expr_bh[rowSums(gtex_expr_bh > 20) >= 20, ]
gene_vars <- apply(gtex_expr_bh_filtered, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:500])
gtex_expr_bh_filtered <- gtex_expr_bh_filtered[top_genes, ]
expr <- gtex_expr_bh_filtered
expr_t <- t(expr)
pca_result <- prcomp(expr_t, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Sample <- rownames(pca_df)
pca_df$Tissue <- tissue_df$tissue
ggplot(pca_df, aes(x = PC1, y = PC2, color=Tissue)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: First Two Principal Components", x = "PC1", y = "PC2")

expr_scaled <- scale(expr_t)
kmeans_result <- kmeans(expr_scaled, centers = 2)
cluster_df <- data.frame(
  Sample = rownames(expr_scaled),
  Cluster = as.factor(kmeans_result$cluster)
)
pca_df$Cluster_Kmean <- cluster_df$Cluster
ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster_Kmean)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "K-means Clustering on Expression Data")