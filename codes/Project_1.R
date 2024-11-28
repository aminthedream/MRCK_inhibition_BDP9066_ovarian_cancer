# Load required libraries
library(ggplot2)
library(pheatmap)
library(cluster)
library(Rtsne)
library(umap)
library(NbClust)
library(clusterProfiler)
library(org.Hs.eg.db)
library(randomForest)
library(ggplot2)
library(tidyr)
library(dplyr)
set.seed(123)


# 1. Data Preprocessing
counts <- read.csv("./data/GSE94304-expression.txt", sep = '\t', row.names = 1)
metadata <- read.csv("./data/GSE94304-metadata.txt", sep = '\t')

# Clean sample titles in metadata
metadata$Sample.Title <- gsub(pattern = "-", replacement = "_", metadata$Sample.Title)
metadata$Sample.Title <- gsub(pattern = "/", replacement = "_", metadata$Sample.Title)

# Update column names in counts to match metadata
colnames(counts) <- metadata$Sample.Title

# Filter low-count genes
keep <- rowSums(counts >= 10) >= ncol(counts)/4  # Keep genes with at least 10 counts in 25% of samples
filtered_counts <- counts[keep,]

# Log transform (add a small constant to avoid log(0))
log_counts <- log2(filtered_counts + 1)


pdf("Project_1_Plots.pdf",width =12, height =12)
# Simple normalization: divide by library size and multiply by 1e6 (CPM)
lib_sizes <- colSums(filtered_counts)
norm_counts <- t(t(filtered_counts) / lib_sizes * 1e6)
log_norm_counts <- log2(norm_counts + 1)

# 2. Exploratory Data Analysis
# K-means clustering with a fixed number of clusters
n_clusters <- 5  # You can adjust this number as needed
kmeans_result <- kmeans(t(log_norm_counts), centers = n_clusters)

# PCA
pca_result <- prcomp(t(log_norm_counts))
pca_df <- as.data.frame(pca_result$x[,1:2])
pca_df$CellSubtype <- metadata$cell.subtype
pca_df$Cluster <- factor(kmeans_result$cluster)

ggplot(pca_df, aes(x = PC1, y = PC2, color = CellSubtype, shape = Cluster)) + 
  geom_point(size = 3) + 
  geom_text(aes(label = rownames(pca_df)), vjust = 2, size = 2) +
  theme_minimal() +
  ggtitle("PCA of Ovarian Cancer Cell Lines") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = 1:n_clusters)
# Heatmap of top variable genes
annotation_col <- data.frame(
  CellSubtype = metadata$cell.subtype,
  Cluster = factor(kmeans_result$cluster)
)
rownames(annotation_col) <- colnames(log_norm_counts)

var_genes <- apply(log_norm_counts, 1, var)
top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:1000]
pheatmap(log_norm_counts[top_var_genes,], 
         scale = "row", 
         show_rownames = FALSE, 
         annotation_col = annotation_col,
         main = "Top 100 Variable Genes")
# 3. Unsupervised Machine Learning Methods


# t-SNE
tsne_result <- Rtsne(t(log_norm_counts), perplexity = 5, check_duplicates = FALSE)
tsne_df <- as.data.frame(tsne_result$Y)
tsne_df$CellSubtype <- metadata$cell.subtype  # Add this line

ggplot(tsne_df, aes(x = V1, y = V2, color = CellSubtype)) + 
  geom_point(size = 3) + 
  geom_text(aes(label = rownames(tsne_df)), vjust = 2) +
  theme_minimal() + 
  labs(color = "Cell Subtype") +
  ggtitle("t-SNE of Ovarian Cancer Cell Lines") +
  scale_color_brewer(palette = "Set1")

# UMAP
umap_result <- umap(t(log_norm_counts))
umap_df <- as.data.frame(umap_result$layout)
umap_df$CellSubtype <- metadata$cell.subtype
umap_df$Cluster <- factor(kmeans_result$cluster)

ggplot(umap_df, aes(x = V1, y = V2, color = CellSubtype, shape = Cluster)) + 
  geom_point(size = 3) + 
  geom_text(aes(label = rownames(umap_df)), vjust = 2, size = 2) +
  theme_minimal() + 
  labs(color = "Cell Subtype", shape = "Cluster") +
  ggtitle("UMAP of Ovarian Cancer Cell Lines") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = 1:n_clusters)
# 4. Functional Analysis

# Get differentially expressed genes for each cluster
cluster_degs <- list()
for(i in 1:n_clusters) {
  cluster_samples <- names(kmeans_result$cluster[kmeans_result$cluster == i])
  cluster_mean <- rowMeans(log_norm_counts[, cluster_samples])
  other_mean <- rowMeans(log_norm_counts[, !colnames(log_norm_counts) %in% cluster_samples])
  fc <- cluster_mean - other_mean
  cluster_degs[[i]] <- names(sort(abs(fc), decreasing = TRUE))[1:100]
}

# Perform GO enrichment analysis for each cluster
go_results <- list()
for(i in 1:n_clusters) {
  go_result <- enrichGO(gene = cluster_degs[[i]],
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
  go_results[[i]] <- go_result
}

# Print top 5 GO terms for each cluster
for(i in 1:n_clusters) {
  cat("Cluster", i, "top GO terms:\n")
  print(head(go_results[[i]]@result$Description, 5))
  cat("\n")
}


###step 2

# Assuming log_norm_counts and kmeans_result are from your previous analysis

# Function to find potential biomarkers for each cluster
find_biomarkers <- function(cluster_number) {
  # Create binary classification: this cluster vs all others
  cluster_indicator <- ifelse(kmeans_result$cluster == cluster_number, 1, 0)
  
  # Calculate fold change and p-value for each gene
  results <- data.frame(
    Gene = rownames(log_norm_counts),
    FoldChange = rep(0, nrow(log_norm_counts)),
    PValue = rep(1, nrow(log_norm_counts))
  )
  
  for (i in 1:nrow(log_norm_counts)) {
    gene_expr <- log_norm_counts[i,]
    fc <- mean(gene_expr[cluster_indicator == 1]) - mean(gene_expr[cluster_indicator == 0])
    p_value <- t.test(gene_expr[cluster_indicator == 1], gene_expr[cluster_indicator == 0])$p.value
    
    results$FoldChange[i] <- fc
    results$PValue[i] <- p_value
  }
  
  # Sort by absolute fold change and filter by p-value
  results <- results[order(abs(results$FoldChange), decreasing = TRUE), ]
  results <- results[results$PValue < 0.05, ]
  
  # Return top 10 genes
  return(head(results, 10))
}

# Find potential biomarkers for each cluster
biomarkers_list <- lapply(1:n_clusters, find_biomarkers)

# Print top 5 potential biomarkers for each cluster
for (i in 1:n_clusters) {
  cat("Top 5 potential biomarkers for Cluster", i, ":\n")
  print(head(biomarkers_list[[i]]$Gene, 5))
  cat("\n")
}

# Visualize expression of top biomarker for each cluster
plot_top_biomarkers <- function() {
  top_biomarkers <- sapply(biomarkers_list, function(x) x$Gene[3])
  plot_data <- log_norm_counts[top_biomarkers, ]
  plot_data <- as.data.frame(t(plot_data))
  plot_data$Cluster <- factor(kmeans_result$cluster)
  plot_data_long <- pivot_longer(plot_data, cols = -Cluster, names_to = "Gene", values_to = "Expression")
  
  ggplot(plot_data_long, aes(x = Cluster, y = Expression, fill = Cluster)) +
    geom_boxplot() +
    facet_wrap(~ Gene, scales = "free_y") +
    theme_minimal() +
    labs(title = "Expression of Top Biomarker Genes Across Clusters",
         y = "Log2 Normalized Expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Generate and display the plot
biomarker_plot <- plot_top_biomarkers()
print(biomarker_plot)






#####another biomarker

plot_top_biomarkers <- function() {
  top_biomarkers <- sapply(biomarkers_list, function(x) x$Gene[1])
  plot_data <- log_norm_counts[top_biomarkers, ]
  plot_data <- as.data.frame(t(plot_data))
  plot_data$Cluster <- factor(kmeans_result$cluster)
  plot_data$CellSubtype <- metadata$cell.subtype
  plot_data_long <- pivot_longer(plot_data, cols = -c(Cluster, CellSubtype), names_to = "Gene", values_to = "Expression")
  
  ggplot(plot_data_long, aes(x = Cluster, y = Expression, fill = CellSubtype)) +
    geom_boxplot() +
    facet_wrap(~ Gene, scales = "free_y") +
    theme_minimal() +
    labs(title = "Expression of Top Biomarker Genes Across Clusters",
         y = "Log2 Normalized Expression",
         x = "Cluster") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set1")
}

biomarker_plot <- plot_top_biomarkers()
print(biomarker_plot)
#################


##################NMF
library(NMF)

# Assuming log_norm_counts is your normalized, log-transformed expression data
# Ensure data is non-negative
nmf_data <- abs(log_norm_counts)

# Determine the optimal number of components (ranks)
estim.r <- nmf(nmf_data, 2:10, nrun = 10, seed = 123)
plot(estim.r)

# Choose the optimal number of components based on the plot
# Let's say we choose 5 components
n_components <- 8

# Run NMF
nmf_result <- nmf(nmf_data, n_components, nrun = 50, seed = 123)

# Extract the basis matrix (metagenes)
metagenes <- basis(nmf_result)

# Extract the coefficient matrix (metasamples)
metasamples <- coef(nmf_result)


basis <- basis(nmf_result)
coeff <- coef(nmf_result)

coeff_df <- as.data.frame(coeff)
coeff_df$Sample <- rownames(coeff_df)

coeff_melt <- melt(coeff_df, id.vars = "Sample", variable.name = "Component", value.name = "Coefficient")




consensusmap(nmf_result)

top_genes <- extractFeatures(nmf_result, 20)  # Top 20 genes per cluster

# Assuming rownames(log_norm_counts) contains gene names
gene_names <- rownames(log_norm_counts)

# Function to convert indices to gene names
indices_to_names <- function(indices) {
  gene_names[indices]
}

# Convert top genes to names
top_genes_named <- lapply(top_genes, indices_to_names)

# Print top genes for each component
for (i in seq_along(top_genes_named)) {
  cat("Component", i, "top genes:\n")
  print(top_genes_named[[i]])
  cat("\n")
}



# Check if CD151 is in the top genes of any component
sapply(top_genes_named, function(x) "CD151" %in% x)

# List of some known cell surface markers (you can expand this list)
cell_surface_markers <- c("CD151", "CD44", "CD24", "EPCAM", "CD133", "CD47")

# Function to find cell surface markers in a gene list
find_markers <- function(gene_list) {
  intersect(gene_list, cell_surface_markers)
}

# Find cell surface markers in each component
markers_in_components <- lapply(top_genes_named, find_markers)

# Print results
for (i in seq_along(markers_in_components)) {
  cat("Component", i, "cell surface markers:\n")
  print(markers_in_components[[i]])
  cat("\n")
}



# Find the index of CD151
cd151_index <- which(gene_names == "CD151")

if (length(cd151_index) > 0) {
  # Get CD151 coefficients across components
  cd151_coeffs <- metagenes[cd151_index,]
  
  # Print CD151 coefficients
  cat("CD151 coefficients across components:\n")
  print(cd151_coeffs)
  
  # Find the component where CD151 is most strongly represented
  strongest_component <- which.max(cd151_coeffs)
  cat("\nCD151 is most strongly represented in component:", strongest_component, "\n")
  
  # Print other top genes in this component
  cat("\nOther top genes in this component:\n")
  print(top_genes_named[[strongest_component]])
} else {
  cat("CD151 not found in the dataset.\n")
}




enrich_results <- lapply(top_genes_named, function(genes) {
  enrichGO(gene = genes,
           OrgDb = org.Hs.eg.db,
           keyType = "SYMBOL",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2)
})

# Print top 5 enriched terms for each component
for (i in seq_along(enrich_results)) {
  cat("Component", i, "top GO terms:\n")
  print(head(enrich_results[[i]]@result$Description, 5))
  cat("\n")
}

barplot(enrich_results[[4]], showCategory=10)






########
library(NMF)

# Assuming log_norm_counts is your normalized, log-transformed expression data
# Ensure data is non-negative
nmf_data1 <- abs(log_norm_counts)

# Determine the optimal number of components (ranks)
estim.r1 <- nmf(nmf_data1, 2:10, nrun = 10, seed = 123)
plot(estim.r1)

# Choose the optimal number of components based on the plot
# Let's say we choose 5 components
n_components1 <- 8

# Run NMF
nmf_result1 <- nmf(nmf_data1, n_components1, nrun = 50, seed = 123)

# Extract the basis matrix (metagenes)
metagenes1 <- basis(nmf_result1)

# Extract the coefficient matrix (metasamples)
metasamples1 <- coef(nmf_result1)


basis1 <- basis(nmf_result1)
coeff1 <- coef(nmf_result1)

coeff_df1 <- as.data.frame(coeff1)
coeff_df1$Sample <- rownames(coeff_df1)

coeff_melt1 <- melt(coeff_df1, id.vars = "Sample", variable.name = "Component", value.name = "Coefficient")

consensusmap(nmf_result1)

top_genes1 <- extractFeatures(nmf_result1, 20)  # Top 20 genes per cluster

# Assuming rownames(log_norm_counts) contains gene names
gene_names1 <- rownames(log_norm_counts)

# Function to convert indices to gene names
indices_to_names1 <- function(indices) {
  gene_names1[indices]
}

# Convert top genes to names
top_genes_named1 <- lapply(top_genes1, indices_to_names1)

# Print top genes for each component
for (i in seq_along(top_genes_named1)) {
  cat("Component", i, "top genes:\n")
  print(top_genes_named1[[i]])
  cat("\n")
}



# Find the index of CD151
cd151_index <- which(gene_names == "CD151")

if (length(cd151_index) > 0) {
  # Get CD151 coefficients across components
  cd151_coeffs <- metagenes1[cd151_index,]
  
  # Print CD151 coefficients
  cat("CD151 coefficients across components:\n")
  print(cd151_coeffs)
  
  # Find the component where CD151 is most strongly represented
  strongest_component <- which.max(cd151_coeffs)
  cat("\nCD151 is most strongly represented in component:", strongest_component, "\n")
  
  # Print other top genes in this component
  cat("\nOther top genes in this component:\n")
  print(top_genes_named[[strongest_component]])
} else {
  cat("CD151 not found in the dataset.\n")
}


enrich_results <- lapply(top_genes_named1, function(genes) {
  enrichGO(gene = genes,
           OrgDb = org.Hs.eg.db,
           keyType = "SYMBOL",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2)
  
})

# Print top 5 enriched terms for each component
for (i in seq_along(enrich_results)) {
  cat("Component", i, "top GO terms:\n")
  print(head(enrich_results[[i]]@result$Description, 5))
  cat("\n")
}

barplot(enrich_results[[8]], showCategory=5)

dev.off()

#######################






# ##### step 3
# 
# # Assuming log_norm_counts is your normalized, log-transformed expression data
# # Transpose the data so that genes are columns and samples are rows
# som_data <- t(log_norm_counts)
# 
# # Scale the data
# som_data_scaled <- scale(som_data)
# 
# 
# # Set up SOM grid (adjust dimensions as needed)
# som_grid <- somgrid(xdim = 5, ydim = 5, topo = "hexagonal")
# 
# # Train the SOM
# set.seed(123)  # for reproducibility
# som_model <- som(som_data_scaled, 
#                  grid = som_grid, 
#                  rlen = 100,  # number of iterations
#                  alpha = c(0.05, 0.01),  # learning rate
#                  keep.data = TRUE)
# 
# 
# 
# 
# # Plot the SOM grid
# plot(som_model, type = "changes")  # learning progress
# plot(som_model, type = "count")    # number of samples in each node
# plot(som_model, type = "quality")  # quality of the mapping
# 
# # Visualize the distribution of a specific gene across the SOM
# gene_of_interest <- "S100P"
# plot(som_model, type = "property", property = som_data_scaled[, gene_of_interest])
# 
# 
# 
# # Assuming som_model is your trained SOM model
# 
# # Extract the codebook vectors
# codebook <- som_model$codes
# 
# # If codebook is a list, we need to convert it to a matrix
# if(is.list(codebook)) {
#   codebook <- do.call(rbind, codebook)
# }
# 
# # Perform hierarchical clustering on the codebook vectors
# som_cluster <- cutree(hclust(dist(codebook)), 5)  # 5 clusters, adjust as needed
# 
# # Add cluster information to the original data
# som_results <- data.frame(
#   sample = rownames(som_data),
#   cluster = som_cluster[som_model$unit.classif]
# )
# 
# # Visualize the clusters
# plot(som_model, type = "mapping", bgcol = rainbow(5)[som_cluster], main = "SOM Clusters")
# 
# 
# 
# 
# # Function to perform GO enrichment for a gene list
# go_enrich <- function(genes) {
#   enrichGO(gene = genes,
#            OrgDb = org.Hs.eg.db,
#            keyType = "SYMBOL",
#            ont = "BP",
#            pAdjustMethod = "BH",
#            pvalueCutoff = 0.05,
#            qvalueCutoff = 0.2)
# }
# 
# # Perform enrichment for each cluster
# cluster_enrichment <- lapply(unique(som_results$cluster), function(c) {
#   cluster_genes <- rownames(log_norm_counts)[som_results$cluster == c]
#   go_enrich(cluster_genes)
# })
# 
# # Print top 5 enriched terms for each cluster
# for (i in seq_along(cluster_enrichment)) {
#   cat("Cluster", i, "top GO terms:\n")
#   print(head(cluster_enrichment[[i]]@result$Description, 5))
#   cat("\n")
# }













