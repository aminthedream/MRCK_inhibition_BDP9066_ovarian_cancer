library(STRINGdb)
library(VennDiagram)
library(biomaRt)
library(ReactomePA)
library(org.Hs.eg.db)
library(gridExtra)
library(clusterProfiler)
library(PharmacoGx)



gdsc <- downloadPSet("GDSC_2020(v2-8.2)")

bdp_sensitivity <- drugSensitivitySig(gdsc, 
                                      drug = "BDP-00009066", 
                                      sensitivity.measure = "ic50_recomputed",
                                      mDataType = "rna")

das_sensitivity <- drugSensitivitySig(gdsc, 
                                      drug = "Dasatinib", 
                                      sensitivity.measure = "ic50_recomputed",
                                      mDataType = "rna")


# Rank cell lines by sensitivity
sensitivity_estimates_bdp <- bdp_sensitivity[,,"estimate"]
sensitivity_estimates_das <- das_sensitivity[,,"estimate"]


sensitivity_df_bdp <- data.frame(
  gene = names(sensitivity_estimates_bdp),
  estimate = as.numeric(sensitivity_estimates_bdp),
  row.names = NULL
)


sensitivity_df_das <- data.frame(
  gene = names(sensitivity_estimates_das),
  estimate = as.numeric(sensitivity_estimates_das),
  row.names = NULL
)


ranked_genes_bdp <- sensitivity_df_bdp[order(sensitivity_df_bdp$estimate), ]
ranked_genes_das <- sensitivity_df_das[order(sensitivity_df_das$estimate), ]




# Assuming lower values indicate higher sensitivity
n_top_genes <- 1000  # You can adjust this number
most_sensitive_genes_bdp <- head(ranked_genes_bdp, n_top_genes)
most_sensitive_genes_das <- head(ranked_genes_das, n_top_genes)

# If you want to filter by a threshold instead
# threshold <- mean(ranked_genes$`BDP-00009066`) - 2*sd(ranked_genes$`BDP-00009066`)
# most_sensitive_genes <- ranked_genes[ranked_genes$`BDP-00009066` < threshold, ]


# Convert ENSEMBL IDs to Entrez IDs
entrez_ids_bdp <- bitr(most_sensitive_genes_bdp$gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids_das <- bitr(most_sensitive_genes_das$gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# GO enrichment analysis
go_enrichment_bdp <- enrichGO(gene = entrez_ids_bdp$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",  # Biological Process
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)


go_enrichment_das <- enrichGO(gene = entrez_ids_das$ENTREZID,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",  # Biological Process
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)

# KEGG pathway analysis
kegg_enrichment_bdp <- enrichKEGG(gene = entrez_ids_bdp$ENTREZID,
                              organism = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05)


kegg_enrichment_das <- enrichKEGG(gene = entrez_ids_das$ENTREZID,
                              organism = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05)


# Visualize results
barplot(go_enrichment_bdp, showCategory=5)
barplot(go_enrichment_bdp, showCategory=5)

barplot(kegg_enrichment_bdp, showCategory=5)
barplot(kegg_enrichment_das, showCategory=5)



string_db <- STRINGdb$new(version="11", species=9606)
string_mapped <- string_db$map(most_sensitive_genes, "gene", removeUnmappedRows = TRUE)
string_db$plot_network(string_mapped$STRING_id)



# Create sets for each drug
bdp_genes <- most_sensitive_genes_bdp$gene
das_genes <- most_sensitive_genes_das$gene

# Find common and unique genes
common_genes <- intersect(bdp_genes, das_genes)
only_bdp <- setdiff(bdp_genes, das_genes)
only_das <- setdiff(das_genes, bdp_genes)

# Create the Venn diagram
venn <- draw.pairwise.venn(
  area1 = length(bdp_genes),
  area2 = length(das_genes),
  cross.area = length(common_genes),
  category = c("BDP-00009066", "Dasatinib"),
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(330, 30),
  cat.dist = 0.07,
  cat.just = list(c(-1, -1), c(1, 1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed"
)

# Create a text grob with gene counts
text_grob <- textGrob(
  paste("Total genes:", length(union(bdp_genes, das_genes)), "\n",
        "BDP-00009066 only:", length(only_bdp), "\n",
        "Dasatinib only:", length(only_das), "\n",
        "Common genes:", length(common_genes)),
  x = 0.5, y = 0.5, just = "center", gp = gpar(fontsize = 12)
)

# Combine the Venn diagram and text
grid.arrange(gTree(children = venn), text_grob, heights = c(4, 1))

# Save the plot
ggsave("informative_drug_sensitivity_overlap.png", width = 10, height = 12, units = "in", dpi = 300)

##############ONLY COMMON GENES##############

# Convert gene names to Entrez IDs
entrez_ids <- bitr(common_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# GO Enrichment Analysis
go_enrichment_common <- enrichGO(gene = entrez_ids$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",  # Biological Process
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)

# KEGG Pathway Analysis
kegg_enrichment_common <- enrichKEGG(gene = entrez_ids$ENTREZID,
                              organism = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05)

# Visualize results
barplot(go_enrichment, showCategory=5)
barplot(kegg_enrichment, showCategory=5)

string_db <- STRINGdb$new(version="11", species=9606)
string_mapped <- string_db$map(data.frame(gene=common_genes), "gene", removeUnmappedRows = TRUE)
string_db$plot_network(string_mapped$STRING_id)



######CDC42#########
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
cdc42_ensembl <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                        filters="hgnc_symbol",
                        values="CD151",
                        mart=ensembl)

cdc42_id <- cdc42_ensembl$ensembl_gene_id[1]
print(paste("cdc42 ENSEMBL ID:", cdc42_id))

# Assuming ranked_genes is your dataframe with 'gene' column as ENSEMBL IDs
# and 'BDP-00009066' column as sensitivity scores

# Get cdc42 sensitivity
cdc42_sensitivity <- ranked_genes_bdp$estimate[ranked_genes_bdp$gene == cdc42_id]

# Calculate difference in sensitivity from cdc42
ranked_genes_bdp$diff_from_cdc42 <- abs(ranked_genes_bdp$estimate - cdc42_sensitivity)

# Get genes with similar sensitivity
similar_sensitivity_genes <- ranked_genes_bdp %>%
  arrange(diff_from_cdc42) %>%
  head(1000)  # Top 100 most similar genes

# Convert ENSEMBL IDs to gene symbols for easier interpretation
ensembl_ids <- similar_sensitivity_genes$gene
gene_symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                      filters="ensembl_gene_id",
                      values=ensembl_ids,
                      mart=ensembl)

# Merge gene symbols with sensitivity data
similar_sensitivity_genes <- merge(similar_sensitivity_genes, gene_symbols, 
                                   by.x="gene", by.y="ensembl_gene_id")

cdc42_entrez <- bitr(similar_sensitivity_genes$gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
cdc42_pathways <- enrichPathway(gene=cdc42_entrez$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)

barplot(cdc42_pathways, showCategory=8)

# Initialize STRING database
string_db <- STRINGdb$new(version="11", species=9606)

# Get cdc42 interactors
cdc42_mapped <- string_db$map(data.frame(gene=cdc42_id), "gene", removeUnmappedRows = TRUE)

cdc42_string_id <- cdc42_mapped$STRING_id
cdc42_interactors <- string_db$get_neighbors(cdc42_string_id)

print("cdc42 interactors:")
print(cdc42_interactors)



###############CDC42 PATHWAY###############

# Convert cdc42 to Entrez ID
cdc42_entrez <- bitr("CD151", fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

# Get pathways involving cdc42
cdc42_pathways <- enrichPathway(gene=cdc42_entrez, pvalueCutoff = 0.05, readable=TRUE)

# Print pathways
print(cdc42_pathways)
barplot(cdc42_pathways, showCategory=10)

# Get genes involved in these pathways
related_genes <- unique(unlist(cdc42_pathways$geneID))

#############################


mical1_string_id <- cdc42_mapped$STRING_id
mical1_interactors <- string_db$get_neighbors(mical1_string_id)

print("MICAL1 interactors:")
print(mical1_interactors)


if(!is.null(mical1_interactors) && nrow(mical1_interactors) > 0) {
  # Convert STRING IDs to ENSEMBL IDs
  interactor_mapping <- string_db$map(mical1_interactors, "STRING_id", removeUnmappedRows = TRUE)
  
  # Get sensitivity data for interactors
  interactor_sensitivity <- ranked_genes %>%
    filter(gene %in% interactor_mapping$ensembl_gene_id)
  
  # Merge with gene symbols for easier interpretation
  gene_symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                        filters="ensembl_gene_id",
                        values=interactor_sensitivity$gene,
                        mart=ensembl)
  
  interactor_sensitivity <- merge(interactor_sensitivity, gene_symbols, 
                                  by.x="gene", by.y="ensembl_gene_id")
  
  # Display results
  print("Sensitivity of MICAL1 interactors:")
  print(interactor_sensitivity[order(interactor_sensitivity$`BDP-00009066`), c("hgnc_symbol", "BDP-00009066")])
  
  # Compare with common genes
  common_interactors <- intersect(interactor_sensitivity$gene, common_genes)
  
  if(length(common_interactors) > 0) {
    print("\nMICAL1 interactors that are also in the common sensitive gene set:")
    print(interactor_sensitivity[interactor_sensitivity$gene %in% common_interactors, c("hgnc_symbol", "BDP-00009066")])
  } else {
    print("\nNo MICAL1 interactors found in the common sensitive gene set")
  }
  
  # Analyze the interaction network
  print("\nNumber of MICAL1 interactors:", nrow(mical1_interactors))
  
  # Plot the interaction network
  png("MICAL1_interaction_network.png", width=800, height=800)
  string_db$plot_network(mical1_interactors$STRING_id)
  dev.off()
  print("Interaction network plot saved as 'MICAL1_interaction_network.png'")
  
  # Functional enrichment of interactors
  enrichment <- string_db$get_enrichment(mical1_interactors$STRING_id)
  
  if(!is.null(enrichment) && nrow(enrichment) > 0) {
    print("\nTop 10 enriched pathways among MICAL1 interactors:")
    print(head(enrichment[order(enrichment$fdr),], 10))
  } else {
    print("\nNo significant pathway enrichment found for MICAL1 interactors")
  }
} else {
  print("No interactors found for MICAL1")
}
#################


cd151_component <- which.max(metagenes1[which(rownames(metagenes1) == "CD151"),])
cd151_coexpressed_genes <- top_genes_named1[[8]]

# Convert gene symbols to ENSEMBL IDs
#ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
d151_genes_ensembl <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                            filters="hgnc_symbol",
                            values=cd151_coexpressed_genes,
                            mart=ensembl)

# Get sensitivity data for these genes
d151_sensitivity <- ranked_genes_bdp %>%
  filter(gene %in% d151_genes_ensembl$ensembl_gene_id)

# Merge with gene symbols
d151_sensitivity <- merge(d151_sensitivity, d151_genes_ensembl, 
                          by.x="gene", by.y="ensembl_gene_id")

# Display results
print("Sensitivity of d151 coexpressed genes to BDP-00009066:")
print(d151_sensitivity[order(d151_sensitivity$estimate), c("hgnc_symbol", "estimate")])




# Calculate mean and standard deviation of overall sensitivity
overall_mean <- mean(ranked_genes_bdp$estimate, na.rm=TRUE)
overall_sd <- sd(ranked_genes_bdp$estimate, na.rm=TRUE)

# Calculate z-scores for d151 genes
d151_sensitivity$z_score <- (d151_sensitivity$estimate - overall_mean) / overall_sd

# Display results with z-scores
print("Sensitivity of d151 coexpressed genes to BDP-00009066 (with z-scores):")
print(d151_sensitivity[order(d151_sensitivity$z_score), c("hgnc_symbol", "estimate", "z_score")])

# Identify significantly sensitive genes (e.g., z-score < -2)
significant_genes <- d151_sensitivity[d151_sensitivity$z_score < -2, ]
if(nrow(significant_genes) > 0) {
  print("Genes significantly more sensitive to BDP-00009066:")
  print(significant_genes[, c("hgnc_symbol", "estimate", "z_score")])
} else {
  print("No genes from the list show significantly higher sensitivity to BDP-00009066.")
}




# Get top 5% sensitive genes
top_5_percent <- quantile(ranked_genes_bdp$estimate, 0.05, na.rm=TRUE)
top_sensitive_genes <- ranked_genes_bdp[ranked_genes_bdp$est <= top_5_percent, ]

# Check overlap with d151 genes
overlap_top_sensitive <- intersect(d151_sensitivity$gene, top_sensitive_genes$gene)

if(length(overlap_top_sensitive) > 0) {
  print("d151 coexpressed genes among top 5% most sensitive genes:")
  print(d151_sensitivity[d151_sensitivity$gene %in% overlap_top_sensitive, c("hgnc_symbol", "estimate")])
} else {
  print("No d151 coexpressed genes are among the top 5% most sensitive genes.")
}





# Check overlap with d151 genes
overlap_top_sensitive <- intersect(d151_sensitivity$gene, top_sensitive_genes$gene)

if(length(overlap_top_sensitive) > 0) {
  print("d151 coexpressed genes among top 5% most sensitive genes:")
  print(d151_sensitivity[d151_sensitivity$gene %in% overlap_top_sensitive, c("hgnc_symbol", "estimate")])
} else {
  print("No d151 coexpressed genes are among the top 5% most sensitive genes.")
}









