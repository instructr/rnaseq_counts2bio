## ----load_packages, eval=TRUE, message=FALSE, warning=FALSE-----------------------------------------------

# # Install the required packages if not already installed
# install.packages(c("pak"))
# 
# pak::pkg_install(c("BiocManager", "remotes", "here", "tidyverse",        
# "DESeq2", "pheatmap", "RColorBrewer", "ggrepel", "clusterProfiler",
# "enrichplot", "org.Mm.eg.db", "patchwork", "ComplexHeatmap"
# ))
# 
# # Install the course data package
# pak::pak("patterninstitute/OSD758")


# Load packages
library("here")            # package to find your current working directory
library("tidyverse")       # packages for data manipulation and visualization
library("DESeq2")          # differential expression analysis
library("pheatmap")        # heatmaps
library("RColorBrewer")    # color palettes
library("ggrepel")         # repel overlapping text labels in ggplot2 plots
library("clusterProfiler") # for enrichment analysis
library("enrichplot")      # to draw functional enrichment results
library("org.Mm.eg.db")    # mouse gene annotation database
library("patchwork")         # combining multiple plots
library("ComplexHeatmap")  # to draw heatmaps

# Install and load package containing the data
library(OSD758)

# Gene expression in Counts
raw_counts <- OSD758::gene_expression(format = "wide", only_expressed_genes = TRUE) 
# View(raw_counts)

# Samples metadata
samples <- OSD758::samples()
# View(samples)



## ----vst, eval=TRUE, message=FALSE, warning=FALSE---------------------------------------------------------

# Create a list to save the QC results
qc <- list()

# You can choose between vst() and rlog() - this tutorial uses vst.
qc$vst <- DESeq2::vst(raw_counts, blind = TRUE)



## ----pca, eval=TRUE---------------------------------------------------------------------------------------

# Run PCA
qc$pca_vst <- prcomp(t(qc$vst)) 

# Extract the components
qc$components <- qc$pca_vst[["x"]]
qc$components <- tibble::as_tibble(qc$components, rownames = "sample_id")

# Add sample annotations to components for plot coloring
qc$components_annot <-
  dplyr::left_join(qc$components, as.data.frame(samples[, c(1,5,6,8)]), by = "sample_id") |>
  dplyr::relocate(spacecraft, acceleration_source, gravity_class, .after = sample_id)

# Calculate the % variance per component
qc$pca_percent_var <- round(qc$pca_vst$sdev^2/sum(qc$pca_vst$sdev^2)*100)

#
# 2D PCA | Using ggplot2
#

# Color by gravity_class
qc$pca_gravity <-
  ggplot(qc$components_annot, aes(x = PC1, y = PC2, color = gravity_class)) +
  geom_point(size = 3) +
  labs(
    title = "PCA gene expression | Colored by gravity_class",
    x = paste0("PC1 (", qc$pca_percent_var[1], "% variance)"),
    y = paste0("PC2 (", qc$pca_percent_var[2], "% variance)")
  ) +
  theme_minimal()

# Color by accelaration_source
qc$pca_acceleration <-
  ggplot(qc$components_annot, aes(x = PC1, y = PC2, color = acceleration_source)) +
  geom_point(size = 3) +
  labs(
    title = "PCA gene expression | Colored by acceleration_source",
    x = paste0("PC1 (", qc$pca_percent_var[1], "% variance)"),
    y = paste0("PC2 (", qc$pca_percent_var[2], "% variance)")
  ) +
  theme_minimal()

# Color by gravity_class and shape by acceleration_source
qc$pca_gravity_acceleration <-
  ggplot(qc$components_annot, aes(x = PC1, y = PC2, 
                               color = gravity_class,
                               shape = acceleration_source)) +
  geom_point(size = 3) +
  labs(
    title = "PCA gene expression | Colored by gravity_class | Shape acceleration_source",
    x = paste0("PC1 (", qc$pca_percent_var[1], "% variance)"),
    y = paste0("PC2 (", qc$pca_percent_var[2], "% variance)")
  ) +
  theme_minimal()


# Assemble pca plots
qc$pca_gravity_acceleration /
(qc$pca_gravity | qc$pca_acceleration)



## ----dist_clust, eval=TRUE--------------------------------------------------------------------------------

# Plot sample to sample distance for hierarchical clustering

# Calculate Euclidean distances between samples (rows) by transposing the matrix with t().
qc$sample_dist_matrix <- as.matrix(dist(t(qc$vst), method = "euclidean"))


# Define a color palette for the heatmap
qc$colors <- colorRampPalette(rev(brewer.pal(9, "Greens")))(255) # function from RColorBrewer package

# Create the heatmap
qc$dist_clustering <- pheatmap::pheatmap(
  qc$sample_dist_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  col = qc$colors,
  fontsize_col = 8,
  fontsize_row = 5
)



## ----corr_clustering, eval=TRUE---------------------------------------------------------------------------

### Compute pairwise correlation values
qc$sample_corr <- cor(qc$vst)

### Plot heatmap using the correlation matrix
qc$corr_clustering <-
  pheatmap::pheatmap(
    qc$sample_corr,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 5,
    fontsize_col = 8
  )



## ----check_ids, eval=TRUE---------------------------------------------------------------------------------

# Create list to save the analysis objects
de_deseq <- list()

# Check that sample ids match between raw_counts and samples 
# Ensure same content
stopifnot(setequal(colnames(raw_counts), samples$sample_id))

# Reorder columns to match sample order
raw_counts <- raw_counts[, samples$sample_id]



## ----de_analysis, eval=TRUE, warning=FALSE, message=FALSE-------------------------------------------------

# Make sure the factor levels are ordered so that the desired baseline comes first.
# DESeq2 uses the first level from factors as the baseline.
samples <-
  samples |>
  dplyr::mutate(gravity_class = factor(
    gravity_class,
    levels = c("1.00_G", "0.33_G", "0.66_G", "micro_G")
  ))


# DE Step 1: Create a DESeqDataSet object (dds)
de_deseq$dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = samples,
  design = ~ gravity_class
)

# DE Step 2: Run the DESeq function to perform the analysis
de_deseq$dds <- DESeq(de_deseq$dds)



## ----dds_view, eval=TRUE, output=FALSE--------------------------------------------------------------------

# Check the design formula
DESeq2::design(de_deseq$dds) 

# Check the sample info
SummarizedExperiment::colData(de_deseq$dds) 

# Display the first rows of the raw counts
head(DESeq2::counts(de_deseq$dds))

# Display the first rows of the normalised counts to compare with raw counts 
head(DESeq2::counts(de_deseq$dds, normalized = TRUE))

# Convert the normalised counts from the DESeq2 object to a tibble
normalised_counts <- tibble::as_tibble(DESeq2::counts(de_deseq$dds, normalized = TRUE),
                                       rownames = "ensembl_gen_id")
head(normalised_counts)



## ----de_res, eval=TRUE, output=FALSE----------------------------------------------------------------------

# Find the names of the estimated effects (coefficients) of the model
DESeq2::resultsNames(de_deseq$dds)

# Extract DE results for each gravity condition vs 1.00 G
    # The results function by default applies the Benjamini-Hochberg method to control FDR
de_deseq$res_033_vs_1G <- DESeq2::results(de_deseq$dds, name = "gravity_class_0.33_G_vs_1.00_G")
de_deseq$res_066_vs_1G <- DESeq2::results(de_deseq$dds, name = "gravity_class_0.66_G_vs_1.00_G")
de_deseq$res_micro_vs_1G <- DESeq2::results(de_deseq$dds, name = "gravity_class_micro_G_vs_1.00_G")


# Summarise the results:
  # Shows the number of tested genes, the number up- and down-regulated (at alpha),
  # and how many were excluded by multiple testing due to low counts.
DESeq2::summary(de_deseq$res_033_vs_1G)
DESeq2::summary(de_deseq$res_066_vs_1G)
DESeq2::summary(de_deseq$res_micro_vs_1G)



## ----sig_res, eval=TRUE, output=FALSE---------------------------------------------------------------------

# Extract significant results (padj < 0.05) and convert to tibble
de_deseq$sig_033_vs_1G <-
  de_deseq$res_033_vs_1G |>
  tibble::as_tibble(rownames = "ensembl_gen_id") |>
  dplyr::filter(!is.na(padj), padj < 0.05) |>
  dplyr::arrange(padj, log2FoldChange)

de_deseq$sig_066_vs_1G <-
  de_deseq$res_066_vs_1G |>
  tibble::as_tibble(rownames = "ensembl_gen_id") |>
  dplyr::filter(!is.na(padj), padj < 0.05) |>
  dplyr::arrange(padj, log2FoldChange)

de_deseq$sig_micro_vs_1G <-
  de_deseq$res_micro_vs_1G |>
  tibble::as_tibble(rownames = "ensembl_gen_id") |>
  dplyr::filter(!is.na(padj), padj < 0.05) |>
  dplyr::arrange(padj, log2FoldChange)

# Look at the top results
head(de_deseq$sig_033_vs_1G)
head(de_deseq$sig_066_vs_1G)
head(de_deseq$sig_micro_vs_1G)



## ----degs_heatmaps, eval=TRUE, fig.height=12, fig.width=12------------------------------------------------

# List to save all the visualization plots
de_plots <- list()

# Extract only gene ids from the significant results
sig_gene_ids <- de_deseq$sig_micro_vs_1G$ensembl_gen_id

# Map between ENSEMBL gene ids and gene symbol
ensembl2symbol <- OSD758::gene_expression("long") |>
  dplyr::select(ensembl_gen_id, gene_symbol) |>
  dplyr::distinct()
  
# Get normalised counts for significant genes 
sig_normalised_counts <- normalised_counts |>
  dplyr::filter(ensembl_gen_id %in% sig_gene_ids) |>
  dplyr::left_join(ensembl2symbol, by = "ensembl_gen_id") |>
  dplyr::select(-ensembl_gen_id) |>
  tibble::column_to_rownames("gene_symbol") |>
  as.matrix()


# Scale each row: subtract mean and divide by SD.
# The 2 transpositions are required because, by default, scale applies to the columns.
sig_normalised_counts_scaled <- t(scale(t(sig_normalised_counts)))   # scale rows, not columns

# Find min and max values to get meaningful colors in heatmaps
range(sig_normalised_counts_scaled)

# Complex heatmap
de_plots$ht <- ComplexHeatmap::Heatmap(sig_normalised_counts_scaled[1:200, ],
                        name = "Exprs (z-score)",
                        column_title = "Microgravity vs Earth's Gravity (1G) | Top 200 DE genes",
                        cluster_columns = TRUE,
                        cluster_rows = TRUE,
                        # number of clusters in K-means to split rows
                        row_km = 2,
                        # add cluster names
                        row_title = c("A", "B"),
                        row_title_rot = 90,
                        row_gap = unit(2, "mm"),
                        # number of clusters in K-means to split columns
                        column_km = 3,
                        column_gap = unit(2, "mm"),
                        border = "grey",
                        na_col = "white",
                        # Color range (min and max values from sig_normalised_counts_scaled)
                        col = circlize::colorRamp2(c(-4, 0, 6), c("skyblue3", "white", "forestgreen")),
                        column_names_gp = grid::gpar(fontsize = 9),
                        row_names_gp = grid::gpar(fontsize = 5),
                        rect_gp = grid::gpar(col = "grey", lwd = 0.5))

# Print the plot
ComplexHeatmap::draw(de_plots$ht, heatmap_legend_side = "right")



## ----volcano, eval=TRUE, warning=FALSE--------------------------------------------------------------------

# Add a column with differential expression status and add gene symbol to the results
sig_res_annot <- 
  de_deseq$sig_micro_vs_1G |>
  dplyr::mutate(diffexpressed = case_when(
    log2FoldChange > 1 & padj < 0.05 ~ 'upregulated',
    log2FoldChange < -1 & padj < 0.05 ~ 'downregulated',
    TRUE ~ 'not_de')) |>
  dplyr::left_join(ensembl2symbol, by = "ensembl_gen_id") |>
  dplyr::select(-ensembl_gen_id) |>
  # add gene symbols
  dplyr::relocate(gene_symbol, .before =  1L) |>
  dplyr::arrange(padj, log2FoldChange)


# Create a volcano plot using ggplot2
de_plots$volcano_plot <-
  ggplot(data = sig_res_annot, aes(
    x = log2FoldChange,
    y = -log10(padj),
    col = diffexpressed))+
  geom_point(size = 0.6) +
  geom_text_repel(data = filter(sig_res_annot, 
                                ((abs(log2FoldChange) > log2(8)) & (padj < -log10(0.05)))), 
                  aes(label = gene_symbol), size = 2.5, max.overlaps = Inf) +
  ggtitle("DE genes micro gravity versus Earth's gravity") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed', linewidth = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed', linewidth = 0.2) +
  theme(plot.title = element_text(size = rel(1.25), hjust = 0.5),
        axis.title = element_text(size = rel(1))) +
  scale_color_manual(values = c("upregulated" = "red",
                                "downregulated" = "blue",
                                "not_de" = "grey")) +
  labs(color = 'DE genes') +
  xlim(-5, 5) +   # Caution: This hides some genes
  ylim(0, 7.5) +  # Caution: This hides some genes
  theme_light()

# Print the volcano plot
de_plots$volcano_plot



## ----func_enrich, eval=TRUE, warning=FALSE, message=FALSE, fig.height=10, fig.width=12--------------------

# Enrichment analysis (ORA)

# Create a list to save the enrichment analysis results
fun_enrich <- list()

# Prepare list of significant DE genes in descending Log2FoldChange
fun_enrich$de_genes_fc <-
  de_deseq$sig_micro_vs_1G |>
  dplyr::select(ensembl_gen_id, log2FoldChange) |>
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Run GO enrichment analysis using the enrichGO function
fun_enrich$ego <- clusterProfiler::enrichGO(
  gene = fun_enrich$de_genes_fc$ensembl_gen_id, # Genes of interest
  universe = ensembl2symbol$ensembl_gen_id,     # Background gene set
  OrgDb = org.Mm.eg.db,                         # Annotation database
  keyType = 'ENSEMBL',                          # Key type for gene identifiers
  readable = TRUE,                              # Convert gene IDs to gene names
  ont = "BP",                                   # Ontology: can be "BP", "MF", "CC", or "ALL"
  pvalueCutoff = 0.05,                          # P-value cutoff for significance
  qvalueCutoff = 0.10                           # Q-value cutoff for significance
)


# Visualize the enriched GO terms
fun_enrich$dotplot <- 
  enrichplot::dotplot(fun_enrich$ego, showCategory = 20, title = "GO BP | Enrichment barplot")

fun_enrich$heatplot <- 
  enrichplot::heatplot(fun_enrich$ego, showCategory = 10, 
                       foldChange = fun_enrich$de_genes_fc$log2FoldChange) +
  ggplot2::ggtitle("GO BP | Enrichment heatplot")

fun_enrich$emapplot <- 
  enrichplot::emapplot(pairwise_termsim(fun_enrich$ego), showCategory = 15, layout = "nicely")

fun_enrich$cnetplot <- 
  enrichplot::cnetplot(fun_enrich$ego, categorySize = "pvalue", showCategory = 5, 
                                 layout = "nicely", foldChange = fun_enrich$de_genes_fc$log2FoldChange)

fun_enrich$treeplot <- 
  enrichplot::treeplot(enrichplot::pairwise_termsim(fun_enrich$ego), 
                       showCategory = 20, nCluster=5, offset = rel(2)) + 
  ggplot2::ggtitle("GO BP | Enrichment treeplot") + 
  ggplot2::theme(text = element_text(size = 8))


# Combine the enrichment plots into panels from a single figure
(fun_enrich$dotplot) |
(fun_enrich$emapplot / fun_enrich$cnetplot)

fun_enrich$treeplot / fun_enrich$heatplot



