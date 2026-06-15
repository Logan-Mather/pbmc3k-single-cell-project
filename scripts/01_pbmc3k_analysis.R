
# 01_pbmc3k_analysis.R

# PBMC3K single-cell RNA-seq analysis using Seurat


# -----------------------------------------
# 0. Load packages and set up project
# -----------------------------------------

library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(mclust)

# Set seed for reproducibility
seed_value <- 6142026

# Create output folders if they do not already exist
dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)


# -----------------------------------------
# 1. Load filtered 10x Genomics matrix
# -----------------------------------------

pbmc.data <- Read10X(
  data.dir = "data/filtered_gene_bc_matrices/hg19/"
)

# Check raw matrix
class(pbmc.data)
dim(pbmc.data)
pbmc.data[1:5, 1:5]


# -----------------------------------------
# 2. Explore raw matrix before creating Seurat object
# -----------------------------------------

# Number of genes detected in each cell
cell_features <- Matrix::colSums(pbmc.data > 0)

# Total RNA counts in each cell
cell_counts <- Matrix::colSums(pbmc.data)

# Number of cells where each gene is detected
gene_cells <- Matrix::rowSums(pbmc.data > 0)

# Summaries
summary(cell_features)
summary(cell_counts)
summary(gene_cells)

# Check how many cells/genes are affected by possible thresholds
sum(cell_features < 200)
sum(gene_cells < 5)


# -----------------------------------------
# 3. Gene-level filter sensitivity check
# -----------------------------------------

min_cells_thresholds <- 1:50

gene_filter_curve <- data.frame(
  min_cells = min_cells_thresholds,
  genes_removed = sapply(min_cells_thresholds, function(x) sum(gene_cells < x)),
  genes_remaining = sapply(min_cells_thresholds, function(x) sum(gene_cells >= x))
)

gene_filter_curve$percent_remaining <- round(
  100 * gene_filter_curve$genes_remaining / length(gene_cells),
  1
)

head(gene_filter_curve, 10)

gene_filter_plot <- ggplot(
  gene_filter_curve,
  aes(x = min_cells, y = percent_remaining)
) +
  geom_line() +
  geom_point(size = 1) +
  scale_x_continuous(breaks = seq(1, 50, by = 5)) +
  labs(
    title = "Gene Retention Across min.cells Thresholds",
    subtitle = "Percent of genes retained as the detection threshold increases",
    x = "Minimum number of cells detecting a gene",
    y = "Genes retained (%)"
  )

gene_filter_plot

ggsave(
  filename = "figures/gene_filter_sensitivity.png",
  plot = gene_filter_plot,
  width = 8,
  height = 5,
  dpi = 300
)

genes_lost_each_step <- data.frame(
  detected_in_cells = 1:50,
  genes_at_detection_level = sapply(1:50, function(x) sum(gene_cells == x))
)

head(genes_lost_each_step, 20)

genes_lost_plot <- ggplot(
  genes_lost_each_step,
  aes(x = detected_in_cells, y = genes_at_detection_level)
) +
  geom_col() +
  coord_cartesian(xlim = c(1, 25)) +
  scale_x_continuous(breaks = 1:25) +
  labs(
    title = "Genes by Number of Detecting Cells",
    subtitle = "Shows how many genes are detected in exactly 1, 2, 3, ... cells",
    x = "Number of cells detecting gene",
    y = "Number of genes"
  )

genes_lost_plot

ggsave(
  filename = "figures/genes_by_detection_level.png",
  plot = genes_lost_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# -----------------------------------------
# 4. Check canonical marker genes before filtering
# -----------------------------------------

marker_genes <- c(
  "IL7R", "CCR7",      # CD4 T cells
  "CD3D", "CD3E",      # T cells
  "CD8A",              # CD8 T cells
  "MS4A1", "CD79A",    # B cells
  "CD14", "LYZ",       # CD14+ monocytes
  "FCGR3A", "MS4A7",   # CD16+ monocytes
  "GNLY", "NKG7",      # NK cells
  "FCER1A", "CST3",    # Dendritic cells
  "PPBP"               # Platelets
)

marker_check <- data.frame(
  gene = marker_genes,
  cells_detected = sapply(marker_genes, function(gene) {
    if (gene %in% names(gene_cells)) {
      as.numeric(gene_cells[gene])
    } else {
      NA
    }
  })
)

marker_check$retained_min_cells_5 <- marker_check$cells_detected >= 5
marker_check$retained_min_cells_10 <- marker_check$cells_detected >= 10
marker_check$retained_min_cells_20 <- marker_check$cells_detected >= 20

marker_check

write.csv(
  marker_check,
  file = "outputs/canonical_marker_filter_check.csv",
  row.names = FALSE
)


# -----------------------------------------
# 5. Create Seurat object
# -----------------------------------------

# Use min.cells = 5 based on the gene-retention sensitivity check.
# This removes genes detected in only 0-4 cells while avoiding an overly aggressive
# global filter that could remove genes useful for rare-cell-type annotation.

pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc3k",
  min.cells = 5,
  min.features = 200
)

pbmc


# -----------------------------------------
# 6. Inspect cell-level metadata
# -----------------------------------------

head(pbmc@meta.data)

summary(pbmc$nFeature_RNA)
summary(pbmc$nCount_RNA)

# nFeature_RNA = number of genes detected in each cell
# nCount_RNA   = total RNA counts detected in each cell


# -----------------------------------------
# 7. Calculate mitochondrial percentage
# -----------------------------------------

head(grep("^MT-", rownames(pbmc), value = TRUE))

pbmc[["percent.mt"]] <- PercentageFeatureSet(
  pbmc,
  pattern = "^MT-"
)

head(pbmc@meta.data)

summary(pbmc$percent.mt)


# -----------------------------------------
# 8. Evaluate QC filtering thresholds
# -----------------------------------------

qc_filter <- pbmc$nFeature_RNA > 200 &
  pbmc$nFeature_RNA < 2500 &
  pbmc$percent.mt < 5

table(qc_filter)

qc_summary <- data.frame(
  total_cells = ncol(pbmc),
  cells_removed = sum(!qc_filter),
  cells_remaining = sum(qc_filter),
  percent_remaining = round(100 * sum(qc_filter) / ncol(pbmc), 1),
  low_feature_cells = sum(pbmc$nFeature_RNA <= 200),
  high_feature_cells = sum(pbmc$nFeature_RNA >= 2500),
  high_mito_cells = sum(pbmc$percent.mt >= 5)
)

qc_summary

write.csv(
  qc_summary,
  file = "outputs/qc_filter_summary.csv",
  row.names = FALSE
)


# -----------------------------------------
# 9. Apply QC filtering
# -----------------------------------------

pbmc_unfiltered <- pbmc

pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    percent.mt < 5
)

pbmc

summary(pbmc$nFeature_RNA)
summary(pbmc$nCount_RNA)
summary(pbmc$percent.mt)


# -----------------------------------------
# 10. QC visualization after filtering
# -----------------------------------------

qc_violin_plot <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

qc_violin_plot

ggsave(
  filename = "figures/qc_violin_after_filtering.png",
  plot = qc_violin_plot,
  width = 10,
  height = 5,
  dpi = 300
)


# -----------------------------------------
# 11. Normalize gene expression values
# -----------------------------------------

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

pbmc


# -----------------------------------------
# 12. Identify highly variable genes
# -----------------------------------------

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

variable_genes <- VariableFeatures(pbmc)

length(variable_genes)
head(variable_genes, 20)

top10_variable_genes <- head(variable_genes, 10)

top10_variable_genes


# -----------------------------------------
# 13. Plot highly variable genes
# -----------------------------------------

variable_feature_plot <- VariableFeaturePlot(pbmc)

labeled_variable_feature_plot <- LabelPoints(
  plot = variable_feature_plot,
  points = top10_variable_genes,
  repel = TRUE
)

labeled_variable_feature_plot

ggsave(
  filename = "figures/variable_features.png",
  plot = labeled_variable_feature_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# -----------------------------------------
# 14. Check relationship between detection frequency and variable feature selection
# -----------------------------------------

counts_current <- LayerData(
  pbmc,
  assay = "RNA",
  layer = "counts"
)

gene_cells_current <- Matrix::rowSums(counts_current > 0)

variable_feature_check <- data.frame(
  gene = rownames(pbmc),
  cells_detected = as.numeric(gene_cells_current[rownames(pbmc)]),
  is_variable_feature = rownames(pbmc) %in% variable_genes
)

# Sanity check: should equal 2000
sum(variable_feature_check$is_variable_feature)

bin_width <- 50

max_detected <- max(variable_feature_check$cells_detected, na.rm = TRUE)

breaks_equal <- seq(
  from = 0,
  to = ceiling(max_detected / bin_width) * bin_width,
  by = bin_width
)

variable_feature_check$detection_bin <- cut(
  variable_feature_check$cells_detected,
  breaks = breaks_equal,
  include.lowest = TRUE,
  right = FALSE
)

genes_per_bin <- aggregate(
  gene ~ detection_bin,
  data = variable_feature_check,
  FUN = length
)

names(genes_per_bin)[2] <- "genes_in_bin"

variable_features_per_bin <- aggregate(
  is_variable_feature ~ detection_bin,
  data = variable_feature_check,
  FUN = sum
)

names(variable_features_per_bin)[2] <- "variable_features"

variable_feature_equal_bins <- merge(
  genes_per_bin,
  variable_features_per_bin,
  by = "detection_bin"
)

variable_feature_equal_bins$percent_variable <- round(
  100 * variable_feature_equal_bins$variable_features /
    variable_feature_equal_bins$genes_in_bin,
  1
)

variable_feature_equal_bins

# Sanity check: should equal 2000
sum(variable_feature_equal_bins$variable_features)

write.csv(
  variable_feature_equal_bins,
  file = "outputs/variable_feature_detection_frequency_summary.csv",
  row.names = FALSE
)

variable_feature_detection_plot <- ggplot(
  variable_feature_equal_bins,
  aes(x = detection_bin, y = percent_variable)
) +
  geom_col() +
  labs(
    title = "Variable Feature Selection by Gene Detection Frequency",
    subtitle = paste0("Genes grouped into equal-width bins of ", bin_width, " cells"),
    x = "Number of cells detecting gene",
    y = "Genes selected as variable features (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

variable_feature_detection_plot

ggsave(
  filename = "figures/variable_feature_detection_frequency.png",
  plot = variable_feature_detection_plot,
  width = 10,
  height = 6,
  dpi = 300
)


# -----------------------------------------
# 15. Save normalized object
# -----------------------------------------

saveRDS(
  pbmc,
  file = "outputs/pbmc_normalized_variable_features.rds"
)


# -----------------------------------------
# 16. Scale data before PCA
# -----------------------------------------

# Scale the highly variable genes so they are comparable for PCA.
# This centers each gene around mean 0 and scales by its standard deviation.

pbmc <- ScaleData(
  pbmc,
  features = variable_genes
)

pbmc


# -----------------------------------------
# 17. Run principal component analysis
# -----------------------------------------

pbmc <- RunPCA(
  pbmc,
  features = variable_genes
)

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


# -----------------------------------------
# 18. PCA diagnostics
# -----------------------------------------

# This plot is used only as a diagnostic check.
# Each point is a cell projected onto the first two principal components.
# The main project visualization will be UMAP after clustering.

pca_plot <- DimPlot(
  pbmc,
  reduction = "pca"
) +
  labs(
    title = "PCA Projection of PBMC3K Cells",
    subtitle = "Each point is a cell projected onto the first two principal components",
    x = "PC1",
    y = "PC2"
  )

pca_plot

ggsave(
  filename = "figures/pca_cells_diagnostic.png",
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# Elbow plot to help decide how many PCs to use downstream
elbow_plot <- ElbowPlot(
  pbmc,
  ndims = 30
) +
  labs(
    title = "PCA Elbow Plot",
    subtitle = "Used to assess how many principal components to carry forward",
    x = "Principal component",
    y = "Standard deviation"
  )

elbow_plot

ggsave(
  filename = "figures/pca_elbow_plot.png",
  plot = elbow_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# Numeric variance summary for PCs
pca_stdev <- pbmc[["pca"]]@stdev
pca_variance <- pca_stdev^2
pca_percent_variance <- 100 * pca_variance / sum(pca_variance)
pca_cumulative_variance <- cumsum(pca_percent_variance)

pca_variance_df <- data.frame(
  PC = 1:length(pca_percent_variance),
  percent_variance = pca_percent_variance,
  cumulative_variance = pca_cumulative_variance
)

head(pca_variance_df, 20)

write.csv(
  pca_variance_df,
  file = "outputs/pca_variance_explained.csv",
  row.names = FALSE
)


# -----------------------------------------
# 19. PC sensitivity check: 6 through 13 PCs
# -----------------------------------------

pc_choices <- 6:13

pc_sensitivity_list <- lapply(pc_choices, function(n_pcs) {
  
  temp_pbmc <- pbmc
  
  temp_pbmc <- FindNeighbors(
    temp_pbmc,
    dims = 1:n_pcs,
    verbose = FALSE
  )
  
  temp_pbmc <- FindClusters(
    temp_pbmc,
    resolution = 0.5,
    verbose = FALSE,
    random.seed = seed_value
  )
  
  cluster_table <- table(Idents(temp_pbmc))
  
  summary_row <- data.frame(
    n_pcs = n_pcs,
    n_clusters = length(cluster_table),
    smallest_cluster = min(cluster_table),
    largest_cluster = max(cluster_table),
    cluster_sizes = paste(as.numeric(cluster_table), collapse = ", ")
  )
  
  cluster_size_rows <- data.frame(
    n_pcs = n_pcs,
    cluster = names(cluster_table),
    cells = as.numeric(cluster_table)
  )
  
  list(
    summary = summary_row,
    cluster_sizes = cluster_size_rows
  )
})

pc_sensitivity_summary <- do.call(
  rbind,
  lapply(pc_sensitivity_list, function(x) x$summary)
)

pc_sensitivity_cluster_sizes <- do.call(
  rbind,
  lapply(pc_sensitivity_list, function(x) x$cluster_sizes)
)

pc_sensitivity_summary
pc_sensitivity_cluster_sizes

write.csv(
  pc_sensitivity_summary,
  file = "outputs/pc_sensitivity_summary.csv",
  row.names = FALSE
)

write.csv(
  pc_sensitivity_cluster_sizes,
  file = "outputs/pc_sensitivity_cluster_sizes.csv",
  row.names = FALSE
)

pc_sensitivity_plot <- ggplot(
  pc_sensitivity_summary,
  aes(x = n_pcs, y = n_clusters)
) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 6:13) +
  scale_y_continuous(
    breaks = seq(0, max(pc_sensitivity_summary$n_clusters), by = 1)
  ) +
  labs(
    title = "Cluster Stability Across Number of PCs",
    subtitle = "Sensitivity check using 6 through 13 principal components",
    x = "Number of principal components used",
    y = "Number of clusters"
  )

pc_sensitivity_plot

ggsave(
  filename = "figures/pc_sensitivity_clusters.png",
  plot = pc_sensitivity_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# -----------------------------------------
# 20. Compare cluster assignments across PC choices
# -----------------------------------------

cluster_assignments <- list()

for (n_pcs in pc_choices) {
  
  temp_pbmc <- pbmc
  
  temp_pbmc <- FindNeighbors(
    temp_pbmc,
    dims = 1:n_pcs,
    verbose = FALSE
  )
  
  temp_pbmc <- FindClusters(
    temp_pbmc,
    resolution = 0.5,
    verbose = FALSE,
    random.seed = seed_value
  )
  
  cluster_assignments[[as.character(n_pcs)]] <- as.character(Idents(temp_pbmc))
}

pc_cluster_summary <- data.frame(
  n_pcs = pc_choices,
  n_clusters = sapply(cluster_assignments, function(x) length(unique(x)))
)

pc_cluster_summary

write.csv(
  pc_cluster_summary,
  file = "outputs/pc_cluster_summary.csv",
  row.names = FALSE
)

ari_vs_8 <- data.frame(
  n_pcs = pc_choices,
  n_clusters = sapply(cluster_assignments, function(x) length(unique(x))),
  ari_compared_to_8_pcs = sapply(pc_choices, function(n_pcs) {
    mclust::adjustedRandIndex(
      cluster_assignments[["8"]],
      cluster_assignments[[as.character(n_pcs)]]
    )
  })
)

ari_vs_8

write.csv(
  ari_vs_8,
  file = "outputs/ari_vs_8_pcs.csv",
  row.names = FALSE
)

pc8_vs_pc10_table <- table(
  pcs_8 = cluster_assignments[["8"]],
  pcs_10 = cluster_assignments[["10"]]
)

pc8_vs_pc10_table

write.csv(
  as.data.frame.matrix(pc8_vs_pc10_table),
  file = "outputs/pc8_vs_pc10_cluster_crosstab.csv"
)


# -----------------------------------------
# 21. Select PCs for downstream analysis
# -----------------------------------------

# The elbow plot suggested that most major structure was captured within the
# first several principal components. A sensitivity check using 6-13 PCs showed
# that 9 clusters were recovered beginning at 8 PCs. Comparing 8 and 10 PCs gave
# a high Adjusted Rand Index, and the cross-tabulation showed that rare/small
# clusters were stable while differences were mostly redistribution among larger
# clusters. Therefore, use the first 8 PCs as a parsimonious choice for neighbor
# graph construction, clustering, and UMAP.

selected_pcs <- 1:8


# -----------------------------------------
# 22. Build nearest-neighbor graph
# -----------------------------------------

pbmc <- FindNeighbors(
  pbmc,
  dims = selected_pcs
)

pbmc


# -----------------------------------------
# 23. Cluster cells
# -----------------------------------------

pbmc <- FindClusters(
  pbmc,
  resolution = 0.5,
  random.seed = seed_value
)

cluster_counts <- table(Idents(pbmc))

cluster_counts

write.csv(
  as.data.frame(cluster_counts),
  file = "outputs/cluster_counts.csv",
  row.names = FALSE
)


# -----------------------------------------
# 24. Run 2D UMAP
# -----------------------------------------

pbmc <- RunUMAP(
  pbmc,
  dims = selected_pcs,
  seed.use = seed_value
)


# -----------------------------------------
# 25. Plot 2D UMAP by Seurat cluster
# -----------------------------------------

umap_cluster_plot <- DimPlot(
  pbmc,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  labs(
    title = "PBMC3K UMAP by Seurat Cluster",
    subtitle = "Each point is a cell; clusters are based on the first 8 principal components",
    x = "UMAP 1",
    y = "UMAP 2"
  )

umap_cluster_plot

ggsave(
  filename = "figures/pbmc3k_umap_clusters.png",
  plot = umap_cluster_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# -----------------------------------------
# 26. Marker sensitivity check: min.pct values
# -----------------------------------------

min_pct_choices <- c(0.10, 0.25)

marker_sensitivity <- lapply(min_pct_choices, function(pct_cutoff) {
  
  markers <- FindAllMarkers(
    pbmc,
    only.pos = TRUE,
    min.pct = pct_cutoff,
    logfc.threshold = 0.25
  )
  
  markers$min_pct_used <- pct_cutoff
  
  markers
})

marker_sensitivity_df <- do.call(rbind, marker_sensitivity)

marker_sensitivity_summary <- marker_sensitivity_df %>%
  group_by(min_pct_used, cluster) %>%
  summarise(
    n_markers = n(),
    .groups = "drop"
  )

marker_sensitivity_summary

write.csv(
  marker_sensitivity_summary,
  file = "outputs/marker_sensitivity_min_pct_summary.csv",
  row.names = FALSE
)


# -----------------------------------------
# 27. Find final marker genes for each cluster
# -----------------------------------------

# Marker genes were identified using only positive markers, min.pct = 0.25,
# and logfc.threshold = 0.25. The min.pct threshold focuses on genes detected
# in a meaningful fraction of cells within a cluster, while the logFC threshold
# removes very small expression differences that are less useful for cell-type
# annotation.

pbmc_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  pbmc_markers,
  file = "outputs/pbmc_cluster_markers.csv",
  row.names = FALSE
)


# -----------------------------------------
# 28. Top marker genes per cluster
# -----------------------------------------

top_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  arrange(cluster, desc(avg_log2FC))

print(top_markers, n = 90)

write.csv(
  top_markers,
  file = "outputs/top_10_markers_per_cluster.csv",
  row.names = FALSE
)


# -----------------------------------------
# 29. Cluster annotation evidence table
# -----------------------------------------

annotation_evidence <- data.frame(
  cluster = as.character(0:8),
  cell_type = c(
    "Naive CD4 T cells",
    "CD14+ Monocytes",
    "Memory CD4 T cells",
    "B cells",
    "CD8 T cells",
    "CD16+ Monocytes",
    "NK cells",
    "Dendritic cells",
    "Platelets"
  ),
  supporting_markers = c(
    "CCR7, LEF1, TCF7, MAL, LDHB",
    "S100A8, S100A9, CD14, LYZ, MS4A6A",
    "IL7R, IL32, CD2, CD40LG, TRAT1",
    "CD79A, MS4A1, TCL1A, FCER2, HLA-DOB",
    "CD8A, GZMK, CCL5, NKG7, CST7",
    "FCGR3A, MS4A7, CDKN1C, LILRA3",
    "GNLY, GZMB, FGFBP2, XCL2, SPON2",
    "FCER1A, CD1C, CLEC10A, SERPINF1",
    "PPBP, PF4, GP9, ITGA2B"
  ),
  annotation_reasoning = c(
    "CCR7, LEF1, and TCF7 support a naive or central-memory CD4 T-cell identity.",
    "CD14, LYZ, S100A8, and S100A9 support a classical CD14+ monocyte identity.",
    "IL7R, IL32, CD2, and CD40LG support a CD4 T-cell identity distinct from the naive CD4 T-cell cluster.",
    "MS4A1 and CD79A are canonical B-cell markers.",
    "CD8A with cytotoxic-associated genes such as GZMK, CCL5, CST7, and NKG7 supports a CD8 T-cell identity.",
    "FCGR3A and MS4A7 support a CD16+ non-classical monocyte identity.",
    "GNLY, GZMB, FGFBP2, XCL2, and SPON2 support an NK-cell identity.",
    "FCER1A, CD1C, and CLEC10A support a dendritic-cell identity.",
    "PPBP, PF4, GP9, and ITGA2B support a platelet identity."
  ),
  stringsAsFactors = FALSE
)

annotation_evidence

write.csv(
  annotation_evidence,
  file = "outputs/cluster_annotation_evidence.csv",
  row.names = FALSE
)


# -----------------------------------------
# 30. Annotate clusters using evidence table
# -----------------------------------------

cluster_annotations <- setNames(
  annotation_evidence$cell_type,
  annotation_evidence$cluster
)

pbmc$cell_type <- unname(cluster_annotations[as.character(Idents(pbmc))])

# Check that every cell received an annotation
table(pbmc$cell_type, useNA = "ifany")

if (any(is.na(pbmc$cell_type))) {
  stop("Some cells were not assigned a cell type. Check cluster IDs in annotation_evidence.")
}

# Confirm metadata column exists
colnames(pbmc@meta.data)

cell_type_counts <- as.data.frame(
  table(cell_type = pbmc$cell_type)
)

cell_type_counts

write.csv(
  cell_type_counts,
  file = "outputs/cell_type_counts.csv",
  row.names = FALSE
)


# -----------------------------------------
# 31. Plot 2D UMAP by annotated cell type
# -----------------------------------------

umap_cell_type_plot <- DimPlot(
  pbmc,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  labs(
    title = "PBMC3K UMAP by Annotated Cell Type",
    subtitle = "Each point is a cell; nearby cells have similar gene-expression profiles",
    x = "UMAP 1",
    y = "UMAP 2"
  )

umap_cell_type_plot

ggsave(
  filename = "figures/pbmc3k_umap_cell_types.png",
  plot = umap_cell_type_plot,
  width = 9,
  height = 6,
  dpi = 300
)


# -----------------------------------------
# 32. Canonical marker dot plot for annotation support
# -----------------------------------------

canonical_markers <- c(
  "CCR7", "LEF1", "TCF7", "IL7R", "CD2", "CD40LG",
  "CD14", "LYZ", "S100A8", "S100A9",
  "MS4A1", "CD79A",
  "CD8A", "GZMK", "CCL5",
  "FCGR3A", "MS4A7",
  "GNLY", "GZMB", "FGFBP2",
  "FCER1A", "CD1C", "CLEC10A",
  "PPBP", "PF4", "GP9", "ITGA2B"
)

canonical_markers_present <- canonical_markers[
  canonical_markers %in% rownames(pbmc)
]

canonical_marker_dotplot <- DotPlot(
  pbmc,
  features = canonical_markers_present
) +
  RotatedAxis() +
  labs(
    title = "Canonical PBMC Marker Expression by Cluster",
    subtitle = "Marker expression used to support manual cluster annotation",
    x = "Marker genes",
    y = "Seurat cluster"
  )

canonical_marker_dotplot

ggsave(
  filename = "figures/canonical_marker_dotplot.png",
  plot = canonical_marker_dotplot,
  width = 12,
  height = 6,
  dpi = 300
)


# -----------------------------------------
# 33. Optional 3D UMAP for Shiny app exploration
# -----------------------------------------

# The 2D UMAP is used for the main static project figures.
# A 3D UMAP is also generated so the Shiny app can provide an
# interactive rotating visualization of cell-type separation.
# UMAP 1, UMAP 2, and UMAP 3 are visualization coordinates, not direct
# biological measurements.

pbmc <- RunUMAP(
  pbmc,
  dims = selected_pcs,
  n.components = 3,
  reduction.name = "umap_3d",
  reduction.key = "UMAP3D_",
  seed.use = seed_value
)

umap_3d_embeddings <- Embeddings(
  pbmc,
  reduction = "umap_3d"
)

umap_3d_df <- data.frame(
  cell = rownames(umap_3d_embeddings),
  UMAP_1 = umap_3d_embeddings[, 1],
  UMAP_2 = umap_3d_embeddings[, 2],
  UMAP_3 = umap_3d_embeddings[, 3],
  seurat_cluster = Idents(pbmc),
  cell_type = pbmc$cell_type
)

head(umap_3d_df)

write.csv(
  umap_3d_df,
  file = "outputs/umap_3d_coordinates.csv",
  row.names = FALSE
)


# -----------------------------------------
# 34. Save annotated Seurat object
# -----------------------------------------

saveRDS(
  pbmc,
  file = "outputs/pbmc_annotated.rds"
)


# -----------------------------------------
# 35. Save session information
# -----------------------------------------

sink("outputs/session_info.txt")
print(sessionInfo())
sink()

sessionInfo()
