# PBMC3K Single-Cell RNA-seq Analysis in R

This project analyzes the 10x Genomics PBMC3K single-cell RNA-seq dataset using R and Seurat. The goal is to build a reproducible single-cell analysis workflow that performs quality control, normalization, dimensionality reduction, clustering, marker gene detection, cell type annotation, and interactive visualization through a Shiny app.

The final Shiny app allows users to explore annotated PBMC cell types, view 2D and 3D UMAP visualizations, examine gene expression patterns, review top marker genes, and understand key single-cell terminology through a built-in glossary.

## Project Summary

This analysis identifies major peripheral blood mononuclear cell populations from the PBMC3K dataset, including:

| Annotated cell type | Number of cells |
| ------------------- | --------------: |
| Naive CD4 T cells   |             734 |
| CD14+ Monocytes     |             484 |
| Memory CD4 T cells  |             451 |
| B cells             |             344 |
| CD8 T cells         |             276 |
| CD16+ Monocytes     |             160 |
| NK cells            |             145 |
| Dendritic cells     |              30 |
| Platelets           |              14 |

Cell type labels were assigned manually using cluster-specific marker genes and canonical PBMC markers. These annotations are evidence-based labels, not experimentally validated ground-truth identities.

## Dataset

The project uses the publicly available 10x Genomics PBMC3K filtered gene-cell matrix.

Expected local data structure:

```text
data/
└── filtered_gene_bc_matrices/
    └── hg19/
        ├── barcodes.tsv
        ├── genes.tsv
        └── matrix.mtx
```

The raw data files are not included in this repository. To reproduce the analysis, download the filtered PBMC3K matrix from 10x Genomics and place the extracted files in the folder shown above.

## Methods

The analysis workflow is implemented in:

```text
scripts/01_pbmc3k_analysis.R
```

Main analysis steps:

1. Load the filtered 10x Genomics PBMC3K matrix
2. Explore raw gene and cell detection patterns
3. Perform gene-level filter sensitivity checks
4. Create a Seurat object using `min.cells = 5`
5. Calculate mitochondrial RNA percentage
6. Apply cell-level QC filtering
7. Normalize expression using `LogNormalize`
8. Identify 2,000 highly variable genes
9. Scale data and run PCA
10. Evaluate PC selection using an elbow plot and PC sensitivity checks
11. Cluster cells using Seurat graph-based clustering
12. Run 2D UMAP and optional 3D UMAP
13. Identify marker genes for each cluster
14. Annotate clusters using marker genes and canonical PBMC markers
15. Save the final annotated Seurat object for the Shiny app

## Key Analysis Decisions

| Step                 | Decision                                                                                                          |
| -------------------- | ----------------------------------------------------------------------------------------------------------------- |
| Gene-level filtering | Genes detected in fewer than 5 cells were removed using `min.cells = 5`.                                          |
| Cell QC              | Cells were retained if `nFeature_RNA > 200`, `nFeature_RNA < 2500`, and `percent.mt < 5`.                         |
| Normalization        | Expression values were normalized using Seurat `LogNormalize` with `scale.factor = 10000`.                        |
| Variable features    | The top 2,000 highly variable genes were selected using the `vst` method.                                         |
| PC selection         | PCs 1–8 were used based on the elbow plot and sensitivity checks across 6–13 PCs.                                 |
| Clustering           | Cells were clustered with Seurat graph-based clustering using `resolution = 0.5`.                                 |
| Marker genes         | Markers were identified with `FindAllMarkers`, `only.pos = TRUE`, `min.pct = 0.25`, and `logfc.threshold = 0.25`. |
| Annotation           | Clusters were manually annotated using top marker genes and canonical PBMC markers.                               |

## Cluster Annotation Evidence

| Cluster | Annotated cell type | Supporting markers                            |
| ------: | ------------------- | --------------------------------------------- |
|       0 | Naive CD4 T cells   | `CCR7`, `LEF1`, `TCF7`, `MAL`, `LDHB`         |
|       1 | CD14+ Monocytes     | `S100A8`, `S100A9`, `CD14`, `LYZ`, `MS4A6A`   |
|       2 | Memory CD4 T cells  | `IL7R`, `IL32`, `CD2`, `CD40LG`, `TRAT1`      |
|       3 | B cells             | `CD79A`, `MS4A1`, `TCL1A`, `FCER2`, `HLA-DOB` |
|       4 | CD8 T cells         | `CD8A`, `GZMK`, `CCL5`, `NKG7`, `CST7`        |
|       5 | CD16+ Monocytes     | `FCGR3A`, `MS4A7`, `CDKN1C`, `LILRA3`         |
|       6 | NK cells            | `GNLY`, `GZMB`, `FGFBP2`, `XCL2`, `SPON2`     |
|       7 | Dendritic cells     | `FCER1A`, `CD1C`, `CLEC10A`, `SERPINF1`       |
|       8 | Platelets           | `PPBP`, `PF4`, `GP9`, `ITGA2B`                |

## Shiny App

The interactive Shiny app is implemented in:

```text
app.R
```

The app loads the completed annotated Seurat object:

```text
outputs/pbmc_annotated.rds
```

The app includes:

* Overview of cell type counts
* Cell type distribution bar chart
* 2D UMAP by annotated cell type
* Interactive 3D UMAP
* Gene expression FeaturePlot
* Gene expression violin plot by annotated cell type
* Suggested marker genes for each cell type
* Top marker gene table with selectable number of markers per cluster
* Methods tab explaining analysis decisions
* Glossary explaining single-cell and immune-cell terminology

The app does not rerun the full analysis pipeline. It uses the saved annotated object produced by the analysis script.

## How to Run

### 1. Install required R packages

```r
install.packages(c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "patchwork",
  "dplyr",
  "mclust",
  "shiny",
  "bslib",
  "plotly",
  "DT"
))
```

### 2. Run the analysis script

From the project root:

```r
source("scripts/01_pbmc3k_analysis.R")
```

This creates the processed outputs, figures, marker tables, annotation evidence table, and final annotated Seurat object.

### 3. Run the Shiny app

After the analysis script completes successfully:

```r
shiny::runApp()
```

## Project Structure

```text
pbmc3k-single-cell-project/
├── app.R
├── README.md
├── data/
│   └── filtered_gene_bc_matrices/
│       └── hg19/
│           ├── barcodes.tsv
│           ├── genes.tsv
│           └── matrix.mtx
├── figures/
│   ├── canonical_marker_dotplot.png
│   ├── gene_filter_sensitivity.png
│   ├── genes_by_detection_level.png
│   ├── pbmc3k_umap_cell_types.png
│   ├── pbmc3k_umap_clusters.png
│   ├── pc_sensitivity_clusters.png
│   ├── pca_cells_diagnostic.png
│   ├── pca_elbow_plot.png
│   ├── qc_violin_after_filtering.png
│   ├── variable_feature_detection_frequency.png
│   └── variable_features.png
├── outputs/
│   ├── ari_vs_8_pcs.csv
│   ├── canonical_marker_filter_check.csv
│   ├── cell_type_counts.csv
│   ├── cluster_annotation_evidence.csv
│   ├── cluster_counts.csv
│   ├── final_sanity_checks.csv
│   ├── marker_sensitivity_min_pct_summary.csv
│   ├── pbmc_annotated.rds
│   ├── pbmc_cluster_markers.csv
│   ├── pc_sensitivity_summary.csv
│   ├── pca_variance_explained.csv
│   ├── qc_filter_summary.csv
│   ├── session_info.txt
│   ├── top_10_markers_per_cluster.csv
│   ├── umap_3d_coordinates.csv
│   └── variable_feature_detection_frequency_summary.csv
└── scripts/
    └── 01_pbmc3k_analysis.R
```

## Main Outputs

| Output                                    | Description                                         |
| ----------------------------------------- | --------------------------------------------------- |
| `outputs/pbmc_annotated.rds`              | Final annotated Seurat object used by the Shiny app |
| `outputs/pbmc_cluster_markers.csv`        | Full marker gene table for all clusters             |
| `outputs/top_10_markers_per_cluster.csv`  | Top marker genes for each cluster                   |
| `outputs/cluster_annotation_evidence.csv` | Manual annotation evidence table                    |
| `outputs/cell_type_counts.csv`            | Number of cells per annotated cell type             |
| `outputs/final_sanity_checks.csv`         | Final reproducibility checks                        |
| `figures/pbmc3k_umap_cell_types.png`      | 2D UMAP by annotated cell type                      |
| `figures/canonical_marker_dotplot.png`    | Canonical PBMC marker expression by cluster         |

## Interpretation Notes

UMAP axes are visualization coordinates and should not be interpreted as direct biological measurements. Cells that appear close together on the UMAP have similar gene-expression profiles, but the exact axis values do not represent specific biological quantities.

Marker gene expression was used to support manual cluster annotation. These labels are intended as evidence-based biological interpretations of the clusters.

## Skills Demonstrated

This project demonstrates:

* Single-cell RNA-seq analysis in R
* Seurat workflow development
* Quality control and filtering
* Dimensionality reduction with PCA and UMAP
* Graph-based clustering
* Marker gene detection and interpretation
* Manual cell type annotation
* Sensitivity analysis for filtering and PC selection
* Shiny app development
* Interactive data visualization
* Reproducible project organization
* Clear documentation of analytical decisions

## Tools Used

* R
* Seurat
* ggplot2
* dplyr
* Shiny
* bslib
* plotly
* DT
* mclust
