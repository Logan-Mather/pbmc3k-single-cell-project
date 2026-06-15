
# app.R

# PBMC3K single-cell RNA-seq Shiny app
# This app loads the completed annotated Seurat object.
# It does not rerun the full analysis pipeline.

library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)


# -----------------------------------------
# Load completed analysis object
# -----------------------------------------

pbmc <- readRDS("outputs/pbmc_annotated.rds")

DefaultAssay(pbmc) <- "RNA"

# Alphabetized gene list for easier search
gene_choices <- sort(rownames(pbmc))

default_gene <- if ("MS4A1" %in% gene_choices) {
  "MS4A1"
} else {
  gene_choices[1]
}


# -----------------------------------------
# App-level summary data
# -----------------------------------------

cell_type_counts <- as.data.frame(
  table(pbmc$cell_type),
  stringsAsFactors = FALSE
)

names(cell_type_counts) <- c("Cell type", "Cells")

cell_type_counts <- cell_type_counts %>%
  arrange(desc(Cells))

# Make cell type order consistent across plots
pbmc$cell_type <- factor(
  as.character(pbmc$cell_type),
  levels = cell_type_counts$`Cell type`
)

n_cells <- ncol(pbmc)
n_cell_types <- length(unique(pbmc$cell_type))
n_clusters <- length(unique(pbmc$seurat_clusters))
n_variable_features <- length(VariableFeatures(pbmc))
n_pcs_used <- 8
cluster_resolution <- 0.5


# -----------------------------------------
# Cell type color palette
# -----------------------------------------

cell_type_colors <- c(
  "Naive CD4 T cells" = "#56B4E9",
  "CD14+ Monocytes" = "#E69F00",
  "Memory CD4 T cells" = "#009E73",
  "B cells" = "#CC79A7",
  "CD8 T cells" = "#F0E442",
  "CD16+ Monocytes" = "#0072B2",
  "NK cells" = "#D55E00",
  "Dendritic cells" = "#00BFC4",
  "Platelets" = "#999999"
)

cell_type_colors <- cell_type_colors[levels(pbmc$cell_type)]


# -----------------------------------------
# 2D UMAP label positions for gene expression plots
# -----------------------------------------

umap_2d_embeddings <- Embeddings(
  pbmc,
  reduction = "umap"
)

umap_label_df <- data.frame(
  UMAP_1 = umap_2d_embeddings[, 1],
  UMAP_2 = umap_2d_embeddings[, 2],
  cell_type = pbmc$cell_type
)

umap_label_positions <- umap_label_df %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )


# -----------------------------------------
# Suggested marker genes for gene expression tab
# -----------------------------------------

marker_examples <- data.frame(
  cell_type = c(
    "Naive CD4 T cells", "Naive CD4 T cells",
    "Memory CD4 T cells", "Memory CD4 T cells",
    "CD14+ Monocytes", "CD14+ Monocytes",
    "CD16+ Monocytes", "CD16+ Monocytes",
    "B cells", "B cells",
    "CD8 T cells", "CD8 T cells",
    "NK cells", "NK cells",
    "Dendritic cells", "Dendritic cells",
    "Platelets", "Platelets"
  ),
  gene = c(
    "CCR7", "LEF1",
    "IL7R", "CD40LG",
    "CD14", "LYZ",
    "FCGR3A", "MS4A7",
    "MS4A1", "CD79A",
    "CD8A", "GZMK",
    "GNLY", "NKG7",
    "FCER1A", "CD1C",
    "PPBP", "PF4"
  ),
  stringsAsFactors = FALSE
)

marker_examples <- marker_examples %>%
  filter(gene %in% gene_choices)


# -----------------------------------------
# Load supporting tables
# -----------------------------------------

annotation_evidence <- read.csv(
  "outputs/cluster_annotation_evidence.csv",
  stringsAsFactors = FALSE
)

annotation_evidence <- annotation_evidence %>%
  mutate(cluster = as.character(cluster))

all_markers <- read.csv(
  "outputs/pbmc_cluster_markers.csv",
  stringsAsFactors = FALSE
)

all_markers <- all_markers %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(
    annotation_evidence %>%
      select(cluster, cell_type),
    by = "cluster"
  ) %>%
  select(
    cluster,
    cell_type,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  )


# -----------------------------------------
# Methods table
# -----------------------------------------

methods_table <- data.frame(
  Step = c(
    "Data source",
    "Gene-level filtering",
    "Cell-level quality control",
    "Normalization",
    "Variable feature selection",
    "PCA",
    "PC selection",
    "Clustering",
    "2D UMAP",
    "Marker gene detection",
    "Cell type annotation",
    "3D UMAP"
  ),
  Decision = c(
    "The analysis uses the 10x Genomics PBMC3K filtered gene-cell matrix.",
    
    "Genes detected in fewer than 5 cells were removed using min.cells = 5. This removes extremely rare genes while avoiding an overly aggressive filter that could remove markers useful for rare cell populations.",
    
    "Cells were retained if nFeature_RNA > 200, nFeature_RNA < 2500, and percent.mt < 5. These thresholds remove very low-feature cells, possible high-feature doublets, and cells with high mitochondrial RNA percentages.",
    
    "Expression values were normalized using Seurat's LogNormalize method with scale.factor = 10000. This adjusts for differences in total RNA counts across cells and applies a log transformation.",
    
    "The top 2,000 highly variable genes were selected using the vst method. These genes capture major expression differences across cells and were used for PCA.",
    
    "PCA was performed using the highly variable genes to reduce the dimensionality of the dataset before clustering and UMAP visualization.",
    
    "The first 8 PCs were used for downstream analysis. This choice was based on the elbow plot and a sensitivity check across 6-13 PCs. Nine clusters were recovered beginning at 8 PCs, and additional PCs did not change the number of clusters.",
    
    "Cells were clustered using Seurat graph-based clustering with resolution = 0.5. In Seurat, resolution controls the granularity of community detection on the shared nearest-neighbor graph; higher values generally produce more clusters.",
    
    "A 2D UMAP was used as the main visualization of cell similarity. UMAP axes are visualization coordinates and are not direct biological measurements.",
    
    "Marker genes were identified using FindAllMarkers with only genes more highly expressed in each cluster, min.pct = 0.25, and logfc.threshold = 0.25. This focuses the marker table on genes useful for identifying each cluster.",
    
    "Clusters were manually annotated using top marker genes and canonical PBMC markers. These labels are evidence-based annotations, not experimentally validated ground-truth identities.",
    
    "A 3D UMAP was generated for optional interactive exploration in the Shiny app. The 2D UMAP remains the primary static visualization."
  ),
  stringsAsFactors = FALSE
)


# -----------------------------------------
# Glossary
# -----------------------------------------

glossary_terms <- data.frame(
  Term = c(
    "Adjusted Rand Index",
    "avg_log2FC",
    "B cells",
    "canonical marker",
    "CD14+ Monocytes",
    "CD16+ Monocytes",
    "CD8 T cells",
    "cell / barcode",
    "cell type annotation",
    "cluster",
    "Dendritic cells",
    "DotPlot",
    "elbow plot",
    "FeaturePlot",
    "gene / feature",
    "LogNormalize",
    "log-normalized expression",
    "marker gene",
    "Memory CD4 T cells",
    "Naive CD4 T cells",
    "nearest-neighbor graph",
    "nCount_RNA",
    "NK cells",
    "normalization",
    "nFeature_RNA",
    "PBMC",
    "PCA",
    "percent.mt",
    "Platelets",
    "principal component",
    "QC filtering",
    "resolution",
    "Seurat",
    "single-cell RNA-seq",
    "UMAP",
    "UMAP 1 / UMAP 2 / UMAP 3",
    "UMI / count",
    "variable feature",
    "VlnPlot"
  ),
  Definition = c(
    "A score used to compare two clustering results. Values closer to 1 indicate more similar cluster assignments.",
    "Average log2 fold-change. Larger positive values mean stronger expression in that cluster compared with other cells.",
    "Antibody-producing immune cells. In this PBMC dataset, B cells are identified by markers such as MS4A1 and CD79A.",
    "A well-established gene used to identify a known cell type.",
    "Classical monocytes involved in innate immune responses, inflammation, and pathogen detection. In this dataset, they are identified by markers such as CD14, LYZ, S100A8, and S100A9.",
    "Non-classical monocytes involved in immune surveillance and inflammatory responses. In this dataset, they are identified by markers such as FCGR3A and MS4A7.",
    "T cells involved in killing infected or abnormal cells. In this dataset, they are identified by markers such as CD8A, GZMK, CCL5, and NKG7.",
    "A single captured cell represented by a barcode.",
    "Assigning biological names to clusters based on marker genes.",
    "A group of cells with similar gene-expression profiles.",
    "Antigen-presenting immune cells that help activate T-cell responses. In this dataset, they are identified by markers such as FCER1A, CD1C, and CLEC10A.",
    "A Seurat plot showing both marker expression level and percent of cells expressing a gene.",
    "A diagnostic plot used to help choose how many principal components to use.",
    "A Seurat plot showing expression of a selected gene across the UMAP.",
    "A gene measured in the dataset. In Seurat output, genes are often called features.",
    "A normalization method that scales each cell by total counts, multiplies by 10,000, and applies a log transformation.",
    "Expression values after library-size scaling and log transformation.",
    "A gene more highly expressed in one cluster compared with other cells.",
    "CD4 T cells with prior immune activation or memory-like features. In this dataset, they are identified by markers such as IL7R, IL32, CD2, and CD40LG.",
    "CD4 T cells that have not yet strongly differentiated into specialized effector states. In this dataset, they are identified by markers such as CCR7, LEF1, TCF7, and MAL.",
    "A graph connecting cells with similar expression profiles. Clustering is based on this graph.",
    "The total RNA count detected in a cell.",
    "Natural killer cells involved in innate immune defense and killing abnormal cells. In this dataset, they are identified by markers such as GNLY, GZMB, NKG7, and FGFBP2.",
    "A process that adjusts expression values so cells can be compared more fairly.",
    "The number of genes detected in a cell.",
    "Peripheral blood mononuclear cells. These include immune cells such as T cells, B cells, monocytes, NK cells, and dendritic cells.",
    "Principal component analysis. A dimensionality reduction method used before clustering and UMAP.",
    "The percentage of a cell's RNA counts coming from mitochondrial genes. High values can indicate stressed or low-quality cells.",
    "Small blood-cell fragments involved in clotting. In this dataset, platelet-like cells are identified by markers such as PPBP, PF4, GP9, and ITGA2B.",
    "A new axis summarizing a major pattern of variation in the gene expression data.",
    "The process of removing low-quality cells before downstream analysis.",
    "A Seurat FindClusters parameter used during graph-based community detection. After Seurat builds a shared nearest-neighbor graph, the resolution value controls the granularity of the modularity optimization step. Higher resolution values generally split the graph into more smaller clusters, while lower values generally produce fewer broader clusters.",
    "An R package used for single-cell RNA-seq processing, clustering, visualization, and marker gene analysis.",
    "A method for measuring gene expression separately in individual cells.",
    "A dimensionality reduction method used to visualize similar cells near each other in 2D or 3D.",
    "Visualization coordinates created by UMAP. They are not direct biological measurements.",
    "A molecule-level count used to estimate gene expression in a cell.",
    "A gene with high variation across cells, used for PCA and clustering.",
    "A plot showing the distribution of expression values across cell types or clusters."
  ),
  stringsAsFactors = FALSE
)

glossary_terms <- glossary_terms %>%
  arrange(tolower(Term))


# -----------------------------------------
# Dark gray theme setup
# -----------------------------------------

app_bg <- "#181818"
panel_bg <- "#242424"
panel_bg_2 <- "#303030"
text_main <- "#f0f0f0"
text_muted <- "#cfcfcf"
grid_col <- "#4a4a4a"
accent_blue <- "#4ea3ff"

app_theme <- bs_theme(
  version = 5,
  bg = app_bg,
  fg = text_main,
  primary = accent_blue
)

dark_plot_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = app_bg, color = NA),
    panel.background = element_rect(fill = app_bg, color = NA),
    panel.grid.major = element_line(color = grid_col),
    panel.grid.minor = element_blank(),
    text = element_text(color = text_main),
    axis.text = element_text(color = text_main),
    axis.title = element_text(color = text_main),
    plot.title = element_text(color = "#ffffff", face = "bold"),
    plot.subtitle = element_text(color = text_muted),
    legend.background = element_rect(fill = app_bg, color = NA),
    legend.key = element_rect(fill = app_bg, color = NA),
    legend.text = element_text(color = text_main),
    legend.title = element_text(color = "#ffffff")
  )

summary_card <- function(label, value) {
  div(
    class = "summary-card",
    div(class = "summary-card-value", value),
    div(class = "summary-card-label", label)
  )
}


# -----------------------------------------
# User interface
# -----------------------------------------

ui <- fluidPage(
  
  theme = app_theme,
  
  tags$head(
    tags$style(
      HTML("
        body {
          background-color: #181818;
          color: #f0f0f0;
        }

        .container-fluid {
          background-color: #181818;
        }

        .tab-content {
          background-color: #181818;
          color: #f0f0f0;
          padding-top: 12px;
        }

        .nav-tabs .nav-link {
          color: #d8d8d8;
          background-color: #202020;
          border-color: #3a3a3a;
        }

        .nav-tabs .nav-link.active {
          background-color: #303030;
          color: #ffffff;
          border-color: #5a5a5a;
        }

        .summary-card {
          background-color: #242424;
          border: 1px solid #4a4a4a;
          border-radius: 10px;
          padding: 16px;
          margin-bottom: 12px;
          text-align: center;
        }

        .summary-card-value {
          font-size: 28px;
          font-weight: 700;
          color: #ffffff;
        }

        .summary-card-label {
          font-size: 13px;
          color: #cfcfcf;
          margin-top: 4px;
        }

        .gene-search-box {
          background-color: #242424;
          border: 1px solid #4a4a4a;
          border-radius: 10px;
          padding: 15px;
          margin-bottom: 12px;
        }

        .marker-button {
          margin: 3px;
        }

        .selectize-input,
        .selectize-dropdown {
          background-color: #303030 !important;
          color: #f0f0f0 !important;
          border-color: #5a5a5a !important;
        }

        .selectize-input input {
          color: #f0f0f0 !important;
        }

        .selectize-dropdown-content .option {
          background-color: #303030 !important;
          color: #f0f0f0 !important;
        }

        .selectize-dropdown-content .active {
          background-color: #4a4a4a !important;
          color: #ffffff !important;
        }

        .form-control {
          background-color: #303030 !important;
          color: #f0f0f0 !important;
          border-color: #5a5a5a !important;
        }

        .dataTables_wrapper {
          color: #f0f0f0 !important;
        }

        table.dataTable {
          color: #f0f0f0 !important;
          background-color: #242424 !important;
        }

        table.dataTable tbody tr {
          background-color: #242424 !important;
          color: #f0f0f0 !important;
        }

        table.dataTable tbody tr:hover {
          background-color: #303030 !important;
        }

        table.dataTable thead th {
          background-color: #303030 !important;
          color: #ffffff !important;
          border-bottom: 1px solid #5a5a5a !important;
        }

        .dataTables_filter input {
          background-color: #303030 !important;
          color: #f0f0f0 !important;
          border: 1px solid #5a5a5a !important;
        }

        .btn {
          margin-top: 4px;
        }

        hr {
          border-top: 1px solid #4a4a4a;
        }
      ")
    )
  ),
  
  titlePanel("PBMC3K Single-Cell RNA-seq Explorer"),
  
  tabsetPanel(
    
    tabPanel(
      "Overview",
      br(),
      
      h3("Project overview"),
      p("This Shiny app explores a PBMC3K single-cell RNA-seq analysis performed with Seurat."),
      p("Cells were filtered, normalized, clustered, visualized with UMAP, and manually annotated using cluster-specific marker genes and canonical PBMC markers."),
      p("Each point in the UMAP plots represents one cell. UMAP shows similarity in gene-expression profiles, not physical distance. Cell-type annotations are evidence-based labels assigned from marker genes."),
      
      br(),
      
      fluidRow(
        column(width = 2, summary_card("Cells analyzed", format(n_cells, big.mark = ","))),
        column(width = 2, summary_card("Cell types", n_cell_types)),
        column(width = 2, summary_card("Clusters", n_clusters)),
        column(width = 2, summary_card("Variable genes", format(n_variable_features, big.mark = ","))),
        column(width = 2, summary_card("PCs used", n_pcs_used)),
        column(width = 2, summary_card("Clustering resolution", cluster_resolution))
      ),
      
      p("Clustering resolution controls how finely Seurat splits cells into graph-based communities. This analysis used resolution = 0.5, which produced 9 clusters."),
      
      br(),
      
      fluidRow(
        column(
          width = 5,
          h4("Annotated cell type counts"),
          DTOutput("cell_type_counts_table"),
          br(),
          downloadButton("download_cell_type_counts", "Download cell type counts")
        ),
        column(
          width = 7,
          h4("Cell type distribution"),
          plotOutput("cell_type_counts_plot", height = "420px")
        )
      ),
      
      br(),
      hr(),
      
      h4("Cluster annotation evidence"),
      p("This table shows how each cluster was manually assigned a biological cell-type label."),
      DTOutput("annotation_evidence_table"),
      br(),
      downloadButton("download_annotation_evidence", "Download annotation evidence")
    ),
    
    tabPanel(
      "2D UMAP",
      br(),
      h3("2D UMAP by annotated cell type"),
      p("Each point is a cell. Nearby cells have similar gene-expression profiles."),
      p("Cell-type labels are shown directly on the plot, so the legend is removed to keep the figure cleaner."),
      plotOutput("umap_2d_plot", height = "700px")
    ),
    
    tabPanel(
      "3D UMAP",
      br(),
      h3("Interactive 3D UMAP"),
      p("The 3D UMAP adds a third visualization coordinate for interactive exploration. UMAP 1, UMAP 2, and UMAP 3 are not direct biological measurements."),
      plotlyOutput("umap_3d_plot", height = "750px")
    ),
    
    tabPanel(
      "Gene Expression",
      br(),
      h3("Gene expression tools"),
      p("Use this tab to explore expression of a selected gene across the annotated PBMC cell types."),
      
      fluidRow(
        column(
          width = 5,
          div(
            class = "gene-search-box",
            h4("Select a gene"),
            p("Choose a gene by typing its symbol or by using the dropdown list."),
            
            textInput(
              inputId = "gene_text",
              label = "Search by typing:",
              value = default_gene,
              placeholder = "Type a gene symbol, e.g. MS4A1, CD14, NKG7"
            ),
            
            actionButton(
              inputId = "use_typed_gene",
              label = "Use typed gene"
            ),
            
            br(),
            br(),
            
            selectizeInput(
              inputId = "gene_dropdown",
              label = "Select from dropdown:",
              choices = gene_choices,
              selected = default_gene,
              multiple = FALSE,
              options = list(
                placeholder = "Select a gene",
                maxOptions = 5000,
                create = FALSE
              )
            ),
            
            br(),
            strong("Current selected gene: "),
            textOutput("current_gene", inline = TRUE)
          )
        ),
        column(
          width = 7,
          h4("Suggested marker genes"),
          p("These are representative marker genes for each annotated cell type. Click a gene to update both gene expression plots."),
          uiOutput("marker_gene_buttons")
        )
      ),
      
      br(),
      hr(),
      
      h3("Gene expression on UMAP"),
      p("The FeaturePlot shows where the selected gene is expressed across cells in the 2D UMAP."),
      p("Color represents normalized log-transformed expression from the RNA assay. Higher values indicate higher expression in that cell."),
      p("The normalized expression is based on log1p(gene counts / total cell counts × 10,000)."),
      plotOutput("feature_plot", height = "650px"),
      
      br(),
      hr(),
      
      h3("Gene expression by annotated cell type"),
      p("The violin plot below shows the distribution of the selected gene's normalized expression across annotated cell types."),
      p("This is useful for checking whether a marker gene is concentrated in the expected cell type."),
      plotOutput("violin_plot", height = "600px")
    ),
    
    tabPanel(
      "Top Markers",
      br(),
      h3("Top marker genes per cluster"),
      p("Markers shown here are genes more highly expressed in each cluster compared with the other cells. Genes lower in a cluster are not included in this marker table."),
      
      tags$ul(
        tags$li(strong("avg_log2FC: "), "Average log2 fold-change. Larger positive values mean the gene is more strongly expressed in that cluster compared with other cells."),
        tags$li(strong("pct.1: "), "Fraction of cells in that cluster where the gene was detected."),
        tags$li(strong("pct.2: "), "Fraction of cells outside that cluster where the gene was detected."),
        tags$li(strong("p_val_adj: "), "Adjusted p-value after correcting for many gene tests. Smaller values provide stronger statistical evidence.")
      ),
      
      selectInput(
        inputId = "n_markers",
        label = "Markers shown per cluster:",
        choices = c(5, 10, 15, 20, 25),
        selected = 10
      ),
      
      DTOutput("top_markers_table"),
      br(),
      downloadButton("download_top_markers", "Download current marker table")
    ),
    
    tabPanel(
      "Methods",
      br(),
      h3("Methods and analysis decisions"),
      p("This tab summarizes the major analysis choices used to generate the annotated PBMC3K object."),
      DTOutput("methods_table")
    ),
    
    tabPanel(
      "Glossary",
      br(),
      h3("Glossary"),
      p("This tab defines common terms used throughout the app and analysis."),
      p("Use the search box to quickly find a term."),
      DTOutput("glossary_table")
    )
  )
)


# -----------------------------------------
# Server logic
# -----------------------------------------

server <- function(input, output, session) {
  
  selected_gene_value <- reactiveVal(default_gene)
  
  match_gene_symbol <- function(gene_input) {
    gene_input <- trimws(gene_input)
    
    matched_gene <- gene_choices[
      toupper(gene_choices) == toupper(gene_input)
    ]
    
    if (length(matched_gene) > 0) {
      matched_gene[1]
    } else {
      gene_input
    }
  }
  
  selected_gene <- reactive({
    selected_gene_value()
  })
  
  output$current_gene <- renderText({
    selected_gene()
  })
  
  observeEvent(input$use_typed_gene, {
    typed_gene <- match_gene_symbol(input$gene_text)
    
    selected_gene_value(typed_gene)
    
    if (typed_gene %in% gene_choices) {
      updateSelectizeInput(
        session,
        inputId = "gene_dropdown",
        selected = typed_gene
      )
    }
  })
  
  observeEvent(input$gene_dropdown, {
    selected_gene_value(input$gene_dropdown)
    
    updateTextInput(
      session,
      inputId = "gene_text",
      value = input$gene_dropdown
    )
  })
  
  
  # -----------------------------------------
  # Suggested marker gene buttons
  # -----------------------------------------
  
  set_gene <- function(gene) {
    if (gene %in% gene_choices) {
      selected_gene_value(gene)
      
      updateTextInput(
        session,
        inputId = "gene_text",
        value = gene
      )
      
      updateSelectizeInput(
        session,
        inputId = "gene_dropdown",
        selected = gene
      )
    }
  }
  
  output$marker_gene_buttons <- renderUI({
    
    marker_groups <- split(marker_examples, marker_examples$cell_type)
    
    tagList(
      lapply(names(marker_groups), function(cell_type_name) {
        
        genes_for_type <- marker_groups[[cell_type_name]]$gene
        
        div(
          style = "margin-bottom: 10px;",
          strong(cell_type_name),
          br(),
          lapply(genes_for_type, function(gene_name) {
            actionButton(
              inputId = paste0("gene_button_", gene_name),
              label = gene_name,
              class = "btn-sm marker-button"
            )
          })
        )
      })
    )
  })
  
  lapply(unique(marker_examples$gene), function(gene_name) {
    observeEvent(input[[paste0("gene_button_", gene_name)]], {
      set_gene(gene_name)
    }, ignoreInit = TRUE)
  })
  
  
  # -----------------------------------------
  # Overview
  # -----------------------------------------
  
  output$cell_type_counts_table <- renderDT({
    datatable(
      cell_type_counts,
      rownames = FALSE,
      options = list(
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        lengthChange = FALSE,
        info = FALSE,
        autoWidth = TRUE
      )
    )
  })
  
  output$cell_type_counts_plot <- renderPlot({
    ggplot(
      cell_type_counts,
      aes(x = reorder(`Cell type`, Cells), y = Cells)
    ) +
      geom_col(fill = accent_blue) +
      geom_text(
        aes(label = Cells),
        hjust = -0.15,
        color = text_main,
        size = 4
      ) +
      coord_flip() +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.15))
      ) +
      labs(
        title = "Cells per annotated cell type",
        x = "Cell type",
        y = "Number of cells"
      ) +
      dark_plot_theme
  })
  
  output$annotation_evidence_table <- renderDT({
    datatable(
      annotation_evidence,
      rownames = FALSE,
      options = list(
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        lengthChange = FALSE,
        info = FALSE,
        autoWidth = TRUE,
        scrollX = TRUE
      )
    )
  })
  
  
  # -----------------------------------------
  # 2D UMAP
  # -----------------------------------------
  
  output$umap_2d_plot <- renderPlot({
    
    DimPlot(
      pbmc,
      reduction = "umap",
      group.by = "cell_type",
      label = TRUE,
      repel = TRUE,
      label.size = 5,
      label.box = TRUE,
      cols = cell_type_colors
    ) +
      NoLegend() +
      labs(
        title = "PBMC3K 2D UMAP by Annotated Cell Type",
        subtitle = "Each point is a cell; nearby cells have similar gene-expression profiles",
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      dark_plot_theme
  })
  
  
  # -----------------------------------------
  # 3D UMAP
  # -----------------------------------------
  
  output$umap_3d_plot <- renderPlotly({
    
    validate(
      need("umap_3d" %in% names(pbmc@reductions), "3D UMAP was not found in the Seurat object.")
    )
    
    umap_3d_embeddings <- Embeddings(
      pbmc,
      reduction = "umap_3d"
    )
    
    umap_3d_meta <- pbmc@meta.data[rownames(umap_3d_embeddings), ]
    
    umap_3d_df <- data.frame(
      cell = rownames(umap_3d_embeddings),
      UMAP_1 = umap_3d_embeddings[, 1],
      UMAP_2 = umap_3d_embeddings[, 2],
      UMAP_3 = umap_3d_embeddings[, 3],
      seurat_cluster = as.character(umap_3d_meta$seurat_clusters),
      cell_type = factor(
        as.character(umap_3d_meta$cell_type),
        levels = levels(pbmc$cell_type)
      ),
      stringsAsFactors = FALSE
    )
    
    umap_3d_df$hover_text <- paste0(
      "Cell: ", umap_3d_df$cell,
      "<br>Cell type: ", umap_3d_df$cell_type,
      "<br>Cluster: ", umap_3d_df$seurat_cluster
    )
    
    p <- plot_ly(
      data = umap_3d_df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      z = ~UMAP_3,
      color = ~cell_type,
      colors = cell_type_colors,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = 2.8,
        opacity = 0.85
      ),
      text = ~hover_text,
      hoverinfo = "text"
    )
    
    p <- layout(
      p,
      title = "PBMC3K 3D UMAP by Annotated Cell Type",
      paper_bgcolor = app_bg,
      plot_bgcolor = app_bg,
      font = list(color = text_main),
      scene = list(
        xaxis = list(
          title = "UMAP 1",
          color = text_main,
          gridcolor = grid_col,
          zerolinecolor = "#666666",
          backgroundcolor = app_bg
        ),
        yaxis = list(
          title = "UMAP 2",
          color = text_main,
          gridcolor = grid_col,
          zerolinecolor = "#666666",
          backgroundcolor = app_bg
        ),
        zaxis = list(
          title = "UMAP 3",
          color = text_main,
          gridcolor = grid_col,
          zerolinecolor = "#666666",
          backgroundcolor = app_bg
        )
      ),
      legend = list(
        title = list(text = "Cell type"),
        itemsizing = "constant",
        font = list(size = 13, color = text_main)
      )
    )
    
    p
  })
  
  
  # -----------------------------------------
  # Gene expression plots
  # -----------------------------------------
  
  output$feature_plot <- renderPlot({
    
    validate(
      need(!is.null(selected_gene()), "Select a gene."),
      need(selected_gene() %in% rownames(pbmc), "Selected gene not found in this dataset.")
    )
    
    FeaturePlot(
      pbmc,
      features = selected_gene(),
      reduction = "umap",
      cols = c("#3a3a3a", accent_blue)
    ) +
      geom_label(
        data = umap_label_positions,
        aes(
          x = UMAP_1,
          y = UMAP_2,
          label = cell_type
        ),
        inherit.aes = FALSE,
        fill = panel_bg,
        color = text_main,
        label.size = 0.25,
        size = 4
      ) +
      labs(
        title = paste("UMAP expression of", selected_gene()),
        subtitle = "Color indicates normalized log-transformed expression",
        color = "Normalized expression",
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      dark_plot_theme
  })
  
  output$violin_plot <- renderPlot({
    
    validate(
      need(!is.null(selected_gene()), "Select a gene."),
      need(selected_gene() %in% rownames(pbmc), "Selected gene not found in this dataset.")
    )
    
    VlnPlot(
      pbmc,
      features = selected_gene(),
      group.by = "cell_type",
      pt.size = 0,
      cols = cell_type_colors
    ) +
      NoLegend() +
      labs(
        title = paste("Expression of", selected_gene(), "by annotated cell type"),
        subtitle = "Expression values are normalized and log-transformed",
        x = "Annotated cell type",
        y = "Normalized expression"
      ) +
      dark_plot_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = text_main)
      )
  })
  
  
  # -----------------------------------------
  # Top markers
  # -----------------------------------------
  
  markers_to_show <- reactive({
    all_markers %>%
      group_by(cluster, cell_type) %>%
      slice_max(
        order_by = avg_log2FC,
        n = as.numeric(input$n_markers),
        with_ties = FALSE
      ) %>%
      ungroup() %>%
      arrange(as.numeric(cluster), desc(avg_log2FC))
  })
  
  output$top_markers_table <- renderDT({
    datatable(
      markers_to_show(),
      rownames = FALSE,
      options = list(
        dom = "ft",
        paging = FALSE,
        searching = TRUE,
        lengthChange = FALSE,
        info = FALSE,
        autoWidth = TRUE,
        scrollX = TRUE
      )
    )
  })
  
  
  # -----------------------------------------
  # Methods and glossary
  # -----------------------------------------
  
  output$methods_table <- renderDT({
    datatable(
      methods_table,
      rownames = FALSE,
      options = list(
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        lengthChange = FALSE,
        info = FALSE,
        autoWidth = TRUE,
        scrollX = TRUE
      )
    )
  })
  
  output$glossary_table <- renderDT({
    datatable(
      glossary_terms,
      rownames = FALSE,
      options = list(
        dom = "ft",
        paging = FALSE,
        searching = TRUE,
        lengthChange = FALSE,
        info = FALSE,
        autoWidth = TRUE,
        scrollX = TRUE
      )
    )
  })
  
  
  # -----------------------------------------
  # Downloads
  # -----------------------------------------
  
  output$download_top_markers <- downloadHandler(
    filename = function() {
      paste0("pbmc3k_top_", input$n_markers, "_markers_per_cluster.csv")
    },
    content = function(file) {
      write.csv(
        markers_to_show(),
        file,
        row.names = FALSE
      )
    }
  )
  
  output$download_annotation_evidence <- downloadHandler(
    filename = function() {
      "pbmc3k_cluster_annotation_evidence.csv"
    },
    content = function(file) {
      write.csv(
        annotation_evidence,
        file,
        row.names = FALSE
      )
    }
  )
  
  output$download_cell_type_counts <- downloadHandler(
    filename = function() {
      "pbmc3k_cell_type_counts.csv"
    },
    content = function(file) {
      write.csv(
        cell_type_counts,
        file,
        row.names = FALSE
      )
    }
  )
}


# -----------------------------------------
# Run app
# -----------------------------------------

shinyApp(ui = ui, server = server)
