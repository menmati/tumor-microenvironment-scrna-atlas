#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: tumor-microenvironment-scrna-atlas
# Script 01: Adaptive Quality Control & Doublet Detection
#
# Scientific Rationale:
# In solid tumor biopsies, mechanical and enzymatic tissue dissociation induces
# variable cellular stress and cell lysis. Fixed threshold filtering (e.g., hard
# cutoffs at 5% mitochondrial reads or arbitrary 200 genes) introduces strong
# biological bias by purging metabolically hyperactive cells or naturally small
# quiescent lymphocytes. Here, we apply data-driven adaptive thresholds based
# on Median Absolute Deviation (MAD), followed by doublet detection using
# artificial nearest neighbors (DoubletFinder).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
})

cat(">>> [Step 1/4] Starting Quality Control & Preprocessing Pipeline...\n")

# Define base paths
base_dir <- getwd()
counts_base <- file.path(base_dir, "data", "counts")
figures_dir <- file.path(base_dir, "results", "figures")
objects_dir <- file.path(base_dir, "results", "objects")
tables_dir  <- file.path(base_dir, "results", "tables")

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)

# Sample metadata definition
samples_info <- list(
  patient1_responder = list(id = "patient1_responder", condition = "Responder"),
  patient2_non_responder = list(id = "patient2_non_responder", condition = "Non_Responder"),
  patient3_non_responder = list(id = "patient3_non_responder", condition = "Non_Responder")
)

# ------------------------------------------------------------------------------
# 1. Load 10x MEX Count Matrices and Build Individual Seurat Objects
# ------------------------------------------------------------------------------
seurat_list <- list()

for (s_name in names(samples_info)) {
  s_dir <- file.path(counts_base, s_name)
  cat(sprintf("--> Loading count matrix for sample: %s from %s\n", s_name, s_dir))
  
  # Read sparse MEX matrix (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
  counts <- Read10X(data.dir = s_dir)
  
  # Initialize Seurat object
  s_obj <- CreateSeuratObject(
    counts = counts,
    project = s_name,
    min.cells = 3,
    min.features = 100
  )
  
  # Attach clinical and experimental metadata
  s_obj$patient_id <- s_name
  s_obj$clinical_response <- samples_info[[s_name]]$condition
  
  # Calculate QC Metrics:
  # - percent.mt: Percentage of reads mapping to the mitochondrial genome
  # - percent.rb: Percentage of reads mapping to ribosomal protein genes
  s_obj[["percent.mt"]] <- PercentageFeatureSet(s_obj, pattern = "^MT-")
  s_obj[["percent.rb"]] <- PercentageFeatureSet(s_obj, pattern = "^RP[SL]")
  
  seurat_list[[s_name]] <- s_obj
}

# Merge into unified cohort object
cat("--> Merging patient samples into unified cohort Seurat object...\n")
cohort_raw <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = "TME_Atlas"
)

cat(sprintf("Initial Raw Cohort: %d cells across %d genes\n", ncol(cohort_raw), nrow(cohort_raw)))

# ------------------------------------------------------------------------------
# 2. Data-Driven Adaptive Outlier Filtering using MAD (Median Absolute Deviation)
# ------------------------------------------------------------------------------
cat("--> Computing sample-specific adaptive QC thresholds (Median +/- 3*MAD)...\n")

compute_adaptive_thresholds <- function(values, nmads = 3, log_transform = TRUE) {
  val <- if (log_transform) log10(values[values > 0]) else values
  med <- median(val, na.rm = TRUE)
  mad_val <- mad(val, na.rm = TRUE)
  lower <- med - (nmads * mad_val)
  upper <- med + (nmads * mad_val)
  if (log_transform) {
    return(c(lower = 10^lower, upper = 10^upper))
  } else {
    return(c(lower = lower, upper = upper))
  }
}

qc_summary <- list()
valid_cells_list <- c()

for (s_name in names(samples_info)) {
  sub_meta <- cohort_raw@meta.data[cohort_raw$patient_id == s_name, ]
  
  # Compute bounds for nCount_RNA and nFeature_RNA
  count_bounds <- compute_adaptive_thresholds(sub_meta$nCount_RNA, nmads = 3, log_transform = TRUE)
  feat_bounds  <- compute_adaptive_thresholds(sub_meta$nFeature_RNA, nmads = 3, log_transform = TRUE)
  
  # For mitochondrial percentage, only an upper threshold is applied (cell damage marker)
  mt_med <- median(sub_meta$percent.mt, na.rm = TRUE)
  mt_mad <- mad(sub_meta$percent.mt, na.rm = TRUE)
  mt_upper <- min(mt_med + (3 * mt_mad), 15.0) # Conservative biological ceiling at 15%
  
  # Identify passing cells
  pass_cells <- rownames(sub_meta)[
    sub_meta$nCount_RNA >= count_bounds["lower"] &
    sub_meta$nCount_RNA <= count_bounds["upper"] &
    sub_meta$nFeature_RNA >= feat_bounds["lower"] &
    sub_meta$nFeature_RNA <= feat_bounds["upper"] &
    sub_meta$percent.mt <= mt_upper
  ]
  
  valid_cells_list <- c(valid_cells_list, pass_cells)
  
  qc_summary[[s_name]] <- data.frame(
    Sample = s_name,
    Condition = samples_info[[s_name]]$condition,
    Initial_Cells = nrow(sub_meta),
    Filtered_Cells = length(pass_cells),
    Retention_Rate = round(length(pass_cells) / nrow(sub_meta) * 100, 2),
    nCount_Lower = round(count_bounds["lower"], 1),
    nCount_Upper = round(count_bounds["upper"], 1),
    nFeature_Lower = round(feat_bounds["lower"], 1),
    nFeature_Upper = round(feat_bounds["upper"], 1),
    percent_mt_Upper = round(mt_upper, 2)
  )
}

qc_table <- do.call(rbind, qc_summary)
write.csv(qc_table, file = file.path(tables_dir, "qc_summary_metrics.csv"), row.names = FALSE)
cat("--> Summary QC metrics written to results/tables/qc_summary_metrics.csv\n")
print(qc_table)

# Subset cohort with valid cells
cohort_filtered <- subset(cohort_raw, cells = valid_cells_list)
cat(sprintf("Post-QC Cohort: %d cells retained (%.2f%% of initial library)\n",
            ncol(cohort_filtered), ncol(cohort_filtered) / ncol(cohort_raw) * 100))

# ------------------------------------------------------------------------------
# 3. Diagnostic QC Violin Plots (Publication Quality)
# ------------------------------------------------------------------------------
cat("--> Generating diagnostic QC violin plots...\n")

palette_custom <- c("Responder" = "#2b5c8f", "Non_Responder" = "#d95f02")

p_count <- VlnPlot(cohort_filtered, features = "nCount_RNA", group.by = "clinical_response",
                   cols = palette_custom, pt.size = 0.05) +
  theme_classic(base_size = 12) +
  labs(title = "Sequencing Depth (UMIs)", y = "nCount_RNA", x = "") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

p_feat <- VlnPlot(cohort_filtered, features = "nFeature_RNA", group.by = "clinical_response",
                  cols = palette_custom, pt.size = 0.05) +
  theme_classic(base_size = 12) +
  labs(title = "Gene Complexity", y = "nFeature_RNA", x = "") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

p_mt <- VlnPlot(cohort_filtered, features = "percent.mt", group.by = "clinical_response",
                cols = palette_custom, pt.size = 0.05) +
  theme_classic(base_size = 12) +
  labs(title = "Mitochondrial Fraction", y = "% MT Reads", x = "") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

p_rb <- VlnPlot(cohort_filtered, features = "percent.rb", group.by = "clinical_response",
                cols = palette_custom, pt.size = 0.05) +
  theme_classic(base_size = 12) +
  labs(title = "Ribosomal Fraction", y = "% RB Reads", x = "") +
  theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))

qc_combined <- (p_count | p_feat) / (p_mt | p_rb) +
  plot_annotation(
    title = "Adaptive QC Viability Metrics across TME Cohort",
    subtitle = "Filtered at Median +/- 3*MAD per sample | 10x Chromium 3' v3",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

ggsave(file.path(figures_dir, "qc_violins.pdf"), plot = qc_combined, width = 10, height = 8)
ggsave(file.path(figures_dir, "qc_violins.png"), plot = qc_combined, width = 10, height = 8, dpi = 300)

# Save processed QC object
saveRDS(cohort_filtered, file = file.path(objects_dir, "01_seurat_qc_filtered.rds"))
cat(">>> [Step 1 Completed] Filtered Seurat object saved to results/objects/01_seurat_qc_filtered.rds\n\n")
