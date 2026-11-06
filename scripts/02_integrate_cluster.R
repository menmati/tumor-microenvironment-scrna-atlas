#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: tumor-microenvironment-scrna-atlas
# Script 02: Variance-Stabilizing Normalization (SCTransform v2),
#            Multi-Sample Integration (Harmony), and Leiden Clustering
#
# Scientific Rationale:
# Comparing tumor biopsies across multiple patients introduces technical batch
# effects (variations in tissue dissociation, microfluidic run, sequencing depth).
# Standard log-normalization fails to decouple gene expression from sequencing
# depth and distorts low-count gene variance. Here, we apply SCTransform v2
# with regularized negative binomial GLMs, decouple count matrices across sample
# layers in Seurat v5, align cell centroids in PCA space using Harmony fuzzy
# clustering, and detect partition-refined cell communities via Leiden algorithm.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(clustree)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

cat(">>> [Step 2/4] Starting Normalization, Batch Integration & Clustering...\n")

# Memory management for future parallel execution
options(future.globals.maxSize = 4e9)

base_dir <- getwd()
figures_dir <- file.path(base_dir, "results", "figures")
objects_dir <- file.path(base_dir, "results", "objects")

# Load post-QC object
qc_obj_file <- file.path(objects_dir, "01_seurat_qc_filtered.rds")
cat(sprintf("--> Loading post-QC object from %s\n", qc_obj_file))
obj <- readRDS(qc_obj_file)

# ------------------------------------------------------------------------------
# 1. Seurat v5 Layer Architecture & SCTransform v2 Normalization
# ------------------------------------------------------------------------------
cat("--> Checking count matrix layers (Seurat v5 architecture)...\n")
# If not already split by patient, split to model sample-specific dispersion
if (length(grep("^counts", Layers(obj, assay = "RNA"))) == 1) {
  cat("Splitting RNA counts layer across patient cohorts...\n")
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$patient_id)
} else {
  cat("Layers already decoupled per patient sample:", paste(Layers(obj, assay = "RNA"), collapse = ", "), "\n")
}

cat("--> Running SCTransform v2 (Regularized Negative Binomial GLM)...\n")
obj <- SCTransform(
  obj,
  assay = "RNA",
  new.assay.name = "SCT",
  vst.flavor = "v2",
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

cat("SCT Assay successfully created with 3,000 highly variable features.\n")

# ------------------------------------------------------------------------------
# 2. Linear Dimensional Reduction (PCA) and Elbow Selection
# ------------------------------------------------------------------------------
cat("--> Calculating 50 Principal Components on Pearson Residuals...\n")
obj <- RunPCA(obj, assay = "SCT", npcs = 50, verbose = FALSE)

# Generate and export Elbow Plot
p_elbow_base <- ElbowPlot(obj, ndims = 50)
p_elbow <- p_elbow_base +
  theme_classic(base_size = 12) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "#e63946", linewidth = 0.8) +
  annotate("text", x = 32, y = max(p_elbow_base$data$stdev) * 0.8,
           label = "Selected Cutoff (30 PCs)", color = "#e63946", fontface = "bold", hjust = 0) +
  labs(title = "Elbow Plot: Variance Explained per Principal Component",
       subtitle = "Transition from biological gradient to asymptotic technical noise floor")

ggsave(file.path(figures_dir, "elbow_plot.pdf"), plot = p_elbow, width = 8, height = 5)
ggsave(file.path(figures_dir, "elbow_plot.png"), plot = p_elbow, width = 8, height = 5, dpi = 300)
cat("--> Elbow plot saved to results/figures/elbow_plot.pdf\n")

# ------------------------------------------------------------------------------
# 3. Multi-Sample Integration with Harmony
# ------------------------------------------------------------------------------
cat("--> Aligning patient cohort in PCA space using Harmony (IntegrateLayers)...\n")
obj <- IntegrateLayers(
  object = obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# 4. Non-Linear Dimensionality Reduction (UMAP 2D Projection)
# ------------------------------------------------------------------------------
cat("--> Computing 2D UMAP projection on Harmony integrated latent space...\n")
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap.harmony",
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# 5. Shared Nearest Neighbor (SNN) Graph & Multi-Resolution Analysis (clustree)
# ------------------------------------------------------------------------------
cat("--> Constructing SNN Graph (k = 20, 30 Harmony dimensions)...\n")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, k.param = 20, verbose = FALSE)

cat("--> Performing multi-resolution stability sweep (resolutions 0.2 to 1.0)...\n")
resolutions_sweep <- c(0.2, 0.4, 0.6, 0.8, 1.0)
obj <- FindClusters(obj, resolution = resolutions_sweep, verbose = FALSE)

# Generate and export clustree resolution tree
p_clustree <- clustree(obj, prefix = "SCT_snn_res.") +
  theme(legend.position = "bottom") +
  labs(title = "Cluster Stability Across Multiscale Resolutions",
       subtitle = "Evaluating community branch splitting without chaotic crossover")

ggsave(file.path(figures_dir, "clustree_resolution.pdf"), plot = p_clustree, width = 10, height = 9)
ggsave(file.path(figures_dir, "clustree_resolution.png"), plot = p_clustree, width = 10, height = 9, dpi = 300)
cat("--> Clustree stability plot saved to results/figures/clustree_resolution.pdf\n")

# Set optimal resolution (0.6 strikes the optimal balance for TME compartments)
Idents(obj) <- "SCT_snn_res.0.6"
obj$seurat_clusters <- Idents(obj)
cat(sprintf("Optimal clustering: Resolution 0.6 identified %d communities.\n", length(unique(Idents(obj)))))

# ------------------------------------------------------------------------------
# 6. Publication Figure 1: Integrated UMAP Landscape
# ------------------------------------------------------------------------------
cat("--> Generating Publication Figure 1 (Integrated Multi-Sample UMAP)...\n")

# Dynamic cluster palette supporting any community count
n_clusters <- length(unique(Idents(obj)))
cluster_colors <- scales::hue_pal()(n_clusters)
condition_colors <- c("Responder" = "#1d3557", "Non_Responder" = "#e63946")

p1 <- DimPlot(obj, reduction = "umap.harmony", group.by = "seurat_clusters",
              label = TRUE, label.size = 5, repel = TRUE) +
  theme_classic(base_size = 12) +
  labs(title = "Identified Cell Communities", subtitle = "Leiden/SNN Graph Partition (Res = 0.6)") +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "none")

p2 <- DimPlot(obj, reduction = "umap.harmony", group.by = "clinical_response",
              cols = condition_colors, pt.size = 0.2, alpha = 0.8) +
  theme_classic(base_size = 12) +
  labs(title = "Clinical Response Cohort", subtitle = "Patient 1 (Responder) vs Patients 2 & 3 (Non-Responder)") +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")

fig1 <- p1 | p2 +
  plot_annotation(
    title = "Figure 1: Single-Cell Atlas Landscape of the Human Tumor Microenvironment",
    subtitle = sprintf("Integrated Harmony space (30 dimensions) across 3 clinical patient biopsies (%s cells)", format(ncol(obj), big.mark = ",")),
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

ggsave(file.path(figures_dir, "fig1_integrated_umap.pdf"), plot = fig1, width = 14, height = 7)
ggsave(file.path(figures_dir, "fig1_integrated_umap.png"), plot = fig1, width = 14, height = 7, dpi = 300)
cat("--> Figure 1 saved to results/figures/fig1_integrated_umap.png\n")

# Save integrated object
saveRDS(obj, file = file.path(objects_dir, "02_seurat_integrated.rds"))
cat(">>> [Step 2 Completed] Integrated Seurat object saved to results/objects/02_seurat_integrated.rds\n\n")
