/*
 * ====================================================================
 *  Module: SEURAT_ANALYSIS (Downstream Systems Biology & Cell Profiling)
 * ====================================================================
 */

process RUN_QC_DOUBLETS {
    tag "QC_Doublets"
    label 'process_medium'
    container 'rocker/r-ver:4.4.2'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts_dir
    path script_path

    output:
    path "results/objects/01_seurat_qc_filtered.rds", emit: qc_rds
    path "results/figures/qc_violins.*",             emit: qc_figures
    path "results/tables/qc_summary_metrics.csv",    emit: qc_table

    script:
    """
    Rscript ${script_path}
    """
}

process RUN_INTEGRATION_CLUSTERING {
    tag "SCT_Harmony_Leiden"
    label 'process_high'
    container 'rocker/r-ver:4.4.2'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path qc_rds
    path script_path

    output:
    path "results/objects/02_seurat_integrated.rds", emit: integrated_rds
    path "results/figures/fig1_integrated_umap.*",   emit: fig1
    path "results/figures/elbow_plot.*",             emit: elbow_fig
    path "results/figures/clustree_resolution.*",    emit: clustree_fig

    script:
    """
    Rscript ${script_path}
    """
}

process RUN_ANNOTATION {
    tag "Marker_Discovery_Annotation"
    label 'process_medium'
    container 'rocker/r-ver:4.4.2'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path integrated_rds
    path script_path

    output:
    path "results/objects/03_seurat_annotated.rds",                 emit: annotated_rds
    path "results/figures/fig2_markers_dotplot.*",                  emit: fig2
    path "results/tables/markers_wilcoxon_all.csv",                 emit: markers_csv
    path "results/tables/confusion_matrix_singler_vs_manual.csv",   emit: confusion_csv

    script:
    """
    Rscript ${script_path}
    """
}

process RUN_SYSTEMS_BIOLOGY {
    tag "Monocle3_CellChat"
    label 'process_medium'
    container 'rocker/r-ver:4.4.2'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path annotated_rds
    path script_path

    output:
    path "results/figures/fig3_tcell_exhaustion_trajectory.*", emit: fig3
    path "results/figures/fig4_cellchat_signaling_networks.*", emit: fig4
    path "results/tables/cellchat_network_centrality.csv",     emit: centrality_csv
    path "results/tables/systems_biology_summary.csv",         emit: summary_csv

    script:
    """
    Rscript ${script_path}
    """
}
