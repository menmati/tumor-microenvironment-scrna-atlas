#!/usr/bin/env nextflow
/*
 * ====================================================================
 *  Pipeline: tumor-microenvironment-scrna-atlas
 *  Deciphering Tumor Microenvironment Heterogeneity via scRNA-seq
 *  Authors: Carmen Macías Portillo (@menmati)
 * ====================================================================
 */
nextflow.enable.dsl = 2

include { STARSOLO } from './modules/starsolo'
include { 
    RUN_QC_DOUBLETS; 
    RUN_INTEGRATION_CLUSTERING; 
    RUN_ANNOTATION; 
    RUN_SYSTEMS_BIOLOGY 
} from './modules/seurat_analysis'

workflow {

    if (params.help) {
        log.info """
        ========================================================================
         T U M O R   M I C R O E N V I R O N M E N T   s c R N A - S E Q   A T L A S
        ========================================================================
         Usage:
           nextflow run main.nf [options]

         Input / Output Options:
           --counts_dir      Directory containing 10x MEX matrix cohorts
                             [default: ${params.counts_dir}]
           --samplesheet     Sample metadata CSV for upstream raw reads
                             [default: ${params.samplesheet}]
           --outdir          Directory to write results
                             [default: ${params.outdir}]

         Execution Toggles:
           --run_upstream    Execute STARsolo FASTQ alignment [default: ${params.run_upstream}]
           --run_downstream  Execute Seurat v5 / Harmony / Systems Biology [default: ${params.run_downstream}]
           --help            Show this help message and exit

         Profiles:
           -profile docker   Containerized execution via Docker Desktop / OCI
           -profile local    Local execution in R / environment
           -profile awsbatch Scalable cloud deployment via AWS Batch
        ========================================================================
        """.stripIndent()
        exit 0
    }

    log.info """
    ========================================================================
     T U M O R   M I C R O E N V I R O N M E N T   s c R N A - S E Q   A T L A S
    ========================================================================
     Samplesheet       : ${params.samplesheet}
     Counts Directory  : ${params.counts_dir}
     Output Directory  : ${params.outdir}
     Run Upstream      : ${params.run_upstream}
     Run Downstream    : ${params.run_downstream}
    ========================================================================
    """.stripIndent()

    // -----------------------------------------------------------------
    // Phase 1: Upstream Raw FASTQ Alignment & Quantification (STARsolo)
    // -----------------------------------------------------------------
    if (params.run_upstream) {
        Channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample, condition: row.condition, chemistry: row.chemistry]
                [meta, [file(row.fastq_1), file(row.fastq_2)]]
            }
            .set { sc_reads_ch }

        def whitelist_file = file(params.whitelist_v3)

        STARSOLO(
            sc_reads_ch,
            file(params.genome_index),
            whitelist_file
        )

        STARSOLO.out.counts_filtered.view { meta, path ->
            "[STARsolo Complete] Count matrix for ${meta.id} generated at: ${path}"
        }
    }

    // -----------------------------------------------------------------
    // Phase 2: Downstream Single-Cell Analytics & Systems Biology
    // -----------------------------------------------------------------
    if (params.run_downstream) {
        // Step 1: Quality Control & Doublet Detection
        def qc_script = file("${projectDir}/scripts/01_qc_and_doublets.R")
        RUN_QC_DOUBLETS(file(params.counts_dir), qc_script)

        // Step 2: Normalization, Integration & Community Detection
        def integ_script = file("${projectDir}/scripts/02_integrate_cluster.R")
        RUN_INTEGRATION_CLUSTERING(RUN_QC_DOUBLETS.out.qc_rds, integ_script)

        // Step 3: Marker Discovery & Phenotypic Annotation
        def annot_script = file("${projectDir}/scripts/03_annotation.R")
        RUN_ANNOTATION(RUN_INTEGRATION_CLUSTERING.out.integrated_rds, annot_script)

        // Step 4: Systems Biology — Monocle3 Trajectory & CellChat v2 Networks
        def sysbio_script = file("${projectDir}/scripts/04_systems_biology.R")
        RUN_SYSTEMS_BIOLOGY(RUN_ANNOTATION.out.annotated_rds, sysbio_script)
    }
}
