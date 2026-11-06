/*
 * STARsolo Process: Alignment, Cell Barcode demultiplexing, and UMI deduplication
 */
process STARSOLO {
    tag { meta.id }
    label 'process_high'

    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_0'

    publishDir path: { "${params.outdir}/starsolo/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path  genome_index
    path  whitelist

    output:
    tuple val(meta), path("Solo.out/Gene/filtered"), emit: counts_filtered
    tuple val(meta), path("Solo.out/Gene/raw"),      emit: counts_raw
    path  "Solo.out/Gene/Summary.csv",               emit: summary
    path  "versions.yml",                            emit: versions

    script:
    def (r1, r2) = reads
    def umi_len  = meta.chemistry == '10xv2' ? 10 : 12
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${genome_index} \\
        --readFilesIn ${r2} ${r1} \\
        --readFilesCommand zcat \\
        --soloType CB_UMI_Simple \\
        --soloCBwhitelist ${whitelist} \\
        --soloCBstart 1 --soloCBlen 16 \\
        --soloUMIstart 17 --soloUMIlen ${umi_len} \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloUMIfiltering MultiGeneUMI_CR \\
        --soloUMIdedup 1MM_CR \\
        --soloCellFilter EmptyDrops_CR \\
        --outSAMtype None \\
        --outFileNamePrefix ${meta.id}_

    mkdir -p Solo.out/Gene
    mv ${meta.id}_Solo.out/Gene/* Solo.out/Gene/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
