// ============================================================
// Módulo: Controle de Qualidade – FastQC e MultiQC
// ============================================================

process FASTQC_PRE {
    tag "${meta.sample}"
    label 'low_mem'
    publishDir "${params.outdir}/qc/pre_trim", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    def (r1, r2) = reads
    """
    fastqc --threads ${task.cpus} --outdir . ${r1} ${r2}
    """
}

process FASTQC_POST {
    tag "${meta.sample}"
    label 'low_mem'
    publishDir "${params.outdir}/qc/post_trim", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    def (r1, r2) = reads
    """
    fastqc --threads ${task.cpus} --outdir . ${r1} ${r2}
    """
}

process MULTIQC {
    label 'low_mem'
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path(reports)
    val(label)

    output:
    path("multiqc_${label}_report.html"),      emit: report
    path("multiqc_${label}_report_data/"),     emit: data

    script:
    """
    multiqc . \
        --filename multiqc_${label}_report \
        --title "MultiQC – ${label}" \
        --force
    """
}
