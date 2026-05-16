// ============================================================
// Módulo: Trimagem – fastp (paired-end)
// ============================================================

process FASTP {
    tag "${meta.sample}"
    label 'low_mem'
    publishDir "${params.outdir}/trimmed", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.json') || fn.endsWith('.html') ? "reports/${fn}" : fn }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample}_R{1,2}_trimmed.fastq.gz"), emit: reads
    path("${meta.sample}_fastp.json"),                                emit: json
    path("${meta.sample}_fastp.html"),                                emit: html

    script:
    def (r1, r2) = reads
    """
    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --out1 ${meta.sample}_R1_trimmed.fastq.gz \\
        --out2 ${meta.sample}_R2_trimmed.fastq.gz \\
        --json ${meta.sample}_fastp.json \\
        --html ${meta.sample}_fastp.html \\
        --length_required ${params.min_length} \\
        --qualified_quality_phred ${params.quality} \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --trim_poly_x \\
        --overrepresentation_analysis \\
        --thread ${task.cpus} \\
        --report_title "${meta.sample}"
    """
}
