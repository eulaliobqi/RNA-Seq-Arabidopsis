// ============================================================
// Módulo: Quantificação por pseudoalinhamento – Salmon
// ============================================================

process SALMON_INDEX {
    label 'medium_mem'
    publishDir "${params.outdir}/genome/salmon_index", mode: 'copy'

    input:
    path(fasta)
    path(gtf)

    output:
    path("salmon_index/"), emit: index

    script:
    """
    # Extrai transcritos do FASTA + GTF
    gffread ${gtf} -g ${fasta} -w transcripts.fa

    salmon index \\
        --transcripts transcripts.fa \\
        --index salmon_index \\
        --threads ${task.cpus} \\
        --gencode
    """
}

process SALMON_QUANT {
    tag "${meta.sample}"
    label 'medium_mem'
    publishDir "${params.outdir}/counts/salmon/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path(index_dir)
    path(gtf)

    output:
    tuple val(meta), path("${meta.sample}/quant.sf"), emit: quant
    path("${meta.sample}/aux_info/"),                 emit: aux

    script:
    def (r1, r2) = reads
    def lib_type = params.strandedness == 1 ? 'SF' :
                   params.strandedness == 2 ? 'SR' : 'A'
    """
    salmon quant \\
        --index ${index_dir} \\
        --libType ${lib_type} \\
        --mates1 ${r1} \\
        --mates2 ${r2} \\
        --threads ${task.cpus} \\
        --validateMappings \\
        --gcBias \\
        --seqBias \\
        --output ${meta.sample}
    """
}

process TXIMPORT {
    label 'low_mem'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(quant_files)
    path(samplesheet)

    output:
    path("salmon_counts.tsv"), emit: counts
    path("salmon_tpm.tsv"),    emit: tpm

    script:
    """
    mamba run -n r-analysis Rscript ${projectDir}/scripts/00_tximport.R \\
        --quant_dir . \\
        --samplesheet ${samplesheet} \\
        --outdir .
    """
}
