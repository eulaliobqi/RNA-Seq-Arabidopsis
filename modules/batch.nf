// ============================================================
// Módulo: Correção de batch effects – ComBat-Seq (SVA)
// Inserido entre PARSE_COUNTS e DESEQ2.
// Aplica correção apenas se PC1 > 40% variância E correlação
// com batch > 0.7; caso contrário passa counts inalterados.
// ============================================================

process COMBAT_SEQ {
    label 'medium_mem'
    publishDir "${params.outdir}/batch_correction", mode: 'copy'

    input:
    path(counts)
    path(metadata)

    output:
    path("counts_corrected.tsv"),   emit: counts
    path("pca_before_batch.pdf"),   emit: pca_before
    path("pca_after_batch.pdf"),    emit: pca_after
    path("batch_report.txt"),       emit: report

    script:
    """
    mamba run -n r-analysis Rscript ${projectDir}/scripts/05_batch_correction.R \\
        --counts   ${counts} \\
        --metadata ${metadata} \\
        --outdir   .
    """
}
