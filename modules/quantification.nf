// ============================================================
// Módulo: Quantificação – featureCounts
// ============================================================

process FEATURECOUNTS {
    label 'medium_mem'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(bams)
    path(gtf)
    val(strandedness)

    output:
    path("counts_matrix.txt"),         emit: counts
    path("counts_matrix.txt.summary"), emit: summary

    script:
    def bam_list  = bams.collect { it.toString() }.join(' ')
    def strand_opt = strandedness == -1 ? 0 : strandedness
    """
    featureCounts \\
        -T ${task.cpus} \\
        -t ${params.feature_type} \\
        -g ${params.gene_attr} \\
        -s ${strand_opt} \\
        -p \\
        --countReadPairs \\
        -B \\
        -C \\
        -a ${gtf} \\
        -o counts_matrix.txt \\
        ${bam_list}
    """
}

process PARSE_COUNTS {
    label 'low_mem'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(raw_counts)
    path(samplesheet)

    output:
    path("counts_clean.tsv"),    emit: counts
    path("sample_metadata.tsv"), emit: metadata

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    df = pd.read_csv("${raw_counts}", sep="\\t", comment="#", index_col=0)

    bam_cols = [c for c in df.columns if c.endswith('.bam') or '_sorted' in c]
    df_counts = df[bam_cols].copy()

    meta = pd.read_csv("${samplesheet}")

    rename_map = {}
    for _, row in meta.iterrows():
        for col in bam_cols:
            if row['sample'] in col:
                rename_map[col] = row['sample']
    df_counts.rename(columns=rename_map, inplace=True)

    ordered = [s for s in meta['sample'].tolist() if s in df_counts.columns]
    df_counts = df_counts[ordered]
    df_counts.index = df_counts.index.str.replace(r'\\.TAIR10$', '', regex=True)
    df_counts.index.name = "gene_id"
    df_counts.to_csv("counts_clean.tsv", sep="\\t")

    meta[['sample','condition','replicate']].to_csv("sample_metadata.tsv", sep="\\t", index=False)
    print(f"Matriz: {df_counts.shape[0]} genes x {df_counts.shape[1]} amostras")
    print(f"Total reads mapeados: {df_counts.values.sum():,}")
    """
}
