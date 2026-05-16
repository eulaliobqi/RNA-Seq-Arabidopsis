// ============================================================
// Módulo: Splicing alternativo – rMATS
// ============================================================

process RMATS {
    label 'high_mem'
    publishDir "${params.outdir}/splicing/rmats_output", mode: 'copy'

    input:
    path(control_bams)
    path(treatment_bams)
    path(gtf)

    output:
    path("rmats_out/"), emit: results_dir

    script:
    def ctrl_list = control_bams.collect { it.toString() }.join(',')
    def trt_list  = treatment_bams.collect { it.toString() }.join(',')
    """
    # Gera arquivos de lista para rMATS
    echo "${ctrl_list}" | tr ',' '\\n' > b1.txt
    echo "${trt_list}"  | tr ',' '\\n' > b2.txt

    python \$(which rmats.py) \\
        --b1 b1.txt \\
        --b2 b2.txt \\
        --gtf ${gtf} \\
        --od rmats_out \\
        --tmp rmats_tmp \\
        -t paired \\
        --readLength ${params.read_length} \\
        --nthread ${task.cpus} \\
        --statoff
    """
}

process PARSE_RMATS {
    label 'low_mem'
    publishDir "${params.outdir}/splicing", mode: 'copy'

    input:
    path(rmats_dir)

    output:
    path("splicing_significant.tsv"), emit: significant
    path("splicing_all.tsv"),         emit: all

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']
    frames = []

    for et in event_types:
        fname = os.path.join("${rmats_dir}", f"{et}.MATS.JC.txt")
        if not os.path.exists(fname):
            fname = os.path.join("${rmats_dir}", f"{et}.MATS.JCEC.txt")
        if not os.path.exists(fname):
            continue
        try:
            df = pd.read_csv(fname, sep="\\t")
            if df.empty:
                continue
            df['event_type'] = et
            frames.append(df)
        except Exception:
            continue

    if not frames:
        pd.DataFrame().to_csv("splicing_all.tsv",         sep="\\t", index=False)
        pd.DataFrame().to_csv("splicing_significant.tsv", sep="\\t", index=False)
        print("Nenhum evento de splicing encontrado.")
        exit(0)

    all_df = pd.concat(frames, ignore_index=True)

    # Colunas de PSI variam por tipo de evento
    psi_cols = [c for c in all_df.columns if 'IncLevel' in c and 'Difference' in c]
    if not psi_cols:
        psi_cols = ['IncLevelDifference'] if 'IncLevelDifference' in all_df.columns else []

    # Filtra: FDR < 0.05 e |ΔPSI| > 0.1
    sig = all_df.copy()
    if 'FDR' in sig.columns:
        sig = sig[sig['FDR'] < 0.05]
    if psi_cols:
        psi_col = psi_cols[0]
        sig = sig[sig[psi_col].abs() > 0.10]

    all_df.to_csv("splicing_all.tsv",         sep="\\t", index=False)
    sig.to_csv(   "splicing_significant.tsv", sep="\\t", index=False)
    print(f"Eventos totais: {len(all_df)} | Significativos: {len(sig)}")
    """
}
