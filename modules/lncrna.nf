// ============================================================
// Módulo: Predição de lncRNAs a partir de transcritos novos
// Requer: run_stringtie=true (GTFs do GFFCOMPARE)
// ============================================================

process LNCRNA_PRED {
    label 'medium_mem'
    publishDir "${params.outdir}/lncrna", mode: 'copy'

    input:
    path(assembled_gtfs)   // GTFs por amostra do STRINGTIE (coletados)
    path(reference_gtf)
    path(genome_fasta)
    path(deseq2_all)

    output:
    path("lncrna_candidates.tsv"), emit: candidates
    path("lncrna_all.tsv"),        emit: all
    path("lncrna_summary.txt"),    emit: summary
    path("figures/"),              emit: figures

    script:
    def gtf_list = assembled_gtfs.collect { it.toString() }.join(' ')
    """
    mkdir -p figures

    # 1. Merge das GTFs por amostra
    stringtie --merge \\
        -G ${reference_gtf} \\
        -o merged.gtf \\
        ${gtf_list}

    # 2. GffCompare para anotar classe de cada transcrito
    gffcompare \\
        -r ${reference_gtf} \\
        -o gffcmp \\
        merged.gtf

    # 3. Extrai transcritos novos (class codes u=intergênico, i=intrónico,
    #    x=antisense, s=shadow) — candidatos a lncRNA
    python3 - <<'PYEOF'
import re, sys

novel_codes = {'u', 'i', 'x', 's', 'o'}
out_lines   = []
code_re     = re.compile(r'class_code "([^"]+)"')

with open("gffcmp.merged.gtf") as fh:
    for line in fh:
        m = code_re.search(line)
        if m and m.group(1) in novel_codes:
            out_lines.append(line)

# Fallback: usa merged.gtf inteiro se gffcmp.merged.gtf não existir
if not out_lines:
    try:
        with open("merged.gtf") as fh:
            for line in fh:
                if 'transcript_id' in line:
                    out_lines.append(line)
    except FileNotFoundError:
        pass

with open("novel_transcripts.gtf", "w") as fh:
    fh.writelines(out_lines)
print(f"Transcritos novos extraídos: {len([l for l in out_lines if '\ttranscript\t' in l])}")
PYEOF

    # 4. Extrai sequências FASTA dos transcritos novos
    if [ -s novel_transcripts.gtf ]; then
        gffread novel_transcripts.gtf -g ${genome_fasta} -w novel_sequences.fa
    else
        touch novel_sequences.fa
        echo "Nenhum transcrito novo encontrado."
    fi

    # 5. Classificação lncRNA em R
    mamba run -n r-analysis Rscript ${projectDir}/scripts/11_lncrna.R \\
        --fasta       novel_sequences.fa \\
        --deseq2      ${deseq2_all} \\
        --outdir      . \\
        --figures_dir figures
    """
}
