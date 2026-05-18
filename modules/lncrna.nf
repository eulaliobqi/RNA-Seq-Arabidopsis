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
import re, os, sys

novel_codes = {'u', 'i', 'x', 's', 'o'}
out_lines   = []
code_re     = re.compile(r'class_code "([^"]+)"')
tid_re      = re.compile(r'transcript_id "([^"]+)"')

# Tentativa 1: gffcmp.annotated.gtf (gffcompare >=0.12)
# Tentativa 2: gffcmp.merged.gtf (gffcompare antigo)
for annotated_gtf in ("gffcmp.annotated.gtf", "gffcmp.merged.gtf"):
    if os.path.exists(annotated_gtf):
        with open(annotated_gtf) as fh:
            for line in fh:
                m = code_re.search(line)
                if m and m.group(1) in novel_codes:
                    out_lines.append(line)
        if out_lines:
            print(f"Usando {annotated_gtf} com class_code", flush=True)
            break

# Tentativa 3: parse gffcmp.tracking → IDs novos → filtra merged.gtf
if not out_lines and os.path.exists("gffcmp.tracking"):
    novel_ids = set()
    with open("gffcmp.tracking") as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4 and parts[3] in novel_codes:
                novel_ids.add(parts[0])
    print(f"Tracking: {len(novel_ids)} transcritos novos encontrados", flush=True)
    if novel_ids and os.path.exists("merged.gtf"):
        with open("merged.gtf") as fh:
            for line in fh:
                m = tid_re.search(line)
                if m and m.group(1) in novel_ids:
                    out_lines.append(line)

# Tentativa 4: usa merged.gtf inteiro como fallback final
if not out_lines and os.path.exists("merged.gtf"):
    print("AVISO: sem class_code disponível; usando merged.gtf completo como candidatos", flush=True)
    with open("merged.gtf") as fh:
        for line in fh:
            if 'transcript_id' in line:
                out_lines.append(line)

with open("novel_transcripts.gtf", "w") as fh:
    fh.writelines(out_lines)
n_tx = len([l for l in out_lines if '\ttranscript\t' in l])
print(f"Transcritos novos extraídos: {n_tx}", flush=True)
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
