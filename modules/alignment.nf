// ============================================================
// Módulo: Alinhamento – HISAT2 + samtools
// ============================================================

process GFFREAD {
    label 'low_mem'
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path(fasta)
    path(gff3)

    output:
    path("annotation.gtf"), emit: gtf

    script:
    """
    gffread ${gff3} -T -o annotation.gtf
    """
}

process HISAT2_BUILD {
    label 'high_mem'
    publishDir "${params.outdir}/genome/hisat2_index", mode: 'copy'

    input:
    path(fasta)
    path(gtf)

    output:
    path("genome_index/"), emit: index

    script:
    """
    mkdir -p genome_index
    # Extrai splice sites e exons para alinhamento ciente de splicing
    hisat2_extract_splice_sites.py ${gtf} > splice_sites.txt
    hisat2_extract_exons.py       ${gtf} > exons.txt
    hisat2-build \\
        -p ${task.cpus} \\
        --ss splice_sites.txt \\
        --exon exons.txt \\
        ${fasta} genome_index/genome
    """
}

process HISAT2_ALIGN {
    tag "${meta.sample}"
    label 'high_mem'
    publishDir "${params.outdir}/aligned", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }

    input:
    tuple val(meta), path(reads)
    path(index_dir)
    path(gtf)

    output:
    tuple val(meta), path("${meta.sample}.sam"), emit: sam
    path("${meta.sample}_hisat2.log"),            emit: log

    script:
    def (r1, r2) = reads
    """
    hisat2 \\
        -p ${task.cpus} \\
        -x ${index_dir}/genome \\
        -1 ${r1} \\
        -2 ${r2} \\
        --dta \\
        --new-summary \\
        --summary-file ${meta.sample}_hisat2.log \\
        -S ${meta.sample}.sam

    # Valida taxa de alinhamento (gate mínimo 40%)
    ALIGN_RATE=\$(grep -i "overall alignment rate" ${meta.sample}_hisat2.log | grep -oP '[0-9.]+' | head -1)
    awk -v rate="\$ALIGN_RATE" 'BEGIN { if (rate+0 < 40) { print "ERRO: taxa de alinhamento " rate "% < 40% para ${meta.sample}"; exit 1 } }'
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "${meta.sample}"
    label 'medium_mem'
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.sample}_sorted.bam"), path("${meta.sample}_sorted.bam.bai"), emit: bam
    path("${meta.sample}_flagstat.txt"), emit: flagstat

    script:
    """
    samtools sort  -@ ${task.cpus} -o ${meta.sample}_sorted.bam ${sam}
    samtools index -@ ${task.cpus}    ${meta.sample}_sorted.bam
    samtools flagstat ${meta.sample}_sorted.bam > ${meta.sample}_flagstat.txt
    rm -f ${sam}
    """
}
