// ============================================================
// Módulo: Alinhamento – STAR (principal) + HISAT2 (legado)
// Selecionar via params.aligner = 'star' | 'hisat2'
// ============================================================

// ── STAR ──────────────────────────────────────────────────────

process STAR_INDEX {
    label 'high_mem'
    publishDir "${params.outdir}/genome/star_index", mode: 'copy'

    input:
    path(fasta)
    path(gtf)

    output:
    path("star_index/"), emit: index

    script:
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus} \\
        --genomeChrBinNbits 16
    """
}

process STAR_ALIGN {
    tag "${meta.sample}"
    label 'high_mem'
    publishDir "${params.outdir}/aligned_star", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.Log.final.out') ? "logs/${fn}" : fn }

    input:
    tuple val(meta), path(reads)
    path(index_dir)
    path(gtf)

    output:
    tuple val(meta), path("${meta.sample}_Aligned.sortedByCoord.out.bam"),
                     path("${meta.sample}_Aligned.sortedByCoord.out.bam.bai"), emit: bam
    path("${meta.sample}_Log.final.out"), emit: log

    script:
    def (r1, r2) = reads
    """
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index_dir} \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${r1} ${r2} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD \\
        --outFileNamePrefix ${meta.sample}_ \\
        --outSAMstrandField intronMotif \\
        --outFilterIntronMotifs RemoveNoncanonical \\
        --alignSoftClipAtReferenceEnds Yes \\
        --quantMode TranscriptomeSAM GeneCounts \\
        --limitBAMsortRAM 10000000000

    samtools index -@ ${task.cpus} ${meta.sample}_Aligned.sortedByCoord.out.bam

    # Valida taxa de alinhamento (gate mínimo 40%)
    ALIGN_RATE=\$(grep "Uniquely mapped reads %" ${meta.sample}_Log.final.out | grep -oP '[0-9.]+' | head -1)
    awk -v rate="\$ALIGN_RATE" 'BEGIN { if (rate+0 < 40) { print "ERRO: taxa de alinhamento " rate "% < 40% para ${meta.sample}"; exit 1 } }'
    """
}

// ── HISAT2 (mantido para transição/comparação) ─────────────────

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
    gffread ${gff3} -T -o annotation_raw.gtf
    # gffread adds 'chr' prefix to Chr* names from TAIR10 GFF3 (Chr5 → chrChr5),
    # which breaks rMATS matching against BAMs built from the FASTA (Chr5).
    sed 's/^chrChr/Chr/' annotation_raw.gtf > annotation.gtf
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
