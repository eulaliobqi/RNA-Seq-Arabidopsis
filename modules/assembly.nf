// ============================================================
// Módulo: Montagem de transcritos – StringTie + GffCompare
// Identifica transcritos novos e isoformas não anotadas.
// ============================================================

process STRINGTIE {
    tag "${meta.sample}"
    label 'medium_mem'
    publishDir "${params.outdir}/stringtie", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(gtf)

    output:
    tuple val(meta), path("${meta.sample}.gtf"), emit: gtf
    path("${meta.sample}_abundance.tsv"),         emit: abundance

    script:
    def strand_opt = params.strandedness == 1 ? '--fr' :
                     params.strandedness == 2 ? '--rf' : ''
    """
    stringtie ${bam} \\
        -G ${gtf} \\
        -o ${meta.sample}.gtf \\
        -A ${meta.sample}_abundance.tsv \\
        -p ${task.cpus} \\
        ${strand_opt} \\
        -v
    """
}

process GFFCOMPARE {
    label 'medium_mem'
    publishDir "${params.outdir}/stringtie", mode: 'copy'

    input:
    path(assembled_gtfs)
    path(reference_gtf)

    output:
    path("gffcmp.annotated.gtf"), emit: annotated_gtf
    path("gffcmp.stats"),          emit: stats
    path("gffcmp.loci"),           emit: loci
    path("gffcmp.tracking"),       emit: tracking

    script:
    def gtf_list = assembled_gtfs.collect { it.toString() }.join(' ')
    """
    # Merge de todos os GTFs por amostra
    stringtie --merge \\
        -G ${reference_gtf} \\
        -o merged.gtf \\
        ${gtf_list}

    # Comparação com anotação de referência
    gffcompare \\
        -r ${reference_gtf} \\
        -o gffcmp \\
        merged.gtf

    cp merged.gtf gffcmp.annotated.gtf
    """
}
