#!/usr/bin/env nextflow
// ============================================================
// RNASeq Insight Platform – Arabidopsis thaliana TAIR10
// Nextflow DSL2 | main.nf
// ============================================================

nextflow.enable.dsl = 2

include { FASTQC as FASTQC_PRE  } from './modules/qc.nf'
include { FASTQC as FASTQC_POST } from './modules/qc.nf'
include { MULTIQC as MULTIQC_PRE  } from './modules/qc.nf'
include { MULTIQC as MULTIQC_POST } from './modules/qc.nf'
include { FASTP                 } from './modules/trimming.nf'
include { GFFREAD; HISAT2_BUILD; HISAT2_ALIGN; SAMTOOLS_SORT_INDEX } from './modules/alignment.nf'
include { FEATURECOUNTS; PARSE_COUNTS } from './modules/quantification.nf'
include { RMATS; PARSE_RMATS        } from './modules/splicing.nf'
include { DESEQ2; ENRICHMENT; WGCNA; INTEGRATION; QUARTO_REPORT } from './modules/analysis.nf'

// ── Validação do samplesheet ──────────────────────────────────
def validate_samplesheet(sheet) {
    def required = ['sample','fastq_1','fastq_2','condition','replicate']
    def header = sheet.first().keySet().toList()
    required.each { col ->
        if (!header.contains(col))
            error "Coluna obrigatória ausente no samplesheet: '${col}'"
    }
}

workflow {

    // ── Lê samplesheet ───────────────────────────────────────
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, strip: true)
        .tap { rows -> validate_samplesheet(rows.toList()) }
        .map { row ->
            def meta = [
                sample:    row.sample,
                condition: row.condition,
                replicate: row.replicate
            ]
            def r1 = file(row.fastq_1, checkIfExists: true)
            def r2 = file(row.fastq_2, checkIfExists: true)
            [ meta, [ r1, r2 ] ]
        }
        .set { reads_ch }

    // ── Genoma / Anotação ────────────────────────────────────
    genome_fasta = file(params.genome_fasta, checkIfExists: true)
    genome_gtf   = params.genome_gtf  ? file(params.genome_gtf,  checkIfExists: true) : null
    genome_gff3  = params.genome_gff3 ? file(params.genome_gff3, checkIfExists: true) : null
    genome_index = params.genome_index ? file(params.genome_index) : null

    // ── Converte GFF3 → GTF se necessário ───────────────────
    if (!genome_gtf) {
        if (!genome_gff3) error "Forneça genome_gtf ou genome_gff3 no params.yaml"
        gtf_ch = GFFREAD(genome_fasta, genome_gff3)
    } else {
        gtf_ch = Channel.value(genome_gtf)
    }

    // ── Constrói índice HISAT2 se necessário ────────────────
    if (!genome_index) {
        index_ch = HISAT2_BUILD(genome_fasta, gtf_ch)
    } else {
        index_ch = Channel.value(genome_index)
    }

    // ── QC pré-trimagem ──────────────────────────────────────
    FASTQC_PRE(reads_ch)
    MULTIQC_PRE(
        FASTQC_PRE.out.zip.collect(),
        Channel.value("pre_trim")
    )

    // ── Trimagem ─────────────────────────────────────────────
    FASTP(reads_ch)
    trimmed_ch = FASTP.out.reads

    // ── QC pós-trimagem ──────────────────────────────────────
    FASTQC_POST(trimmed_ch)
    MULTIQC_POST(
        FASTQC_POST.out.zip.collect().mix(FASTP.out.json.collect()),
        Channel.value("post_trim")
    )

    // ── Alinhamento ──────────────────────────────────────────
    HISAT2_ALIGN(trimmed_ch, index_ch, gtf_ch)
    SAMTOOLS_SORT_INDEX(HISAT2_ALIGN.out.sam)
    bam_ch = SAMTOOLS_SORT_INDEX.out.bam

    // ── Quantificação ────────────────────────────────────────
    FEATURECOUNTS(bam_ch.map { m, b, i -> b }.collect(), gtf_ch, params.strandedness)
    PARSE_COUNTS(FEATURECOUNTS.out.counts, file(params.samplesheet))

    counts_ch   = PARSE_COUNTS.out.counts
    metadata_ch = PARSE_COUNTS.out.metadata

    // ── Expressão diferencial (DESeq2) ───────────────────────
    DESEQ2(counts_ch, metadata_ch)

    // ── Enriquecimento GO/KEGG/GSEA ──────────────────────────
    ENRICHMENT(DESEQ2.out.results_all, DESEQ2.out.norm_counts)

    // ── WGCNA ────────────────────────────────────────────────
    WGCNA(DESEQ2.out.norm_counts, metadata_ch)

    // ── Splicing alternativo (rMATS) ─────────────────────────
    control_bams = bam_ch
        .filter { meta, bam, bai -> meta.condition == params.control_group }
        .map    { meta, bam, bai -> bam }
        .collect()

    treatment_bams = bam_ch
        .filter { meta, bam, bai -> meta.condition == params.treatment_group }
        .map    { meta, bam, bai -> bam }
        .collect()

    RMATS(control_bams, treatment_bams, gtf_ch)
    PARSE_RMATS(RMATS.out.results_dir)
    splicing_ch = PARSE_RMATS.out.significant

    // ── Integração multi-ômica ───────────────────────────────
    INTEGRATION(
        DESEQ2.out.results_all,
        DESEQ2.out.norm_counts,
        splicing_ch.ifEmpty(file("${params.outdir}/splicing/splicing_significant.tsv")),
        ENRICHMENT.out.go_bp,
        ENRICHMENT.out.kegg,
        WGCNA.out.modules
    )

    // ── Relatório Quarto ─────────────────────────────────────
    QUARTO_REPORT(
        DESEQ2.out.results_all,
        DESEQ2.out.results_sig,
        DESEQ2.out.norm_counts,
        ENRICHMENT.out.go_bp,
        ENRICHMENT.out.kegg,
        splicing_ch.ifEmpty(file("${params.outdir}/splicing/splicing_significant.tsv")),
        WGCNA.out.modules,
        INTEGRATION.out.ranking,
        INTEGRATION.out.candidates
    )
}
