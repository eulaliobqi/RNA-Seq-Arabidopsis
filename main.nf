#!/usr/bin/env nextflow
// ============================================================
// RNASeq Insight Platform – Arabidopsis thaliana TAIR10
// Nextflow DSL2 | main.nf | v2.0
// ============================================================

nextflow.enable.dsl = 2

include { FASTQC_PRE; FASTQC_POST }                              from './modules/qc.nf'
include { MULTIQC as MULTIQC_PRE  }                              from './modules/qc.nf'
include { MULTIQC as MULTIQC_POST }                              from './modules/qc.nf'
include { FASTP                   }                              from './modules/trimming.nf'
include { GFFREAD;
          HISAT2_BUILD; HISAT2_ALIGN;
          STAR_INDEX;   STAR_ALIGN;
          SAMTOOLS_SORT_INDEX      }                              from './modules/alignment.nf'
include { FEATURECOUNTS; PARSE_COUNTS }                          from './modules/quantification.nf'
include { SALMON_INDEX; SALMON_QUANT; TXIMPORT }                 from './modules/salmon.nf'
include { COMBAT_SEQ              }                              from './modules/batch.nf'
include { RMATS; PARSE_RMATS; RMATS_FILTER }                     from './modules/splicing.nf'
include { STRINGTIE; GFFCOMPARE   }                              from './modules/assembly.nf'
include { DESEQ2; ENRICHMENT; GOSEQ; WGCNA; INTEGRATION;
          METANALYSIS; PLANTFDB; GENIE3; MACHINE_LEARNING; PPI_NETWORK;
          QUARTO_REPORT            }                              from './modules/analysis.nf'
include { LNCRNA_PRED             }                              from './modules/lncrna.nf'

workflow {

    // ── Lê e valida samplesheet ──────────────────────────────
    def sheet_rows = file(params.samplesheet).splitCsv(header: true, strip: true)
    def required   = ['sample','fastq_1','fastq_2','condition','replicate']
    def header     = sheet_rows[0].keySet().toList()
    required.each { col ->
        if (!header.contains(col))
            error "Coluna obrigatória ausente no samplesheet: '${col}'"
    }

    Channel
        .from(sheet_rows)
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

    // ── QC pré-trimagem ──────────────────────────────────────
    FASTQC_PRE(reads_ch)
    MULTIQC_PRE(
        FASTQC_PRE.out.zip.flatMap { meta, zips -> zips }.collect(),
        Channel.value("pre_trim")
    )

    // ── Trimagem ─────────────────────────────────────────────
    FASTP(reads_ch)
    trimmed_ch = FASTP.out.reads

    // ── QC pós-trimagem ──────────────────────────────────────
    FASTQC_POST(trimmed_ch)
    MULTIQC_POST(
        FASTQC_POST.out.zip.flatMap { meta, zips -> zips }.collect().mix(FASTP.out.json.collect()),
        Channel.value("post_trim")
    )

    // ── Alinhamento (STAR ou HISAT2 via params.aligner) ──────
    if (params.aligner == 'star') {
        if (!genome_index) {
            index_ch = STAR_INDEX(genome_fasta, gtf_ch)
        } else {
            index_ch = Channel.value(genome_index)
        }
        STAR_ALIGN(trimmed_ch, index_ch, gtf_ch)
        bam_ch = STAR_ALIGN.out.bam
    } else {
        if (!genome_index) {
            index_ch = HISAT2_BUILD(genome_fasta, gtf_ch)
        } else {
            index_ch = Channel.value(genome_index)
        }
        HISAT2_ALIGN(trimmed_ch, index_ch, gtf_ch)
        SAMTOOLS_SORT_INDEX(HISAT2_ALIGN.out.sam)
        bam_ch = SAMTOOLS_SORT_INDEX.out.bam
    }

    // ── Salmon (pseudoalinhamento paralelo) ──────────────────
    if (params.run_salmon) {
        salmon_idx_ch = SALMON_INDEX(genome_fasta, gtf_ch)
        SALMON_QUANT(trimmed_ch, salmon_idx_ch, gtf_ch)
        // Coleta os diretórios por amostra (sample1/, sample2/ ...)
        // Nextflow preserva o nome do diretório no work dir do TXIMPORT
        quant_dirs_ch = SALMON_QUANT.out.quant_dir.map { meta, dir -> dir }.collect()
        TXIMPORT(quant_dirs_ch, file(params.samplesheet))
    }

    // ── Quantificação (featureCounts) ────────────────────────
    FEATURECOUNTS(bam_ch.map { m, b, i -> b }.collect(), gtf_ch, params.strandedness)
    PARSE_COUNTS(FEATURECOUNTS.out.counts, file(params.samplesheet))

    raw_counts_ch = PARSE_COUNTS.out.counts
    metadata_ch   = PARSE_COUNTS.out.metadata

    // ── Correção de batch (ComBat-Seq) ───────────────────────
    if (params.run_combat_seq) {
        COMBAT_SEQ(raw_counts_ch, metadata_ch)
        counts_ch = COMBAT_SEQ.out.counts
    } else {
        counts_ch = raw_counts_ch
    }

    // ── Expressão diferencial (DESeq2) ───────────────────────
    DESEQ2(counts_ch, metadata_ch)

    // ── Enriquecimento GO/KEGG/GSEA ──────────────────────────
    ENRICHMENT(DESEQ2.out.results_all, DESEQ2.out.norm_counts)

    // ── GOseq (correção de viés de tamanho) ──────────────────
    GOSEQ(DESEQ2.out.results_all, gtf_ch)

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
    RMATS_FILTER(PARSE_RMATS.out.all)
    splicing_ch = PARSE_RMATS.out.significant

    // ── Montagem de transcritos – StringTie (opcional) ───────
    if (params.run_stringtie) {
        STRINGTIE(bam_ch, gtf_ch)
        assembled_gtfs_ch = STRINGTIE.out.gtf.map { meta, gtf -> gtf }.collect()
        GFFCOMPARE(assembled_gtfs_ch, gtf_ch)

        // ── lncRNA prediction ────────────────────────────────
        LNCRNA_PRED(assembled_gtfs_ch, gtf_ch, genome_fasta, DESEQ2.out.results_all)
    }

    // ── Integração multi-ômica ───────────────────────────────
    INTEGRATION(
        DESEQ2.out.results_all,
        DESEQ2.out.norm_counts,
        splicing_ch.ifEmpty(file("${params.outdir}/splicing/splicing_significant.tsv")),
        ENRICHMENT.out.go_bp,
        ENRICHMENT.out.kegg,
        WGCNA.out.modules
    )

    // ── Metanálise GEO/SRA ───────────────────────────────────
    METANALYSIS(DESEQ2.out.results_all)

    // ── PlantTFDB ─────────────────────────────────────────────
    tf_file_ch = params.plantfdb_file ?
        Channel.value(file(params.plantfdb_file)) :
        Channel.value(file("NO_FILE"))
    PLANTFDB(DESEQ2.out.results_all, tf_file_ch)

    // ── GENIE3 ────────────────────────────────────────────────
    GENIE3(
        DESEQ2.out.norm_counts,
        DESEQ2.out.results_all,
        PLANTFDB.out.classified,
        tf_file_ch
    )

    // ── Machine Learning ─────────────────────────────────────
    if (params.run_ml) {
        MACHINE_LEARNING(DESEQ2.out.norm_counts, DESEQ2.out.results_all, metadata_ch)
    }

    // ── Rede PPI (STRING) ─────────────────────────────────────
    if (params.run_ppi) {
        PPI_NETWORK(DESEQ2.out.results_all, INTEGRATION.out.ranking)
    }

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
