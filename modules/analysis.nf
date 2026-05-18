// ============================================================
// Módulo: Análise R – DESeq2, Enriquecimento, WGCNA, Integração, Relatório
// ============================================================

process DESEQ2 {
    label 'medium_mem'
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path(counts)
    path(metadata)

    output:
    path("deseq2_results_all.tsv"),  emit: results_all
    path("deseq2_results_sig.tsv"),  emit: results_sig
    path("normalized_counts.tsv"),   emit: norm_counts
    path("deseq2_summary.txt"),      emit: summary
    path("figures/"),                emit: figures

    script:
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/01_deseq2.R \\
        --counts      ${counts} \\
        --metadata    ${metadata} \\
        --control     ${params.control_group} \\
        --treatment   ${params.treatment_group} \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process ENRICHMENT {
    label 'medium_mem'
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path(deseq2_all)
    path(norm_counts)

    output:
    path("go_bp_results.tsv"),   emit: go_bp
    path("go_mf_results.tsv"),   emit: go_mf
    path("go_cc_results.tsv"),   emit: go_cc
    path("kegg_results.tsv"),    emit: kegg
    path("gsea_go_results.tsv"), emit: gsea_go
    path("gsea_kegg_results.tsv"),emit: gsea_kegg
    path("figures/"),             emit: figures

    script:
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/02_enrichment.R \\
        --deseq2      ${deseq2_all} \\
        --norm_counts ${norm_counts} \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --organism    ${params.kegg_organism} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process WGCNA {
    label 'high_mem'
    publishDir "${params.outdir}/wgcna", mode: 'copy'

    input:
    path(norm_counts)
    path(metadata)

    output:
    path("wgcna_modules.tsv"),       emit: modules
    path("wgcna_hub_genes.tsv"),     emit: hub_genes
    path("wgcna_eigengenes.tsv"),    emit: eigengenes
    path("wgcna_module_summary.tsv"),emit: summary
    path("figures/"),                 emit: figures

    script:
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/03_wgcna.R \\
        --norm_counts ${norm_counts} \\
        --metadata    ${metadata} \\
        --min_genes   5000 \\
        --soft_power  0 \\
        --outdir      . \\
        --figures_dir figures
    """
}

process INTEGRATION {
    label 'medium_mem'
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path(deseq2_all)
    path(norm_counts)
    path(splicing)
    path(go_bp)
    path(kegg)
    path(wgcna_modules)

    output:
    path("integrated_genes.tsv"),  emit: integrated
    path("gene_ranking.tsv"),      emit: ranking
    path("key_candidates.tsv"),    emit: candidates
    path("candidates_table.tsv"),  emit: table
    path("figures/"),               emit: figures

    script:
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/04_integration.R \\
        --deseq2      ${deseq2_all} \\
        --norm_counts ${norm_counts} \\
        --splicing    ${splicing} \\
        --go_bp       ${go_bp} \\
        --kegg        ${kegg} \\
        --wgcna       ${wgcna_modules} \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process GOSEQ {
    label 'medium_mem'
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path(deseq2_all)
    path(gtf)

    output:
    path("goseq_bp_results.tsv"), emit: go_bp
    path("goseq_mf_results.tsv"), emit: go_mf
    path("goseq_cc_results.tsv"), emit: go_cc
    path("figures/"),              emit: figures

    script:
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/02b_goseq.R \\
        --deseq2      ${deseq2_all} \\
        --gtf         ${gtf} \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process METANALYSIS {
    label 'medium_mem'
    publishDir "${params.outdir}/metanalysis", mode: 'copy'

    input:
    path(deseq2_all)

    output:
    path("metanalysis_overlap.tsv"),         emit: overlap
    path("metanalysis_validated_genes.tsv"), emit: validated
    path("metanalysis_summary.tsv"),         emit: summary
    path("metanalysis_report.txt"),          emit: report
    path("figures/"),                        emit: figures

    script:
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/10_metanalysis.R \\
        --deseq2      ${deseq2_all} \\
        --geo_ids     "${params.geo_accessions}" \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process PLANTFDB {
    label 'low_mem'
    publishDir "${params.outdir}/plantfdb", mode: 'copy'

    input:
    path(deseq2_all)
    path(tf_file)       // arquivo PlantTFDB ou arquivo vazio (optional)

    output:
    path("tf_deg_classified.tsv"),   emit: classified
    path("tf_family_summary.tsv"),   emit: summary
    path("tf_family_enrichment.tsv"),emit: enrichment
    path("plantfdb_summary.txt"),    emit: report
    path("figures/"),                emit: figures

    script:
    def tf_arg = tf_file.name != 'NO_FILE' ? "--tf_file ${tf_file}" : ""
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/08_plantfdb.R \\
        --deseq2      ${deseq2_all} \\
        ${tf_arg} \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process GENIE3 {
    label 'high_mem'
    publishDir "${params.outdir}/genie3", mode: 'copy'

    input:
    path(norm_counts)
    path(deseq2_all)
    path(tf_classified)   // PLANTFDB.out.classified (pode ser vazio)
    path(plantfdb_file)   // arquivo PlantTFDB completo ou NO_FILE

    output:
    path("genie3_network.tsv"),  emit: network
    path("genie3_hub_tfs.tsv"),  emit: hubs
    path("genie3_summary.txt"),  emit: summary
    path("figures/"),            emit: figures

    script:
    def tf_cl_arg  = tf_classified.name  != 'NO_FILE' ? "--tf_classified  ${tf_classified}"  : ""
    def tf_db_arg  = plantfdb_file.name  != 'NO_FILE' ? "--plantfdb_file  ${plantfdb_file}"  : ""
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/09_genie3.R \\
        --norm_counts  ${norm_counts} \\
        --deseq2       ${deseq2_all} \\
        ${tf_cl_arg} \\
        ${tf_db_arg} \\
        --padj         ${params.padj_cutoff} \\
        --lfc          ${params.lfc_cutoff} \\
        --n_trees      ${params.genie3_trees} \\
        --n_links      ${params.genie3_links} \\
        --ncores       ${task.cpus} \\
        --outdir       . \\
        --figures_dir  figures
    """
}

process MACHINE_LEARNING {
    label 'medium_mem'
    publishDir "${params.outdir}/ml", mode: 'copy'

    input:
    path(norm_counts)
    path(deseq2_all)
    path(metadata)

    output:
    path("ml_results.tsv"),         emit: results
    path("feature_importance.tsv"), emit: features
    path("figures/"),               emit: figures

    script:
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/06_machine_learning.R \\
        --norm_counts ${norm_counts} \\
        --deseq2      ${deseq2_all} \\
        --metadata    ${metadata} \\
        --n_features  500 \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --outdir      . \\
        --figures_dir figures
    """
}

process PPI_NETWORK {
    label 'low_mem'
    publishDir "${params.outdir}/network", mode: 'copy'

    input:
    path(deseq2_all)
    path(integration)

    output:
    path("ppi_edges.tsv"),       emit: edges
    path("ppi_nodes.tsv"),       emit: nodes
    path("hub_genes.tsv"),       emit: hubs
    path("network_summary.txt"), emit: summary
    path("figures/"),            emit: figures

    script:
    """
    mkdir -p figures
    mamba run -n r-analysis Rscript ${projectDir}/scripts/07_ppi_network.R \\
        --deseq2      ${deseq2_all} \\
        --integration ${integration} \\
        --padj        ${params.padj_cutoff} \\
        --lfc         ${params.lfc_cutoff} \\
        --score       400 \\
        --outdir      . \\
        --figures_dir figures
    """
}

process QUARTO_REPORT {
    label 'medium_mem'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path(deseq2_all)
    path(deseq2_sig)
    path(norm_counts)
    path(go_bp)
    path(kegg)
    path(splicing)
    path(wgcna_modules)
    path(gene_ranking)
    path(key_candidates)

    output:
    path("rnaseq_report.html"), emit: report

    script:
    """
    # Copia dados para subpastas esperadas pelo template Quarto
    mkdir -p deseq2_data enrichment_data splicing_data wgcna_data integration_data
    cp ${deseq2_all}     deseq2_data/
    cp ${deseq2_sig}     deseq2_data/
    cp ${norm_counts}    deseq2_data/
    cp ${go_bp}          enrichment_data/
    cp ${kegg}           enrichment_data/
    cp ${splicing}       splicing_data/
    cp ${wgcna_modules}  wgcna_data/
    cp ${gene_ranking}   integration_data/
    cp ${key_candidates} integration_data/

    cp ${projectDir}/report/rnaseq_report.qmd .
    mamba run -n r-analysis quarto render rnaseq_report.qmd \\
        -P report_title:"${params.report_title}" \\
        -P report_author:"${params.report_author}" \\
        --output rnaseq_report.html
    """
}
