#!/usr/bin/env Rscript
# ============================================================
# 02b_goseq.R – Enriquecimento GO com correção de viés de tamanho
# Usa goseq (Young et al. 2010) com bias.data = comprimentos de genes
# Corre em paralelo ao 02_enrichment.R (não o substitui).
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(goseq)
  library(GenomicFeatures)
  library(org.At.tair.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
})

# txdbmaker é necessário em GenomicFeatures >= 1.61 (makeTxDbFromGFF foi movido)
make_txdb <- function(gtf) {
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdbmaker::makeTxDbFromGFF(gtf, format = "GTF")
  } else {
    GenomicFeatures::makeTxDbFromGFF(gtf, format = "GTF")
  }
}

opt_list <- list(
  make_option("--deseq2",      type = "character"),
  make_option("--gtf",         type = "character"),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1.0),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("deseq2", "gtf"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Leitura DESeq2 ────────────────────────────────────────────
res <- read_tsv(opt$deseq2, show_col_types = FALSE) |>
  mutate(gene_id = gsub("\\.TAIR10$", "", gene_id)) |>
  filter(!is.na(padj))

all_genes <- res$gene_id

# ── Vetor binário DEG ─────────────────────────────────────────
is_deg <- as.integer(
  res$padj < opt$padj & abs(res$log2FoldChange) > opt$lfc
)
names(is_deg) <- all_genes

cat(sprintf("Total genes testados: %d | DEGs: %d\n", length(is_deg), sum(is_deg)))

# ── Comprimentos de genes a partir do GTF ────────────────────
cat("Extraindo comprimentos de genes do GTF...\n")
txdb <- tryCatch(
  make_txdb(opt$gtf),
  error = function(e) {
    message("Aviso: makeTxDbFromGFF falhou, usando comprimentos uniformes. ", e$message)
    NULL
  }
)

if (!is.null(txdb)) {
  exons_by_gene <- exonsBy(txdb, by = "gene")
  gene_lengths   <- sapply(exons_by_gene, function(ex) {
    sum(width(reduce(ex)))
  })
  # Remove sufixo .TAIR10 dos nomes
  names(gene_lengths) <- gsub("\\.TAIR10$", "", names(gene_lengths))
  # Alinha com vetor de genes
  gene_lengths <- gene_lengths[names(is_deg)]
  gene_lengths[is.na(gene_lengths)] <- median(gene_lengths, na.rm = TRUE)
} else {
  # Fallback: comprimentos uniformes (sem correção efetiva)
  gene_lengths <- rep(1000L, length(is_deg))
  names(gene_lengths) <- names(is_deg)
}

# ── GOseq: viés de probabilidade ─────────────────────────────
# Se comprimentos uniformes, nullp não consegue ajustar spline (< 6 valores únicos)
# → cria pwf manual e usa Hypergeometric (equivalente a Fisher sem correção de viés)
n_unique_lengths <- length(unique(gene_lengths))
if (n_unique_lengths >= 6) {
  pwf          <- nullp(is_deg, bias.data = gene_lengths, plot.fit = FALSE)
  goseq_method <- "Wallenius"
  cat("PWF ajustado com comprimentos reais de genes (método Wallenius).\n")
} else {
  cat(sprintf("Aviso: apenas %d comprimentos únicos — usando método Hypergeometric (sem correção de viés).\n",
              n_unique_lengths))
  pwf <- data.frame(
    DEgenes   = is_deg,
    bias.data = gene_lengths,
    pwf       = sum(is_deg) / length(is_deg),
    stringsAsFactors = FALSE
  )
  rownames(pwf) <- names(is_deg)
  goseq_method <- "Hypergeometric"
}

# ── Anotações GO de org.At.tair.db (gene2cat customizado) ────
# genome="tair10" / id="tair" não está no banco interno do goseq ≥1.58;
# usa AnnotationDbi para construir o mapeamento gene→GO por ontologia.
get_gene2go <- function(gene_ids, ontology) {
  tryCatch({
    go_df <- AnnotationDbi::select(
      org.At.tair.db,
      keys    = gene_ids,
      columns = c("GO", "ONTOLOGY"),
      keytype = "TAIR"
    )
    go_df <- go_df[!is.na(go_df$GO) & go_df$ONTOLOGY == ontology, , drop = FALSE]
    if (nrow(go_df) == 0) return(list())
    split(go_df$GO, go_df$TAIR)
  }, error = function(e) {
    message("Erro ao carregar anotações GO (", ontology, "): ", e$message)
    list()
  })
}

run_goseq <- function(pwf, ontology, method = "Wallenius") {
  gene2cat <- get_gene2go(rownames(pwf), ontology)
  if (length(gene2cat) == 0) {
    message(sprintf("Nenhuma anotação GO:%s em org.At.tair.db.", ontology))
    return(data.frame())
  }
  tryCatch({
    res_go <- goseq(pwf, gene2cat = gene2cat, method = method)
    res_go |>
      filter(numDEInCat > 0) |>
      mutate(p.adjust = p.adjust(over_represented_pvalue, method = "BH")) |>
      arrange(p.adjust)
  }, error = function(e) {
    message(sprintf("GOseq %s falhou: %s", ontology, e$message))
    data.frame()
  })
}

go_bp <- run_goseq(pwf, "BP", method = goseq_method)
go_mf <- run_goseq(pwf, "MF", method = goseq_method)
go_cc <- run_goseq(pwf, "CC", method = goseq_method)

# ── Salva resultados ──────────────────────────────────────────
write_tsv(go_bp, file.path(opt$outdir, "goseq_bp_results.tsv"))
write_tsv(go_mf, file.path(opt$outdir, "goseq_mf_results.tsv"))
write_tsv(go_cc, file.path(opt$outdir, "goseq_cc_results.tsv"))

# ── Gráfico lollipop top 20 GO-BP ────────────────────────────
if (nrow(go_bp) > 0) {
  top20 <- head(go_bp[go_bp$p.adjust < 0.05, ], 20)
  if (nrow(top20) > 0) {
    p <- ggplot(top20, aes(x = -log10(p.adjust),
                           y = reorder(term, -log10(p.adjust)),
                           size = numDEInCat)) +
      geom_segment(aes(xend = 0, yend = reorder(term, -log10(p.adjust))),
                   color = "grey70") +
      geom_point(color = "#377EB8") +
      labs(x = "-log10(FDR)", y = NULL,
           title = "GOseq GO-BP (com correção de viés de tamanho)",
           size = "DEGs") +
      theme_bw()
    ggsave(file.path(opt$figures_dir, "goseq_bp_top20.pdf"), p, width = 9, height = 6)
    ggsave(file.path(opt$figures_dir, "goseq_bp_top20.png"), p, width = 9, height = 6, dpi = 300)
  }
}

cat(sprintf("GOseq concluído: %d termos BP | %d MF | %d CC\n",
            nrow(go_bp), nrow(go_mf), nrow(go_cc)))
