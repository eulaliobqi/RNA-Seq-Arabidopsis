#!/usr/bin/env Rscript
# ============================================================
# 02_enrichment.R – GO, KEGG e GSEA para Arabidopsis thaliana
# Usa org.At.tair.db (disponível Bioconductor 3.20+)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(org.At.tair.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ── CLI ──────────────────────────────────────────────────────
opt_list <- list(
  make_option("--deseq2",      type = "character"),
  make_option("--norm_counts", type = "character"),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1.0),
  make_option("--organism",    type = "character", default = "ath"),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))
for (p in c("deseq2","norm_counts"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Leitura ──────────────────────────────────────────────────
res <- read_tsv(opt$deseq2, show_col_types = FALSE) |>
  mutate(gene_id = gsub("\\.TAIR10$", "", gene_id))

gene_ids <- res$gene_id

sig_genes <- res |>
  filter(padj < opt$padj, abs(log2FoldChange) > opt$lfc, !is.na(padj)) |>
  pull(gene_id)

gsea_list <- res |>
  filter(!is.na(log2FoldChange), !is.na(padj)) |>
  arrange(desc(log2FoldChange)) |>
  (\(df) setNames(df$log2FoldChange, df$gene_id))()

# ── Conversão TAIR → ENTREZID (necessária para enrichGO/gseGO) ─
# clusterProfiler 4.x não mapeia keyType="TAIR" diretamente para
# termos GO em org.At.tair.db; converter para ENTREZID é obrigatório
entrez_map <- suppressMessages(suppressWarnings(
  bitr(gene_ids, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)
))
cat(sprintf("IDs mapeados: %d/%d TAIR → Entrez\n", nrow(entrez_map), length(gene_ids)))
tair2ent <- setNames(entrez_map$ENTREZID, entrez_map$TAIR)

sig_entrez <- na.omit(tair2ent[sig_genes])
cat(sprintf("DEGs com Entrez ID: %d/%d\n", length(sig_entrez), length(sig_genes)))

# Gene ranking para GSEA com Entrez IDs
gsea_idx    <- intersect(names(gsea_list), names(tair2ent))
gsea_entrez <- sort(setNames(gsea_list[gsea_idx], tair2ent[gsea_idx]), decreasing = TRUE)
gsea_entrez <- gsea_entrez[!duplicated(names(gsea_entrez))]

# ── Função de fallback para outputs vazios ────────────────────
empty_tsv <- function(path, cols = c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count")) {
  df <- as.data.frame(matrix(ncol = length(cols), nrow = 0))
  colnames(df) <- cols
  write_tsv(df, path)
  invisible(df)
}

# ── GO (BP, MF, CC) ───────────────────────────────────────────
# Usa ENTREZID: único keyType que mapeia genes → termos GO de forma
# confiável no org.At.tair.db com clusterProfiler 4.x
run_go <- function(genes_entrez, ontology, label) {
  out_path <- file.path(opt$outdir, sprintf("go_%s_results.tsv", tolower(label)))
  if (length(genes_entrez) < 5) { empty_tsv(out_path); return(invisible(NULL)) }

  tryCatch({
    ego <- enrichGO(
      gene          = genes_entrez,
      OrgDb         = org.At.tair.db,
      keyType       = "ENTREZID",
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write_tsv(as.data.frame(ego), out_path)
      p <- dotplot(ego, showCategory = 20, title = sprintf("GO %s – A. thaliana", label))
      ggsave(file.path(opt$figures_dir, sprintf("go_%s_dotplot.pdf", tolower(label))), p, width = 8, height = 7)
      ggsave(file.path(opt$figures_dir, sprintf("go_%s_dotplot.png", tolower(label))), p, width = 8, height = 7, dpi = 300)
    } else {
      empty_tsv(out_path)
    }
    return(ego)
  }, error = function(e) {
    message(sprintf("Aviso GO %s: %s", label, e$message))
    empty_tsv(out_path)
    return(NULL)
  })
}

ego_bp <- run_go(sig_entrez, "BP", "BP")
ego_mf <- run_go(sig_entrez, "MF", "MF")
ego_cc <- run_go(sig_entrez, "CC", "CC")

# Emap para BP (se ≥ 2 termos)
if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) >= 2) {
  tryCatch({
    ego_bp2 <- pairwise_termsim(ego_bp)
    p_emap  <- emapplot(ego_bp2, showCategory = 30)
    ggsave(file.path(opt$figures_dir, "go_bp_emap.pdf"), p_emap, width = 10, height = 9)
    ggsave(file.path(opt$figures_dir, "go_bp_emap.png"), p_emap, width = 10, height = 9, dpi = 300)
  }, error = function(e) message("emapplot: ", e$message))
}

# ── KEGG ──────────────────────────────────────────────────────
# KEGG Arabidopsis (ath) usa TAIR IDs no banco de dados REST.
# keyType = "kegg" instrui o clusterProfiler a não tentar converter IDs
kegg_path <- file.path(opt$outdir, "kegg_results.tsv")
if (length(sig_genes) >= 5) {
  tryCatch({
    ekegg <- enrichKEGG(
      gene          = sig_genes,
      organism      = opt$organism,
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      use_internal_data = FALSE
    )
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      write_tsv(as.data.frame(ekegg), kegg_path)
      p_kegg <- dotplot(ekegg, showCategory = 20, title = "KEGG – A. thaliana")
      ggsave(file.path(opt$figures_dir, "kegg_dotplot.pdf"), p_kegg, width = 8, height = 7)
      ggsave(file.path(opt$figures_dir, "kegg_dotplot.png"), p_kegg, width = 8, height = 7, dpi = 300)
      p_kbar <- barplot(ekegg, showCategory = 20, title = "KEGG Barplot")
      ggsave(file.path(opt$figures_dir, "kegg_barplot.pdf"), p_kbar, width = 8, height = 7)
      ggsave(file.path(opt$figures_dir, "kegg_barplot.png"), p_kbar, width = 8, height = 7, dpi = 300)
    } else {
      empty_tsv(kegg_path)
    }
  }, error = function(e) {
    message("Aviso KEGG: ", e$message)
    empty_tsv(kegg_path)
  })
} else {
  empty_tsv(kegg_path)
}

# ── GSEA – GO ─────────────────────────────────────────────────
gsea_go_path <- file.path(opt$outdir, "gsea_go_results.tsv")
if (length(gsea_entrez) >= 10) {
  tryCatch({
    gsea_go <- gseGO(
      geneList      = gsea_entrez,
      OrgDb         = org.At.tair.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      verbose       = FALSE
    )
    if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
      write_tsv(as.data.frame(gsea_go), gsea_go_path)
      p_gsea <- dotplot(gsea_go, showCategory = 20, split = ".sign", title = "GSEA GO-BP") +
        facet_grid(. ~ .sign)
      ggsave(file.path(opt$figures_dir, "gsea_go_dotplot.pdf"), p_gsea, width = 10, height = 8)
      ggsave(file.path(opt$figures_dir, "gsea_go_dotplot.png"), p_gsea, width = 10, height = 8, dpi = 300)
    } else {
      empty_tsv(gsea_go_path)
    }
  }, error = function(e) {
    message("Aviso GSEA GO: ", e$message)
    empty_tsv(gsea_go_path)
  })
} else {
  empty_tsv(gsea_go_path)
}

# ── GSEA – KEGG ───────────────────────────────────────────────
gsea_kegg_path <- file.path(opt$outdir, "gsea_kegg_results.tsv")
if (length(gsea_list) >= 10) {
  tryCatch({
    gsea_kegg <- gseKEGG(
      geneList          = gsea_list,
      organism          = opt$organism,
      keyType           = "kegg",
      pAdjustMethod     = "BH",
      pvalueCutoff      = 0.05,
      verbose           = FALSE,
      use_internal_data = FALSE
    )
    if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
      write_tsv(as.data.frame(gsea_kegg), gsea_kegg_path)
      p_gsea_k <- dotplot(gsea_kegg, showCategory = 20, split = ".sign", title = "GSEA KEGG") +
        facet_grid(. ~ .sign)
      ggsave(file.path(opt$figures_dir, "gsea_kegg_dotplot.pdf"), p_gsea_k, width = 10, height = 8)
      ggsave(file.path(opt$figures_dir, "gsea_kegg_dotplot.png"), p_gsea_k, width = 10, height = 8, dpi = 300)
    } else {
      empty_tsv(gsea_kegg_path)
    }
  }, error = function(e) {
    message("Aviso GSEA KEGG: ", e$message)
    empty_tsv(gsea_kegg_path)
  })
} else {
  empty_tsv(gsea_kegg_path)
}

cat("Enriquecimento funcional concluído.\n")
