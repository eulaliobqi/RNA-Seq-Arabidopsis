#!/usr/bin/env Rscript
# ============================================================
# 05_batch_correction.R – Detecção e correção de batch effects
# Lógica:
#   1. PCA dos counts normalizados
#   2. Se PC1 > 40% variância E cor(PC1, batch) > 0.7 → ComBat_seq
#   3. Senão → passa counts originais sem alteração
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(sva)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(dplyr)
  library(tibble)
})

opt_list <- list(
  make_option("--counts",   type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--outdir",   type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("counts", "metadata"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Leitura ──────────────────────────────────────────────────
counts_df <- read_tsv(opt$counts, show_col_types = FALSE)
meta      <- read_tsv(opt$metadata, show_col_types = FALSE) |> as.data.frame()
rownames(meta) <- meta$sample

gene_ids <- counts_df$gene_id
counts_mat <- counts_df |> select(-gene_id) |> as.matrix()
rownames(counts_mat) <- gene_ids

# Alinha ordem de amostras
common <- intersect(colnames(counts_mat), meta$sample)
counts_mat <- counts_mat[, common]
meta <- meta[common, ]

# ── PCA antes da correção ─────────────────────────────────────
plot_pca <- function(mat, meta, title, file) {
  vst_mat <- tryCatch({
    dds_tmp <- DESeqDataSetFromMatrix(round(mat), meta, ~condition)
    assay(vst(dds_tmp, blind = TRUE))
  }, error = function(e) log1p(mat))

  pca    <- prcomp(t(vst_mat), scale. = TRUE)
  pvar   <- summary(pca)$importance[2, ]
  df_pca <- data.frame(
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    condition = meta$condition,
    sample    = meta$sample,
    batch     = if ("batch" %in% colnames(meta)) meta$batch else "unknown"
  )
  p <- ggplot(df_pca, aes(PC1, PC2, color = condition, shape = batch, label = sample)) +
    geom_point(size = 3) +
    geom_text_repel(size = 2.5) +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", pvar[1] * 100),
         y = sprintf("PC2 (%.1f%%)", pvar[2] * 100)) +
    theme_bw()
  ggsave(file, p, width = 7, height = 5)
  list(pca = pca, pvar = pvar)
}

res_before <- plot_pca(counts_mat, meta,
                       "PCA – antes da correção de batch",
                       file.path(opt$outdir, "pca_before_batch.pdf"))

# ── Detecção automática de batch ──────────────────────────────
batch_detected <- FALSE
correction_msg <- "Batch não detectado. Counts originais mantidos."

if ("batch" %in% colnames(meta)) {
  pc1     <- res_before$pca$x[, 1]
  pvar_pc1 <- res_before$pvar[1]
  batch_num <- as.numeric(factor(meta$batch))

  if (length(unique(batch_num)) > 1) {
    cor_val <- abs(cor(pc1, batch_num))
    cat(sprintf("PC1 variância: %.1f%% | cor(PC1, batch): %.2f\n",
                pvar_pc1 * 100, cor_val))

    if (pvar_pc1 > 0.40 && cor_val > 0.7) {
      batch_detected <- TRUE
      correction_msg <- sprintf(
        "Batch detectado (PC1=%.1f%%, cor=%.2f). ComBat-seq aplicado.",
        pvar_pc1 * 100, cor_val
      )
      cat(correction_msg, "\n")

      # Aplica ComBat_seq
      counts_corrected <- ComBat_seq(
        counts   = counts_mat,
        batch    = meta$batch,
        group    = meta$condition,
        full_mod = TRUE
      )
      counts_mat <- counts_corrected

      # PCA após correção
      plot_pca(counts_mat, meta,
               "PCA – após ComBat-seq",
               file.path(opt$outdir, "pca_after_batch.pdf"))
    }
  }
}

# Se não houve correção, gera PCA "after" idêntico ao before para manter outputs
if (!batch_detected) {
  file.copy(
    file.path(opt$outdir, "pca_before_batch.pdf"),
    file.path(opt$outdir, "pca_after_batch.pdf"),
    overwrite = TRUE
  )
}

# ── Salva relatório ───────────────────────────────────────────
writeLines(c(
  sprintf("Amostras: %d", ncol(counts_mat)),
  sprintf("Genes:    %d", nrow(counts_mat)),
  sprintf("Coluna batch disponível: %s", "batch" %in% colnames(meta)),
  correction_msg
), file.path(opt$outdir, "batch_report.txt"))

# ── Salva counts (corrigidos ou originais) ────────────────────
out_df <- as.data.frame(counts_mat) |> rownames_to_column("gene_id")
write_tsv(out_df, file.path(opt$outdir, "counts_corrected.tsv"))

cat(correction_msg, "\n")
cat("Batch correction concluído.\n")
