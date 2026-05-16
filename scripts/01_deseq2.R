#!/usr/bin/env Rscript
# ============================================================
# 01_deseq2.R – Expressão diferencial com DESeq2
# Uso: Rscript 01_deseq2.R --counts counts_clean.tsv \
#                           --metadata sample_metadata.tsv \
#                           --control WT --treatment mutant
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ── CLI ──────────────────────────────────────────────────────
opt_list <- list(
  make_option("--counts",      type = "character"),
  make_option("--metadata",    type = "character"),
  make_option("--control",     type = "character"),
  make_option("--treatment",   type = "character"),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1.0),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("counts","metadata","control","treatment"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório ausente: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Leitura ──────────────────────────────────────────────────
counts <- read_tsv(opt$counts,   show_col_types = FALSE) |>
          column_to_rownames("gene_id") |>
          as.matrix()
meta   <- read_tsv(opt$metadata, show_col_types = FALSE) |>
          as.data.frame()
rownames(meta) <- meta$sample

# Garante mesma ordem
common <- intersect(colnames(counts), meta$sample)
counts <- counts[, common]
meta   <- meta[common, ]
meta$condition <- factor(meta$condition, levels = c(opt$control, opt$treatment))

# Filtro: remove genes com < 10 reads totais
counts <- counts[rowSums(counts) >= 10, ]
cat(sprintf("Genes após filtro: %d\n", nrow(counts)))

# ── DESeq2 ───────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(counts, meta, ~condition)
dds <- DESeq(dds)

# LFC shrinkage com fallback automático
coef_name <- resultsNames(dds)[2]
res_shrunk <- tryCatch(
  lfcShrink(dds, coef = coef_name, type = "apeglm"),
  error = function(e) tryCatch(
    lfcShrink(dds, coef = coef_name, type = "ashr"),
    error = function(e2) results(dds)
  )
)

res_df <- as.data.frame(res_shrunk) |>
  tibble::rownames_to_column("gene_id") |>
  arrange(padj)

res_df$regulation <- case_when(
  res_df$padj < opt$padj & res_df$log2FoldChange >  opt$lfc ~ "up",
  res_df$padj < opt$padj & res_df$log2FoldChange < -opt$lfc ~ "down",
  TRUE ~ "ns"
)

# ── Contagens normalizadas (VST) ──────────────────────────────
vst <- vst(dds, blind = FALSE)
norm_counts <- as.data.frame(assay(vst)) |>
  tibble::rownames_to_column("gene_id")

# ── Outputs ──────────────────────────────────────────────────
write_tsv(res_df,                   file.path(opt$outdir, "deseq2_results_all.tsv"))
write_tsv(filter(res_df, padj < opt$padj), file.path(opt$outdir, "deseq2_results_sig.tsv"))
write_tsv(norm_counts,              file.path(opt$outdir, "normalized_counts.tsv"))

n_up   <- sum(res_df$regulation == "up",   na.rm = TRUE)
n_down <- sum(res_df$regulation == "down", na.rm = TRUE)
writeLines(c(
  sprintf("Contraste: %s vs %s", opt$treatment, opt$control),
  sprintf("Genes testados:    %d", nrow(res_df)),
  sprintf("DEGs (padj<%.2f, |lfc|>%.1f):", opt$padj, opt$lfc),
  sprintf("  Up-regulated:   %d", n_up),
  sprintf("  Down-regulated: %d", n_down),
  sprintf("  Total:          %d", n_up + n_down)
), file.path(opt$outdir, "deseq2_summary.txt"))
cat(readLines(file.path(opt$outdir, "deseq2_summary.txt")), sep = "\n")

# ── PCA ──────────────────────────────────────────────────────
pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
pvar     <- attr(pca_data, "percentVar")
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  labs(
    title = "PCA – Amostras",
    x = sprintf("PC1 (%.1f%%)", pvar[1] * 100),
    y = sprintf("PC2 (%.1f%%)", pvar[2] * 100)
  ) +
  theme_bw()
ggsave(file.path(opt$figures_dir, "pca_samples.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "pca_samples.png"), p_pca, width = 7, height = 5, dpi = 300)

# ── Volcano ───────────────────────────────────────────────────
colors <- c("up" = "#E41A1C", "down" = "#377EB8", "ns" = "grey70")
p_vol <- ggplot(res_df |> filter(!is.na(padj)),
                aes(log2FoldChange, -log10(padj), color = regulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = c(-opt$lfc, opt$lfc), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(opt$padj),      linetype = "dashed", color = "grey40") +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(padj)") +
  theme_bw()
ggsave(file.path(opt$figures_dir, "volcano_plot.pdf"), p_vol, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "volcano_plot.png"), p_vol, width = 7, height = 5, dpi = 300)

# ── MA plot ───────────────────────────────────────────────────
p_ma <- ggplot(res_df |> filter(!is.na(padj)),
               aes(log10(baseMean + 1), log2FoldChange, color = regulation)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "MA Plot", x = "log10(baseMean + 1)", y = "log2 Fold Change") +
  theme_bw()
ggsave(file.path(opt$figures_dir, "ma_plot.pdf"), p_ma, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "ma_plot.png"), p_ma, width = 7, height = 5, dpi = 300)

# ── Heatmap top 50 genes ──────────────────────────────────────
top50 <- res_df |> filter(regulation != "ns") |>
         arrange(padj) |> head(50) |> pull(gene_id)
if (length(top50) >= 2) {
  mat <- assay(vst)[top50, , drop = FALSE]
  mat_scaled <- t(scale(t(mat)))
  ann <- data.frame(condition = meta$condition, row.names = meta$sample)

  pdf(file.path(opt$figures_dir, "heatmap_top_genes.pdf"), width = 9, height = 10)
  pheatmap(mat_scaled, annotation_col = ann, show_rownames = length(top50) <= 50,
           fontsize_row = 7, main = "Top 50 DEGs – Z-score")
  dev.off()

  png(file.path(opt$figures_dir, "heatmap_top_genes.png"), width = 900, height = 1000, res = 120)
  pheatmap(mat_scaled, annotation_col = ann, show_rownames = length(top50) <= 50,
           fontsize_row = 7, main = "Top 50 DEGs – Z-score")
  dev.off()
}

# ── Barplot DE ────────────────────────────────────────────────
de_bar <- data.frame(
  regulation = c("Up-regulated", "Down-regulated"),
  count      = c(n_up, n_down)
)
p_bar <- ggplot(de_bar, aes(regulation, count, fill = regulation)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8")) +
  geom_text(aes(label = count), vjust = -0.3, size = 5) +
  labs(title = "DEGs por Regulação", x = NULL, y = "Número de genes") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(opt$figures_dir, "de_barplot.pdf"), p_bar, width = 5, height = 4)
ggsave(file.path(opt$figures_dir, "de_barplot.png"), p_bar, width = 5, height = 4, dpi = 300)

cat("DESeq2 concluído.\n")
