#!/usr/bin/env Rscript
# ============================================================
# 04_integration.R – Integração multi-ômica e ranking de genes
# Integration score (0–10):
#   lfc_score×3 + mean_score×2 + sig_score×2
#   + has_splicing×1.5 + in_pathway×1.0 + is_hub×0.5
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
})

opt_list <- list(
  make_option("--deseq2",      type = "character"),
  make_option("--norm_counts", type = "character"),
  make_option("--splicing",    type = "character"),
  make_option("--go_bp",       type = "character"),
  make_option("--kegg",        type = "character"),
  make_option("--wgcna",       type = "character"),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1.0),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))
for (p in c("deseq2","norm_counts"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Leitura ──────────────────────────────────────────────────
res  <- read_tsv(opt$deseq2,      show_col_types = FALSE)
expr <- read_tsv(opt$norm_counts, show_col_types = FALSE)

safe_read <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  tryCatch(read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
}

splicing <- safe_read(opt$splicing)
go_bp    <- safe_read(opt$go_bp)
kegg     <- safe_read(opt$kegg)
wgcna    <- safe_read(opt$wgcna)

# ── Genes spliced ─────────────────────────────────────────────
spliced_genes <- character(0)
if (!is.null(splicing) && nrow(splicing) > 0) {
  gene_col <- intersect(c("geneSymbol","GeneID","gene_id"), colnames(splicing))[1]
  if (!is.na(gene_col)) spliced_genes <- unique(splicing[[gene_col]])
}

# ── Genes em vias enriquecidas ────────────────────────────────
pathway_genes <- character(0)
for (df in list(go_bp, kegg)) {
  if (!is.null(df) && nrow(df) > 0 && "geneID" %in% colnames(df)) {
    new_g <- unlist(strsplit(df$geneID, "/"))
    pathway_genes <- union(pathway_genes, new_g)
  }
}

# ── Hub genes WGCNA ──────────────────────────────────────────
hub_genes_list <- character(0)
if (!is.null(wgcna) && nrow(wgcna) > 0 && "gene_id" %in% colnames(wgcna)) {
  hub_genes_list <- unique(wgcna$gene_id)
}

# ── Calcula baseMean por gene (média das contagens normalizadas) ──
mean_expr <- expr |>
  select(-gene_id) |>
  rowMeans()
expr_means <- tibble(gene_id = expr$gene_id, mean_expr = mean_expr)

# ── Monta tabela integrada ────────────────────────────────────
integrated <- res |>
  filter(!is.na(padj)) |>
  left_join(expr_means, by = "gene_id") |>
  mutate(
    is_de          = padj < opt$padj & abs(log2FoldChange) > opt$lfc,
    has_splicing   = gene_id %in% spliced_genes,
    in_pathway     = gene_id %in% pathway_genes,
    is_hub         = gene_id %in% hub_genes_list,
    # Scores normalizados 0–1
    lfc_score  = pmin(abs(log2FoldChange) / 5, 1),
    mean_score = pmin(log10(mean_expr + 1) / 5, 1),
    sig_score  = ifelse(padj > 0, pmin(-log10(padj) / 10, 1), 1),
    # Integration score 0–10
    integration_score = lfc_score * 3 +
                        mean_score * 2 +
                        sig_score  * 2 +
                        as.numeric(has_splicing) * 1.5 +
                        as.numeric(in_pathway)   * 1.0 +
                        as.numeric(is_hub)       * 0.5,
    evidence_layers = as.integer(is_de) +
                      as.integer(has_splicing) +
                      as.integer(in_pathway) +
                      as.integer(is_hub)
  ) |>
  arrange(desc(integration_score))

write_tsv(integrated, file.path(opt$outdir, "integrated_genes.tsv"))

# Ranking: apenas DEGs
ranking <- integrated |> filter(is_de) |>
           select(gene_id, log2FoldChange, padj, regulation,
                  integration_score, evidence_layers,
                  has_splicing, in_pathway, is_hub)
write_tsv(ranking, file.path(opt$outdir, "gene_ranking.tsv"))

# Candidatos chave: evidência em ≥ 2 camadas
candidates <- integrated |> filter(evidence_layers >= 2) |>
              arrange(desc(integration_score))
write_tsv(candidates, file.path(opt$outdir, "key_candidates.tsv"))

# Tabela resumida de candidatos
candidates_table <- candidates |>
  select(gene_id, log2FoldChange, padj, integration_score,
         evidence_layers, has_splicing, in_pathway, is_hub) |>
  head(100)
write_tsv(candidates_table, file.path(opt$outdir, "candidates_table.tsv"))

cat(sprintf("Genes integrados: %d | Candidatos chave: %d\n",
            nrow(integrated), nrow(candidates)))

# ── Figura: camadas de evidência ─────────────────────────────
ev_df <- integrated |> filter(is_de) |>
  summarise(
    DEG         = sum(is_de),
    Splicing    = sum(has_splicing),
    Pathways    = sum(in_pathway),
    `Hub WGCNA` = sum(is_hub)
  ) |>
  pivot_longer(everything(), names_to = "layer", values_to = "count")

p_ev <- ggplot(ev_df, aes(reorder(layer, -count), count, fill = layer)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = count), vjust = -0.3, size = 4.5) +
  labs(title = "Genes por Camada de Evidência",
       x = NULL, y = "Número de genes") +
  theme_bw()

# Figura: top 20 genes por integration score
top20 <- ranking |> head(20)
p_top <- ggplot(top20, aes(reorder(gene_id, integration_score),
                             integration_score, fill = regulation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8")) +
  labs(title = "Top 20 Genes – Integration Score",
       x = NULL, y = "Integration Score (0–10)") +
  theme_bw()

p_comb <- p_ev / p_top
ggsave(file.path(opt$figures_dir, "evidence_layers.pdf"), p_comb, width = 8, height = 10)
ggsave(file.path(opt$figures_dir, "evidence_layers.png"), p_comb, width = 8, height = 10, dpi = 300)
ggsave(file.path(opt$figures_dir, "top_genes_integration.pdf"), p_top, width = 8, height = 7)
ggsave(file.path(opt$figures_dir, "top_genes_integration.png"), p_top, width = 8, height = 7, dpi = 300)

cat("Integração multi-ômica concluída.\n")
