#!/usr/bin/env Rscript
# ============================================================
# 09_genie3.R – Inferência de rede regulatória via GENIE3
# Método: Random Forest (Huynh-Thu et al. 2010)
# Reguladores: TFs do PlantTFDB presentes na matriz de expressão
# Alvos: DEGs
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GENIE3)
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
})

opt_list <- list(
  make_option("--norm_counts",  type = "character"),
  make_option("--deseq2",       type = "character"),
  make_option("--tf_classified",type = "character", default = NULL),
  make_option("--plantfdb_file",type = "character", default = NULL),
  make_option("--padj",         type = "double",    default = 0.05),
  make_option("--lfc",          type = "double",    default = 1.0),
  make_option("--n_trees",      type = "integer",   default = 500L),
  make_option("--n_links",      type = "integer",   default = 5000L),
  make_option("--ncores",       type = "integer",   default = 1L),
  make_option("--outdir",       type = "character", default = "."),
  make_option("--figures_dir",  type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("norm_counts", "deseq2"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ─────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "genie3_network.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "genie3_hub_tfs.tsv"))
writeLines("GENIE3: inicializando...", file.path(opt$outdir, "genie3_summary.txt"))
writeLines("", file.path(opt$figures_dir, ".genie3_placeholder"))

# ── Leitura ────────────────────────────────────────────────────
counts_df <- read_tsv(opt$norm_counts, show_col_types = FALSE)
deseq2    <- read_tsv(opt$deseq2, show_col_types = FALSE) |>
  mutate(gene_id = gsub("\\.TAIR10$", "", gene_id)) |>
  filter(!is.na(padj))

deg_genes <- deseq2 |>
  filter(padj < opt$padj, abs(log2FoldChange) > opt$lfc) |>
  pull(gene_id)

cat(sprintf("DEGs (alvos): %d\n", length(deg_genes)))

if (length(deg_genes) < 2) {
  writeLines("Menos de 2 DEGs – GENIE3 não executado.",
             file.path(opt$outdir, "genie3_summary.txt"))
  quit(save = "no", status = 0)
}

# ── Monta matriz de expressão (genes × amostras) ───────────────
expr_mat <- counts_df |>
  column_to_rownames("gene_id") |>
  as.matrix()

# Remove genes com variância zero
expr_mat <- expr_mat[apply(expr_mat, 1, var) > 0, , drop = FALSE]
cat(sprintf("Genes na matriz: %d | Amostras: %d\n", nrow(expr_mat), ncol(expr_mat)))

# ── Carrega lista de TFs reguladores ──────────────────────────
# Prioridade: (1) PlantTFDB full, (2) tf_classified (só DEG TFs),
# (3) fallback: todos os genes como reguladores
load_tf_list <- function(plantfdb_file, tf_classified) {
  # Tenta PlantTFDB completo primeiro (todos os TFs background)
  if (!is.null(plantfdb_file) && file.exists(plantfdb_file) &&
      file.size(plantfdb_file) > 0) {
    df <- tryCatch(read_tsv(plantfdb_file, show_col_types = FALSE),
                   error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      cols_lower <- tolower(colnames(df))
      col_gene <- colnames(df)[which(cols_lower %in% c("gene_id","geneid","tair_id","locus"))[1]]
      if (!is.na(col_gene)) {
        ids <- gsub("\\.TAIR10$|\\.[0-9]+$", "", as.character(df[[col_gene]]))
        ids <- unique(ids[!is.na(ids) & nchar(ids) > 0])
        cat(sprintf("TFs carregados do PlantTFDB completo: %d\n", length(ids)))
        return(ids)
      }
    }
  }

  # Tenta tf_classified (apenas TFs DEGs)
  if (!is.null(tf_classified) && file.exists(tf_classified) &&
      file.size(tf_classified) > 0) {
    df <- tryCatch(read_tsv(tf_classified, show_col_types = FALSE),
                   error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0 && "gene_id" %in% colnames(df)) {
      ids <- unique(df$gene_id)
      cat(sprintf("TFs carregados de tf_classified (só DEG TFs): %d\n", length(ids)))
      return(ids)
    }
  }

  cat("AVISO: nenhuma lista de TFs carregada. Usando todos os genes como reguladores.\n")
  NULL
}

tf_arg_classified <- if (!is.null(opt$tf_classified) &&
                          opt$tf_classified != "NO_FILE" &&
                          file.exists(opt$tf_classified)) opt$tf_classified else NULL
tf_arg_plantfdb   <- if (!is.null(opt$plantfdb_file) &&
                          opt$plantfdb_file != "NO_FILE" &&
                          file.exists(opt$plantfdb_file)) opt$plantfdb_file else NULL

all_tf_ids <- load_tf_list(tf_arg_plantfdb, tf_arg_classified)

# Restringe reguladores aos que estão na matriz de expressão
if (!is.null(all_tf_ids)) {
  regulators <- intersect(all_tf_ids, rownames(expr_mat))
  cat(sprintf("TFs na matriz de expressão: %d\n", length(regulators)))
  if (length(regulators) == 0) {
    cat("AVISO: nenhum TF encontrado na matriz. Usando todos os genes.\n")
    regulators <- NULL
  }
} else {
  regulators <- NULL
}

# Restringe alvos (DEGs) aos que estão na matriz
targets <- intersect(deg_genes, rownames(expr_mat))
cat(sprintf("Alvos (DEGs na matriz): %d\n", length(targets)))

if (length(targets) < 2) {
  writeLines("DEGs insuficientes na matriz de expressão.",
             file.path(opt$outdir, "genie3_summary.txt"))
  quit(save = "no", status = 0)
}

# ── GENIE3 ────────────────────────────────────────────────────
cat(sprintf("Iniciando GENIE3: %d reguladores × %d alvos | %d árvores | %d cores\n",
            if (is.null(regulators)) nrow(expr_mat) else length(regulators),
            length(targets), opt$n_trees, opt$ncores))

run_genie3 <- function(expr_mat, regulators, targets, n_trees, ncores) {
  tryCatch(
    GENIE3(exprMatrix = expr_mat, regulators = regulators, targets = targets,
           treeMethod = "RF", K = "sqrt", nTrees = n_trees, nCores = ncores),
    error = function(e) {
      message(sprintf("GENIE3 (nCores=%d) falhou: %s", ncores, e$message))
      if (ncores > 1L) {
        message("Repetindo com nCores=1 (modo serial)...")
        tryCatch(
          GENIE3(exprMatrix = expr_mat, regulators = regulators, targets = targets,
                 treeMethod = "RF", K = "sqrt", nTrees = n_trees, nCores = 1L),
          error = function(e2) {
            message("GENIE3 (nCores=1) também falhou: ", e2$message)
            NULL
          }
        )
      } else NULL
    }
  )
}

set.seed(42)
weight_mat <- run_genie3(expr_mat, regulators, targets, opt$n_trees, opt$ncores)

if (is.null(weight_mat)) {
  writeLines("GENIE3 falhou durante execução.", file.path(opt$outdir, "genie3_summary.txt"))
  quit(save = "no", status = 0)
}

cat("GENIE3 concluído. Extraindo link list...\n")

# ── Link list ──────────────────────────────────────────────────
links <- getLinkList(weight_mat, reportMax = opt$n_links) |>
  as_tibble() |>
  rename(tf = regulatoryGene, target = targetGene, weight = weight) |>
  filter(weight > 0) |>
  left_join(deseq2 |> select(gene_id, lfc_tf = log2FoldChange,
                              padj_tf = padj),
            by = c("tf" = "gene_id")) |>
  left_join(deseq2 |> select(gene_id, lfc_target = log2FoldChange,
                              padj_target = padj),
            by = c("target" = "gene_id"))

write_tsv(links, file.path(opt$outdir, "genie3_network.tsv"))
cat(sprintf("Links exportados: %d\n", nrow(links)))

# ── Hub TFs (por out-degree ponderado) ────────────────────────
hub_tfs <- links |>
  group_by(tf) |>
  summarise(
    n_targets    = n(),
    mean_weight  = mean(weight),
    total_weight = sum(weight),
    .groups = "drop"
  ) |>
  left_join(
    deseq2 |> select(gene_id, log2FoldChange, padj),
    by = c("tf" = "gene_id")
  ) |>
  mutate(regulation = case_when(
    !is.na(log2FoldChange) & log2FoldChange >  opt$lfc ~ "up",
    !is.na(log2FoldChange) & log2FoldChange < -opt$lfc ~ "down",
    TRUE ~ "not_deg"
  )) |>
  arrange(desc(n_targets))

write_tsv(hub_tfs, file.path(opt$outdir, "genie3_hub_tfs.tsv"))

# ── Sumário ────────────────────────────────────────────────────
top10_hubs <- head(hub_tfs, 10)
summary_lines <- c(
  sprintf("GENIE3 – Rede Regulatória Arabidopsis thaliana"),
  sprintf("Reguladores usados: %d",
          if (is.null(regulators)) nrow(expr_mat) else length(regulators)),
  sprintf("Alvos (DEGs): %d", length(targets)),
  sprintf("Árvores por gene: %d", opt$n_trees),
  sprintf("Links na rede (top %d): %d", opt$n_links, nrow(links)),
  sprintf("TFs reguladores identificados: %d", nrow(hub_tfs)),
  "",
  "Top 10 TFs hub (por número de alvos):",
  sprintf("  %-15s alvos=%d peso_total=%.4f [%s]",
          top10_hubs$tf, top10_hubs$n_targets,
          top10_hubs$total_weight, top10_hubs$regulation)
)
writeLines(summary_lines, file.path(opt$outdir, "genie3_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# ── Figura 1: top hub TFs ─────────────────────────────────────
top30 <- head(hub_tfs, 30)
reg_colors <- c("up" = "#E41A1C", "down" = "#377EB8", "not_deg" = "#888888")

p_hub <- ggplot(top30, aes(reorder(tf, n_targets), n_targets,
                            fill = regulation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = reg_colors) +
  labs(title = "Top 30 TFs Hub – GENIE3",
       subtitle = sprintf("Total de links: %d | Alvos: %d DEGs",
                          nrow(links), length(targets)),
       x = NULL, y = "Número de alvos DEG", fill = "Status") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$figures_dir, "genie3_hub_tfs.pdf"), p_hub, width = 8, height = 8)
ggsave(file.path(opt$figures_dir, "genie3_hub_tfs.png"), p_hub, width = 8, height = 8, dpi = 300)

# ── Figura 2: distribuição de pesos ───────────────────────────
p_weights <- ggplot(links, aes(weight)) +
  geom_histogram(bins = 50, fill = "#2166AC", color = "white") +
  scale_x_log10() +
  labs(title = "Distribuição de Pesos GENIE3",
       x = "Peso (escala log10)", y = "Frequência") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$figures_dir, "genie3_weight_dist.pdf"), p_weights, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "genie3_weight_dist.png"), p_weights, width = 7, height = 5, dpi = 300)

# ── Figura 3: sub-rede dos top 10 TFs ─────────────────────────
top10_tf_names <- top10_hubs$tf
sub_links <- links |>
  filter(tf %in% top10_tf_names) |>
  arrange(desc(weight)) |>
  head(200)   # limita para visualização

if (nrow(sub_links) >= 2) {
  g_sub <- graph_from_data_frame(sub_links[, c("tf","target","weight")],
                                  directed = TRUE)

  # Cor dos nós: TF = círculo colorido por regulação; target = cinza
  is_tf     <- V(g_sub)$name %in% top10_tf_names
  tf_reg    <- setNames(hub_tfs$regulation, hub_tfs$tf)
  node_col  <- ifelse(is_tf,
                      ifelse(tf_reg[V(g_sub)$name] == "up",   "#E41A1C",
                             ifelse(tf_reg[V(g_sub)$name] == "down", "#377EB8",
                                    "#F39C12")),
                      "#DDDDDD")
  node_size <- ifelse(is_tf, 12, 5)

  pdf(file.path(opt$figures_dir, "genie3_subnetwork.pdf"), width = 12, height = 10)
  par(mar = c(1, 1, 2, 1))
  set.seed(42)
  plot(g_sub,
       vertex.color = node_col,
       vertex.size  = node_size,
       vertex.label = ifelse(is_tf, V(g_sub)$name, NA),
       vertex.label.cex   = 0.7,
       vertex.label.color = "black",
       edge.arrow.size    = 0.3,
       edge.width         = E(g_sub)$weight * 20,
       edge.color         = "grey70",
       layout             = layout_with_fr(g_sub),
       main = "Sub-rede Top 10 TFs Reguladores – GENIE3")
  legend("bottomleft",
         legend = c("TF up", "TF down", "TF não-DEG", "Gene alvo"),
         fill   = c("#E41A1C","#377EB8","#F39C12","#DDDDDD"),
         bty = "n", cex = 0.8)
  dev.off()
}

cat("GENIE3 concluído.\n")
