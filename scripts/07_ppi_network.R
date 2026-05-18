#!/usr/bin/env Rscript
# ============================================================
# 07_ppi_network.R – Rede PPI via API STRING (species=3702)
# Identifica hubs por grau e centralidade de betweenness
# Integra com integration_score do 04_integration.R
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(httr)
  library(jsonlite)
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
})

opt_list <- list(
  make_option("--deseq2",      type = "character"),
  make_option("--integration", type = "character"),
  make_option("--padj",        type = "double",  default = 0.05),
  make_option("--lfc",         type = "double",  default = 1.0),
  make_option("--score",       type = "integer", default = 400L),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("deseq2"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Leitura ────────────────────────────────────────────────────
deseq2 <- read_tsv(opt$deseq2, show_col_types = FALSE) |>
  mutate(gene_id = gsub("\\.TAIR10$", "", gene_id))

# Lê integration opcionalmente
integration <- NULL
if (!is.null(opt$integration) && file.exists(opt$integration)) {
  integration <- tryCatch(
    read_tsv(opt$integration, show_col_types = FALSE),
    error = function(e) NULL
  )
}

# ── Filtra DEGs ────────────────────────────────────────────────
degs <- deseq2 |>
  filter(!is.na(padj), padj < opt$padj, abs(log2FoldChange) > opt$lfc)

cat(sprintf("DEGs para consulta STRING: %d\n", nrow(degs)))

if (nrow(degs) < 2) {
  cat("AVISO: Menos de 2 DEGs – rede PPI não gerada.\n")
  write_tsv(data.frame(), file.path(opt$outdir, "ppi_edges.tsv"))
  write_tsv(data.frame(), file.path(opt$outdir, "hub_genes.tsv"))
  writeLines("Menos de 2 DEGs disponíveis.", file.path(opt$outdir, "network_summary.txt"))
  quit(save = "no", status = 0)
}

# ── Consulta API STRING ────────────────────────────────────────
query_string <- function(genes, species = 3702, min_score = 400) {
  cat(sprintf("Consultando STRING API para %d genes (score≥%d)...\n",
              length(genes), min_score))
  resp <- tryCatch(
    POST(
      url    = "https://string-db.org/api/json/network",
      body   = list(
        identifiers     = paste(genes, collapse = "\r"),
        species         = as.character(species),
        required_score  = as.character(min_score),
        caller_identity = "rnaseq-arabidopsis-pipeline"
      ),
      encode  = "form",
      timeout(60)
    ),
    error = function(e) { message("STRING API erro: ", e$message); NULL }
  )
  if (is.null(resp) || http_error(resp)) {
    message("STRING API indisponível ou erro HTTP.")
    return(NULL)
  }
  parsed <- tryCatch(content(resp, as = "parsed", simplifyVector = TRUE),
                     error = function(e) NULL)
  if (is.null(parsed) || length(parsed) == 0) return(NULL)
  as.data.frame(parsed)
}

edges_raw <- query_string(degs$gene_id, min_score = opt$score)

if (is.null(edges_raw) || nrow(edges_raw) == 0) {
  cat("AVISO: STRING não retornou interações. Verificar conectividade ou IDs.\n")
  write_tsv(data.frame(), file.path(opt$outdir, "ppi_edges.tsv"))
  write_tsv(data.frame(), file.path(opt$outdir, "hub_genes.tsv"))
  writeLines("Nenhuma interação retornada pela API STRING.",
             file.path(opt$outdir, "network_summary.txt"))
  quit(save = "no", status = 0)
}

cat(sprintf("Interações retornadas: %d\n", nrow(edges_raw)))

# ── Constrói edge list padronizado ─────────────────────────────
# STRING retorna preferredName_A / preferredName_B (e.g. AT1G01010)
edges <- edges_raw |>
  select(
    gene_a = preferredName_A,
    gene_b = preferredName_B,
    score  = score
  ) |>
  mutate(score = as.numeric(score)) |>
  filter(!is.na(score), gene_a != gene_b) |>
  distinct()

write_tsv(edges, file.path(opt$outdir, "ppi_edges.tsv"))

# ── Constrói grafo ─────────────────────────────────────────────
g <- graph_from_data_frame(edges, directed = FALSE)
E(g)$weight <- edges$score / 1000   # normaliza para 0-1

# ── Métricas de centralidade ───────────────────────────────────
deg_vals <- degree(g)
btw_vals <- betweenness(g, normalized = TRUE)
cls_vals <- closeness(g, normalized = TRUE)

node_df <- tibble(
  gene_id     = V(g)$name,
  degree      = deg_vals,
  betweenness = btw_vals,
  closeness   = cls_vals
) |>
  left_join(
    deseq2 |> select(gene_id, log2FoldChange, padj,
                     regulation = if ("regulation" %in% colnames(deseq2)) "regulation" else "log2FoldChange"),
    by = "gene_id"
  )

# Resolve coluna regulation se não existir
if (!"regulation" %in% colnames(node_df)) {
  node_df <- node_df |>
    mutate(regulation = case_when(
      log2FoldChange >  opt$lfc ~ "up",
      log2FoldChange < -opt$lfc ~ "down",
      TRUE                      ~ "ns"
    ))
}

# Adiciona integration_score se disponível
if (!is.null(integration) && "integration_score" %in% colnames(integration)) {
  node_df <- node_df |>
    left_join(integration |> select(gene_id, integration_score), by = "gene_id")
}

# ── Define hub genes: top 10% por grau ─────────────────────────
hub_threshold   <- quantile(deg_vals, 0.90)
node_df$is_hub  <- node_df$degree >= hub_threshold

hub_genes <- node_df |>
  filter(is_hub) |>
  arrange(desc(degree), desc(betweenness))

write_tsv(node_df,   file.path(opt$outdir, "ppi_nodes.tsv"))
write_tsv(hub_genes, file.path(opt$outdir, "hub_genes.tsv"))

cat(sprintf("Nós na rede: %d | Arestas: %d | Hubs (top 10%%): %d\n",
            vcount(g), ecount(g), sum(node_df$is_hub)))

# ── Sumário ────────────────────────────────────────────────────
summary_lines <- c(
  sprintf("Nós (proteínas): %d",       vcount(g)),
  sprintf("Arestas (interações): %d",  ecount(g)),
  sprintf("Score mínimo STRING: %d",   opt$score),
  sprintf("Hub genes (top 10%%): %d",  sum(node_df$is_hub)),
  "",
  "Top 10 hub genes:",
  sprintf("  %s (grau=%d, btw=%.4f)",
          hub_genes$gene_id[seq_len(min(10, nrow(hub_genes)))],
          hub_genes$degree[seq_len(min(10, nrow(hub_genes)))],
          hub_genes$betweenness[seq_len(min(10, nrow(hub_genes)))])
)
writeLines(summary_lines, file.path(opt$outdir, "network_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# ── Figura 1: distribuição de grau ────────────────────────────
p_deg <- ggplot(node_df, aes(degree, fill = regulation)) +
  geom_histogram(bins = 20, color = "white") +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8", "ns" = "#AAAAAA"),
                    na.value = "#AAAAAA") +
  labs(title = "Distribuição de Grau – Rede PPI (STRING)",
       x = "Grau (número de interações)", y = "Frequência", fill = "DEG") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$figures_dir, "ppi_degree_dist.pdf"), p_deg, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "ppi_degree_dist.png"), p_deg, width = 7, height = 5, dpi = 300)

# ── Figura 2: hub genes lollipop ──────────────────────────────
if (nrow(hub_genes) > 0) {
  top_hubs <- head(hub_genes, 30)
  reg_colors <- c("up" = "#E41A1C", "down" = "#377EB8", "ns" = "#999999")
  p_hub <- ggplot(top_hubs, aes(reorder(gene_id, degree), degree,
                                 color = regulation)) +
    geom_segment(aes(xend = gene_id, yend = 0), linewidth = 0.8) +
    geom_point(aes(size = betweenness)) +
    coord_flip() +
    scale_color_manual(values = reg_colors, na.value = "#999999") +
    labs(title = "Hub Genes – Rede PPI STRING",
         x = NULL, y = "Grau (degree)", size = "Betweenness", color = "DEG") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$figures_dir, "ppi_hub_genes.pdf"), p_hub, width = 8, height = 7)
  ggsave(file.path(opt$figures_dir, "ppi_hub_genes.png"), p_hub, width = 8, height = 7, dpi = 300)
}

# ── Figura 3: visualização da rede (subgrafo dos hubs) ────────
hub_names <- hub_genes$gene_id
g_hub     <- induced_subgraph(g, vids = V(g)[V(g)$name %in% hub_names])

if (vcount(g_hub) >= 2) {
  node_colors <- setNames(
    c("#E41A1C", "#377EB8", "#AAAAAA"),
    c("up", "down", "ns")
  )
  vcolors <- node_colors[node_df$regulation[match(V(g_hub)$name, node_df$gene_id)]]
  vcolors[is.na(vcolors)] <- "#AAAAAA"
  vsizes  <- log1p(degree(g_hub)) * 4 + 3

  pdf(file.path(opt$figures_dir, "ppi_hub_network.pdf"), width = 10, height = 10)
  par(mar = c(1, 1, 2, 1))
  set.seed(42)
  plot(g_hub,
       vertex.color = vcolors,
       vertex.size  = vsizes,
       vertex.label = V(g_hub)$name,
       vertex.label.cex   = 0.65,
       vertex.label.color = "black",
       edge.width   = E(g_hub)$weight * 2,
       edge.color   = "grey70",
       layout       = layout_with_fr(g_hub),
       main         = "Subgraph Hub Genes – STRING PPI (Arabidopsis)")
  legend("bottomleft", legend = c("Up-regulated","Down-regulated","Não-DEG"),
         fill = c("#E41A1C","#377EB8","#AAAAAA"), bty = "n", cex = 0.8)
  dev.off()
}

cat("Rede PPI concluída.\n")
