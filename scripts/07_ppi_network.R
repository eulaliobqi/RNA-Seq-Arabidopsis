#!/usr/bin/env Rscript
# ============================================================
# 07_ppi_network.R – Rede PPI via STRING (species=3702)
# Tenta STRINGdb (pacote R), depois httr (API REST).
# Todos os outputs são escritos no início (vazios) para que
# Nextflow não falhe por arquivos ausentes em saída antecipada.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
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

# ── Inicializa todos os outputs vazios ─────────────────────────
# Garante que Nextflow não falhe por arquivo ausente caso qualquer
# etapa termine antes da hora.
write_tsv(data.frame(), file.path(opt$outdir, "ppi_edges.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "ppi_nodes.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "hub_genes.tsv"))
writeLines("Rede PPI: inicializando...", file.path(opt$outdir, "network_summary.txt"))
# Placeholder para o diretório de figuras
writeLines("", file.path(opt$figures_dir, ".ppi_placeholder"))

# ── Leitura ────────────────────────────────────────────────────
deseq2 <- read_tsv(opt$deseq2, show_col_types = FALSE) |>
  mutate(gene_id = gsub("\\.TAIR10$", "", gene_id))

integration <- NULL
if (!is.null(opt$integration) && file.exists(opt$integration)) {
  integration <- tryCatch(read_tsv(opt$integration, show_col_types = FALSE),
                          error = function(e) NULL)
}

degs <- deseq2 |>
  filter(!is.na(padj), padj < opt$padj, abs(log2FoldChange) > opt$lfc)

cat(sprintf("DEGs para consulta STRING: %d\n", nrow(degs)))

if (nrow(degs) < 2) {
  writeLines("Menos de 2 DEGs – rede PPI não gerada.",
             file.path(opt$outdir, "network_summary.txt"))
  cat("AVISO: Menos de 2 DEGs.\n")
  quit(save = "no", status = 0)
}

# ── Tentativa 1: STRINGdb (pacote R, mais robusto) ────────────
query_via_stringdb <- function(genes, species = 3702, score = 400) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) return(NULL)
  cat("Usando STRINGdb R package...\n")
  db <- tryCatch(
    STRINGdb::STRINGdb$new(version = "12.0", species = species,
                            score_threshold = score, network_type = "full"),
    error = function(e) { message("STRINGdb$new falhou: ", e$message); NULL }
  )
  if (is.null(db)) return(NULL)

  mapped <- tryCatch(
    db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE),
    error = function(e) { message("STRINGdb map falhou: ", e$message); NULL }
  )
  if (is.null(mapped) || nrow(mapped) == 0) {
    message("STRINGdb: nenhum gene mapeado."); return(NULL)
  }
  cat(sprintf("Genes mapeados STRING: %d / %d\n", nrow(mapped), length(genes)))

  ints <- tryCatch(db$get_interactions(mapped$STRING_id),
                   error = function(e) NULL)
  if (is.null(ints) || nrow(ints) == 0) return(NULL)

  id_map <- setNames(mapped$gene, mapped$STRING_id)
  ints$gene_a <- id_map[ints$from]
  ints$gene_b <- id_map[ints$to]
  ints <- ints[!is.na(ints$gene_a) & !is.na(ints$gene_b), ]
  data.frame(gene_a = ints$gene_a, gene_b = ints$gene_b,
             score  = as.numeric(ints$combined_score))
}

# ── Tentativa 2: httr REST API ────────────────────────────────
query_via_httr <- function(genes, species = 3702, score = 400) {
  if (!requireNamespace("httr", quietly = TRUE)) return(NULL)
  cat("Usando STRING REST API...\n")
  resp <- tryCatch(
    httr::POST(
      url    = "https://string-db.org/api/json/network",
      body   = list(
        identifiers     = paste(genes, collapse = "\r"),
        species         = as.character(species),
        required_score  = as.character(score),
        caller_identity = "rnaseq-arabidopsis-pipeline"
      ),
      encode  = "form",
      httr::timeout(60)
    ),
    error = function(e) { message("httr falhou: ", e$message); NULL }
  )
  if (is.null(resp) || httr::http_error(resp)) return(NULL)

  parsed <- tryCatch(
    httr::content(resp, as = "parsed", simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(parsed) || length(parsed) == 0 || !is.data.frame(parsed)) return(NULL)
  if (!all(c("preferredName_A", "preferredName_B", "score") %in% colnames(parsed)))
    return(NULL)

  data.frame(
    gene_a = parsed$preferredName_A,
    gene_b = parsed$preferredName_B,
    score  = as.numeric(parsed$score)
  )
}

# ── Executa consulta ───────────────────────────────────────────
edges <- query_via_stringdb(degs$gene_id, score = opt$score)
if (is.null(edges)) edges <- query_via_httr(degs$gene_id, score = opt$score)

if (is.null(edges) || nrow(edges) == 0) {
  msg <- paste(
    "STRING não retornou interações.",
    "Possíveis causas: sem acesso à internet no servidor,",
    "IDs não reconhecidos, ou nenhuma interação acima do score mínimo.",
    "Instale STRINGdb: mamba install -n r-analysis bioconductor-stringdb",
    sep = "\n"
  )
  writeLines(msg, file.path(opt$outdir, "network_summary.txt"))
  cat(sprintf("AVISO: %s\n", gsub("\n", " | ", msg)))
  quit(save = "no", status = 0)
}

cat(sprintf("Interações retornadas: %d\n", nrow(edges)))

edges <- edges |>
  filter(!is.na(score), gene_a != gene_b) |>
  distinct()

write_tsv(edges, file.path(opt$outdir, "ppi_edges.tsv"))

# ── Constrói grafo ─────────────────────────────────────────────
g <- graph_from_data_frame(edges, directed = FALSE)
E(g)$weight <- edges$score / 1000

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
    deseq2 |> select(gene_id, log2FoldChange, padj),
    by = "gene_id"
  ) |>
  mutate(regulation = case_when(
    !is.na(log2FoldChange) & log2FoldChange >  opt$lfc ~ "up",
    !is.na(log2FoldChange) & log2FoldChange < -opt$lfc ~ "down",
    TRUE                                                ~ "ns"
  ))

if (!is.null(integration) && "integration_score" %in% colnames(integration)) {
  node_df <- node_df |>
    left_join(integration |> select(gene_id, integration_score), by = "gene_id")
}

hub_threshold  <- quantile(deg_vals, 0.90)
node_df$is_hub <- node_df$degree >= hub_threshold

hub_genes <- node_df |> filter(is_hub) |> arrange(desc(degree), desc(betweenness))

write_tsv(node_df,   file.path(opt$outdir, "ppi_nodes.tsv"))
write_tsv(hub_genes, file.path(opt$outdir, "hub_genes.tsv"))

cat(sprintf("Nós: %d | Arestas: %d | Hubs (top 10%%): %d\n",
            vcount(g), ecount(g), sum(node_df$is_hub)))

# ── Sumário ────────────────────────────────────────────────────
top10 <- head(hub_genes, 10)
summary_lines <- c(
  sprintf("Nós (proteínas): %d",       vcount(g)),
  sprintf("Arestas (interações): %d",  ecount(g)),
  sprintf("Score mínimo STRING: %d",   opt$score),
  sprintf("Hub genes (top 10%%): %d",  sum(node_df$is_hub)),
  "",
  "Top 10 hub genes:",
  if (nrow(top10) > 0)
    sprintf("  %s (grau=%d, btw=%.4f)", top10$gene_id, top10$degree, top10$betweenness)
  else "  (nenhum)"
)
writeLines(summary_lines, file.path(opt$outdir, "network_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# ── Figura 1: distribuição de grau ────────────────────────────
p_deg <- ggplot(node_df, aes(degree, fill = regulation)) +
  geom_histogram(bins = 20, color = "white") +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8", "ns" = "#AAAAAA"),
                    na.value = "#AAAAAA") +
  labs(title = "Distribuição de Grau – Rede PPI (STRING)",
       x = "Grau", y = "Frequência", fill = "DEG") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$figures_dir, "ppi_degree_dist.pdf"), p_deg, width = 7, height = 5)
ggsave(file.path(opt$figures_dir, "ppi_degree_dist.png"), p_deg, width = 7, height = 5, dpi = 300)

# ── Figura 2: hub genes ────────────────────────────────────────
if (nrow(hub_genes) > 0) {
  top_hubs <- head(hub_genes, 30)
  p_hub <- ggplot(top_hubs, aes(reorder(gene_id, degree), degree, color = regulation)) +
    geom_segment(aes(xend = gene_id, yend = 0), linewidth = 0.8) +
    geom_point(aes(size = betweenness)) +
    coord_flip() +
    scale_color_manual(values = c("up" = "#E41A1C", "down" = "#377EB8", "ns" = "#999999"),
                       na.value = "#999999") +
    labs(title = "Hub Genes – STRING PPI", x = NULL, y = "Grau",
         size = "Betweenness", color = "DEG") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$figures_dir, "ppi_hub_genes.pdf"), p_hub, width = 8, height = 7)
  ggsave(file.path(opt$figures_dir, "ppi_hub_genes.png"), p_hub, width = 8, height = 7, dpi = 300)
}

# ── Figura 3: rede dos hubs ────────────────────────────────────
hub_names <- hub_genes$gene_id
g_hub     <- induced_subgraph(g, vids = V(g)[V(g)$name %in% hub_names])

if (vcount(g_hub) >= 2) {
  reg_map  <- setNames(node_df$regulation, node_df$gene_id)
  vcolors  <- c("up" = "#E41A1C", "down" = "#377EB8", "ns" = "#AAAAAA")
  vcol     <- vcolors[reg_map[V(g_hub)$name]]
  vcol[is.na(vcol)] <- "#AAAAAA"

  pdf(file.path(opt$figures_dir, "ppi_hub_network.pdf"), width = 10, height = 10)
  par(mar = c(1, 1, 2, 1))
  set.seed(42)
  plot(g_hub,
       vertex.color = vcol,
       vertex.size  = log1p(degree(g_hub)) * 4 + 3,
       vertex.label = V(g_hub)$name,
       vertex.label.cex   = 0.65,
       vertex.label.color = "black",
       edge.width   = E(g_hub)$weight * 2,
       edge.color   = "grey70",
       layout       = layout_with_fr(g_hub),
       main         = "Hub Genes – STRING PPI (Arabidopsis thaliana)")
  legend("bottomleft", legend = c("Up","Down","Não-DEG"),
         fill = c("#E41A1C","#377EB8","#AAAAAA"), bty = "n", cex = 0.8)
  dev.off()
}

cat("Rede PPI concluída.\n")
