#!/usr/bin/env Rscript
# ============================================================
# 08_plantfdb.R – Classificação de TFs via PlantTFDB v5
# Identifica TFs entre os DEGs, classifica por família e testa
# enriquecimento (Fisher) de famílias entre DEGs vs background.
#
# Fonte: http://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz
# Se --tf_file não for fornecido, tenta download automático.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

opt_list <- list(
  make_option("--deseq2",      type = "character"),
  make_option("--tf_file",     type = "character", default = NULL),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1.0),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("deseq2"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs vazios ──────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "tf_deg_classified.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "tf_family_summary.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "tf_family_enrichment.tsv"))
writeLines("PlantTFDB: inicializando...", file.path(opt$outdir, "plantfdb_summary.txt"))
writeLines("", file.path(opt$figures_dir, ".plantfdb_placeholder"))

# ── Leitura DESeq2 ─────────────────────────────────────────────
deseq2 <- read_tsv(opt$deseq2, show_col_types = FALSE) |>
  mutate(gene_id = gsub("\\.TAIR10$", "", gene_id)) |>
  filter(!is.na(padj))

all_genes <- deseq2$gene_id
degs      <- deseq2 |> filter(padj < opt$padj, abs(log2FoldChange) > opt$lfc)
cat(sprintf("Genes background: %d | DEGs: %d\n", length(all_genes), nrow(degs)))

# ── Carrega PlantTFDB ──────────────────────────────────────────
load_plantfdb <- function(tf_file) {
  # 1) arquivo fornecido pelo usuário
  if (!is.null(tf_file) && file.exists(tf_file)) {
    cat(sprintf("Usando arquivo PlantTFDB: %s\n", tf_file))
    return(read_tsv(tf_file, show_col_types = FALSE))
  }

  # 2) download automático
  urls <- c(
    "http://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz",
    "https://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz"
  )
  tmp_gz <- tempfile(fileext = ".txt.gz")
  tmp    <- sub("\\.gz$", "", tmp_gz)

  for (url in urls) {
    cat(sprintf("Tentando download: %s\n", url))
    ok <- tryCatch({
      download.file(url, tmp_gz, quiet = TRUE, timeout = 60)
      TRUE
    }, error = function(e) { message("  falhou: ", e$message); FALSE })

    if (ok) {
      tryCatch({
        # R.utils::gunzip ou conexão gzip nativa
        con <- gzcon(file(tmp_gz, "rb"))
        raw <- readLines(con)
        close(con)
        writeLines(raw, tmp)
        df <- read_tsv(tmp, show_col_types = FALSE, comment = "#")
        unlink(c(tmp_gz, tmp))
        cat(sprintf("PlantTFDB baixado: %d linhas, colunas: %s\n",
                    nrow(df), paste(colnames(df), collapse=", ")))
        return(df)
      }, error = function(e) message("  parse falhou: ", e$message))
    }
  }
  NULL
}

tf_db <- load_plantfdb(opt$tf_file)

if (is.null(tf_db) || nrow(tf_db) == 0) {
  msg <- paste(
    "PlantTFDB não carregado.",
    "Forneça --tf_file /path/Ath_TF_list.txt",
    "ou garanta acesso à internet no servidor.",
    "Download manual: http://planttfdb.gao-lab.org/download.php",
    sep = "\n"
  )
  writeLines(msg, file.path(opt$outdir, "plantfdb_summary.txt"))
  cat(sprintf("AVISO: %s\n", gsub("\n", " | ", msg)))
  quit(save = "no", status = 0)
}

# ── Normaliza colunas ──────────────────────────────────────────
# PlantTFDB v5 usa: Gene_ID, Family (possíveis variações de nome)
cols_lower <- tolower(colnames(tf_db))

find_col <- function(candidates) {
  for (c in candidates) {
    hit <- which(cols_lower == c)
    if (length(hit) > 0) return(colnames(tf_db)[hit[1]])
  }
  NULL
}

col_gene   <- find_col(c("gene_id", "geneid", "tair_id", "locus"))
col_family <- find_col(c("family", "tf_family", "tf family"))

if (is.null(col_gene) || is.null(col_family)) {
  cat(sprintf("AVISO: colunas Gene_ID/Family não encontradas. Colunas disponíveis: %s\n",
              paste(colnames(tf_db), collapse=", ")))
  # Tentativa: usa primeiras duas colunas
  if (ncol(tf_db) >= 2) {
    col_gene   <- colnames(tf_db)[1]
    col_family <- colnames(tf_db)[2]
    cat(sprintf("Usando colunas: %s (gene), %s (família)\n", col_gene, col_family))
  } else {
    writeLines("Formato PlantTFDB não reconhecido.",
               file.path(opt$outdir, "plantfdb_summary.txt"))
    quit(save = "no", status = 0)
  }
}

tf_clean <- tf_db |>
  select(gene_id = all_of(col_gene), family = all_of(col_family)) |>
  mutate(gene_id = gsub("\\.TAIR10$|\\.[0-9]+$", "", as.character(gene_id)),
         family  = as.character(family)) |>
  filter(!is.na(gene_id), !is.na(family), nchar(gene_id) > 0) |>
  distinct(gene_id, family)

cat(sprintf("TFs no PlantTFDB (Arabidopsis): %d em %d famílias\n",
            nrow(tf_clean), n_distinct(tf_clean$family)))

# ── Classifica DEGs ────────────────────────────────────────────
tf_degs <- degs |>
  inner_join(tf_clean, by = "gene_id") |>
  arrange(family, padj)

tf_all_bg <- all_genes[all_genes %in% tf_clean$gene_id]

cat(sprintf("TFs entre DEGs: %d | TFs no background: %d\n",
            nrow(tf_degs), length(tf_all_bg)))

if (nrow(tf_degs) == 0) {
  writeLines("Nenhum TF encontrado entre os DEGs.", file.path(opt$outdir, "plantfdb_summary.txt"))
  quit(save = "no", status = 0)
}

write_tsv(tf_degs, file.path(opt$outdir, "tf_deg_classified.tsv"))

# ── Contagens por família ──────────────────────────────────────
n_deg   <- nrow(degs)
n_bg    <- length(all_genes)

family_summary <- tf_clean |>
  group_by(family) |>
  summarise(
    n_tf_bg  = n_distinct(gene_id[gene_id %in% all_genes]),
    n_tf_deg = n_distinct(gene_id[gene_id %in% degs$gene_id]),
    .groups = "drop"
  ) |>
  filter(n_tf_bg > 0) |>
  mutate(
    pct_deg_in_family = round(n_tf_deg / n_tf_bg * 100, 1),
    pct_family_in_degs = round(n_tf_deg / n_deg * 100, 2)
  ) |>
  arrange(desc(n_tf_deg))

write_tsv(family_summary, file.path(opt$outdir, "tf_family_summary.tsv"))

# ── Teste de enriquecimento (Fisher por família) ───────────────
fisher_family <- function(fam, tf_db, all_genes, deg_genes) {
  tf_fam  <- tf_db$gene_id[tf_db$family == fam]
  a <- sum(deg_genes %in% tf_fam)               # DEG & TF família
  b <- sum(!(deg_genes %in% tf_fam))            # DEG & não-TF família
  c <- sum(all_genes %in% tf_fam) - a           # não-DEG & TF família
  d <- length(all_genes) - a - b - c            # não-DEG & não-TF família
  m <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(m, alternative = "greater")
  data.frame(family = fam, n_tf_deg = a, n_tf_bg = a + c,
             OR = round(ft$estimate, 3), pvalue = ft$p.value)
}

enrichment <- lapply(unique(tf_clean$family), fisher_family,
                     tf_db    = tf_clean,
                     all_genes = all_genes,
                     deg_genes = degs$gene_id) |>
  bind_rows() |>
  mutate(padj = p.adjust(pvalue, method = "BH")) |>
  filter(n_tf_deg > 0) |>
  arrange(padj)

write_tsv(enrichment, file.path(opt$outdir, "tf_family_enrichment.tsv"))

sig_families <- enrichment |> filter(padj < 0.05)
cat(sprintf("Famílias de TF enriquecidas (FDR<0.05): %d\n", nrow(sig_families)))

# ── Sumário ────────────────────────────────────────────────────
top_fams <- head(tf_degs |> count(family, sort = TRUE), 10)
summary_lines <- c(
  sprintf("TFs no PlantTFDB (background): %d", length(tf_all_bg)),
  sprintf("TFs entre DEGs: %d (%d famílias)",
          nrow(tf_degs), n_distinct(tf_degs$family)),
  sprintf("Famílias enriquecidas (FDR<0.05): %d", nrow(sig_families)),
  "",
  "Top 10 famílias nos DEGs:",
  sprintf("  %-20s %d TFs", top_fams$family, top_fams$n)
)
writeLines(summary_lines, file.path(opt$outdir, "plantfdb_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# ── Figura 1: barplot de famílias nos DEGs ─────────────────────
top20_fams <- head(family_summary |> filter(n_tf_deg > 0), 20)

p_bar <- ggplot(top20_fams, aes(reorder(family, n_tf_deg), n_tf_deg,
                                 fill = n_tf_deg)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n_tf_deg), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Famílias de TF entre os DEGs (PlantTFDB v5)",
       subtitle = sprintf("Total: %d TFs em %d famílias",
                          nrow(tf_degs), n_distinct(tf_degs$family)),
       x = NULL, y = "Número de TFs DEG") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$figures_dir, "tf_families_barplot.pdf"), p_bar, width = 8, height = 7)
ggsave(file.path(opt$figures_dir, "tf_families_barplot.png"), p_bar, width = 8, height = 7, dpi = 300)

# ── Figura 2: up vs down por família ──────────────────────────
tf_reg <- tf_degs |>
  mutate(regulation = case_when(
    log2FoldChange >  opt$lfc ~ "up",
    log2FoldChange < -opt$lfc ~ "down",
    TRUE ~ "ns"
  )) |>
  count(family, regulation) |>
  filter(family %in% top20_fams$family)

p_stack <- ggplot(tf_reg, aes(reorder(family, n, sum), n, fill = regulation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8")) +
  labs(title = "TFs DEG por família e direção",
       x = NULL, y = "Número de TFs", fill = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(opt$figures_dir, "tf_families_updown.pdf"), p_stack, width = 8, height = 7)
ggsave(file.path(opt$figures_dir, "tf_families_updown.png"), p_stack, width = 8, height = 7, dpi = 300)

# ── Figura 3: bubble plot de enriquecimento ────────────────────
enr_plot <- enrichment |>
  filter(n_tf_deg >= 2, padj < 0.2) |>
  head(25)

if (nrow(enr_plot) > 0) {
  p_enr <- ggplot(enr_plot,
                  aes(x = OR, y = reorder(family, -padj),
                      size = n_tf_deg, color = -log10(padj))) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_gradient(low = "#AED6F1", high = "#1B4F72") +
    labs(title = "Enriquecimento de famílias TF (Fisher, FDR<0.2)",
         x = "Odds Ratio", y = NULL,
         size = "TFs DEG", color = "-log10(FDR)") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$figures_dir, "tf_enrichment_bubble.pdf"), p_enr, width = 8, height = 6)
  ggsave(file.path(opt$figures_dir, "tf_enrichment_bubble.png"), p_enr, width = 8, height = 6, dpi = 300)
}

cat("PlantTFDB concluído.\n")
