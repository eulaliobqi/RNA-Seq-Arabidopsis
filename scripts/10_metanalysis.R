#!/usr/bin/env Rscript
# ============================================================
# 10_metanalysis.R – Metanálise GEO/SRA para validação externa
# Baixa datasets públicos via GEOquery, calcula DEGs com limma
# e mede sobreposição com nossos DEGs (Fisher + Jaccard).
#
# Uso: fornecer --geo_ids "GSE12345,GSE67890" em params.yaml
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
  make_option("--geo_ids",     type = "character", default = ""),
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

# ── Inicializa outputs ─────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "metanalysis_overlap.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "metanalysis_validated_genes.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "metanalysis_summary.tsv"))
writeLines("Metanálise: inicializando...", file.path(opt$outdir, "metanalysis_report.txt"))
writeLines("", file.path(opt$figures_dir, ".meta_placeholder"))

# ── Leitura dos nossos DEGs ────────────────────────────────────
our_res <- read_tsv(opt$deseq2, show_col_types = FALSE) |>
  mutate(gene_id = gsub("\\.TAIR10$", "", gene_id)) |>
  filter(!is.na(padj))

our_degs <- our_res |>
  filter(padj < opt$padj, abs(log2FoldChange) > opt$lfc) |>
  pull(gene_id)

our_all  <- our_res$gene_id
cat(sprintf("Nossos DEGs: %d | Background: %d\n", length(our_degs), length(our_all)))

# ── Verifica se IDs GEO foram fornecidos ──────────────────────
geo_ids <- trimws(unlist(strsplit(opt$geo_ids, ",")))
geo_ids <- geo_ids[nchar(geo_ids) > 0]

if (length(geo_ids) == 0) {
  msg <- c(
    "Nenhum accession GEO fornecido.",
    "",
    "Para ativar a metanálise:",
    "  1. Identifique datasets relevantes em: https://www.ncbi.nlm.nih.gov/geo/",
    "     Sugestão de busca: 'Arabidopsis thaliana[organism] stress[title] RNA-seq'",
    "  2. Adicione em params.yaml:",
    "     geo_accessions: \"GSE12345,GSE67890\"",
    "  3. Re-execute: bash run_pipeline.sh local -resume",
    "",
    "Datasets de referência para Arabidopsis:",
    "  - GSE52979: resposta a estresse hídrico (RNA-seq)",
    "  - GSE114065: respostas a patógenos (RNA-seq)",
    "  - GSE80460: resposta a calor/frio (RNA-seq)"
  )
  writeLines(msg, file.path(opt$outdir, "metanalysis_report.txt"))
  cat(paste(msg, collapse="\n"), "\n")
  quit(save = "no", status = 0)
}

# ── Carrega pacotes para download GEO ─────────────────────────
if (!requireNamespace("GEOquery", quietly = TRUE) ||
    !requireNamespace("limma",    quietly = TRUE)) {
  msg <- "GEOquery ou limma não instalados. Execute: mamba install -n r-analysis bioconductor-geoquery bioconductor-limma"
  writeLines(msg, file.path(opt$outdir, "metanalysis_report.txt"))
  cat("AVISO:", msg, "\n")
  quit(save = "no", status = 0)
}
library(GEOquery)
library(limma)

# ── Processa cada dataset GEO ──────────────────────────────────
process_geo <- function(gse_id) {
  cat(sprintf("\nProcessando %s...\n", gse_id))

  gse <- tryCatch(
    getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE),
    error = function(e) { message(gse_id, " download falhou: ", e$message); NULL }
  )
  if (is.null(gse)) return(NULL)

  eset <- if (is.list(gse)) gse[[1]] else gse

  # Matriz de expressão
  expr <- exprs(eset)
  if (nrow(expr) == 0) { message(gse_id, ": matriz vazia"); return(NULL) }

  # Fenótipo
  pheno <- pData(eset)
  cat(sprintf("  Amostras: %d | Features: %d\n", ncol(expr), nrow(expr)))

  # Tenta identificar coluna de condição
  cond_col <- intersect(c("treatment","condition","source_name_ch1",
                           "characteristics_ch1","title"), colnames(pheno))[1]
  if (is.na(cond_col)) {
    message(gse_id, ": coluna de condição não identificada"); return(NULL)
  }

  conditions <- as.character(pheno[[cond_col]])
  unique_conds <- unique(conditions)

  if (length(unique_conds) < 2) {
    message(gse_id, ": apenas 1 condição detectada"); return(NULL)
  }

  cat(sprintf("  Condições: %s\n", paste(unique_conds, collapse=" | ")))

  # Limma DEG (contraste binário entre primeiras 2 condições)
  cond_factor  <- factor(ifelse(conditions == unique_conds[1], "A", "B"))
  design       <- model.matrix(~ 0 + cond_factor)
  colnames(design) <- c("A","B")
  contrast_mat <- makeContrasts(A - B, levels = design)

  fit <- tryCatch({
    fit1 <- lmFit(expr, design)
    fit2 <- contrasts.fit(fit1, contrast_mat)
    eBayes(fit2)
  }, error = function(e) { message("limma falhou: ", e$message); NULL })
  if (is.null(fit)) return(NULL)

  top <- topTable(fit, coef = 1, number = Inf, adjust.method = "BH") |>
    rownames_to_column("probe_id") |>
    filter(adj.P.Val < opt$padj, abs(logFC) > opt$lfc)

  cat(sprintf("  DEGs limma: %d\n", nrow(top)))

  list(
    gse_id   = gse_id,
    deg_ids  = top$probe_id,
    n_degs   = nrow(top),
    cond_col = cond_col,
    conds    = paste(unique_conds[1:2], collapse=" vs ")
  )
}

results <- lapply(geo_ids, process_geo)
results <- Filter(Negate(is.null), results)

if (length(results) == 0) {
  msg <- c("Nenhum dataset processado com sucesso.",
           "Verifique conectividade ou formato dos datasets GEO.")
  writeLines(msg, file.path(opt$outdir, "metanalysis_report.txt"))
  cat(paste(msg, collapse="\n"), "\n")
  quit(save = "no", status = 0)
}

# ── Sobreposição com nossos DEGs ───────────────────────────────
# Nota: IDs GEO podem ser probe IDs, não TAIR; sobreposição direta pode ser zero.
# Loga situação claramente.

overlap_df <- lapply(results, function(r) {
  geo_degs  <- r$deg_ids
  overlap   <- intersect(our_degs, geo_degs)
  n_overlap <- length(overlap)

  # Fisher exact test
  a <- n_overlap
  b <- length(our_degs)  - a
  c <- length(geo_degs)  - a
  d <- length(union(our_all, geo_degs)) - a - b - c
  ft <- tryCatch(
    fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater"),
    error = function(e) list(p.value=NA, estimate=NA)
  )

  # Jaccard
  jaccard <- ifelse(length(union(our_degs, geo_degs)) > 0,
                    n_overlap / length(union(our_degs, geo_degs)), 0)

  data.frame(
    gse_id    = r$gse_id,
    contraste = r$conds,
    n_geo_degs = r$n_degs,
    n_our_degs = length(our_degs),
    n_overlap  = n_overlap,
    jaccard    = round(jaccard, 4),
    OR         = round(as.numeric(ft$estimate), 3),
    pvalue     = ft$p.value,
    note       = if (n_overlap == 0)
                   "IDs GEO possivelmente probe IDs (não TAIR) – sem sobreposição direta"
                 else ""
  )
}) |> bind_rows() |>
  mutate(padj = p.adjust(pvalue, method = "BH"))

write_tsv(overlap_df, file.path(opt$outdir, "metanalysis_overlap.tsv"))

# Genes validados em ≥ 1 dataset
validated <- Reduce(union, lapply(results, function(r) intersect(our_degs, r$deg_ids)))
validated_df <- our_res |>
  filter(gene_id %in% validated) |>
  mutate(n_datasets_validated = sapply(gene_id, function(g)
    sum(sapply(results, function(r) g %in% r$deg_ids))))

write_tsv(validated_df, file.path(opt$outdir, "metanalysis_validated_genes.tsv"))
write_tsv(overlap_df, file.path(opt$outdir, "metanalysis_summary.tsv"))

# ── Relatório ──────────────────────────────────────────────────
report_lines <- c(
  sprintf("Metanálise GEO – Arabidopsis thaliana"),
  sprintf("Datasets analisados: %d", length(results)),
  sprintf("Nossos DEGs: %d", length(our_degs)),
  "",
  "Sobreposição por dataset:",
  sprintf("  %-12s alvos=%d overlap=%d Jaccard=%.4f OR=%.2f padj=%.4f",
          overlap_df$gse_id, overlap_df$n_geo_degs,
          overlap_df$n_overlap, overlap_df$jaccard,
          overlap_df$OR, replace(overlap_df$padj, is.na(overlap_df$padj), 1)),
  "",
  sprintf("Genes validados em ≥1 dataset: %d", length(validated))
)
writeLines(report_lines, file.path(opt$outdir, "metanalysis_report.txt"))
cat(paste(report_lines, collapse="\n"), "\n")

# ── Figura: sobreposição ───────────────────────────────────────
if (nrow(overlap_df) > 0 && any(!is.na(overlap_df$OR))) {
  p <- ggplot(overlap_df, aes(gse_id, n_overlap, fill = -log10(pvalue + 1e-300))) +
    geom_col() +
    geom_text(aes(label = n_overlap), vjust = -0.3, size = 4) +
    scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
    labs(title = "Sobreposição de DEGs com Datasets GEO",
         x = "Dataset GEO", y = "Genes em comum",
         fill = "-log10(p-value)") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$figures_dir, "metanalysis_overlap.pdf"), p, width = 7, height = 5)
  ggsave(file.path(opt$figures_dir, "metanalysis_overlap.png"), p, width = 7, height = 5, dpi = 300)
}

cat("Metanálise concluída.\n")
