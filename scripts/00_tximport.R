#!/usr/bin/env Rscript
# ============================================================
# 00_tximport.R – Importação de quantificação Salmon via tximport
# Uso: Rscript 00_tximport.R --quant_dir . \
#                             --samplesheet samplesheet.csv \
#                             --outdir .
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(readr)
  library(dplyr)
  library(tibble)
})

opt_list <- list(
  make_option("--quant_dir",   type = "character"),
  make_option("--samplesheet", type = "character"),
  make_option("--outdir",      type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("quant_dir", "samplesheet"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Lê samplesheet – suporta CSV e TSV ───────────────────────
ext  <- tolower(tools::file_ext(opt$samplesheet))
meta <- if (ext == "csv") {
  read_csv(opt$samplesheet, show_col_types = FALSE)
} else {
  read_tsv(opt$samplesheet, show_col_types = FALSE)
}

if (!"sample" %in% colnames(meta))
  stop(sprintf("Coluna 'sample' não encontrada. Colunas: %s",
               paste(colnames(meta), collapse = ", ")))

samples <- meta$sample

# ── Localiza quant.sf: procura pelo path direto, depois por glob ──
# O Nextflow pode renomear diretórios staged (ex: ColWT_rep1_1/) em
# caso de colisão, então buscamos qualquer subdir com quant.sf
find_quant_sf <- function(base_dir, sample_name) {
  # Primeiro: caminho esperado
  direct <- file.path(base_dir, sample_name, "quant.sf")
  if (file.exists(direct)) return(direct)

  # Fallback: busca em qualquer subdir que contenha o nome da amostra
  candidates <- list.files(
    base_dir, pattern = "quant.sf",
    recursive = TRUE, full.names = TRUE
  )
  # Filtra pelo nome da amostra no path
  match <- grep(sample_name, candidates, value = TRUE, fixed = TRUE)
  if (length(match) > 0) return(match[1])

  return(NA_character_)
}

cat(sprintf("Buscando quant.sf em: %s\n", normalizePath(opt$quant_dir)))
cat(sprintf("Subdirs disponíveis: %s\n",
            paste(list.dirs(opt$quant_dir, recursive = FALSE, full.names = FALSE),
                  collapse = ", ")))

quant_files <- vapply(samples, find_quant_sf,
                      base_dir = opt$quant_dir,
                      FUN.VALUE = character(1))
names(quant_files) <- samples

for (i in seq_along(quant_files)) {
  status <- if (!is.na(quant_files[i]) && file.exists(quant_files[i])) "OK"
            else "NAO ENCONTRADO"
  cat(sprintf("  %-20s → %s [%s]\n", names(quant_files)[i],
              quant_files[i], status))
}

missing <- is.na(quant_files) | !file.exists(quant_files)
if (any(missing))
  stop(sprintf("quant.sf não encontrado para: %s",
               paste(samples[missing], collapse = ", ")))

# ── tximport: agrega transcritos → genes ─────────────────────
txi <- tximport(
  quant_files,
  type                = "salmon",
  txOut               = FALSE,
  ignoreTxVersion     = TRUE,
  countsFromAbundance = "scaledTPM"
)

# ── Limpa IDs (.TAIR10) ───────────────────────────────────────
rownames(txi$counts)    <- gsub("\\.TAIR10$", "", rownames(txi$counts))
rownames(txi$abundance) <- gsub("\\.TAIR10$", "", rownames(txi$abundance))

# ── Salva outputs ─────────────────────────────────────────────
counts_df <- as.data.frame(round(txi$counts)) |> rownames_to_column("gene_id")
tpm_df    <- as.data.frame(txi$abundance)     |> rownames_to_column("gene_id")

write_tsv(counts_df, file.path(opt$outdir, "salmon_counts.tsv"))
write_tsv(tpm_df,    file.path(opt$outdir, "salmon_tpm.tsv"))

cat(sprintf("tximport concluído: %d genes x %d amostras\n",
            nrow(counts_df), ncol(counts_df) - 1))
