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

# ── Lê samplesheet (CSV ou TSV por extensão) ──────────────────
ext  <- tolower(tools::file_ext(opt$samplesheet))
meta <- if (ext == "csv") read_csv(opt$samplesheet, show_col_types = FALSE) else
                          read_tsv(opt$samplesheet, show_col_types = FALSE)

if (!"sample" %in% colnames(meta))
  stop(sprintf("Coluna 'sample' não encontrada. Colunas: %s",
               paste(colnames(meta), collapse = ", ")))

samples <- meta$sample

# ── Localiza quant.sf por amostra (com fallback por glob) ────
find_quant_sf <- function(base_dir, sample_name) {
  direct <- file.path(base_dir, sample_name, "quant.sf")
  if (file.exists(direct)) return(direct)
  candidates <- list.files(base_dir, pattern = "quant.sf",
                            recursive = TRUE, full.names = TRUE)
  match <- grep(sample_name, candidates, value = TRUE, fixed = TRUE)
  if (length(match) > 0) return(match[1])
  return(NA_character_)
}

cat(sprintf("Buscando quant.sf em: %s\n", normalizePath(opt$quant_dir)))
cat(sprintf("Subdirs: %s\n",
    paste(list.dirs(opt$quant_dir, recursive = FALSE, full.names = FALSE),
          collapse = ", ")))

quant_files <- vapply(samples, find_quant_sf,
                      base_dir = opt$quant_dir, FUN.VALUE = character(1))
names(quant_files) <- samples

for (i in seq_along(quant_files)) {
  status <- if (!is.na(quant_files[i]) && file.exists(quant_files[i])) "OK"
            else "NAO ENCONTRADO"
  cat(sprintf("  %-20s → %s [%s]\n", names(quant_files)[i],
              quant_files[i], status))
}

missing <- is.na(quant_files) | !file.exists(quant_files)
if (any(missing))
  stop(sprintf("quant.sf não encontrado: %s",
               paste(samples[missing], collapse = ", ")))

# ── Constrói tx2gene a partir dos IDs do primeiro quant.sf ───
# tximport com txOut=FALSE exige tx2gene para agregar tx → gene.
# IDs TAIR10: AT1G01010.1 → gene AT1G01010
#              AT1G01010.TAIR10.1 → gene AT1G01010 (remove .N e .TAIR10)
first_quant <- read_tsv(quant_files[1], show_col_types = FALSE)
tx_ids      <- first_quant$Name
# Remove sufixo de versão (.1, .2, ...)
gene_ids    <- sub("\\.[0-9]+$", "", tx_ids)
# Remove sufixo .TAIR10 que gffread pode adicionar
gene_ids    <- gsub("\\.TAIR10$", "", gene_ids)
tx2gene     <- data.frame(tx = tx_ids, gene = gene_ids)

cat(sprintf("tx2gene: %d transcritos → %d genes únicos\n",
            nrow(tx2gene), length(unique(tx2gene$gene))))
cat(sprintf("Exemplo: %s → %s\n", tx2gene$tx[1], tx2gene$gene[1]))

# ── tximport: agrega transcritos → genes ─────────────────────
txi <- tximport(
  quant_files,
  type                = "salmon",
  tx2gene             = tx2gene,
  txOut               = FALSE,
  ignoreTxVersion     = FALSE,  # gerenciado manualmente no tx2gene
  countsFromAbundance = "scaledTPM"
)

# Remove sufixo .TAIR10 residual nos rownames (segurança extra)
rownames(txi$counts)    <- gsub("\\.TAIR10$", "", rownames(txi$counts))
rownames(txi$abundance) <- gsub("\\.TAIR10$", "", rownames(txi$abundance))

# ── Salva outputs ─────────────────────────────────────────────
counts_df <- as.data.frame(round(txi$counts)) |> rownames_to_column("gene_id")
tpm_df    <- as.data.frame(txi$abundance)     |> rownames_to_column("gene_id")

write_tsv(counts_df, file.path(opt$outdir, "salmon_counts.tsv"))
write_tsv(tpm_df,    file.path(opt$outdir, "salmon_tpm.tsv"))

cat(sprintf("tximport concluído: %d genes x %d amostras\n",
            nrow(counts_df), ncol(counts_df) - 1))
