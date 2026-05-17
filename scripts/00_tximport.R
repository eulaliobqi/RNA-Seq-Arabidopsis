#!/usr/bin/env Rscript
# ============================================================
# 00_tximport.R – Importação de quantificação Salmon via tximport
# Uso: Rscript 00_tximport.R --quant_dir . \
#                             --samplesheet sample_metadata.tsv \
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

# ── Lê metadados ─────────────────────────────────────────────
meta <- read_tsv(opt$samplesheet, show_col_types = FALSE)

# ── Localiza arquivos quant.sf por amostra ───────────────────
quant_files <- file.path(opt$quant_dir, meta$sample, "quant.sf")
names(quant_files) <- meta$sample

missing <- !file.exists(quant_files)
if (any(missing)) {
  stop(sprintf("quant.sf não encontrado para: %s",
               paste(meta$sample[missing], collapse = ", ")))
}

# ── tximport: agrega transcritos → genes ─────────────────────
# tx2gene: lê a coluna Name (transcript) dos quant.sf e extrai gene_id
# Usa gene-level summarization via countsFromAbundance="scaledTPM"
txi <- tximport(
  quant_files,
  type            = "salmon",
  txOut           = FALSE,               # agrega para gene-level
  ignoreTxVersion = TRUE,
  countsFromAbundance = "scaledTPM"
)

# ── Limpa IDs (remove sufixo .TAIR10 se presente) ────────────
rownames(txi$counts) <- gsub("\\.TAIR10$", "", rownames(txi$counts))
rownames(txi$abundance) <- gsub("\\.TAIR10$", "", rownames(txi$abundance))

# ── Salva counts (escalados) ──────────────────────────────────
counts_df <- as.data.frame(round(txi$counts)) |>
  rownames_to_column("gene_id")
write_tsv(counts_df, file.path(opt$outdir, "salmon_counts.tsv"))

# ── Salva TPM ─────────────────────────────────────────────────
tpm_df <- as.data.frame(txi$abundance) |>
  rownames_to_column("gene_id")
write_tsv(tpm_df, file.path(opt$outdir, "salmon_tpm.tsv"))

cat(sprintf(
  "tximport concluído: %d genes x %d amostras\n",
  nrow(counts_df), ncol(counts_df) - 1
))
