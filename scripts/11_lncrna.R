#!/usr/bin/env Rscript
# ============================================================
# 11_lncrna.R – Classificação de lncRNAs via análise de ORF
# Input: sequências FASTA de transcritos novos (StringTie/GffCompare)
# Critérios: comprimento >= 200 nt E maior ORF < 100 aa
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
})

opt_list <- list(
  make_option("--fasta",       type = "character"),
  make_option("--deseq2",      type = "character", default = NULL),
  make_option("--min_len",     type = "integer",   default = 200L),
  make_option("--max_orf_aa",  type = "integer",   default = 100L),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("fasta"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Inicializa outputs ─────────────────────────────────────────
write_tsv(data.frame(), file.path(opt$outdir, "lncrna_candidates.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "lncrna_all.tsv"))
writeLines("lncRNA: inicializando...", file.path(opt$outdir, "lncrna_summary.txt"))
writeLines("", file.path(opt$figures_dir, ".lncrna_placeholder"))

# ── Verifica arquivo FASTA ─────────────────────────────────────
if (!file.exists(opt$fasta) || file.size(opt$fasta) == 0) {
  msg <- c(
    "FASTA de transcritos novos não encontrado ou vazio.",
    "Certifique-se de que run_stringtie=true e o pipeline foi executado com StringTie."
  )
  writeLines(msg, file.path(opt$outdir, "lncrna_summary.txt"))
  cat(paste(msg, collapse = "\n"), "\n")
  quit(save = "no", status = 0)
}

# ── Leitura das sequências ─────────────────────────────────────
cat(sprintf("Lendo sequências: %s\n", opt$fasta))
seqs <- tryCatch(
  readDNAStringSet(opt$fasta),
  error = function(e) { message("Erro ao ler FASTA: ", e$message); NULL }
)
if (is.null(seqs) || length(seqs) == 0) {
  writeLines("Nenhuma sequência no FASTA.", file.path(opt$outdir, "lncrna_summary.txt"))
  quit(save = "no", status = 0)
}
cat(sprintf("Transcritos lidos: %d\n", length(seqs)))

# ── Detecta maior ORF em todas as 6 fases ─────────────────────
stop_codons <- c("TAA","TAG","TGA")

find_max_orf_aa <- function(dna_seq) {
  fwd  <- as.character(dna_seq)
  rev  <- as.character(reverseComplement(dna_seq))
  max_orf <- 0L
  for (seq_str in c(fwd, rev)) {
    for (frame in 0:2) {
      s     <- substring(seq_str, frame + 1)
      n_cod <- nchar(s) %/% 3
      if (n_cod < 2) next
      codons <- substring(s, seq(1, n_cod * 3 - 2, 3), seq(3, n_cod * 3, 3))
      start  <- which(codons == "ATG")
      stops  <- which(codons %in% stop_codons)
      for (st in start) {
        end_pos <- stops[stops > st][1]
        if (!is.na(end_pos)) {
          orf_len <- end_pos - st   # in codons (aa)
          if (orf_len > max_orf) max_orf <- orf_len
        }
      }
    }
  }
  max_orf
}

cat("Calculando ORFs (pode levar alguns minutos)...\n")
tx_names  <- names(seqs)
tx_lengths <- width(seqs)

# Filtra primeiro por comprimento mínimo para economizar tempo
long_idx <- which(tx_lengths >= opt$min_len)
cat(sprintf("Transcritos >= %d nt: %d / %d\n",
            opt$min_len, length(long_idx), length(seqs)))

max_orfs <- integer(length(seqs))
max_orfs[long_idx] <- vapply(
  seqs[long_idx], find_max_orf_aa, integer(1)
)

# ── Classifica transcritos ─────────────────────────────────────
tx_df <- tibble(
  transcript_id = tx_names,
  length_nt     = tx_lengths,
  max_orf_aa    = max_orfs,
  is_lncrna     = length_nt >= opt$min_len & max_orf_aa < opt$max_orf_aa,
  class         = case_when(
    length_nt < opt$min_len                              ~ "short (<200nt)",
    max_orf_aa >= opt$max_orf_aa                         ~ "putative_coding",
    TRUE                                                 ~ "lncrna_candidate"
  )
)

write_tsv(tx_df, file.path(opt$outdir, "lncrna_all.tsv"))

lncrna_df <- tx_df |> filter(is_lncrna) |> arrange(desc(length_nt))
write_tsv(lncrna_df, file.path(opt$outdir, "lncrna_candidates.tsv"))

cat(sprintf("Candidatos lncRNA (>=%dnt, ORF<%daa): %d / %d\n",
            opt$min_len, opt$max_orf_aa, nrow(lncrna_df), length(seqs)))

# ── Integra com DEGs (se disponível) ──────────────────────────
deseq2_info <- ""
if (!is.null(opt$deseq2) && file.exists(opt$deseq2) &&
    file.size(opt$deseq2) > 0) {
  deseq2 <- tryCatch(
    read_tsv(opt$deseq2, show_col_types = FALSE) |>
      mutate(gene_id = gsub("\\.TAIR10$", "", gene_id)),
    error = function(e) NULL
  )
  if (!is.null(deseq2)) {
    # Tenta match por gene ID (MSTRG IDs não vão coincidir com TAIR IDs diretamente)
    # Extrai possível gene base do transcript ID (MSTRG.X → MSTRG)
    lncrna_df <- lncrna_df |>
      mutate(gene_base = gsub("\\.[0-9]+$", "", transcript_id))
    deseq2_info <- sprintf(
      "\nNota: IDs de transcritos novos (MSTRG.X) não coincidem diretamente com TAIR IDs.\n%d lncRNA candidatos identificados para análise manual / anotação downstream.",
      nrow(lncrna_df)
    )
  }
}

# ── Sumário ────────────────────────────────────────────────────
class_counts <- table(tx_df$class)
summary_lines <- c(
  sprintf("lncRNA Prediction – StringTie Novel Transcripts"),
  sprintf("Comprimento mínimo: %d nt | ORF máximo: %d aa", opt$min_len, opt$max_orf_aa),
  sprintf("Total de transcritos novos: %d", nrow(tx_df)),
  "",
  "Classificação:",
  sprintf("  %-25s %d", names(class_counts), as.integer(class_counts)),
  "",
  sprintf("Candidatos lncRNA: %d (%.1f%%)",
          nrow(lncrna_df), nrow(lncrna_df) / nrow(tx_df) * 100),
  deseq2_info
)
writeLines(summary_lines, file.path(opt$outdir, "lncrna_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# ── Figura 1: comprimento vs maior ORF ────────────────────────
p_scatter <- ggplot(tx_df |> filter(length_nt >= opt$min_len),
                    aes(length_nt, max_orf_aa, color = class)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = opt$max_orf_aa, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("lncrna_candidate" = "#2ECC71",
                                 "putative_coding"  = "#E74C3C",
                                 "short (<200nt)"   = "#95A5A6")) +
  scale_x_log10() +
  labs(title = "Classificação de Transcritos Novos (StringTie)",
       x = "Comprimento do transcrito (nt, log10)",
       y = "Maior ORF (aa)",
       color = "Classe") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$figures_dir, "lncrna_scatter.pdf"), p_scatter, width = 8, height = 6)
ggsave(file.path(opt$figures_dir, "lncrna_scatter.png"), p_scatter, width = 8, height = 6, dpi = 300)

# ── Figura 2: distribuição de comprimentos dos candidatos ──────
if (nrow(lncrna_df) > 0) {
  p_len <- ggplot(lncrna_df, aes(length_nt)) +
    geom_histogram(bins = 40, fill = "#2ECC71", color = "white") +
    labs(title = sprintf("Distribuição de Comprimento – %d lncRNA Candidatos",
                         nrow(lncrna_df)),
         x = "Comprimento (nt)", y = "Frequência") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$figures_dir, "lncrna_length_dist.pdf"), p_len, width = 7, height = 5)
  ggsave(file.path(opt$figures_dir, "lncrna_length_dist.png"), p_len, width = 7, height = 5, dpi = 300)
}

# ── Figura 3: barplot de classes ───────────────────────────────
class_df <- tx_df |> count(class)
p_class  <- ggplot(class_df, aes(reorder(class, n), n, fill = class)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.2, size = 4) +
  coord_flip() +
  scale_fill_manual(values = c("lncrna_candidate" = "#2ECC71",
                                "putative_coding"  = "#E74C3C",
                                "short (<200nt)"   = "#95A5A6")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Classificação de Transcritos Novos", x = NULL, y = "n") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$figures_dir, "lncrna_class_barplot.pdf"), p_class, width = 7, height = 4)
ggsave(file.path(opt$figures_dir, "lncrna_class_barplot.png"), p_class, width = 7, height = 4, dpi = 300)

cat("lncRNA prediction concluída.\n")
