#!/usr/bin/env Rscript
# ============================================================
# 06_machine_learning.R – Classificação e seleção de biomarcadores
# Modelos: Random Forest | SVM-RBF | ElasticNet
# CV: repeatedcv (k≤5, 3 repeats) | Métrica: AUC (ROC)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(caret)
  library(randomForest)
  library(kernlab)
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(patchwork)
})

opt_list <- list(
  make_option("--norm_counts", type = "character"),
  make_option("--deseq2",      type = "character"),
  make_option("--metadata",    type = "character"),
  make_option("--n_features",  type = "integer",   default = 500L),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1.0),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))

for (p in c("norm_counts", "deseq2", "metadata"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Leitura ────────────────────────────────────────────────────
counts_df <- read_tsv(opt$norm_counts, show_col_types = FALSE)
deseq2    <- read_tsv(opt$deseq2,      show_col_types = FALSE)
meta      <- read_tsv(opt$metadata,    show_col_types = FALSE) |> as.data.frame()
rownames(meta) <- meta$sample

# ── Verifica mínimo de amostras ────────────────────────────────
if (nrow(meta) < 4) {
  cat("AVISO: Menos de 4 amostras – ML não executado.\n")
  write_tsv(data.frame(), file.path(opt$outdir, "ml_results.tsv"))
  write_tsv(data.frame(), file.path(opt$outdir, "feature_importance.tsv"))
  quit(save = "no", status = 0)
}

# ── Seleciona features (top DEGs por padj) ─────────────────────
degs <- deseq2 |>
  filter(!is.na(padj), padj < opt$padj, abs(log2FoldChange) > opt$lfc) |>
  arrange(padj)

if (nrow(degs) == 0) {
  cat("AVISO: Sem DEGs. Usando top 200 por padj.\n")
  degs <- deseq2 |> filter(!is.na(padj)) |> arrange(padj) |> head(200)
}

feature_genes <- head(degs$gene_id, opt$n_features)
cat(sprintf("Features: %d genes | DEGs totais: %d\n", length(feature_genes), nrow(degs)))

# ── Monta matriz ────────────────────────────────────────────────
counts_mat <- counts_df |>
  filter(gene_id %in% feature_genes) |>
  column_to_rownames("gene_id") |>
  as.matrix()

common <- intersect(colnames(counts_mat), meta$sample)
counts_mat <- counts_mat[, common, drop = FALSE]
meta       <- meta[common, ]

X <- t(counts_mat) |> as.data.frame()
X <- X[, apply(X, 2, var) > 0, drop = FALSE]   # remove variância zero

y       <- factor(meta$condition)
classes <- levels(y)
cat(sprintf("Classes: %s (n=%d) vs %s (n=%d)\n",
            classes[1], sum(y == classes[1]),
            classes[2], sum(y == classes[2])))

# ── Cross-validation ───────────────────────────────────────────
k_folds <- min(5L, min(table(y)))
if (k_folds < 2) {
  cat("AVISO: Amostras insuficientes para CV.\n")
  quit(save = "no", status = 0)
}

ctrl <- trainControl(
  method          = "repeatedcv",
  number          = k_folds,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(42)

# ── Treinamento ────────────────────────────────────────────────
safe_train <- function(method, grid = NULL) {
  args <- list(x = X, y = y, method = method, metric = "ROC",
               trControl = ctrl, preProcess = c("center", "scale"))
  if (!is.null(grid)) args$tuneGrid <- grid
  tryCatch(do.call(train, args), error = function(e) {
    message(sprintf("%s falhou: %s", method, e$message)); NULL
  })
}

cat("Treinando Random Forest...\n")
m_rf  <- safe_train("rf",
           expand.grid(mtry = unique(c(2L, floor(sqrt(ncol(X))), floor(ncol(X)/3L)))))

cat("Treinando SVM-RBF...\n")
m_svm <- safe_train("svmRadial")

cat("Treinando ElasticNet...\n")
m_en  <- safe_train("glmnet",
           expand.grid(alpha  = c(0, 0.5, 1),
                       lambda = 10^seq(-4, 0, length.out = 8)))

models <- Filter(Negate(is.null),
                 list(RandomForest = m_rf, SVM_RBF = m_svm, ElasticNet = m_en))

if (length(models) == 0) stop("Todos os modelos falharam.")

# ── AUC de treino ──────────────────────────────────────────────
auc_df <- lapply(names(models), function(nm) {
  perf <- getTrainPerf(models[[nm]])
  data.frame(
    model = nm,
    AUC   = round(perf$TrainROC,  4),
    Sens  = round(perf$TrainSens, 4),
    Spec  = round(perf$TrainSpec, 4),
    best_params = paste(names(models[[nm]]$bestTune),
                        unlist(models[[nm]]$bestTune),
                        sep = "=", collapse = "; ")
  )
}) |> bind_rows()

write_tsv(auc_df, file.path(opt$outdir, "ml_results.tsv"))
cat("\n── Resultados (CV) ──\n")
print(auc_df, row.names = FALSE)

# ── Importância de variáveis ───────────────────────────────────
get_imp <- function(m, nm) {
  tryCatch({
    vi  <- varImp(m, scale = TRUE)$importance
    col <- if ("Overall" %in% colnames(vi)) "Overall" else colnames(vi)[1]
    data.frame(gene_id = rownames(vi), importance = vi[[col]], model = nm) |>
      arrange(desc(importance)) |>
      head(50)
  }, error = function(e) NULL)
}

imp_list <- Filter(Negate(is.null),
                   mapply(get_imp, models, names(models), SIMPLIFY = FALSE))

if (length(imp_list) > 0) {
  write_tsv(bind_rows(imp_list), file.path(opt$outdir, "feature_importance.tsv"))

  best_nm  <- auc_df$model[which.max(auc_df$AUC)]
  top20    <- imp_list[[best_nm]]
  if (!is.null(top20) && nrow(top20) > 0) {
    top20 <- head(top20, 20)
    p_imp <- ggplot(top20, aes(reorder(gene_id, importance), importance)) +
      geom_col(fill = "#2166AC") +
      coord_flip() +
      labs(title = sprintf("Top 20 Genes Preditivos – %s (AUC=%.3f)",
                           best_nm, auc_df$AUC[auc_df$model == best_nm]),
           x = NULL, y = "Importância (normalizada)") +
      theme_bw(base_size = 11)
    ggsave(file.path(opt$figures_dir, "feature_importance.pdf"), p_imp, width = 8, height = 6)
    ggsave(file.path(opt$figures_dir, "feature_importance.png"), p_imp, width = 8, height = 6, dpi = 300)
  }
}

# ── Curvas ROC ────────────────────────────────────────────────
filter_best <- function(m) {
  pred <- m$pred
  if (is.null(pred)) return(NULL)
  for (col in colnames(m$bestTune))
    if (col %in% colnames(pred))
      pred <- pred[abs(as.numeric(pred[[col]]) - as.numeric(m$bestTune[[col]])) < 1e-9, ]
  pred
}

pos_class <- classes[1]
roc_list  <- lapply(names(models), function(nm) {
  pred <- filter_best(models[[nm]])
  if (is.null(pred) || !pos_class %in% colnames(pred)) return(NULL)
  tryCatch({
    r <- roc(pred$obs, pred[[pos_class]], quiet = TRUE)
    data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities,
               model = sprintf("%s (AUC=%.3f)", nm, as.numeric(auc(r))))
  }, error = function(e) NULL)
})
roc_list <- Filter(Negate(is.null), roc_list)

if (length(roc_list) > 0) {
  roc_df <- bind_rows(roc_list)
  p_roc <- ggplot(roc_df, aes(FPR, TPR, color = model)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Curvas ROC – Classificação de Condição",
         x = "1 - Especificidade", y = "Sensibilidade", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(opt$figures_dir, "roc_curves.pdf"), p_roc, width = 7, height = 6)
  ggsave(file.path(opt$figures_dir, "roc_curves.png"), p_roc, width = 7, height = 6, dpi = 300)
}

cat(sprintf("\nML concluído. Melhor modelo: %s (AUC=%.3f)\n",
            auc_df$model[which.max(auc_df$AUC)],
            max(auc_df$AUC, na.rm = TRUE)))
