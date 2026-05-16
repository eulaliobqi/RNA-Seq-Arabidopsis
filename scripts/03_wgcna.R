#!/usr/bin/env Rscript
# ============================================================
# 03_wgcna.R – Co-expressão gênica (WGCNA)
# Detecção automática de soft power + módulos + hub genes
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(readr)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

opt_list <- list(
  make_option("--norm_counts", type = "character"),
  make_option("--metadata",    type = "character"),
  make_option("--min_genes",   type = "integer",   default = 5000L),
  make_option("--soft_power",  type = "integer",   default = 0L),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--figures_dir", type = "character", default = "figures")
)
opt <- parse_args(OptionParser(option_list = opt_list))
for (p in c("norm_counts","metadata"))
  if (is.null(opt[[p]])) stop(sprintf("Parâmetro obrigatório: --%s", p))

dir.create(opt$outdir,      showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── Leitura ──────────────────────────────────────────────────
expr <- read_tsv(opt$norm_counts, show_col_types = FALSE) |>
        column_to_rownames("gene_id") |>
        as.matrix()
meta <- read_tsv(opt$metadata, show_col_types = FALSE) |> as.data.frame()

# Transpõe: amostras nas linhas, genes nas colunas
datExpr <- t(expr)

# ── Seleciona genes com maior variância ───────────────────────
n_genes <- min(opt$min_genes, ncol(datExpr))
vars     <- apply(datExpr, 2, var)
datExpr  <- datExpr[, order(vars, decreasing = TRUE)[1:n_genes]]
cat(sprintf("Genes para WGCNA: %d\n", ncol(datExpr)))

# ── Verifica amostras boas ────────────────────────────────────
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("Após goodSamplesGenes: %d amostras, %d genes\n",
              nrow(datExpr), ncol(datExpr)))
}

# ── Clustering hierárquico de amostras ───────────────────────
sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file.path(opt$figures_dir, "sample_clustering_tree.pdf"), width = 10, height = 5)
par(cex = 0.8, mar = c(0,4,2,0))
plot(sampleTree, main = "Clustering Hierárquico de Amostras",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1.2)
dev.off()

# ── Seleção do soft-thresholding power ───────────────────────
if (opt$soft_power == 0) {
  powers <- c(1:10, seq(12, 20, 2))
  sft    <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

  pdf(file.path(opt$figures_dir, "soft_threshold.pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "R² Scale Free Topology",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = 0.9, col = "red")
  abline(h = 0.8, col = "red", lty = 2)
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = "Mean connectivity")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
  dev.off()

  # Auto-seleciona power com R² > 0.8
  chosen <- sft$fitIndices |>
    filter(-sign(V3) * V2 >= 0.8) |>
    slice(1) |>
    pull(V1)
  soft_power <- if (length(chosen) == 0) 6 else chosen
  cat(sprintf("Soft power selecionado: %d\n", soft_power))
} else {
  soft_power <- opt$soft_power
}

# ── Constrói rede ─────────────────────────────────────────────
net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  TOMType           = "unsigned",
  minModuleSize     = 30,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  verbose           = 0
)
cat(sprintf("Módulos detectados: %s\n", paste(table(net$colors), collapse = ", ")))

# ── Dendrograma de genes + módulos ───────────────────────────
pdf(file.path(opt$figures_dir, "gene_dendrogram_modules.pdf"), width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Dendrograma de genes e módulos")
dev.off()

# ── Eigengenes e correlação com traits ────────────────────────
MEs0   <- moduleEigengenes(datExpr, net$colors)$eigengenes
MEs    <- orderMEs(MEs0)
traits <- as.data.frame(as.numeric(meta$condition == unique(meta$condition)[2]))
rownames(traits) <- meta$sample
colnames(traits) <- "treatment"

moduleTraitCor    <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

pdf(file.path(opt$figures_dir, "module_trait_heatmap.pdf"), width = 8, height = 10)
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")")
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix       = moduleTraitCor,
               xLabels      = colnames(traits),
               yLabels      = rownames(moduleTraitCor),
               ySymbols     = rownames(moduleTraitCor),
               colorLabels  = FALSE,
               colors       = blueWhiteRed(50),
               textMatrix   = textMatrix,
               setStdMargins = FALSE,
               cex.text     = 0.7,
               zlim         = c(-1, 1),
               main         = "Correlação Módulo–Trait")
dev.off()

# ── Outputs TSV ──────────────────────────────────────────────
modules_df <- data.frame(gene_id = colnames(datExpr), module = net$colors)
write_tsv(modules_df, file.path(opt$outdir, "wgcna_modules.tsv"))

# Hub genes: top 10 por módulo (maior kWithin)
adj     <- adjacency(datExpr, power = soft_power)
kIM     <- intramodularConnectivity(adj, net$colors)
hub_df  <- kIM |>
  tibble::rownames_to_column("gene_id") |>
  mutate(module = net$colors[gene_id]) |>
  group_by(module) |>
  slice_max(kWithin, n = 10) |>
  ungroup()
write_tsv(hub_df, file.path(opt$outdir, "wgcna_hub_genes.tsv"))

eigen_df <- MEs |> tibble::rownames_to_column("sample")
write_tsv(eigen_df, file.path(opt$outdir, "wgcna_eigengenes.tsv"))

summary_df <- modules_df |>
  count(module, name = "n_genes") |>
  left_join(
    data.frame(module = sub("ME","",colnames(MEs)),
               trait_cor = moduleTraitCor[,1],
               trait_pval = moduleTraitPvalue[,1]),
    by = "module"
  )
write_tsv(summary_df, file.path(opt$outdir, "wgcna_module_summary.tsv"))

cat("WGCNA concluído.\n")
