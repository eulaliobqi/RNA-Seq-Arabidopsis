#!/usr/bin/env Rscript
# ============================================================
# install_r_bioc_packages.R – Arabidopsis thaliana
# Execute após criar o ambiente r-analysis:
#   mamba run -n r-analysis Rscript scripts/install_r_bioc_packages.R
# ============================================================

cat("Instalando pacotes Bioconductor para Arabidopsis thaliana...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c(
  "WGCNA",
  "org.At.tair.db",   # Arabidopsis annotation – disponível em Bioconductor 3.20+
  "GO.db",
  "impute",
  "preprocessCore"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Instalando %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("  ✓ %s já instalado\n", pkg))
  }
}

cat("\nValidação:\n")
for (pkg in bioc_pkgs) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %s %s\n", ifelse(ok, "✓", "✗"), pkg))
}

cat("\nConcluído.\n")
