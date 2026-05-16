# RNASeq Insight – Arabidopsis thaliana TAIR10

Pipeline RNA-Seq completo para *Arabidopsis thaliana*, adaptado do projeto *Glycine max* (soja).

## Stack tecnológico
- **Orquestração**: Nextflow DSL2 ≥ 24.04 (env `rnaseq-tools`)
- **QC**: FastQC 0.12.1 + MultiQC 1.21
- **Trimagem**: fastp 0.23.4
- **Alinhamento**: HISAT2 2.2.1 (splice-aware)
- **Quantificação**: featureCounts (subread 2.0.6)
- **Splicing**: rMATS ≥ 4.1.0
- **Análise R**: DESeq2, clusterProfiler, WGCNA, Shiny, Quarto (env `r-analysis`)

## Parâmetros organismo-específicos (A. thaliana)
- **KEGG**: `ath`
- **Org.db**: `org.At.tair.db` (Bioconductor 3.20+)
- **biomaRt dataset**: `"athaliana_eg_gene"` (host: plants.ensembl.org)
- **IDs TAIR10**: formato `AT[1-5MC]G[0-9]{5}` — **não precisam de limpeza de sufixos**
- **Cache GO**: `at_go_cache.rds` (se implementado)

## Critérios de análise
- FDR < 0.05, |log2FC| > 1.0
- Splicing: FDR < 0.05, |ΔPSI| > 0.1
- Integration score: `lfc×3 + mean×2 + sig×2 + splicing×1.5 + pathway×1.0 + hub×0.5`
- Key candidates: genes com evidência em ≥ 2 camadas

## Estrutura de arquivos
```
main.nf                    # Pipeline principal
nextflow.config            # Configuração técnica (recursos, profiles)
params.yaml                # Parâmetros biológicos (edite para cada experimento)
samplesheet.csv            # Amostras (edite com paths reais)
modules/                   # Módulos Nextflow DSL2
scripts/                   # Scripts R (01–04)
dashboard/app.R            # Dashboard Shiny interativo
envs/                      # Ambientes Conda
report/rnaseq_report.qmd   # Template Quarto
```

## Agentes especializados
- `/validate-qc` — valida QC das amostras
- `/interpret-results` — interpretação biológica dos resultados
- `/debug-pipeline` — diagnóstico de falhas no Nextflow

## Executar pipeline
```bash
# 1. Setup (uma vez)
bash setup.sh

# 2. Execute
bash run_pipeline.sh local

# 3. Retomar após falha
bash run_pipeline.sh local --resume

# 4. Dashboard
RESULTS_DIR=results Rscript -e "shiny::runApp('dashboard/app.R', port=3838)"
```

## Genoma TAIR10
- Baixar de: https://www.arabidopsis.org/download
- `TAIR10_chr_all.fas` + `TAIR10_GFF3_genes.gff`
- Atualizar paths em `params.yaml`

## Bugs conhecidos e soluções (herdados do projeto soja)
- `org.Gmax.eg.db` removido Bioconductor 3.20+ → para Arabidopsis usar `org.At.tair.db` (disponível)
- rMATS exige ≥ 2 BAMs por grupo → validado no `main.nf`
- featureCounts nomeia colunas com path BAM → `PARSE_COUNTS` renomeia para sample names
- LFC shrinkage: fallback automático `apeglm → ashr → normal`

## GitHub
https://github.com/eulaliobqi/RNA-Seq-Arabidopsis
