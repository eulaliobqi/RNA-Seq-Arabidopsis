# RNASeq Insight Platform – *Arabidopsis thaliana* TAIR10

Pipeline automatizado e reprodutível para análise completa de RNA-Seq em *Arabidopsis thaliana*, incluindo expressão diferencial, splicing alternativo, enriquecimento funcional, co-expressão gênica, integração multi-ômica e dashboard interativo.

---

## Visão Geral

```
FASTQ (paired-end)
  │
  ├─ FastQC (pré) ─────────────────────────────────── MultiQC
  │
  ├─ fastp (trimagem)
  │
  ├─ FastQC (pós) ─────────────────────────────────── MultiQC
  │
  ├─ HISAT2 (alinhamento splice-aware) + samtools
  │
  ├─ featureCounts ──── DESeq2 ──── GO/KEGG/GSEA (org.At.tair.db)
  │                        │
  │                      WGCNA (co-expressão)
  │
  └─ rMATS (splicing) ────────── Integração Multi-ômica
                                        │
                               Dashboard Shiny + Relatório Quarto
```

---

## Funcionalidades

| Módulo | Ferramenta | Output |
|---|---|---|
| QC | FastQC + MultiQC | Relatórios HTML |
| Trimagem | fastp | Reads limpos + métricas |
| Alinhamento | HISAT2 | BAM ordenados e indexados |
| Quantificação | featureCounts | Matriz de contagens |
| Expressão diferencial | DESeq2 | Tabelas + figuras |
| Enriquecimento | clusterProfiler + org.At.tair.db | GO-BP/MF/CC, KEGG, GSEA |
| Splicing alternativo | rMATS | SE, A5SS, A3SS, MXE, RI |
| Co-expressão | WGCNA | Módulos + hub genes |
| Integração | R | Gene ranking + candidatos |
| Dashboard | Shiny + Plotly | App interativo |
| Relatório | Quarto | HTML completo |

---

## Requisitos

- Linux (Ubuntu 20.04+ ou servidor HPC)
- [Miniforge / Mamba](https://github.com/conda-forge/miniforge)
- Java ≥ 11 (para Nextflow)

---

## Instalação

```bash
# Clone o repositório
git clone https://github.com/eulaliobqi/RNA-Seq-Arabidopsis.git
cd RNA-Seq-Arabidopsis

# Instala ambientes conda e valida ferramentas (≈ 20–40 min)
bash setup.sh
```

---

## Uso

### 1. Configure o genoma

Baixe o genoma TAIR10 de [arabidopsis.org](https://www.arabidopsis.org/download):
- `TAIR10_chr_all.fas`
- `TAIR10_GFF3_genes.gff`

### 2. Configure suas amostras

Edite [samplesheet.csv](samplesheet.csv):

```csv
sample,fastq_1,fastq_2,condition,replicate
WT_rep1,/data/WT_rep1_R1.fastq.gz,/data/WT_rep1_R2.fastq.gz,WT,1
WT_rep2,/data/WT_rep2_R1.fastq.gz,/data/WT_rep2_R2.fastq.gz,WT,2
WT_rep3,/data/WT_rep3_R1.fastq.gz,/data/WT_rep3_R2.fastq.gz,WT,3
mutant_rep1,/data/mutant_rep1_R1.fastq.gz,/data/mutant_rep1_R2.fastq.gz,mutant,1
mutant_rep2,/data/mutant_rep2_R1.fastq.gz,/data/mutant_rep2_R2.fastq.gz,mutant,2
mutant_rep3,/data/mutant_rep3_R1.fastq.gz,/data/mutant_rep3_R2.fastq.gz,mutant,3
```

### 3. Configure os parâmetros

Edite [params.yaml](params.yaml):

```yaml
genome_fasta: "/data/genomes/arabidopsis/TAIR10_chr_all.fas"
genome_gff3:  "/data/genomes/arabidopsis/TAIR10_GFF3_genes.gff"
control_group:   "WT"
treatment_group: "mutant"
```

### 4. Execute o pipeline

```bash
# Execução local
bash run_pipeline.sh local

# Cluster SLURM
bash run_pipeline.sh slurm

# Retomar execução interrompida
bash run_pipeline.sh local --resume
```

---

## Estrutura de Resultados

```
results/
├── qc/
│   ├── pre_trim/          FastQC pré-trimagem
│   ├── post_trim/         FastQC pós-trimagem
│   └── multiqc/           Relatórios MultiQC consolidados
├── trimmed/               Reads trimados + relatórios fastp
├── aligned/               BAMs + logs HISAT2 + flagstat
├── counts/                Matriz de contagens (featureCounts)
├── genome/                annotation.gtf + índice HISAT2
├── deseq2/
│   ├── deseq2_results_all.tsv
│   ├── deseq2_results_sig.tsv
│   ├── normalized_counts.tsv
│   └── figures/           PCA, volcano, heatmap, MA plot
├── enrichment/
│   ├── go_bp_results.tsv
│   ├── kegg_results.tsv
│   ├── gsea_go_results.tsv
│   ├── gsea_kegg_results.tsv
│   └── figures/
├── splicing/
│   ├── splicing_significant.tsv
│   ├── splicing_all.tsv
│   └── rmats_output/
├── wgcna/
│   ├── wgcna_modules.tsv
│   ├── wgcna_hub_genes.tsv
│   └── figures/
├── integration/
│   ├── gene_ranking.tsv
│   ├── key_candidates.tsv
│   └── figures/
├── report/
│   └── rnaseq_report.html
└── nextflow_*.html         Logs e timeline Nextflow
```

---

## Dashboard Interativo

```bash
RESULTS_DIR=results Rscript -e "shiny::runApp('dashboard/app.R', port=3838)"
```

Acesse em `http://localhost:3838`

**Funcionalidades:**
- PCA interativo das amostras
- Volcano plot com filtros dinâmicos (FDR, |log2FC|)
- Tabelas filtráveis e exportáveis (CSV/Excel)
- Dotplots GO e KEGG
- Eventos de splicing por tipo
- Módulos WGCNA e hub genes
- Ranking de candidatos por integration score

---

## Agentes e Skills (Claude Code)

Este projeto inclui agentes especializados para uso com Claude Code:

| Comando | Função |
|---|---|
| `/validate-qc` | Valida QC das amostras e taxa de alinhamento |
| `/interpret-results` | Interpretação biológica dos resultados |
| `/debug-pipeline` | Diagnóstico de falhas no Nextflow |

---

## Parâmetros de Análise

| Parâmetro | Valor padrão | Arquivo |
|---|---|---|
| FDR (padj) | < 0.05 | `params.yaml` |
| \|log2FC\| | > 1.0 | `params.yaml` |
| rMATS FDR | < 0.05 | `scripts/04_integration.R` |
| rMATS \|ΔPSI\| | > 0.10 | `scripts/04_integration.R` |
| WGCNA min genes | 5000 | `modules/analysis.nf` |
| Alinhamento mínimo | 40% | `modules/alignment.nf` |

---

## Citações

Se utilizar este pipeline, cite:

- **DESeq2**: Love MI, Huber W, Anders S (2014). *Genome Biology*, 15:550.
- **HISAT2**: Kim D et al. (2019). *Nature Methods*, 16:3.
- **clusterProfiler**: Wu T et al. (2021). *The Innovation*, 2(3):100141.
- **rMATS**: Shen S et al. (2014). *PNAS*, 111:E5593.
- **WGCNA**: Langfelder P, Horvath S (2008). *BMC Bioinformatics*, 9:559.
- **fastp**: Chen S et al. (2018). *Bioinformatics*, 34:i884.

---

## Autor

**Eulalio Santos** | Universidade Federal de Viçosa  
GitHub: [@eulaliobqi](https://github.com/eulaliobqi)  
Email: eulalio.santos@ufv.br

---

*Desenvolvido com [Claude Code](https://claude.ai/claude-code) – Anthropic*
