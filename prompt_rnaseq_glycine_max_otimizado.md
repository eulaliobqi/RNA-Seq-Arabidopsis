# RNASeq Insight Platform – *Glycine max* Wm82.a4.v1
## Guia Mestre do Pipeline | Adaptado de *Arabidopsis thaliana* TAIR10

> **Base:** Pipeline Arabidopsis TAIR10 — 31 processos, todas as 3 fases concluídas com sucesso em 17-May-2026.
> Este documento consolida todas as lições aprendidas, bugs corrigidos e boas práticas
> para aplicação imediata em soja (*Glycine max* Wm82.a4.v1).

---

## Índice

1. [Visão Geral da Arquitetura](#1-visão-geral-da-arquitetura)
2. [Stack Tecnológico](#2-stack-tecnológico)
3. [Parâmetros Organismo-Específicos](#3-parâmetros-organismo-específicos-glycine-max)
4. [Estrutura de Arquivos](#4-estrutura-de-arquivos)
5. [Fases do Pipeline](#5-fases-do-pipeline)
6. [Agentes Especializados](#6-agentes-especializados)
7. [Skills (Slash Commands)](#7-skills-slash-commands)
8. [Lições Críticas — Bugs do Arabidopsis → Prevenção](#8-lições-críticas--bugs-do-arabidopsis--prevenção)
9. [Instruções de Análise por Módulo](#9-instruções-de-análise-por-módulo)
10. [Checklist de Revisão de Código](#10-checklist-de-revisão-de-código)
11. [Configuração do Servidor](#11-configuração-do-servidor)
12. [Referências e Recursos](#12-referências-e-recursos)

---

## 1. Visão Geral da Arquitetura

```
┌─────────────────────────────────────────────────────────────────────┐
│                   GLYCINE MAX RNA-SEQ PIPELINE                      │
│                    Nextflow DSL2 | v1.0 | 31+ processos             │
└─────────────────────────────────────────────────────────────────────┘

FASE 0 – PRÉ-PROCESSAMENTO
  FASTQ → [FASTQC_PRE] → [FASTP] → [FASTQC_POST] → [MULTIQC]

FASE 1 – ALINHAMENTO E QUANTIFICAÇÃO
  [GFFREAD] → GTF
  GTF + FASTA → [STAR_INDEX] → [STAR_ALIGN] → BAM
  FASTQ → [SALMON_INDEX] → [SALMON_QUANT] → [TXIMPORT]
  BAM → [FEATURECOUNTS] → [PARSE_COUNTS]
  Counts → [COMBAT_SEQ] (opcional)

FASE 1 – EXPRESSÃO DIFERENCIAL
  Counts → [DESEQ2] → DEGs
  DEGs → [ENRICHMENT] (GO/KEGG/GSEA via biomaRt)
  DEGs + GTF → [GOSEQ] (correção viés tamanho)
  Counts → [WGCNA] → Módulos
  BAM → [RMATS] → [PARSE_RMATS] → [RMATS_FILTER] → Splicing
  BAM → [STRINGTIE] → [GFFCOMPARE] → [LNCRNA_PRED]

FASE 2 – INTEGRAÇÃO
  DEGs + Counts + Splicing + GO + KEGG + WGCNA → [INTEGRATION]

FASE 3 – ALTO IMPACTO
  [PLANTFDB] → TFs classificados
  [GENIE3]   → Rede regulatória TF→gene
  [MACHINE_LEARNING] → Biomarcadores (RF + SVM + ElasticNet)
  [PPI_NETWORK] → Rede PPI via STRING (espécie 3847)
  [METANALYSIS] → Validação cruzada GEO/SRA
  [QUARTO_REPORT] → Relatório HTML interativo
```

### Diagrama de Dependências (ordem de execução)

```
GFFREAD ──────────────────────────────────────────────────────┐
                                                               ↓
FASTQ → FASTQC_PRE → MULTIQC_PRE                         STAR_INDEX
FASTQ → FASTP → FASTQC_POST → MULTIQC_POST → STAR_ALIGN ←────┘
                                                    ↓
                               FEATURECOUNTS ← BAM ─┤
                                     ↓              ↓
                               PARSE_COUNTS    RMATS → PARSE_RMATS → RMATS_FILTER
                                     ↓              ↓
                               COMBAT_SEQ    STRINGTIE → GFFCOMPARE → LNCRNA_PRED
                                     ↓
SALMON_INDEX ← GTF          DESEQ2 (results_all, results_sig, norm_counts)
SALMON_QUANT → TXIMPORT          ↓              ↓              ↓
                             ENRICHMENT      GOSEQ         WGCNA
                                     ↓                         ↓
                               INTEGRATION ←─────────────────────
                                     ↓
                    PLANTFDB → GENIE3 → hub TFs
                    MACHINE_LEARNING → features
                    PPI_NETWORK → rede proteica
                    METANALYSIS → validação GEO
                                     ↓
                              QUARTO_REPORT → HTML
```

---

## 2. Stack Tecnológico

### Ambiente `rnaseq-tools` (Nextflow + bioinformática)

| Ferramenta | Versão | Propósito |
|-----------|--------|-----------|
| Nextflow | ≥24.04 (DSL2) | Orquestração |
| FastQC | 0.12.1 | QC das reads |
| MultiQC | 1.21 | Relatório QC agregado |
| fastp | 0.23.4 | Trimagem adaptativa |
| STAR | **2.7.10b** | Alinhamento splice-aware |
| featureCounts | 2.0.6 (subread) | Quantificação por gene |
| Salmon | 1.10.x | Pseudoalinhamento |
| rMATS | ≥4.1.0 | Splicing alternativo |
| StringTie | 2.2.x | Montagem de transcritos |
| gffcompare | 0.12.x | Comparação GTF |
| gffread | 0.12.x | Conversão GFF3→GTF |
| samtools | 1.18 | Processamento BAM |

> **ATENÇÃO:** Usar STAR **2.7.10b**, NÃO 2.7.11b. A versão 2.7.11b via bioconda instala um
> wrapper que seleciona STAR-avx2 causando SIGSEGV (exit 139) em CPUs sem AVX2.
> 2.7.10b tem binário único sem wrapper. Ver [Lição 1](#lição-1-star-271xb-avx2-segfault).

### Ambiente `r-analysis` (R/Bioconductor)

| Pacote | Versão | Propósito |
|--------|--------|-----------|
| DESeq2 | ≥1.46 | Expressão diferencial |
| edgeR | ≥4.4 | filterByExpr |
| tximport | ≥1.30 | Importa Salmon |
| sva/ComBat-Seq | ≥3.52 | Correção de batch |
| clusterProfiler | ≥4.14 | GO/KEGG/GSEA |
| enrichplot | ≥1.24 | Visualização enriquecimento |
| fgsea | ≥1.32 | GSEA rápido |
| biomaRt | ≥2.62 | Anotação via Ensembl Plants |
| AnnotationHub | ≥3.14 | OrgDb customizado Glyma |
| goseq | ≥1.58 | GO com viés de tamanho |
| GenomicFeatures | ≥1.58 | Comprimentos de genes |
| txdbmaker | ≥1.2 | makeTxDbFromGFF (movido de GenomicFeatures) |
| Biostrings | ≥2.72 | Detecção ORF (lncRNA) |
| GEOquery | ≥2.72 | Download GEO |
| limma | ≥3.62 | DE em datasets GEO |
| WGCNA | ≥1.72 | Rede co-expressão |
| GENIE3 | ≥1.28 | Rede regulatória |
| STRINGdb | ≥2.8 | PPI (espécie 3847) |
| caret | ≥6.0 | Machine Learning |
| randomForest | ≥4.7 | RF classificador |
| kernlab | ≥0.9 | SVM-RBF |
| glmnet | ≥4.1 | ElasticNet |
| pROC | ≥1.18 | AUC/ROC |
| igraph | ≥2.0 | Análise de grafos |
| ggplot2 | ≥3.5 | Visualizações |
| pheatmap | ≥1.0.12 | Heatmaps |
| EnhancedVolcano | ≥1.20 | Volcano plot |
| patchwork | ≥1.2 | Layout de figuras |
| Quarto | latest | Relatório HTML |

---

## 3. Parâmetros Organismo-Específicos: *Glycine max*

### Identificadores e Bancos

```yaml
# params.yaml – Glycine max Wm82.a4.v1
organism:         "Glycine max"
assembly:         "Wm82.a4.v1"
kegg_organism:    "gma"          # código KEGG Glycine max
string_species:   3847           # taxon ID STRING para Glycine max
plantfdb_prefix:  "Gma"         # PlantTFDB: Gma_TF_list.txt.gz
biomart_dataset:  "gmax_eg_gene" # Ensembl Plants dataset
```

### IDs de Genes Wm82.a4.v1

```
Formato padrão:  Glyma.01G000100
Transcritos:     Glyma.01G000100.1, Glyma.01G000100.2
GFF3 Phytozome:  gene_id "Glyma.01G000100"
                 transcript_id "Glyma.01G000100.1"
```

**Regra crítica:** Sempre strip do sufixo de versão do transcrito antes de qualquer busca em banco:
```r
gene_id <- gsub("\\.[0-9]+$", "", gene_id)  # Glyma.01G000100.1 → Glyma.01G000100
```

Diferente de Arabidopsis (que tinha `.TAIR10`), *Glycine max* Wm82.a4.v1 usa sufixos numéricos
de versão (`.1`, `.2`) que devem ser removidos apenas nos transcritos, não nos genes.

### OrgDb — Problema Crítico

`org.Gmax.eg.db` foi **removido do Bioconductor 3.20+** e não está disponível.
Estratégia substituta (usar em todos os scripts R):

```r
# Opção A: AnnotationHub (recomendada — funciona sem internet em cache)
library(AnnotationHub)
ah  <- AnnotationHub()
# Procurar: query(ah, c("Glycine max", "OrgDb"))
# Usar o mais recente disponível
orgdb_gmax <- ah[["AH100000"]]   # substituir pelo ID correto após query

# Opção B: biomaRt via Ensembl Plants (requer internet)
library(biomaRt)
mart <- useMart("plants_mart",
                dataset  = "gmax_eg_gene",
                host     = "https://plants.ensembl.org")
# Mapping Glyma → GO:
go_map <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
                mart = mart)

# Opção C: Arquivo de anotação GO do Phytozome/SoyBase (offline)
# Download: https://www.soybase.org/resources/data_downloads.php
# Arquivo: Gmax_508_Wm82.a4.v1.annotation_info.txt
# Coluna gene_id → go_terms (pipe-separated)
```

**Para clusterProfiler sem OrgDb:**
```r
# Constrói TERM2GENE e TERM2NAME manualmente a partir do biomaRt/SoyBase
go_df <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
               mart = mart) |>
  filter(!is.na(go_id), go_id != "")

TERM2GENE <- go_df |> select(go_id, ensembl_gene_id)
TERM2NAME <- go_df |> select(go_id, name_1006) |> distinct()

# enricher() aceita TERM2GENE/TERM2NAME sem OrgDb:
ego <- enricher(gene = deg_genes,
                TERM2GENE = TERM2GENE,
                TERM2NAME = TERM2NAME,
                pAdjustMethod = "BH")
```

### Genoma de Referência

```
Fonte: JGI Phytozome v13 (https://phytozome-next.jgi.doe.gov/)
  Glyma.Wm82.a4.v1.fa           (cromossomos Chr01–Chr20 + scaffolds)
  Glyma.Wm82.a4.v1.gene.gff3   (anotação GFF3)

Alternativa SoyBase (https://www.soybase.org/):
  Gmax_508_Wm82.a4.v1.fa.gz
  Gmax_508_Wm82.a4.v1.gene_exons.gff3.gz

Cromossomos: Chr01 a Chr20 (+ scaffolds Scaffold_*)
  ⚠ Mesmos nomes no FASTA e no GFF3/GTF → sem problemas de chrChr como em Arabidopsis
  ⚠ Verificar que gffread NÃO adiciona prefixo indesejado
```

### Recursos Computacionais

Glycine max (~978 Mb) é ~7× maior que Arabidopsis (135 Mb):

| Processo | RAM mínima | Recomendado |
|---------|-----------|-------------|
| STAR_INDEX | 32 GB | 64 GB |
| STAR_ALIGN | 16 GB | 32 GB |
| SALMON_INDEX | 8 GB | 16 GB |
| DESEQ2 | 16 GB | 32 GB |
| WGCNA | 32 GB | 64 GB |
| GENIE3 | 32 GB | 64 GB |

---

## 4. Estrutura de Arquivos

```
glycine-max-rnaseq/
├── CLAUDE.md                    # Instruções para Claude (este documento adaptado)
├── main.nf                      # Pipeline principal DSL2
├── nextflow.config              # Recursos (CPUs, memória) por processo
├── params.yaml                  # Parâmetros biológicos (editar por experimento)
├── samplesheet.csv              # Amostras com paths reais
├── setup.sh                     # Script de setup inicial
├── run_pipeline.sh              # Script de execução
│
├── modules/                     # Módulos Nextflow DSL2
│   ├── qc.nf                   # FASTQC_PRE, FASTQC_POST, MULTIQC
│   ├── trimming.nf              # FASTP
│   ├── alignment.nf             # GFFREAD, STAR_INDEX, STAR_ALIGN, SAMTOOLS
│   ├── quantification.nf        # FEATURECOUNTS, PARSE_COUNTS
│   ├── salmon.nf                # SALMON_INDEX, SALMON_QUANT, TXIMPORT
│   ├── batch.nf                 # COMBAT_SEQ
│   ├── splicing.nf              # RMATS, PARSE_RMATS, RMATS_FILTER
│   ├── assembly.nf              # STRINGTIE, GFFCOMPARE
│   ├── lncrna.nf               # LNCRNA_PRED
│   └── analysis.nf              # DESEQ2, ENRICHMENT, GOSEQ, WGCNA,
│                               #   INTEGRATION, METANALYSIS, PLANTFDB,
│                               #   GENIE3, MACHINE_LEARNING, PPI_NETWORK,
│                               #   QUARTO_REPORT
│
├── scripts/                     # Scripts R (numerados por ordem de fase)
│   ├── 00_tximport.R           # Importa Salmon
│   ├── 01_deseq2.R             # DESeq2 + visualizações
│   ├── 02_enrichment.R          # GO/KEGG/GSEA via biomaRt (SEM OrgDb)
│   ├── 02b_goseq.R             # GOseq com gene2cat customizado
│   ├── 03_wgcna.R              # WGCNA modular
│   ├── 04_integration.R         # Score multi-camada
│   ├── 05_batch_correction.R    # ComBat-Seq + PCA diagnóstico
│   ├── 06_machine_learning.R    # RF + SVM + ElasticNet (caret)
│   ├── 07_ppi_network.R        # STRING PPI (STRINGdb + httr fallback)
│   ├── 08_plantfdb.R           # PlantTFDB TF classification
│   ├── 09_genie3.R             # GENIE3 inferência de rede
│   ├── 10_metanalysis.R        # Metanálise GEO
│   └── 11_lncrna.R             # Predição lncRNA via Biostrings
│
├── envs/
│   ├── rnaseq-tools.yml        # Conda env: STAR, featureCounts, etc.
│   └── r-analysis.yml          # Conda env: R/Bioconductor
│
├── report/
│   └── rnaseq_report.qmd       # Template Quarto
│
├── dashboard/
│   └── app.R                    # Dashboard Shiny interativo
│
└── results/                     # Saídas do pipeline (gitignore)
    ├── qc/
    ├── counts/
    ├── deseq2/
    ├── enrichment/
    ├── splicing/
    ├── wgcna/
    ├── integration/
    ├── ml/
    ├── plantfdb/
    ├── genie3/
    ├── network/
    ├── lncrna/
    ├── metanalysis/
    └── report/
```

---

## 5. Fases do Pipeline

### Fase 0 – Setup (uma vez)

```bash
# 1. Clonar repositório
git clone https://github.com/eulaliobqi/glycine-max-rnaseq.git
cd glycine-max-rnaseq

# 2. Criar ambientes conda
mamba env create -f envs/rnaseq-tools.yml
mamba env create -f envs/r-analysis.yml

# 3. Baixar genoma Wm82.a4.v1 do Phytozome (requer conta JGI)
# https://phytozome-next.jgi.doe.gov/
# Arquivos: Glyma.Wm82.a4.v1.fa + Glyma.Wm82.a4.v1.gene.gff3

# 4. Atualizar params.yaml com paths reais
```

### Fase 1 – Core (obrigatório)

**Ordem estrita de execução:**

```
1. GFFREAD        → GFF3 → GTF (valida nomes de cromossomo!)
2. FASTQC_PRE     → QC das reads brutas
3. FASTP          → Trimagem (min_length=36, quality=20)
4. FASTQC_POST    → QC pós-trimagem
5. MULTIQC        → Relatório agregado (pré e pós)
6. STAR_INDEX     → Índice genômico (uma vez, cache)
7. STAR_ALIGN     → Alinhamento splice-aware
8. SALMON_INDEX   → Índice Salmon (paralelo)
9. SALMON_QUANT   → Quantificação Salmon
10. TXIMPORT      → Importa Salmon (tx→gene)
11. FEATURECOUNTS → Contagens por gene
12. PARSE_COUNTS  → Normaliza nomes de colunas + strip versão
13. COMBAT_SEQ    → Correção de batch (se detectado)
14. DESEQ2        → DEGs (padj<0.05, |lfc|>1.0)
15. ENRICHMENT    → GO/KEGG/GSEA
16. GOSEQ         → GO com correção viés tamanho
17. WGCNA         → Módulos co-expressão
18. RMATS         → Splicing alternativo
19. PARSE_RMATS   → Parse resultados rMATS
20. RMATS_FILTER  → Filtra eventos significativos
21. INTEGRATION   → Score multi-camada
22. QUARTO_REPORT → Relatório HTML
```

### Fase 2 – Avançado (run_stringtie=true)

```
23. STRINGTIE     → Montagem de transcritos por amostra
24. GFFCOMPARE    → Compara com anotação de referência
25. LNCRNA_PRED   → Predição de lncRNAs novos
```

### Fase 3 – Alto Impacto (run_ml=true, run_ppi=true)

```
26. PLANTFDB      → Classifica TFs DEGs por família
27. GENIE3        → Rede regulatória TF→gene (RF)
28. MACHINE_LEARNING → Biomarcadores (RF+SVM+ElasticNet)
29. PPI_NETWORK   → Rede PPI (STRINGdb espécie 3847)
30. METANALYSIS   → Validação cruzada GEO/SRA
```

### Critérios de Qualidade (cutoffs)

```yaml
padj_cutoff:       0.05      # FDR DESeq2
lfc_cutoff:        1.0       # |log2FC| mínimo
splicing_fdr:      0.05      # FDR rMATS
splicing_dpsi:     0.1       # |ΔPSI| mínimo
wgcna_r2:          0.85      # R² mínimo soft power
integration_score: "lfc×3 + mean×2 + sig×2 + splicing×1.5 + pathway×1.0 + hub×0.5"
key_candidates:    ≥2        # evidências em ≥2 camadas
```

---

## 6. Agentes Especializados

### `/validate-qc` — Agente de Validação QC

**Descrição:** Analisa saídas do FastQC/MultiQC e valida qualidade antes de prosseguir.

**Trigger:** Após MULTIQC_POST completar.

**O que verificar:**
- Per-base sequence quality ≥ Q20 em >90% das posições
- %GC dentro de ±5% do esperado para Glycine max (~44%)
- Taxa de duplicação < 50% (RNA-Seq)
- Tamanho de reads após trimagem ≥ 36 bp
- % reads alinhadas pelo STAR ≥ 70% (`Overall alignment rate`)
- % reads mapeadas pelo Salmon ≥ 60% (Mapping rate)
- Correlação inter-réplicas (Pearson no PCA) ≥ 0.95

**Comando de inspeção:**
```bash
# Log de alinhamento STAR (grep -i para capturar maiúsculas/minúsculas)
grep -i "overall alignment rate" logs/star_align_*.log

# Taxa de mapeamento Salmon
grep "Mapping rate" results/counts/salmon/*/logs/salmon_quant.log
```

**Ação:** Se qualidade insatisfatória → ajustar `quality` e `min_length` em `params.yaml`
e re-executar a partir de FASTP com `-resume`.

---

### `/interpret-results` — Agente de Interpretação Biológica

**Descrição:** Interpreta resultados DESeq2, enriquecimento GO/KEGG, WGCNA e integração
no contexto biológico de *Glycine max* e do experimento.

**Trigger:** Após INTEGRATION + ENRICHMENT completarem.

**Framework de interpretação:**

1. **DEGs:**
   - Quantos up vs down? Qual proporção?
   - Genes com |lfc| > 2 são os mais responsivos
   - Verificar genes conhecidos da literatura de soja (e.g., GmWRKY, GmERF, GmNAC)

2. **Enriquecimento GO (BP prioritário):**
   - Termos de resposta a estresse aiótico e biótico?
   - Processos de desenvolvimento (floração, fixação de N₂)?
   - Via de síntese de isoflavonoides (importante em soja)?

3. **KEGG:**
   - Vias de metabolismo de nitrogênio (gma00910)?
   - Biossíntese de flavonoides (gma00941)?
   - Fotossíntese, resposta a hormônios?

4. **WGCNA:**
   - Módulos com |r| > 0.7 e p < 0.05 são biologicamente relevantes
   - Hub genes: conectores centrais da rede

5. **Integração:**
   - Genes com score > 4.0 e ≥3 camadas de evidência são candidatos fortes
   - Priorizar: DEG + módulo WGCNA relevante + anotação funcional conhecida

6. **TFs (PlantTFDB):**
   - Famílias WRKY, ERF/AP2, MYB, NAC são chave em respostas a estresse em soja
   - GmWRKY = reguladores centrais de defesa

---

### `/debug-pipeline` — Agente de Diagnóstico

**Descrição:** Diagnóstico de falhas Nextflow, logs de erro e sugestões de correção.

**Trigger:** Qualquer processo com exit status ≠ 0.

**Protocolo de diagnóstico:**

```bash
# 1. Identificar o processo que falhou
cat .nextflow.log | grep "ERROR\|failed\|exit"

# 2. Ir ao work directory do processo
cd work/XX/YYYYYY*

# 3. Verificar os 3 arquivos de diagnóstico
cat .command.sh   # script executado
cat .command.out  # stdout
cat .command.err  # stderr (geralmente onde está o erro R)

# 4. Reproduzir manualmente (substitui variáveis NF)
bash .command.sh 2>&1 | head -100
```

**Erros comuns e soluções:**

| Erro | Causa | Solução |
|------|-------|---------|
| `exit 139 (SIGSEGV)` STAR | STAR-avx2 em CPU sem AVX2 | `star=2.7.10b` em rnaseq-tools.yml |
| `Missing output file: X.tsv` | Script R falhou antes de gerar | Ler `.command.err`; o script inicializa outputs vazios no início |
| `cannot rescale constant column` | prcomp em gene de variância zero | `vst_mat[apply(vst_mat,1,var)>0,]` |
| `makeTxDbFromGFF deprecated` | GenomicFeatures ≥1.61 | helper `make_txdb()` com txdbmaker |
| `x has insufficient unique values` | GOseq com comprimentos uniformes | Usar `method="Hypergeometric"` |
| `SyntaxError: unterminated string` Python | `\n` em heredoc Groovy | `rstrip()` sem argumento |
| `FileNotFoundError: gffcmp.merged.gtf` | gffcompare ≥0.12 usa `.annotated.gtf` | Cadeia de fallback; ver lncrna.nf |
| `GENIE3 falhou durante execução` | nCores > 1 sem doParallel | fallback `nCores=1` em `run_genie3()` |
| `GOseq genome=tair10 não encontrado` | banco interno goseq desatualizado | `gene2cat` de org.db via AnnotationDbi |

---

### `/review-code` — Agente de Revisão de Código

**Descrição:** Revisão sistemática de scripts R e módulos Nextflow antes de commit.

**Trigger:** Antes de qualquer `git push` com novos scripts.

**Checklist de revisão automática:**
- [ ] Todos os outputs inicializados vazios no início do script?
- [ ] Escapes `\n` e `\t` ausentes em heredocs Python dentro de Nextflow?
- [ ] Strip de sufixo de versão nos gene IDs (`.1`, `.2`)?
- [ ] `tryCatch` em todas as chamadas a funções externas?
- [ ] `suppressPackageStartupMessages` em todos os scripts R?
- [ ] Parâmetro `show_col_types = FALSE` em todos os `read_tsv`?
- [ ] `nCores` com fallback serial nos processos paralelos?
- [ ] Nenhum `genome=X, id=Y` sem verificar disponibilidade no banco?

---

### `/check-ids` — Agente de Validação de IDs

**Descrição:** Verifica consistência dos gene IDs entre GTF, counts, DEGs e bancos de anotação.

**Trigger:** Após PARSE_COUNTS e DESEQ2.

**Verificações:**

```r
# 1. Formato correto dos IDs?
head(gene_ids)
# Esperado: "Glyma.01G000100" (sem sufixo de versão)
# Problema: "Glyma.01G000100.1" → aplicar gsub("\\.[0-9]+$", "", .)

# 2. IDs reconhecidos pelo banco KEGG?
kegg_map <- clusterProfiler::bitr_kegg(gene_ids[1:10], fromType="ncbi-geneid",
                                        toType="kegg", organism="gma")

# 3. IDs reconhecidos pelo biomaRt?
test_map <- getBM(attributes=c("ensembl_gene_id"),
                  filters="ensembl_gene_id",
                  values=gene_ids[1:10], mart=mart)
cat("IDs reconhecidos:", nrow(test_map), "/", 10, "\n")

# 4. Cobertura das anotações GO
go_coverage <- sum(gene_ids %in% go_map$ensembl_gene_id) / length(gene_ids)
cat("Cobertura GO:", round(go_coverage*100, 1), "%\n")
# Esperado: >60% para Wm82.a4.v1 (genoma bem anotado)
```

---

## 7. Skills (Slash Commands)

### Uso no Claude Code

```bash
# Validar QC após multiqc
/validate-qc

# Interpretar resultados biológicos
/interpret-results

# Debugar processo falho
/debug-pipeline

# Revisar código antes de commit
/review-code

# Verificar IDs de genes
/check-ids
```

### Definição no CLAUDE.md

```markdown
## Agentes especializados
- `/validate-qc`        — Valida QC das amostras (FastQC/MultiQC + STAR + Salmon)
- `/interpret-results`  — Interpretação biológica em contexto de Glycine max
- `/debug-pipeline`     — Diagnóstico de falhas Nextflow (logs + work dir)
- `/review-code`        — Revisão de código R/Nextflow antes de commit
- `/check-ids`          — Validação de consistência dos gene IDs Glyma
```

---

## 8. Lições Críticas — Bugs do Arabidopsis → Prevenção

### Lição 1: STAR 2.7.11b AVX2 Segfault

**Arabidopsis:** STAR 2.7.11b instalado via bioconda inclui wrapper que seleciona
automaticamente o binário STAR-avx2 em sistemas x86_64. Em CPUs sem AVX2 (hardware antigo
ou VMs), causa SIGSEGV (exit 139) com o message "Killed" no log.

**Prevenção em Glycine max:**
```yaml
# envs/rnaseq-tools.yml — fixar versão!
- star=2.7.10b    # NÃO usar 2.7.11b; 2.7.10b tem binário único sem wrapper
```

**Teste de verificação:**
```bash
mamba run -n rnaseq-tools STAR --version
# Deve mostrar: 2.7.10b
# Se mostrar 2.7.11b com "*-avx2" no path → reinstalar
```

---

### Lição 2: Sufixo de Versão nos Gene IDs

**Arabidopsis:** gffread adicionava `.TAIR10` aos IDs (e.g., `AT3G30775.TAIR10`).
`org.At.tair.db` não reconhecia esse sufixo → falhas silenciosas em todas as buscas.

**Glycine max equivalente:** gffread pode adicionar `.1`, `.2` (versão do transcrito)
quando o GFF3 do Phytozome usa `mRNA` como feature primária em vez de `gene`.

**Prevenção:**
```r
# Em TODOS os scripts R que leem DESeq2/counts:
gene_id <- gsub("\\.[0-9]+$", "", gene_id)
# Também remover sufixo Wm82.a4.v1 se aparecer:
gene_id <- gsub("\\.Wm82\\.a4\\.v1.*$", "", gene_id)
```

**No PARSE_COUNTS (módulo Nextflow):**
```bash
# Após featureCounts, strip de sufixo nas linhas do counts:
sed 's/\.[0-9]*$//' counts_raw.txt > counts_clean.txt
```

---

### Lição 3: chrChr nos Cromossomos (GFFREAD)

**Arabidopsis:** gffread adicionava prefixo `chr` a cromossomos já capitalizados (`Chr5` → `chrChr5`),
causando incompatibilidade entre GTF e BAM.

**Glycine max (Wm82.a4.v1):** Os cromossomos são `Chr01`...`Chr20` no FASTA do Phytozome.
gffread pode adicionar `chr` → `chrChr01`.

**Prevenção no módulo GFFREAD:**
```bash
# Corrige automaticamente após conversão GFF3→GTF:
gffread ${gff3} -T -o - | sed 's/^chrChr/Chr/' > annotation.gtf
# Alternativa segura (cobre Chr e chr):
gffread ${gff3} -T -o - | \
  sed 's/^chrChr/Chr/; s/^chrchr/Chr/' > annotation.gtf
```

**Validação:**
```bash
# Confirmar que nomes batem entre GTF e BAM
cut -f1 annotation.gtf | sort -u | head -25
samtools view -H sample.bam | grep "^@SQ" | cut -f2 | sed 's/SN://' | head -25
# Devem ser idênticos (ex: Chr01, Chr02, ..., Chr20)
```

---

### Lição 4: keyType="TAIR" vs Glycine max

**Arabidopsis:** `keyType="TAIR"` falha silenciosamente em clusterProfiler 4.x com
`org.At.tair.db` → solução: converter para ENTREZID via `bitr()`.

**Glycine max:** `org.Gmax.eg.db` foi removido do Bioconductor 3.20+.
A solução correta é usar `biomaRt` + `enricher()` com TERM2GENE/TERM2NAME manuais.

**Prevenção:**
```r
# NÃO usar (vai falhar):
# enrichGO(gene = glyma_ids, OrgDb = org.Gmax.eg.db, keyType = "GID")

# USAR (funciona):
# 1. Obter GO via biomaRt (Ensembl Plants)
go_df <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
               mart = mart)

# 2. Separar por ontologia
make_term2gene <- function(df, ont) {
  df |>
    filter(namespace_1003 == ont, go_id != "") |>
    select(go_id, ensembl_gene_id)
}

TERM2GENE_BP <- make_term2gene(go_df, "biological_process")
TERM2GENE_MF <- make_term2gene(go_df, "molecular_function")
TERM2GENE_CC <- make_term2gene(go_df, "cellular_component")

# 3. enricher() sem OrgDb
ego_bp <- enricher(gene = deg_genes,
                   TERM2GENE = TERM2GENE_BP,
                   TERM2NAME = go_df |> filter(namespace_1003=="biological_process") |>
                               select(go_id, name_1006) |> distinct(),
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
```

---

### Lição 5: ComBat-Seq prcomp Variância Zero

**Arabidopsis:** Após VST, genes com expressão constante entre amostras têm variância=0.
`prcomp(scale.=TRUE)` falha com "cannot rescale a constant/zero column".

**Prevenção:**
```r
# Em 05_batch_correction.R, ANTES do prcomp:
vst_mat <- vst_mat[apply(vst_mat, 1, var) > 0, , drop = FALSE]
# Também aplicar antes de qualquer prcomp no pipeline
```

---

### Lição 6: GOseq — genome/id não disponível

**Arabidopsis:** `genome="tair10", id="tair"` não estava no banco interno do goseq ≥1.58
→ erro silencioso capturado pelo tryCatch → outputs com 0 bytes.

**Prevenção (Glycine max — nenhum código gma no goseq):**
```r
# NUNCA usar genome="gma" ou similar — não existe no banco do goseq
# SEMPRE usar gene2cat customizado:

get_gene2go <- function(gene_ids, ontology, mart) {
  go_df <- getBM(
    attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
    filters    = "ensembl_gene_id",
    values     = gene_ids,
    mart       = mart
  )
  ont_map <- c("BP"="biological_process", "MF"="molecular_function",
               "CC"="cellular_component")
  go_df <- go_df[go_df$namespace_1003 == ont_map[ontology] & !is.na(go_df$go_id), ]
  if (nrow(go_df) == 0) return(list())
  split(go_df$go_id, go_df$ensembl_gene_id)
}

# Usar:
gene2cat_bp <- get_gene2go(names(is_deg), "BP", mart)
res_bp <- goseq(pwf, gene2cat = gene2cat_bp, method = goseq_method)
```

---

### Lição 7: makeTxDbFromGFF — Migração de Pacote

**Arabidopsis:** `GenomicFeatures::makeTxDbFromGFF()` foi movido para `txdbmaker::makeTxDbFromGFF()`
em GenomicFeatures ≥1.61.1 → erro de função não encontrada.

**Prevenção:**
```r
# Helper em TODOS os scripts que usam makeTxDbFromGFF:
make_txdb <- function(gtf) {
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdbmaker::makeTxDbFromGFF(gtf, format = "GTF")
  } else {
    GenomicFeatures::makeTxDbFromGFF(gtf, format = "GTF")
  }
}
```

---

### Lição 8: GOseq nullp — Comprimentos Uniformes

**Arabidopsis:** Quando `txdb` falha e fallback usa comprimentos uniformes (1000L para todos),
`nullp()` não consegue ajustar spline com <6 valores únicos → erro de knots.

**Prevenção:**
```r
n_unique <- length(unique(gene_lengths))
if (n_unique >= 6) {
  pwf          <- nullp(is_deg, bias.data = gene_lengths, plot.fit = FALSE)
  goseq_method <- "Wallenius"
} else {
  cat("Aviso: comprimentos uniformes → usando Hypergeometric\n")
  pwf <- data.frame(DEgenes=is_deg, bias.data=gene_lengths,
                    pwf=sum(is_deg)/length(is_deg))
  rownames(pwf) <- names(is_deg)
  goseq_method <- "Hypergeometric"
}
```

---

### Lição 9: rMATS — Formato b1.txt / b2.txt

**Arabidopsis:** rMATS exige arquivos `b1.txt` e `b2.txt` com **todos os BAMs numa só linha,
separados por vírgula** (sem newline). Formato errado causa falha de parsing.

**Prevenção:**
```bash
# CORRETO (uma linha, vírgulas):
echo "ctrl_rep1.bam,ctrl_rep2.bam,ctrl_rep3.bam" > b1.txt
echo "trt_rep1.bam,trt_rep2.bam,trt_rep3.bam"  > b2.txt

# ERRADO (uma por linha):
# ctrl_rep1.bam
# ctrl_rep2.bam
```

**No módulo Nextflow:**
```groovy
script:
def b1_bams = control_bams.collect { it.toString() }.join(',')
def b2_bams = treatment_bams.collect { it.toString() }.join(',')
"""
echo "${b1_bams}" > b1.txt
echo "${b2_bams}" > b2.txt
rmats.py --b1 b1.txt --b2 b2.txt ...
"""
```

---

### Lição 10: LNCRNA_PRED — gffcompare.annotated.gtf

**Arabidopsis:** Script Python abria `gffcmp.merged.gtf` diretamente. gffcompare ≥0.12
cria `gffcmp.annotated.gtf` (não `merged.gtf`) → FileNotFoundError.

**Prevenção — cadeia de fallback obrigatória:**
```python
# 1. Tenta gffcmp.annotated.gtf (gffcompare >=0.12)
# 2. Tenta gffcmp.merged.gtf (versão antiga)
# 3. Parse gffcmp.tracking → filtra merged.gtf por TCONS IDs novos
# 4. Usa merged.gtf completo como último recurso

novel_codes = {'u', 'i', 'x', 's', 'o'}
for annotated_gtf in ("gffcmp.annotated.gtf", "gffcmp.merged.gtf"):
    if os.path.exists(annotated_gtf):
        # ... processar

if not out_lines and os.path.exists("gffcmp.tracking"):
    novel_ids = set()
    with open("gffcmp.tracking") as fh:
        for line in fh:
            parts = line.rstrip().split('\t')  # rstrip() SEM argumento!
            if len(parts) >= 4 and parts[3] in novel_codes:
                novel_ids.add(parts[0])
```

---

### Lição 11: \n e \t em Heredocs Python dentro do Nextflow

**Arabidopsis:** Nextflow Groovy interpola `\n` como newline literal e `\t` como tab literal
dentro do bloco `script:` (string triple-quoted Groovy) **antes** de escrever o `.command.sh`.
Isso quebra strings Python: `rstrip('\n')` → `rstrip('[newline]')` → SyntaxError.

**Prevenção:**
- **Nunca usar `\n` em strings Python** dentro de scripts Nextflow
- `\t` dentro de strings Python é válido (tab literal) mas dificulta leitura
- Substituições seguras:

```python
# PROIBIDO em heredoc Nextflow:
line.rstrip('\n')           # → SyntaxError
line.split('\t')            # → tab literal (funciona mas evitar)
if '\ttranscript\t' in l   # → tab literal (funciona mas evitar)

# USAR:
line.rstrip()               # strip qualquer whitespace
line.split()                # split em qualquer whitespace (cuidado se campos têm espaços)
# Ou para tab específico:
sep = chr(9)                # tab via chr() — não usa escape
line.split(sep)
```

---

### Lição 12: GENIE3 — Paralelismo

**Arabidopsis:** `GENIE3(nCores=8)` falhou silenciosamente (provável falta de doParallel
registrado) → outputs com 0 bytes + script exit 0 (Nextflow cacheou o resultado vazio).

**Prevenção — sempre usar fallback serial:**
```r
run_genie3 <- function(expr_mat, regulators, targets, n_trees, ncores) {
  tryCatch(
    GENIE3(exprMatrix=expr_mat, regulators=regulators, targets=targets,
           treeMethod="RF", K="sqrt", nTrees=n_trees, nCores=ncores),
    error = function(e) {
      message(sprintf("GENIE3 (nCores=%d) falhou: %s", ncores, e$message))
      if (ncores > 1L) {
        message("Repetindo com nCores=1...")
        tryCatch(
          GENIE3(exprMatrix=expr_mat, regulators=regulators, targets=targets,
                 treeMethod="RF", K="sqrt", nTrees=n_trees, nCores=1L),
          error = function(e2) { message("GENIE3 serial falhou: ", e2$message); NULL }
        )
      } else NULL
    }
  )
}
```

---

### Lição 13: Inicialização de Outputs no Início dos Scripts R

**Arabidopsis:** Nextflow falha com "Missing output file" quando o script R termina antes
de criar todos os arquivos declarados no bloco `output:` do processo.

**Regra obrigatória:** Todo script R deve inicializar **todos** os outputs vazios na primeira
seção, antes de qualquer lógica de negócio:

```r
# ── Inicializa outputs (evita "Missing output file" em saídas antecipadas) ──
write_tsv(data.frame(), file.path(opt$outdir, "results.tsv"))
write_tsv(data.frame(), file.path(opt$outdir, "summary.tsv"))
writeLines("Script inicializando...", file.path(opt$outdir, "report.txt"))
writeLines("", file.path(opt$figures_dir, ".placeholder"))
```

---

### Lição 14: STRINGdb — Acesso à Internet no Servidor

**Arabidopsis:** O servidor não tinha acesso à internet → STRINGdb R package falhou
na primeira conexão ao STRINGdb server.

**Prevenção:**
```r
# 07_ppi_network.R — sempre inicializar outputs antes de tentar STRINGdb:
# Outputs vazios no início do script (ver Lição 13)

# Tentar STRINGdb (pode falhar sem internet):
ppi_result <- tryCatch({
  string_db <- STRINGdb$new(version="11.5", species=3847, # 3847 = Glycine max
                             score_threshold=400,
                             input_directory="/tmp/stringdb_cache")
  # ... lógica PPI
}, error = function(e) {
  message("STRINGdb falhou (sem internet?): ", e$message)
  NULL
})

# Se falhou, outputs vazios já foram criados → pipeline continua
```

---

## 9. Instruções de Análise por Módulo

### 01_deseq2.R

**Objetivo:** Identificar genes diferencialmente expressos.

**Parâmetros:**
```r
# LFC shrinkage: tentar apeglm → ashr → normal (fallback automático)
tryCatch(
  res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm"),
  error = function(e) {
    tryCatch(
      res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="ashr"),
      error = function(e2) results(dds)
    )
  }
)
```

**filterByExpr (edgeR) para pré-filtragem:**
```r
library(edgeR)
keep <- filterByExpr(counts_mat, group = metadata$condition)
dds  <- DESeqDataSetFromMatrix(counts_mat[keep, ], colData = metadata, design = ~condition)
```

**Saídas esperadas:**
- `deseq2_results_all.tsv` — todos os genes testados
- `deseq2_results_sig.tsv` — apenas DEGs (padj<0.05, |lfc|>1.0)
- `normalized_counts.tsv` — rlog ou VST normalizados

---

### 02_enrichment.R (Glycine max — sem OrgDb)

**Objetivo:** GO/KEGG/GSEA sem `org.Gmax.eg.db`.

**Estratégia obrigatória:**
```r
# 1. Conectar ao Ensembl Plants
mart <- useMart("plants_mart",
                dataset = "gmax_eg_gene",
                host    = "https://plants.ensembl.org")

# 2. Mapear Glyma IDs → NCBI Entrez (para KEGG)
id_map <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                filters    = "ensembl_gene_id",
                values     = all_genes,
                mart       = mart)

# 3. KEGG com Entrez IDs
entrez_degs <- id_map$entrezgene_id[id_map$ensembl_gene_id %in% deg_genes]
kegg_res    <- enrichKEGG(gene = na.omit(entrez_degs),
                           universe   = na.omit(id_map$entrezgene_id),
                           organism   = "gma",
                           pvalueCutoff = 0.05)

# 4. GO via enricher() com TERM2GENE manual (ver Lição 4)
```

**GSEA (lista ordenada por stat):**
```r
gsea_list <- deseq2_all |>
  filter(!is.na(stat)) |>
  arrange(desc(stat)) |>
  { x <- x$stat; names(x) <- x$gene_id; x }()
```

---

### 02b_goseq.R (Glycine max)

**Objetivo:** GO com correção de viés de comprimento do gene.

**Fonte de comprimentos:** TxDb do GTF Phytozome (via txdbmaker).

**gene2cat:** Construir via biomaRt (ver Lição 6 acima).

**Fluxo:**
```
GTF → txdbmaker → comprimentos por gene
Glyma IDs → biomaRt → GO mapping por ontologia
PWF → goseq (Wallenius ou Hypergeometric)
```

---

### 03_wgcna.R

**Objetivo:** Identificar módulos de co-expressão correlacionados com o tratamento.

**Ajustes para Glycine max (genoma maior):**
```r
# maxBlockSize maior para comportar mais genes
WGCNA::blockwiseModules(
  datExpr    = expr_mat,
  power      = soft_power,
  maxBlockSize = 30000,    # Arabidopsis usava 20000; soja tem mais genes
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  TOMType    = "unsigned",
  saveTOMs   = FALSE,
  verbose    = 3
)
```

**Soft power:** R² ≥ 0.85 (mesmo critério do Arabidopsis).

---

### 05_batch_correction.R

**Objetivo:** Detectar e corrigir efeito de batch usando PCA.

**Critério de aplicação:** PC1 > 40% da variância AND correlação PC1 × batch > 0.7.

**OBRIGATÓRIO antes do prcomp:**
```r
vst_mat <- vst_mat[apply(vst_mat, 1, var) > 0, , drop = FALSE]
```

---

### 06_machine_learning.R

**Importante sobre AUC=1.0:**
Com poucos replicatas (n=6), cross-validation retorna AUC=1.0 trivialmente
(cada fold de teste tem apenas 1 amostra). **Usar `feature_importance.tsv`**
em vez dos resultados de AUC para seleção de biomarcadores.

Para experimentos maiores (n≥12):
```r
# Aumentar k-fold para 5 folds se n>=12
k_fold <- min(5, floor(nrow(meta)/2))
# Usar LOOCV se n < 10:
if (nrow(meta) < 10) {
  ctrl <- trainControl(method="LOOCV", classProbs=TRUE, summaryFunction=twoClassSummary)
}
```

---

### 07_ppi_network.R

**Espécie STRING para Glycine max:** `3847` (não 3702 do Arabidopsis).

```r
string_db <- STRINGdb$new(version    = "11.5",
                           species    = 3847,       # Glycine max
                           score_threshold = 400,
                           input_directory = "/tmp/stringdb_gmax")
```

---

### 08_plantfdb.R

**PlantTFDB para Glycine max:**
- URL: `http://planttfdb.gao-lab.org/download.php`
- Arquivo: `Gma_TF_list.txt.gz` (prefixo `Gma`, não `Ath`)

```r
# Atualizar em load_plantfdb():
plantfdb_url <- "http://planttfdb.gao-lab.org/download/TF_list/Gma_TF_list.txt.gz"
```

**Famílias TF importantes em Glycine max:**
- **WRKY:** Resistência a patógenos (GmWRKY)
- **MYB:** Biossíntese de flavonoides/isoflavonoides
- **ERF/AP2:** Resposta a etileno, estresse
- **NAC:** Senescência, resposta a seca
- **bHLH:** Regulação de ferro, antocianinas
- **ARF:** Sinalização de auxina (desenvolvimento de vagem)

---

### 09_genie3.R

**Usar sempre** `run_genie3()` com fallback serial (ver Lição 12).

**Reguladores:** TFs do PlantTFDB Gma presentes na matriz de expressão.

```r
# Intersecção: TFs Gma × genes na matriz
regulators <- intersect(gma_tf_ids, rownames(expr_mat))
cat(sprintf("TFs na matriz: %d\n", length(regulators)))
```

---

### 10_metanalysis.R

**IDs GEO relevantes para Glycine max:**
```yaml
# Datasets de soja relacionados a estresse/expressão diferencial
geo_accessions: "GSE99698,GSE107900,GSE143156"
# GSE99698: Glycine max response to Phytophthora sojae
# GSE107900: Soybean drought stress
# GSE143156: Soybean seed development
```

---

### 11_lncrna.R

**Parâmetros para Glycine max:**
```r
# Mesmo critério: comprimento ≥200 nt + ORF máximo <100 aa
# Genome size maior = mais lncRNA candidatos esperados
min_len_nt  <- 200L
max_orf_aa  <- 100L
```

---

## 10. Checklist de Revisão de Código

### Antes de cada `git push`

#### Scripts R

```
□ Todos os outputs inicializados vazios no início?
□ gene_id strip de versão (.1, .2, .Wm82.a4.v1)?
□ prcomp precedido de filtro de variância zero?
□ makeTxDbFromGFF com helper make_txdb() (txdbmaker)?
□ GOseq usando gene2cat customizado (não genome/id)?
□ clusterProfiler sem OrgDb → enricher() com TERM2GENE?
□ GENIE3 com run_genie3() e fallback nCores=1?
□ STRINGdb com species=3847 (Glycine max)?
□ PlantTFDB com prefixo Gma_ (não Ath_)?
□ KEGG com organism="gma" (não "ath")?
□ tryCatch em todas as chamadas biomaRt (pode ter timeout)?
□ suppressPackageStartupMessages nos library()?
□ show_col_types=FALSE nos read_tsv()?
□ Figuras salvas em PDF + PNG (300 dpi)?
```

#### Módulos Nextflow

```
□ Nenhum \n em strings Python dentro de heredoc?
□ rstrip() sem argumento (não rstrip('\n'))?
□ b1.txt/b2.txt: todos os BAMs em uma linha com vírgulas?
□ gffread → sed 's/^chrChr/Chr/' após conversão?
□ LNCRNA_PRED com cadeia de fallback para gffcmp.annotated.gtf?
□ Todos os processos com label de recursos (low_mem/medium_mem/high_mem)?
□ publishDir configurado para todos os processos?
□ Output de figures/ declarado como path("figures/"), emit: figures?
□ Inputs opcionais com sentinel NO_FILE?
```

#### params.yaml

```
□ genome_fasta: path real no servidor?
□ genome_gff3: path real no servidor?
□ control_group: sem hífens (quebra fórmula R ~condition)?
□ treatment_group: sem hífens?
□ kegg_organism: "gma"?
□ read_length: valor correto para o experimento?
□ strandedness: validado pelo Salmon?
```

---

## 11. Configuração do Servidor

### Setup Inicial (uma vez)

```bash
# 1. Clonar repositório
git clone https://github.com/eulaliobqi/glycine-max-rnaseq.git
cd glycine-max-rnaseq

# 2. Instalar Mamba (se não instalado)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh

# 3. Criar ambientes conda
mamba env create -f envs/rnaseq-tools.yml
mamba env create -f envs/r-analysis.yml

# 4. Verificar STAR versão (crítico!)
mamba run -n rnaseq-tools STAR --version   # deve ser 2.7.10b

# 5. Verificar espaço em disco
# STAR index Glycine max: ~30–40 GB
# BAMs (6 amostras × 4–6 GB): ~30 GB
# Pipeline completo: ~100–150 GB total
df -h /home/eulalio/glycine-max-rnaseq/
```

### Execução do Pipeline

```bash
# Execução local completa
bash run_pipeline.sh local

# Retomar após falha (NÃO use --resume; use -resume com um hífen)
bash run_pipeline.sh local -resume

# Execução em SLURM (cluster)
bash run_pipeline.sh slurm

# Forçar re-execução de um processo específico (apagar cache):
rm -rf work/XX/YYYYYY*    # hash do processo do .nextflow.log
bash run_pipeline.sh local -resume
```

### Atualização de Ambientes

```bash
# Após mudanças nos envs/*.yml:
mamba env update -n rnaseq-tools -f envs/rnaseq-tools.yml --prune
mamba env update -n r-analysis -f envs/r-analysis.yml --prune
```

### Dashboard Shiny

```bash
RESULTS_DIR=results Rscript -e "shiny::runApp('dashboard/app.R', port=3838)"
# Acessar: http://servidor:3838
```

---

## 12. Referências e Recursos

### Glycine max Wm82.a4.v1

| Recurso | URL |
|---------|-----|
| Phytozome (genoma + GFF3) | https://phytozome-next.jgi.doe.gov/ |
| SoyBase (alternativa livre) | https://www.soybase.org/ |
| PlantTFDB Gma | http://planttfdb.gao-lab.org/ |
| KEGG gma | https://www.genome.jp/kegg-bin/show_organism?org=gma |
| STRING 3847 | https://string-db.org/ |
| Ensembl Plants biomaRt | https://plants.ensembl.org |
| SoyKB (knowledge base) | https://soykb.org/ |
| GEO – Soybean | https://www.ncbi.nlm.nih.gov/geo/ |

### Publicações de Referência do Pipeline

| Ferramenta | Referência |
|-----------|-----------|
| STAR | Dobin et al. 2013, Bioinformatics |
| featureCounts | Liao et al. 2014, Bioinformatics |
| Salmon | Patro et al. 2017, Nature Methods |
| DESeq2 | Love et al. 2014, Genome Biology |
| clusterProfiler | Yu et al. 2012, OMICS |
| WGCNA | Langfelder & Horvath 2008, BMC Bioinformatics |
| GENIE3 | Huynh-Thu et al. 2010, PLoS ONE |
| rMATS | Shen et al. 2014, PNAS |
| ComBat-Seq | Zhang et al. 2020, NAR Genomics |
| GOseq | Young et al. 2010, Genome Biology |

---

## Apêndice A: Template params.yaml para Glycine max

```yaml
# ============================================================
# params.yaml – Parâmetros biológicos – Glycine max Wm82.a4.v1
# ============================================================

# ── Genoma de referência ─────────────────────────────────────
genome_fasta: "/home/eulalio/glycine-max-rnaseq/Glyma.Wm82.a4.v1.fa"
genome_gff3:  "/home/eulalio/glycine-max-rnaseq/Glyma.Wm82.a4.v1.gene.gff3"
genome_gtf:   ""            # vazio → converte GFF3 automaticamente
genome_index: ""            # vazio → constrói automaticamente

# ── Amostras ──────────────────────────────────────────────────
samplesheet: "samplesheet.csv"

# ── Contraste ─────────────────────────────────────────────────
# IMPORTANTE: sem hífens nos nomes de grupo
control_group:   "WT"
treatment_group: "Tratamento"

# ── Qualidade / Trimagem ──────────────────────────────────────
min_length: 36
quality:    20

# ── Alinhador principal ───────────────────────────────────────
aligner: "star"             # star (2.7.10b!) | hisat2

# ── Quantificação ─────────────────────────────────────────────
feature_type: "exon"
gene_attr:    "gene_id"
strandedness: 0             # 0=unstranded | 1=stranded | 2=reversely-stranded

# ── Salmon (pseudoalinhamento paralelo) ───────────────────────
run_salmon: true

# ── Correção de batch (ComBat-Seq) ────────────────────────────
run_combat_seq: true

# ── Montagem de transcritos (StringTie) ──────────────────────
run_stringtie: false        # true para predição de lncRNAs

# ── Fase 3 ────────────────────────────────────────────────────
run_ml:  true
run_ppi: true
plantfdb_file: ""           # vazio → download automático (Gma_TF_list.txt.gz)
geo_accessions: ""          # ex: "GSE99698,GSE107900"
genie3_trees: 500
genie3_links: 5000

# ── Análise estatística ───────────────────────────────────────
padj_cutoff:   0.05
lfc_cutoff:    1.0

# ── Anotação funcional – Glycine max ─────────────────────────
kegg_organism: "gma"       # código KEGG para Glycine max

# ── rMATS (splicing alternativo) ──────────────────────────────
read_length: 150            # ajustar para o experimento

# ── Relatório ─────────────────────────────────────────────────
report_title:  "RNASeq Insight – Glycine max Wm82.a4.v1"
report_author: "Eulalio Santos – UFV"

# ── Output ────────────────────────────────────────────────────
outdir: "results"
```

---

## Apêndice B: Validação Rápida Pós-Pipeline

```bash
# Verificar tamanhos dos arquivos principais (0 bytes = erro)
find results/ -name "*.tsv" -size 0 | sort
# Esperado: nenhum arquivo 0 bytes nos resultados principais

# Verificar DEGs
wc -l results/deseq2/deseq2_results_sig.tsv
# Esperado: >1 linha (header + genes)

# Verificar enriquecimento GO-BP
wc -l results/enrichment/go_bp_results.tsv
# Esperado: >1 linha

# Verificar GENIE3
wc -l results/genie3/genie3_network.tsv
# Esperado: >1 linha (se reguladores encontrados)

# Verificar GOseq
wc -l results/enrichment/goseq_bp_results.tsv
# Esperado: >1 linha

# Verificar relatório
ls -lh results/report/rnaseq_report.html
# Esperado: >1 MB
```

---

*Documento gerado em 18-May-2026 | Pipeline Arabidopsis TAIR10 v2.0 → Glycine max Wm82.a4.v1 v1.0*
*Autores: Eulalio Santos – UFV | Claude Sonnet 4.6 – Anthropic*
