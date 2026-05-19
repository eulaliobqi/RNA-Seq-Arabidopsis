# 🌿 RNASeq Insight – Arabidopsis thaliana TAIR10
# Prompt de Otimização para Claude Code
# Pipeline atual: executando corretamente com HISAT2 + featureCounts + DESeq2 + rMATS + WGCNA + GO/KEGG/GSEA + Integração
# Objetivo: elevar para alto impacto científico com mínimo risco de regressão

---

## 📋 ESTADO ATUAL DO PIPELINE (BASE SÓLIDA)

**Pipeline funcional e validado:**
- Nextflow DSL2 ≥ 24.04 com profiles local/slurm
- fastqc → fastp → fastqc+multiqc (QC robusto)
- HISAT2 2.2.1 + samtools (alinhamento splice-aware)
- featureCounts (subread 2.0.6) → DESeq2 (DE)
- rMATS ≥ 4.1.0 (splicing alternativo: SE, A5SS, A3SS, MXE, RI)
- WGCNA (co-expressão + hub genes + eigengenes)
- clusterProfiler + org.At.tair.db (GO BP/MF/CC + KEGG + GSEA)
- Integração multi-ômica (integration score: lfc×3 + mean×2 + sig×2 + splicing×1.5 + pathway×1.0 + hub×0.5)
- Dashboard Shiny + Relatório Quarto
- Ambientes Mamba: `rnaseq-tools` + `r-analysis`

**Parâmetros atuais:**
- FDR < 0.05, |log2FC| > 1.0
- rMATS: FDR < 0.05, |ΔPSI| > 0.1
- KEGG organism: `ath`, Org.db: `org.At.tair.db`
- Genoma: TAIR10 (`Athaliana_genome.fa` + `Athaliana_annot.gff3`)

**Bugs conhecidos já mitigados:**
- Sufixo `.TAIR10` em IDs: `gsub("\.TAIR10$", "", gene_id)` no `02_enrichment.R`
- `keyType="TAIR"` falha: converter para ENTREZID via `bitr()` antes de `enrichGO`
- rMATS FDR=1.0/contagens zero: fix `sed 's/^chrChr/Chr/'` no processo `GFFREAD`
- LFC shrinkage fallback: `apeglm → ashr → normal`

---

## 🎯 MIGRAÇÃO INCREMENTAL: 3 FASES

### FASE 1 — Alto Impacto / Baixo Risco (~1 semana)
**Princípio: Adicionar sem quebrar o que funciona. Cada item é independente.**

#### 1.1 SUBSTITUIR HISAT2 → STAR
**Justificativa:** STAR é 2-3x mais rápido, mapeia 95.9-99.5% das reads em Arabidopsis, e produz BAMs mais compatíveis com rMATS. HISAT2 está obsoleto para genomas pequenos.

**Implementação segura:**

Criar processo STAR_ALIGN paralelo ao HISAT2 inicialmente, depois remove HISAT2 após validação.

Configuração Nextflow:
- cpus: 16
- memory: 16.GB (Arabidopsis = genoma pequeno, 135 Mb)
- time: 2.h
- conda: bioconda::star=2.7.11b + bioconda::samtools=1.21

Indexação TAIR10:
```bash
STAR --runMode genomeGenerate     --genomeDir star_index_tair10/     --genomeFastaFiles Athaliana_genome.fa     --sjdbGTFfile Athaliana_annot.gtf     --runThreadN 16     --genomeChrBinNbits 16
```

**Validação:** Comparar % mapeamento HISAT2 vs STAR nas mesmas amostras. Esperado: STAR >= HISAT2.

---

#### 1.2 ADICIONAR SALMON (pseudoalinhamento de isoformas)
**Justificativa:** 10-100x mais rápido que alinhamento completo; quantifica isoformas; economiza storage. Kallisto e Salmon têm 98% de overlap em DEGs em Arabidopsis.

**Implementação:**

Processos Nextflow:
- SALMON_INDEX: cpus 8, memory 8.GB
- SALMON_QUANT: cpus 8, memory 8.GB, time 30.min
- conda: bioconda::salmon=1.10.3

**Integração com DESeq2 via tximport (novo script `00_tximport.R`):**

Executar ANTES do 01_deseq2.R. Importar counts de Salmon via tximport, salvar counts e TPM para downstream.

**Decisão de design:** Manter featureCounts (STAR) E Salmon em paralelo inicialmente. Comparar concordância de DEGs. Se >95% overlap, migrar completamente para Salmon como principal.

---

#### 1.3 FILTRAGEM DE GENES DE BAIXA EXPRESSÃO (filterByExpr)
**Justificativa:** Remove ruído antes do DE; evita FDR inflacionado. O pipeline atual usa filtro simples (`rowSums(counts) >= 10`).

**Modificação em `01_deseq2.R`:**

Substituir filtro simples por edgeR filterByExpr:
```r
library(edgeR)
dge <- DGEList(counts = counts, group = meta$condition)
keep <- filterByExpr(dge, group = meta$condition, min.count = 10, min.total.count = 15)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
counts <- dge_filtered$counts
```

**Impacto:** Reduz ~15-30% de genes de baixa expressão, melhorando power estatístico.

---

#### 1.4 CORREÇÃO DE BATCH EFFECTS (ComBat-Seq)
**Justificativa:** Crítico para experimentos com múltiplas corridas de sequenciamento, técnicos, ou locais de coleta. O pipeline atual não detecta batch automaticamente.

**Novo processo `COMBAT_SEQ` no Nextflow:**

- conda: bioconductor-sva=3.52.0 + bioconductor-deseq2=1.46.0
- cpus 2, memory 8.GB

Lógica:
1. PCA antes da correção
2. Detectar batch: se PC1 > 40% variância E correlaciona com batch (cor > 0.7)
3. Se detectado: aplicar ComBat_seq, gerar PCA pós-correção
4. Se não detectado: passar counts originais adiante

**Integração no workflow:** Inserir entre PARSE_COUNTS e DESEQ2.

---

#### 1.5 GOSEQ (correção de viés de tamanho de gene)
**Justificativa:** Genes longos têm mais chance de serem detectados em RNA-Seq. GOseq corrige esse viés. O pipeline atual usa clusterProfiler simples.

**Novo script `02b_goseq.R` (paralelo ao `02_enrichment.R`):**

Usar goseq com bias.data = comprimentos de genes do TAIR10. Rodar em paralelo ao clusterProfiler; comparar resultados.

**Comprimentos de genes:** Extrair do GTF com script Python ou GenomicFeatures (R).

---

### FASE 2 — Alto Impacto / Médio Esforço (~2 semanas)

#### 2.1 WGCNA MELHORADO (MAD filter + blockwiseModules)
**Justificativa:** Pipeline atual usa filtro simples por variância. MAD é mais robusto a outliers. Arabidopsis tem ~27k genes, então blockwise pode ser opcional.

**Modificação em `03_wgcna.R`:**

- Substituir filtro por variância (var) por MAD (Median Absolute Deviation)
- Adicionar detecção automática de soft-threshold power (pickSoftThreshold)
- Selecionar power onde R² > 0.85 na topologia scale-free
- maxBlockSize: 20000 (Arabidopsis permite maior que soja devido ao genoma menor)

---

#### 2.2 STRINGTIE (montagem de transcritos de novo)
**Justificativa:** Identifica transcritos novos e isoformas não anotadas. Complementa Salmon (que usa transcriptoma de referência).

**Novo processo no Nextflow:**

- STRINGTIE: montagem por amostra, cpus 8, memory 8.GB
- GFFCOMPARE: merge e comparação com anotação TAIR10, cpus 4, memory 8.GB
- conda: stringtie=2.2.3 + gffcompare=0.12.6

---

#### 2.3 rMATS COM FILTROS RIGOROSOS
**Justificativa:** Pipeline atual gera todos os eventos. Filtros pós-processamento reduzem falsos positivos.

**Modificação no processo de pós-processamento rMATS:**

Criar processo RMATS_FILTER que:
1. Lê todos os eventos (SE, MXE, A5SS, A3SS, RI)
2. Filtra: FDR < 0.05 E |IncLevelDifference| > 0.1
3. Gera summary por tipo de evento
4. Output: splicing_significant_filtered.tsv + splicing_summary.txt

---

### FASE 3 — Muito Alto Impacto / Alto Esforço (~3-4 semanas)

#### 3.1 GENIE3 (Redes Regulatórias Gene-Specific)
**Justificativa:** Inferência de redes TF→gene. Arabidopsis tem ~1.500 TFs bem caracterizados e database ConnecTF para validação.

**Container Docker (não disponível no Mamba):**

Criar imagem rocker/tidyverse:4.4.0 com GENIE3 instalado via BiocManager.

**Processo Nextflow:**
- container: rnaseq-genie3:latest
- cpus 16, memory 32.GB, time 8.h
- Input: norm_counts_tsv + tf_list_txt (~1.500 TFs de PlantTFDB)
- Output: genie3_weight_matrix.rds + genie3_edges.tsv

**Validação:** Cruzar edges preditos com ConnecTF (database curada de GRNs de Arabidopsis). Reportar % de validação.

---

#### 3.2 PLANTTFDB + CLASSIFICAÇÃO DE TFs
**Justificativa:** Arabidopsis tem ~1.500 TFs em ~60 famílias. Identificar famílias enriquecidas nos DEGs é diferencial.

**Script R `05_tfs.R`:**

1. Download lista de TFs do PlantTFDB (1.534 TFs em 58 famílias)
2. Cruzar com DEGs
3. Enriquecimento de famílias: Fisher's exact test com correção FDR
4. Output: tf_family_counts.tsv + tf_family_enrichment.tsv + plot

---

#### 3.3 PPI + HUB PROTEINS (STRING + Cytoscape headless)
**Justificativa:** Identificar proteínas centrais nas redes de resposta. STRING tem dados de Arabidopsis (organism 3702).

**Processo Nextflow:**
- conda: r-analysis (RCy3)
- Input: lista de DEGs
- Baixar STRING data previamente (3702.protein.links.full.v12.0)
- Filtrar: DEGs em ambos os nós + combined_score >= 400
- Identificar hubs: degree > mean + 2SD
- Exportar para Cytoscape via RCy3 (se disponível)
- Output: ppi_network.tsv + hub_proteins.tsv + ppi_cytoscape.xml

---

#### 3.4 LNCRNAS (CPC2 + CPAT + CNCI)
**Justificativa:** Arabidopsis tem 1.359 ncRNAs anotados no TAIR10. Pipeline atual não explora lncRNAs.

**Container Docker:**

Criar imagem python:3.12-slim com:
- CPC2 (git clone + setup.py install)
- CPAT (pip install)
- CNCI (git clone)
- HMMER (apt-get install)

**Pipeline:**
1. Extrair sequências de transcritos (gffread)
2. CPC2 + CPAT + CNCI (3 preditores)
3. HMMER (Pfam) — verificar domínios conhecidos
4. Consenso: 2/3 preditores + ORF < 300 aa + sem Pfam
5. Output: lncrna_candidates.tsv + lncrna_pipeline_summary.txt

---

#### 3.5 MACHINE LEARNING (caret)
**Justificativa:** Classificar condições/genótipos a partir de perfis de expressão. Validar se o transcriptoma é preditivo.

**Script `06_machine_learning.R`:**

1. Input: TPM matrix (Salmon) + metadata
2. Feature selection: top 1000 genes por variância
3. Split 80/20 (createDataPartition)
4. Modelos: Random Forest + SVM (radial) + Elastic Net
5. Validação: 10-fold cross-validation, metric = ROC
6. Avaliação: AUC no teste
7. Output: model_comparison.txt + top_features_rf.tsv + roc_curve.pdf

**Target:** AUC > 0.85 para considerar transcriptoma preditivo.

---

#### 3.6 METANÁLISE COM DADOS PÚBLICOS (GEO/SRA)
**Justificativa:** Validar DEGs contra milhares de estudos públicos de Arabidopsis.

**Script `07_meta_analysis.R`:**

1. Baixar série GEO relevante (mesmo tratamento/fenótipo)
2. Harmonização com ComBat-seq (batch = estudo)
3. Calcular correlação de log2FC entre estudos
4. Identificar meta-DEGs (consistentes em ambos)
5. Output: meta_degs.tsv + correlation_report.txt

---

## 🔧 MODIFICAÇÕES NO NEXTFLOW CONFIG

### Adições ao nextflow.config existente:

**FASE 1 — Processos de baixo risco:**
- STAR_ALIGN: cpus 16, memory 16.GB, time 2.h
- STAR_INDEX: cpus 16, memory 32.GB, time 1.h
- SALMON_INDEX: cpus 8, memory 8.GB, time 30.min
- SALMON_QUANT: cpus 8, memory 8.GB, time 30.min
- COMBAT_SEQ: cpus 2, memory 8.GB, time 30.min
- GOSEQ: cpus 2, memory 4.GB, time 20.min

**FASE 2 — Processos de médio risco:**
- STRINGTIE: cpus 8, memory 8.GB, time 1.h
- GFFCOMPARE: cpus 4, memory 8.GB, time 30.min
- RMATS_FILTER: cpus 2, memory 4.GB, time 15.min

**FASE 3 — Containers Docker:**
- GENIE3: container rnaseq-genie3:latest, cpus 16, memory 32.GB, time 8.h
- LNCRNA_IDENTIFICATION: container rnaseq-lncrna:latest, cpus 8, memory 16.GB, time 4.h
- STRING_PPI: container rnaseq-ppi:latest, cpus 4, memory 8.GB, time 1.h

---

## 📦 AMBIENTES MAMBA ATUALIZADOS

### `envs/rnaseq-tools.yml` — ADIÇÕES:

**Já existem (manter):**
- fastqc=0.12.1, fastp=0.23.4, multiqc=1.21
- hisat2=2.2.1 (MANTER durante transição)
- samtools=1.21, subread=2.0.6, rmats=4.1.0
- nextflow=24.04

**Novos — Fase 1:**
- star=2.7.11b
- salmon=1.10.3
- kallisto=0.50.1 (alternativa rápida)
- stringtie=2.2.3
- gffcompare=0.12.6

### `envs/r-analysis.yml` — ADIÇÕES:

**Já existem (manter):**
- r-base=4.4.0, r-deseq2=1.46.0, r-wgcna=1.73
- r-ggplot2, r-pheatmap, r-readr, r-dplyr, r-tidyr, r-stringr, r-tibble
- r-ggrepel, r-patchwork, r-optparse
- bioconductor-clusterprofiler=4.14.0, bioconductor-enrichplot=1.24.0
- bioconductor-org.at.tair.db=3.19.1
- quarto=1.5.0, r-shiny, r-plotly, r-dt

**Novos — Fase 1:**
- bioconductor-tximport=1.30.0
- bioconductor-edger=4.4.0
- bioconductor-sva=3.52.0
- bioconductor-goseq=1.58.0
- bioconductor-genomicfeatures=1.58.0
- bioconductor-fgsea=1.32.0
- r-factoextra=1.0.7

**Novos — Fase 3:**
- r-caret=6.0.94
- r-randomforest=4.7.1
- r-kernlab=0.9.33
- r-glmnet=4.1.8
- r-pROC=1.18.5
- r-igraph=2.1.1
- bioconductor-rcy3=2.24.0

---

## 🗂️ ESTRUTURA DE OUTPUT EXPANDIDA

```
results/
├── qc/                          [EXISTE]
│   ├── pre_trim/
│   ├── post_trim/
│   └── multiqc/
├── trimmed/                     [EXISTE]
├── aligned/                     [EXISTE — HISAT2]
├── aligned_star/                [NOVO — STAR]
│   ├── *.bam
│   ├── *.bai
│   └── *.ReadsPerGene.out.tab
├── counts/                      [EXISTE — featureCounts]
│   ├── counts_clean.tsv
│   ├── salmon_counts.tsv        [NOVO]
│   └── salmon_tpm.tsv           [NOVO]
├── genome/                      [EXISTE]
│   ├── star_index/              [NOVO]
│   ├── salmon_index/            [NOVO]
│   └── gene_lengths.tsv         [NOVO — para GOseq]
├── deseq2/                      [EXISTE]
│   ├── deseq2_results_all.tsv
│   ├── deseq2_results_sig.tsv
│   ├── normalized_counts.tsv
│   └── figures/
├── enrichment/                  [EXISTE]
│   ├── go_bp_results.tsv
│   ├── go_mf_results.tsv
│   ├── go_cc_results.tsv
│   ├── kegg_results.tsv
│   ├── gsea_go_results.tsv      [NOVO]
│   ├── gsea_kegg_results.tsv    [NOVO]
│   ├── goseq_results.tsv        [NOVO]
│   └── figures/
├── splicing/                    [EXISTE]
│   ├── splicing_significant.tsv
│   ├── splicing_all.tsv
│   ├── splicing_filtered.tsv    [NOVO — pós-rMATS]
│   └── rmats_output/
├── wgcna/                       [EXISTE]
│   ├── wgcna_modules.tsv
│   ├── wgcna_hub_genes.tsv
│   ├── wgcna_eigengenes.tsv
│   └── figures/
├── stringtie/                   [NOVO]
│   ├── *.gtf
│   ├── gffcmp.annotated.gtf
│   └── gffcmp.stats
├── regulation/                  [NOVO]
│   ├── genie3_weight_matrix.rds
│   ├── genie3_edges.tsv
│   ├── genie3_validation.txt
│   └── tf_family_enrichment.tsv
├── lncrna/                      [NOVO]
│   ├── lncrna_candidates.tsv
│   └── lncrna_pipeline_summary.txt
├── ppi/                         [NOVO]
│   ├── ppi_network.tsv
│   ├── hub_proteins.tsv
│   └── ppi_cytoscape.xml
├── machine_learning/            [NOVO]
│   ├── model_comparison.txt
│   ├── top_features_rf.tsv
│   └── roc_curve.pdf
├── meta_analysis/               [NOVO]
│   ├── meta_degs.tsv
│   └── correlation_report.txt
├── integration/                 [EXISTE]
│   ├── integrated_genes.tsv
│   ├── gene_ranking.tsv
│   ├── key_candidates.tsv
│   └── figures/
├── report/                      [EXISTE]
│   └── rnaseq_report.html
├── dashboard/                   [EXISTE]
│   └── app.R
├── batch_correction/            [NOVO]
│   ├── pca_before.pdf
│   ├── pca_after.pdf
│   └── batch_report.txt
└── nextflow_*.html              [EXISTE]
```

---

## ⚠️ REGRAS DE MIGRAÇÃO (CRÍTICO)

1. **Nunca remover HISAT2 antes de validar STAR:** Manter ambos em paralelo por pelo menos 3 amostras piloto. Comparar:
   - % mapeamento total
   - % uniquely mapped
   - Concordância de DEGs (esperado: >95%)
   - Tempo de execução

2. **featureCounts continua como principal:** Salmon/tximport é secundário inicialmente. Só migrar se concordância >95%.

3. **Manter `02_enrichment.R` intacto:** Adicionar `02b_goseq.R` como paralelo, não substituto. Comparar resultados.

4. **Batch correction só se detectado:** PCA automático decide se aplica ComBat-Seq. Se não detectado, passa adiante sem alterar counts.

5. **Containers Docker para Fase 3:** Criar imagens versionadas (`:v1.0`, `:v1.1`) para reprodutibilidade.

6. **Todos os scripts novos em `scripts/`:** Manter convenção de nomenclatura (`00_`, `01_`, `02_`, etc.).

7. **Parâmetros em `params.yaml`:** Adicionar novos parâmetros com defaults conservadores.

---

## 🎯 CHECKLIST DE IMPLEMENTAÇÃO

### Fase 1 (Semana 1):
- [ ] Criar STAR_INDEX e STAR_ALIGN processos
- [ ] Criar SALMON_INDEX e SALMON_QUANT processos
- [ ] Criar `00_tximport.R`
- [ ] Modificar `01_deseq2.R` (adicionar filterByExpr)
- [ ] Criar COMBAT_SEQ processo
- [ ] Criar `02b_goseq.R`
- [ ] Atualizar `envs/rnaseq-tools.yml` e `envs/r-analysis.yml`
- [ ] Testar com 3 amostras piloto
- [ ] Comparar HISAT2 vs STAR vs Salmon

### Fase 2 (Semanas 2-3):
- [ ] Melhorar `03_wgcna.R` (MAD filter + auto soft-threshold)
- [ ] Criar STRINGTIE e GFFCOMPARE processos
- [ ] Criar RMATS_FILTER processo
- [ ] Validar outputs contra pipeline atual

### Fase 3 (Semanas 4-7):
- [ ] Criar containers Docker (GENIE3, lncRNA, PPI)
- [ ] Criar `05_tfs.R` (PlantTFDB)
- [ ] Criar `06_machine_learning.R` (caret)
- [ ] Criar `07_meta_analysis.R` (GEO/SRA)
- [ ] Criar processos Nextflow para Fase 3
- [ ] Validar redes contra ConnecTF
- [ ] Testar ML com AUC > 0.85

---

## 📊 MÉTRICAS DE SUCESSO

| Métrica | Pipeline Atual | Target Fase 1 | Target Fase 3 |
|---------|---------------|---------------|---------------|
| % mapeamento | HISAT2 | STAR >= HISAT2 | STAR >= 95% |
| DEGs concordantes | 100% (baseline) | >95% vs atual | >95% vs atual |
| Isoformas quantificadas | 0 (featureCounts) | Salmon TPM | Salmon + StringTie |
| Batch correction | Manual | Automático (PCA) | Automático + report |
| GO/KEGG enrichment | clusterProfiler | + GOseq | + GSEA + GOseq |
| Redes regulatórias | 0 | 0 | GENIE3 + ConnecTF validação |
| lncRNAs | 0 | 0 | CPC2/CPAT/CNCI pipeline |
| ML classificação | 0 | 0 | AUC > 0.85 |
| Meta-análise | 0 | 0 | GEO/SRA integração |
| Tempo total (30 amostras) | ~8h | ~6h | ~6h |

---

## 🚀 COMANDO ÚNICO DE INSTALAÇÃO (Fase 1)

```bash
# Atualizar ambiente existente
mamba env update -n rnaseq-tools -f envs/rnaseq-tools.yml --prune
mamba env update -n r-analysis -f envs/r-analysis.yml --prune

# Ou criar novo ambiente completo
mamba create -n rnaseq-arabidopsis-v2 -c conda-forge -c bioconda     star=2.7.11b salmon=1.10.3 kallisto=0.50.1     fastqc=0.12.1 fastp=0.23.4 multiqc=1.21     samtools=1.21 subread=2.0.6 stringtie=2.2.3 gffcompare=0.12.6     rmats=4.1.0 nextflow=24.04     r-base=4.4.0 bioconductor-deseq2=1.46.0 bioconductor-tximport=1.30.0     bioconductor-edger=4.4.0 bioconductor-sva=3.52.0     bioconductor-clusterprofiler=4.14.0 bioconductor-fgsea=1.32.0     bioconductor-goseq=1.58.0 bioconductor-org.at.tair.db=3.19.1     r-wgcna=1.73 r-caret=6.0.94 r-ggplot2=3.5.1     -y
```

---

## 📚 REFERÊNCIAS

- **STAR:** Dobin et al. (2013). *Bioinformatics*, 29(1):15-21.
- **Salmon:** Patro et al. (2017). *Nature Methods*, 14(4):417-419.
- **tximport:** Soneson et al. (2015). *F1000Research*, 4:1521.
- **ComBat-seq:** Zhang et al. (2020). *NAR*, 48(8):e42.
- **GOseq:** Young et al. (2010). *Genome Biology*, 11:R14.
- **GENIE3:** Huynh-Thu et al. (2010). *PLoS ONE*, 5(9):e12776.
- **ConnecTF:** O'Malley et al. (2016). *Cell*, 165(5):1280-1292.
- **PlantTFDB:** Jin et al. (2017). *NAR*, 45(D1):D1042-D1048.
- **WGCNA:** Langfelder & Horvath (2008). *BMC Bioinformatics*, 9:559.
- **rMATS:** Shen et al. (2014). *PNAS*, 111(45):E5593-E5601.
- **nf-core/rnaseq:** https://github.com/nf-core/rnaseq

---

*Prompt gerado para otimização incremental do pipeline RNASeq Insight – Arabidopsis thaliana TAIR10*
*Pipeline atual funcional -> Target: alto impacto científico com mínimo risco de regressão*
