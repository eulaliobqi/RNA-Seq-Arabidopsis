# RNASeq Insight Platform — Guia de Implementação e Modelo para Projetos Futuros

> **Projeto referência:** *Glycine max* (soja IAC100 vs BR16)  
> **Próximo projeto:** *Arabidopsis thaliana* TAIR10 — `C:\Users\eulal\.claude\RNA-Seq-Arabidopsis`  
> **Autor:** Eulalio Santos · UFV  
> **Data de referência:** Maio 2026

---

## 1. Contexto e Visão do Produto

Este projeto implementou um pipeline completo de RNA-Seq vegetal como **produto analítico de alto valor agregado**, não apenas como script de análise. A visão central é:

> **"RNASeq Insight Platform (Plant Edition)"** — transforma dados brutos FASTQ em produto analítico completo: expressão diferencial + splicing alternativo + enriquecimento funcional + co-expressão + integração multi-ômica + dashboard interativo + relatório automatizado.

### Comparativo cultivares analisados (soja)
| Item | Valor |
|------|-------|
| Organismo | *Glycine max* (soja) |
| Contraste | IAC100 vs BR16 |
| Anotação | Phytozome GFF3 |
| KEGG code | `gmx` (corrigido de `gma`) |
| Reads | Paired-end 150 bp |
| Réplicas | 3–4 por grupo |
| Alinhador | HISAT2 2.2.1 (splice-aware) |

---

## 2. Arquitetura do Sistema

```
[INPUT]
FASTQ (R1+R2) + samplesheet.csv + genome FASTA + GFF3/GTF
        ↓
[PIPELINE CORE — Nextflow DSL2]
GFF3→GTF → HISAT2-index → FastQC → fastp → FastQC → HISAT2-align
→ samtools sort/index → featureCounts → rMATS
        ↓
[DATA LAYER]
counts_clean.tsv + sample_metadata.tsv + splicing_significant.tsv
        ↓
[ANÁLISE R]
DESeq2 → GO/KEGG/GSEA → WGCNA → Integração multi-ômica
        ↓
[INTERFACE]
Dashboard Shiny (interativo) + Relatório Quarto (HTML)
        ↓
[OUTPUT FINAL]
gene_ranking.tsv + key_candidates.tsv + rnaseq_report.html
```

### Módulos Nextflow (`modules/`)

| Módulo | Processos | Saída principal |
|--------|-----------|-----------------|
| `qc.nf` | FASTQC, MULTIQC | Relatórios HTML/ZIP |
| `trimming.nf` | FASTP | FASTQ trimados + JSON/HTML |
| `alignment.nf` | GFFREAD, HISAT2_BUILD, HISAT2_ALIGN, SAMTOOLS_SORT_INDEX | BAM + BAI + flagstat |
| `quantification.nf` | FEATURECOUNTS, PARSE_COUNTS | counts_clean.tsv |
| `splicing.nf` | RMATS, PARSE_RMATS | splicing_significant.tsv |
| `analysis.nf` | DESEQ2, ENRICHMENT, WGCNA, INTEGRATION, QUARTO_REPORT | Todos os resultados finais |

---

## 3. Scripts R — Descrição e Parâmetros

### `scripts/01_deseq2.R`

**Função:** Expressão diferencial com DESeq2 + figuras

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `--counts` | — | Matriz de contagens (TSV) |
| `--metadata` | — | Metadados das amostras (TSV) |
| `--control` | — | Nome do grupo controle |
| `--treatment` | — | Nome do grupo tratamento |
| `--padj` | 0.05 | Cutoff FDR |
| `--lfc` | 1.0 | Cutoff log2FoldChange |
| `--outdir` | `.` | Diretório de saída (tabelas) |
| `--figures_dir` | `figures` | Diretório de figuras |

**Figuras geradas:** `pca_samples`, `volcano_plot`, `ma_plot`, `heatmap_top_genes`, `de_barplot` (PDF + PNG)

**Outputs:** `deseq2_results_all.tsv`, `deseq2_results_sig.tsv`, `normalized_counts.tsv`, `deseq2_summary.txt`

**Detalhes importantes:**
- Filtro inicial: remove genes com < 10 reads totais
- LFC shrinkage: tenta `apeglm` → `ashr` → `normal` (fallback automático)
- Normalização VST para visualização
- Classificação: `up` (padj<0.05, lfc>1), `down` (padj<0.05, lfc<−1), `ns`

---

### `scripts/02_enrichment.R`

**Função:** Enriquecimento GO, KEGG e GSEA

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `--deseq2` | — | Resultados DESeq2 (todos os genes) |
| `--norm_counts` | — | Contagens normalizadas |
| `--padj` | 0.05 | Cutoff FDR |
| `--lfc` | 1.0 | Cutoff log2FC |
| `--organism` | `gma` | Código KEGG do organismo |
| `--go_annot` | NULL | Arquivo annotation_info Phytozome (opcional) |
| `--outdir` | `.` | Diretório de saída |
| `--figures_dir` | `figures` | Diretório de figuras |

**Estratégia de anotação (prioridade):**
1. Arquivo Phytozome (`--go_annot`) — mais completo para plantas
2. Cache local `gmax_go_cache.rds` — reutiliza entre execuções
3. biomaRt via Ensembl Plants — com múltiplos hosts de fallback

**Figuras:** `go_bp_dotplot`, `go_bp_emap`, `kegg_dotplot`, `kegg_barplot`, `gsea_go_dotplot`, `gsea_kegg_dotplot`

**Outputs:** `go_bp_results.tsv`, `go_mf_results.tsv`, `go_cc_results.tsv`, `kegg_results.tsv`, `gsea_go_results.tsv`, `gsea_kegg_results.tsv`

**Correções aplicadas no projeto soja:**
- Código KEGG correto para *Glycine max*: `gmx` (não `gma`)
- Limpeza de IDs Phytozome: `Glyma.08G200100.Wm82.a2.v1` → `Glyma.08G200100`
- `org.Gmax.eg.db` removido do Bioconductor 3.20+ → substituído por AnnotationHub
- Fallback gracioso: produz arquivos vazios se sem anotação

---

### `scripts/03_wgcna.R`

**Função:** Análise de co-expressão + detecção de módulos

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `--norm_counts` | — | Contagens VST normalizadas |
| `--metadata` | — | Metadados das amostras |
| `--min_genes` | 5000 | Mínimo genes retidos (por variância) |
| `--soft_power` | 0 | Power manual (0 = auto-detecção) |
| `--outdir` | `.` | Saída tabelas |
| `--figures_dir` | `figures` | Saída figuras |

**Detecção automática de soft power:** R² > 0.8, topologia scale-free

**Outputs:** `wgcna_modules.tsv`, `wgcna_hub_genes.tsv`, `wgcna_eigengenes.tsv`, `wgcna_module_summary.tsv`

**Figuras:** `soft_threshold.pdf`, `sample_clustering_tree.pdf`, `gene_dendrogram_modules.pdf`, `module_trait_heatmap.pdf`

---

### `scripts/04_integration.R`

**Função:** Integração multi-ômica e ranking de genes candidatos

**Fórmula do score de integração (0–10):**
```
integration_score = lfc_score × 3 + mean_score × 2 + sig_score × 2
                  + has_splicing × 1.5 + in_pathway × 1.0 + is_hub × 0.5
```

**Definição de "key candidates":** genes com evidência em ≥ 2 camadas (DE, splicing, vias enriquecidas, hub WGCNA)

**Outputs:** `integrated_genes.tsv`, `gene_ranking.tsv`, `key_candidates.tsv`, `candidates_table.tsv`

**Figuras:** `evidence_layers.{pdf,png}`, `top_genes_integration.{pdf,png}`

---

## 4. Agentes Claude Code Utilizados

O projeto foi construído inteiramente com **Claude Code** (claude-sonnet-4-6), usando os seguintes agentes especializados:

### 4.1 `general-purpose` — Agente de implementação principal

**Usado para:**
- Criar e refinar todos os scripts R (`01_deseq2.R`, `02_enrichment.R`, `03_wgcna.R`, `04_integration.R`)
- Criar módulos Nextflow (`qc.nf`, `trimming.nf`, `alignment.nf`, `quantification.nf`, `splicing.nf`, `analysis.nf`)
- Debugging de erros de execução (leitura de logs + correção)
- Refatoração iterativa baseada em outputs reais

**Estratégia eficaz:** Passar logs de erro completos + contexto do script para diagnóstico preciso

---

### 4.2 `Explore` — Agente de exploração de codebase

**Usado para:**
- Localizar definições de funções/variáveis em múltiplos arquivos
- Verificar consistência de parâmetros entre `main.nf` e módulos
- Identificar onde IDs de genes eram manipulados (bug de normalização)
- Auditar estrutura de diretórios de resultados

**Estratégia eficaz:** Usar `search breadth: very thorough` para bugs cross-file

---

### 4.3 `Plan` — Agente de planejamento arquitetural

**Usado para:**
- Planejar estrutura modular Nextflow antes da implementação
- Definir interfaces entre módulos (channels, outputs)
- Desenhar estratégia de fallback para anotações GO/KEGG
- Planejar fórmula de integration score

---

## 5. Skills Utilizadas

### 5.1 `simplify` — Revisão e simplificação de código

**Quando usar:** Após implementação funcional de um script, antes de integrar ao pipeline

**O que detectou no projeto soja:**
- Código de limpeza de IDs Phytozome duplicado em funções diferentes
- Imports R desnecessários aumentando tempo de carregamento
- Lógica de fallback do biomaRt que podia ser simplificada

---

### 5.2 `fewer-permission-prompts` — Redução de prompts de permissão

**Quando usar:** Após as primeiras sessões de execução do pipeline

**O que adicionou:** Allowlist de comandos Bash read-only frequentes (ls, cat de logs, grep em resultados)

---

### 5.3 `init` — Inicialização de CLAUDE.md

**Quando usar:** No início do projeto Arabidopsis, antes de implementar qualquer código

**Gera:** Documentação do codebase, convenções, stack tecnológico — contexto automático para todas as sessões futuras

---

## 6. O Que Funcionou — Lições Validadas

### 6.1 Bugs corrigidos e soluções aplicadas

| Bug | Causa | Solução |
|-----|-------|---------|
| GO/KEGG vazio para *Glycine max* | `org.Gmax.eg.db` removido Bioconductor 3.20+ | Substituição por AnnotationHub + biomaRt Ensembl Plants |
| KEGG sem resultados | Código `gma` inválido | Código correto: `gmx` |
| IDs não matchavam entre datasets | Sufixos Phytozome (.Wm82.a2.v1) | Limpeza com `sub("\\.Wm82.*", "", id)` antes de qualquer join |
| featureCounts: nomes de colunas = paths BAM | Comportamento padrão featureCounts | Script Python `PARSE_COUNTS` renomeia para sample names do samplesheet |
| rMATS falha com menos de 2 replicatas | rMATS exige ≥2 BAMs por grupo | Validação no `main.nf` antes de chamar RMATS |
| Quarto: arquivos não encontrados | Report recebia dirs de figuras, não arquivos TSV | Processo `QUARTO_REPORT` refatorado: copia TSV para subpastas `*_data/` |

### 6.2 Decisões arquiteturais que funcionaram bem

**Conda/Mamba com dois ambientes separados:**
- `rnaseq-tools`: ferramentas bioinformáticas (Nextflow, HISAT2, featureCounts, rMATS, fastp, samtools)
- `r-analysis`: R + Bioconductor (DESeq2, clusterProfiler, WGCNA, Shiny, Quarto)
- Evita conflitos de dependências entre Python e R

**`params.yaml` separado do `nextflow.config`:**
- `nextflow.config`: configuração técnica (recursos, profiles, processos)
- `params.yaml`: parâmetros biológicos (genome paths, contraste, cutoffs, organismo)
- Permite reutilizar o mesmo config técnico para diferentes experimentos

**LFC shrinkage com fallback automático:**
```r
tryCatch(
  lfcShrink(dds, coef=2, type="apeglm"),
  error = function(e) tryCatch(
    lfcShrink(dds, coef=2, type="ashr"),
    error = function(e2) results(dds)
  )
)
```

**Cache de anotações GO:**
- `gmax_go_cache.rds` salvo após primeira consulta biomaRt
- Economiza 5–15 min por re-execução do pipeline
- Aplicar no Arabidopsis com `at_go_cache.rds`

**Integration score ponderado:**
- Prioriza DE forte (peso 3) + expressão base (peso 2) + significância (peso 2)
- Bônus por splicing (1.5) + vias (1.0) + hub WGCNA (0.5)
- Escala 0–10 normalizada facilita interpretação

### 6.3 Validações que evitaram erros silenciosos

- Taxa de alinhamento HISAT2: gate mínimo 40% (script para com erro se abaixo)
- `goodSamplesGenes()` do WGCNA remove outliers automáticamente
- `parse_args()` com `stop()` em parâmetros obrigatórios faltando
- Verificação de FASTQ antes de criar canais no `main.nf`

### 6.4 Dashboard Shiny — O que funcionou

- Carregamento lazy (lê TSV só quando aba é acessada)
- `tryCatch` em cada `read_tsv()` — não quebra se arquivo não existe
- Plotly para interatividade: zoom, hover, toggle de legenda
- DT com `dom = 'Bfrtip'` para exportar CSV/Excel direto do browser

---

## 7. Stack Tecnológico Completo

| Camada | Ferramenta | Versão | Conda env |
|--------|-----------|--------|-----------|
| Orquestração | Nextflow DSL2 | ≥ 24.04 | rnaseq-tools |
| QC | FastQC | 0.12.1 | rnaseq-tools |
| QC agregado | MultiQC | 1.21 | rnaseq-tools |
| Trimagem | fastp | 0.23.4 | rnaseq-tools |
| Conversão anotação | gffread | 0.12.7 | rnaseq-tools |
| Alinhamento | HISAT2 | 2.2.1 | rnaseq-tools |
| BAM | samtools | 1.20 | rnaseq-tools |
| Quantificação | featureCounts (subread) | 2.0.6 | rnaseq-tools |
| Splicing alternativo | rMATS | ≥ 4.1.0 | rnaseq-tools |
| Expressão diferencial | DESeq2 | ≥ 1.40 | r-analysis |
| Enriquecimento | clusterProfiler + enrichplot | ≥ 4.8 + 1.20 | r-analysis |
| GSEA | fgsea | ≥ 1.24 | r-analysis |
| Anotação | biomaRt + AnnotationHub | ≥ 2.54 + 3.8 | r-analysis |
| Co-expressão | WGCNA | Bioconductor | r-analysis |
| Visualização | ggplot2 + pheatmap + ggrepel | 3.4+ | r-analysis |
| Dashboard | Shiny + shinydashboard + plotly + DT | 1.7+ | r-analysis |
| Relatório | Quarto | Latest | r-analysis |

---

## 8. Estrutura de Resultados Produzidos

```
results/
├── qc/pre_trim/          # FastQC HTML/ZIP por amostra
├── qc/post_trim/         # FastQC pós-trimagem
├── qc/multiqc/           # Relatórios MultiQC consolidados
├── trimmed/              # FASTQ trimados + relatórios fastp
├── aligned/              # BAM sorted+indexed + flagstat + logs HISAT2
├── counts/               # counts_matrix.txt + counts_clean.tsv + sample_metadata.tsv
├── splicing/rmats_output/ # Arquivos brutos rMATS (SE, A5SS, A3SS, MXE, RI)
├── splicing/             # splicing_significant.tsv + splicing_all.tsv
├── deseq2/               # Tabelas DE + normalized_counts.tsv + figures/
├── enrichment/           # GO/KEGG/GSEA TSVs + figures/ + cache RDS
├── wgcna/                # Módulos + hub genes + eigengenes + figures/
├── integration/          # gene_ranking.tsv + key_candidates.tsv + figures/
├── report/               # rnaseq_report.html (Quarto)
├── genome/               # annotation.gtf + hisat2_index/ (se gerado)
├── nextflow_report.html  # Métricas de execução Nextflow
├── nextflow_trace.txt    # CPU/memória/tempo por processo
├── nextflow_timeline.html
└── nextflow_dag.html
```

---

## 9. Adaptação para *Arabidopsis thaliana* TAIR10

### 9.1 Mudanças obrigatórias em `params.yaml`

```yaml
# Trocar
samplesheet: "samplesheet.csv"
genome_fasta: "/path/to/TAIR10_genome.fa"
genome_gff3: "/path/to/TAIR10_GFF3.gff3"   # ou genome_gtf: TAIR10.gtf
genome_index: ""                              # deixar vazio para gerar

contrast: "seu_tratamento_vs_controle"
control_group: "WT"
treatment_group: "mutante"

kegg_organism: "ath"            # Arabidopsis thaliana (substituir "gmx")
report_title: "RNASeq – Arabidopsis thaliana TAIR10"
report_author: "Eulalio Santos"
```

### 9.2 Mudanças em `scripts/02_enrichment.R`

```r
# Organismo KEGG
# gmx → ath

# biomaRt dataset
# "soybean_gene_ensembl" → "athaliana_eg_gene"
# host: "https://plants.ensembl.org"

# Cache local
# "gmax_go_cache.rds" → "at_go_cache.rds"

# Limpeza de IDs
# sub("\\.Wm82.*", "", id)  →  remover apenas se IDs TAIR tiverem sufixos
# IDs TAIR10 já limpos: AT1G01010, AT3G55280, etc.
```

### 9.3 Package de anotação Arabidopsis

O `org.At.tair.db` ainda está disponível no Bioconductor 3.20+:

```r
# Em install_r_bioc_packages.R
BiocManager::install("org.At.tair.db")

# Em 02_enrichment.R
library(org.At.tair.db)
# keytype = "TAIR" para IDs AT*G*
```

### 9.4 Checklist de adaptação

- [ ] Baixar genoma TAIR10: `TAIR10_chr_all.fas` + `TAIR10_GFF3_genes.gff`
- [ ] Atualizar `params.yaml`: paths genoma + `kegg_organism: "ath"`
- [ ] Atualizar `scripts/02_enrichment.R`: dataset biomaRt + org.db + cache name
- [ ] Atualizar `scripts/install_r_bioc_packages.R`: trocar `org.Gmax.eg.db` → `org.At.tair.db`
- [ ] Criar `samplesheet.csv` com amostras Arabidopsis
- [ ] Ajustar `contrast`, `control_group`, `treatment_group` no `params.yaml`
- [ ] Atualizar título/autor no `params.yaml`
- [ ] Rodar `bash setup.sh` para criar ambientes conda no novo diretório
- [ ] Testar com `bash run_pipeline.sh test` antes da execução completa

### 9.5 O que NÃO muda (reutilizar integralmente)

- `main.nf` — lógica de orquestração é organismo-agnóstica
- `nextflow.config` — profiles local/slurm/test permanecem iguais
- `modules/qc.nf`, `modules/trimming.nf`, `modules/alignment.nf` — genéricos
- `modules/quantification.nf` — featureCounts é genérico
- `modules/splicing.nf` — rMATS é genérico (só GTF muda)
- `scripts/01_deseq2.R` — DESeq2 é totalmente genérico
- `scripts/03_wgcna.R` — WGCNA é totalmente genérico
- `scripts/04_integration.R` — fórmula de integration score é genérica
- `dashboard/app.R` — Shiny lê TSVs padronizados, sem código organismo-específico
- `setup.sh`, `run_pipeline.sh` — scripts de execução são genéricos
- `envs/rnaseq-tools.yml` — ferramentas bioinformáticas são genéricas
- `envs/r-analysis.yml` — apenas adicionar `r-org.at.tair.db`

---

## 10. Referências de Genoma TAIR10

| Recurso | URL |
|---------|-----|
| Genoma FASTA | TAIR: `https://www.arabidopsis.org/download` |
| GFF3 anotação | TAIR: mesmo site, seção "Genome annotation" |
| Ensembl Plants | `https://plants.ensembl.org/Arabidopsis_thaliana` |
| biomaRt dataset | `"athaliana_eg_gene"` (Ensembl Plants) |
| KEGG organism | `ath` |
| org.db | `org.At.tair.db` (Bioconductor) |

---

## 11. Fluxo Recomendado para Iniciar o Projeto Arabidopsis

```bash
# 1. Criar diretório
mkdir "C:\Users\eulal\.claude\RNA-Seq-Arabidopsis"
cd "C:\Users\eulal\.claude\RNA-Seq-Arabidopsis"

# 2. Copiar estrutura do projeto soja como template
cp -r /mnt/c/Users/eulal/.claude/soja-iac/* .

# 3. Editar params.yaml com paths Arabidopsis
# 4. Atualizar scripts/02_enrichment.R (organismo)
# 5. Atualizar scripts/install_r_bioc_packages.R (org.At.tair.db)
# 6. Criar samplesheet.csv com amostras do experimento
# 7. Executar setup
bash setup.sh

# 8. Validar com dados de teste
bash run_pipeline.sh test

# 9. Executar pipeline completo
bash run_pipeline.sh local
```

---

## 12. Prompt de Instrução para Claude Code — Projeto Arabidopsis

Ao iniciar o projeto Arabidopsis no Claude Code, usar este contexto:

```
Você está construindo um pipeline RNA-Seq para Arabidopsis thaliana TAIR10.
O projeto é uma adaptação do pipeline de Glycine max (soja) já implementado e validado.

Referência: /mnt/c/Users/eulal/.claude/soja-iac/
Destino: C:\Users\eulal\.claude\RNA-Seq-Arabidopsis\

Mudanças necessárias em relação ao projeto soja:
1. kegg_organism: "ath" (era "gmx")
2. biomaRt dataset: "athaliana_eg_gene" (era "soybean_gene_ensembl")
3. org.db: org.At.tair.db (era AnnotationHub para Gmax)
4. IDs TAIR: AT[1-5MC]G[0-9]{5} — não precisam de limpeza de sufixos
5. Cache GO: "at_go_cache.rds"
6. Genoma: TAIR10_chr_all.fas + TAIR10_GFF3_genes.gff

O restante da arquitetura (Nextflow, R scripts estruturais, dashboard, Quarto) é reutilizado sem modificação.

Critérios de análise mantidos:
- FDR < 0.05, |log2FC| > 1.0
- Splicing: FDR < 0.05, |ΔPSI| > 0.1
- Integration score: lfc×3 + mean×2 + sig×2 + splicing×1.5 + pathway×1.0 + hub×0.5
- Key candidates: genes com evidência em ≥ 2 camadas
```

---

*Documento gerado em Maio 2026 — serve como modelo permanente para projetos RNA-Seq de plantas baseados nesta plataforma.*
