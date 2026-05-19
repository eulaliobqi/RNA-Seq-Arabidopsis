# Metodologia Computacional — RNA-Seq *Arabidopsis thaliana* TAIR10

**Projeto:** RNASeq Insight Platform — *Arabidopsis thaliana* TAIR10  
**Contraste:** BP12 vs. ColWT  
**Autores:** Eulalio Santos — UFV  
**Versão do pipeline:** v2.0 (Nextflow DSL2 ≥ 24.04)

---

## Sumário

1. [Genoma de Referência](#1-genoma-de-referência)
2. [Controle de Qualidade das Leituras](#2-controle-de-qualidade-das-leituras)
3. [Trimagem de Adaptadores e Filtro de Qualidade](#3-trimagem-de-adaptadores-e-filtro-de-qualidade)
4. [Alinhamento Genômico — STAR](#4-alinhamento-genômico--star)
5. [Pseudoalinhamento — Salmon](#5-pseudoalinhamento--salmon)
6. [Importação de Abundâncias — tximport](#6-importação-de-abundâncias--tximport)
7. [Quantificação por Gene — featureCounts](#7-quantificação-por-gene--featurecounts)
8. [Pré-filtragem de Genes — filterByExpr](#8-pré-filtragem-de-genes--filterbyexpr)
9. [Correção de Efeito de Lote — ComBat-Seq](#9-correção-de-efeito-de-lote--combat-seq)
10. [Expressão Diferencial — DESeq2](#10-expressão-diferencial--deseq2)
11. [Enriquecimento Funcional — GO e KEGG](#11-enriquecimento-funcional--go-e-kegg)
12. [GSEA — Gene Set Enrichment Analysis](#12-gsea--gene-set-enrichment-analysis)
13. [GOseq — Correção de Viés de Comprimento Gênico](#13-goseq--correção-de-viés-de-comprimento-gênico)
14. [WGCNA — Análise de Co-expressão](#14-wgcna--análise-de-co-expressão)
15. [rMATS — Splicing Alternativo](#15-rmats--splicing-alternativo)
16. [StringTie e GffCompare — Montagem de Transcritos](#16-stringtie-e-gffcompare--montagem-de-transcritos)
17. [Predição de lncRNA](#17-predição-de-lncrna)
18. [Integração Multi-ômica](#18-integração-multi-ômica)
19. [PlantTFDB — Classificação de Fatores de Transcrição](#19-planttfdb--classificação-de-fatores-de-transcrição)
20. [GENIE3 — Inferência de Rede Regulatória](#20-genie3--inferência-de-rede-regulatória)
21. [Machine Learning — Seleção de Biomarcadores](#21-machine-learning--seleção-de-biomarcadores)
22. [Rede de Interação Proteína–Proteína — STRINGdb](#22-rede-de-interação-proteína-proteína--stringdb)
23. [Metanálise — Validação Cruzada GEO](#23-metanálise--validação-cruzada-geo)
24. [Relatório e Visualização](#24-relatório-e-visualização)
25. [Parâmetros Consolidados](#25-parâmetros-consolidados)
26. [Softwares e Versões](#26-softwares-e-versões)

---

## 1. Genoma de Referência

O genoma de *Arabidopsis thaliana* ecótipo Columbia (Col-0) na versão **TAIR10** foi utilizado como referência. O arquivo FASTA do genoma (`TAIR10_chr_all.fas`) e a anotação no formato GFF3 (`TAIR10_GFF3_genes.gff3`) foram obtidos do banco de dados oficial TAIR (*The Arabidopsis Information Resource*, www.arabidopsis.org).

A anotação GFF3 foi convertida para o formato GTF utilizando a ferramenta **gffread** (v0.12.x), com a seguinte correção de nomes cromossômicos aplicada via `sed`: prefixos `chrChr` erroneamente inseridos pelo gffread — decorrentes da capitalização `Chr*` dos cromossomos do TAIR10 — foram substituídos pelo prefixo correto `Chr` (i.e., `chrChr5` → `Chr5`), garantindo consistência entre o GTF e os arquivos BAM gerados pelo alinhamento.

```
Cromossomos: Chr1, Chr2, Chr3, Chr4, Chr5, ChrC, ChrM
Total de genes anotados: ~27.416 loci codificantes
```

---

## 2. Controle de Qualidade das Leituras

A qualidade das leituras *paired-end* brutas e pós-trimagem foi avaliada com **FastQC** (v0.12.1), que reporta:
- Distribuição de qualidade Phred por posição da leitura
- Conteúdo de GC e viés de nucleotídeo
- Sequências superrepresentadas e contaminações por adaptadores
- Comprimento de leituras e nível de duplicação

Os relatórios individuais foram agregados com **MultiQC** (v1.21) em um único relatório HTML, gerado separadamente para leituras pré-trimagem (MULTIQC_PRE) e pós-trimagem (MULTIQC_POST).

**Critério de aprovação de qualidade:**
- Qualidade Phred médio ≥ Q20 em >90% das posições
- Taxa de leituras alinhadas pelo STAR ≥ 40% (gate automatizado no pipeline)

---

## 3. Trimagem de Adaptadores e Filtro de Qualidade

A trimagem de adaptadores e o controle de qualidade das leituras foram realizados com **fastp** (v0.23.4) em modo *paired-end*.

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `--length_required` | 36 bp | Comprimento mínimo de leitura após trimagem |
| `--qualified_quality_phred` | 20 | Qualidade Phred mínima por base (Q20) |
| `--detect_adapter_for_pe` | ativado | Detecção automática de adaptadores *paired-end* |
| `--trim_poly_g` | ativado | Remove caudas poli-G (artefato de sequenciadores NextSeq/NovaSeq) |
| `--trim_poly_x` | ativado | Remove caudas de qualquer base homopolimérica |
| `--overrepresentation_analysis` | ativado | Analisa sequências superrepresentadas |

O relatório de qualidade do fastp foi exportado nos formatos JSON e HTML por amostra.

---

## 4. Alinhamento Genômico — STAR

O alinhamento das leituras ao genoma de referência foi realizado com **STAR** (v2.7.10b, *Spliced Transcripts Alignment to a Reference*) em modo *paired-end*, utilizando o algoritmo de alinhamento *splice-aware* baseado em grafos de splicing.

> **Nota sobre versão:** A versão 2.7.11b foi excluída pois a distribuição via bioconda instala um
> binário wrapper que seleciona automaticamente `STAR-avx2`, causando falha de segmentação
> (SIGSEGV, exit 139) em processadores sem suporte à instrução AVX2. A versão 2.7.10b
> distribui um binário único, sem seleção automática de instrução.

### 4.1 Construção do Índice Genômico

```
STAR --runMode genomeGenerate
     --genomeDir star_index/
     --genomeFastaFiles TAIR10_chr_all.fa
     --sjdbGTFfile annotation.gtf
     --genomeChrBinNbits 16
```

O parâmetro `--genomeChrBinNbits 16` é recomendado para genomas menores que o humano, evitando o consumo excessivo de memória RAM.

### 4.2 Alinhamento por Amostra

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `--outSAMtype` | BAM SortedByCoordinate | Saída em BAM ordenado por coordenada |
| `--outSAMattributes` | NH HI AS NM MD | Atributos SAM para análises downstream |
| `--outSAMstrandField` | intronMotif | Inferência de orientação de fita via motivos de íntron |
| `--outFilterIntronMotifs` | RemoveNoncanonical | Remove junções com motivos de splicing não-canônicos |
| `--alignSoftClipAtReferenceEnds` | Yes | Permite *soft-clipping* nas extremidades do genoma |
| `--quantMode` | TranscriptomeSAM GeneCounts | Gera BAM para Salmon e contagens por gene |
| `--limitBAMsortRAM` | 10 GB | Limite de RAM para ordenação do BAM |
| `--readFilesCommand` | zcat | Descompressão automática de arquivos `.fastq.gz` |

Os arquivos BAM foram indexados com **samtools index** (v1.18) para acesso aleatório eficiente.

**Gate de qualidade automatizado:** Amostras com taxa de alinhamento único < 40% são rejeitadas automaticamente pelo pipeline com mensagem de erro explícita.

---

## 5. Pseudoalinhamento — Salmon

Em paralelo ao alinhamento genômico pelo STAR, a quantificação por pseudoalinhamento foi realizada com **Salmon** (v1.10.x), que estima a abundância de transcritos sem a necessidade de alinhamento posicional explícito. Essa análise paralela serve como validação cruzada dos resultados do featureCounts.

### 5.1 Construção do Índice de Transcritos

O arquivo FASTA de transcritos foi gerado a partir do GTF e do genoma via `gffread`, e o índice foi construído com `salmon index`.

### 5.2 Quantificação

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `--libType` | A (auto-detect) | Detecção automática do tipo de biblioteca; SF para stranded-forward, SR para stranded-reverse |
| `--validateMappings` | ativado | Validação de mapeamento para maior precisão |
| `--gcBias` | ativado | Correção de viés de conteúdo GC |
| `--seqBias` | ativado | Correção de viés de sequência nos fragmentos |

A correção de viés de GC e de sequência (*sequence bias correction*) permite estimativas de abundância mais precisas ao modelar as não-uniformidades sistemáticas da captura e sequenciamento de fragmentos.

---

## 6. Importação de Abundâncias — tximport

As estimativas de abundância do Salmon (TPM e contagens) foram importadas para o R com o pacote **tximport** (v1.30+), realizando a agregação de nível transcrito para nível gênico (*transcript-level to gene-level summarization*). O método `"lengthScaledTPM"` é utilizado por padrão, que escala as contagens pelo comprimento efetivo do transcrito normalizado, produzindo valores adequados para análises de expressão diferencial com DESeq2 ou edgeR.

---

## 7. Quantificação por Gene — featureCounts

A quantificação das leituras alinhadas pelo STAR por locus gênico foi realizada com **featureCounts** do pacote **subread** (v2.0.6).

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `-t` | exon | Tipo de *feature* para quantificação |
| `-g` | gene_id | Atributo GTF para agrupamento de *features* |
| `-s` | 0 | *Strandedness* (0 = não-stranded) |
| `-p` | ativado | Modo *paired-end* |
| `--countReadPairs` | ativado | Conta pares de leitura, não leituras individuais |
| `-B` | ativado | Conta apenas pares com ambas as leituras mapeadas |
| `-C` | ativado | Não conta pares com leituras mapeadas em cromossomos diferentes (*chimeric pairs*) |

Os nomes de colunas gerados pelo featureCounts (que incluem o caminho completo do BAM) foram renomeados para os identificadores de amostra correspondentes via script Python (PARSE_COUNTS). O sufixo `.TAIR10` adicionado pelo gffread aos identificadores de gene foi removido por `str.replace('.TAIR10', '')`.

---

## 8. Pré-filtragem de Genes — filterByExpr

Antes da análise de expressão diferencial, genes com contagens muito baixas foram removidos utilizando a função `filterByExpr()` do pacote **edgeR** (v4.4+). Esta abordagem é estatisticamente superior ao simples limiar de soma de leituras (e.g., `rowSums >= 10`), pois considera o tamanho de biblioteca de cada amostra e o delineamento experimental.

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `min.count` | 10 | Contagem mínima em pelo menos um grupo |
| `min.total.count` | 15 | Contagem total mínima entre todas as amostras |
| `group` | `meta$condition` | Variável de grupo para o delineamento experimental |

Apenas genes que satisfazem simultaneamente os critérios de `min.count` e `min.total.count` são retidos para análise. Isso preserva genes com expressão biológica real enquanto remove ruído técnico de baixa contagem.

---

## 9. Correção de Efeito de Lote — ComBat-Seq

O efeito de lote (*batch effect*) foi avaliado e corrigido automaticamente pelo módulo COMBAT_SEQ. A detecção baseia-se em análise de componentes principais (PCA) aplicada aos dados normalizados por VST (*Variance Stabilizing Transformation*):

**Critério de aplicação do ComBat-Seq:**
- PC1 explica > 40% da variância total **E**
- Correlação de Pearson entre PC1 e variável de lote > 0,7

Quando ambas as condições são satisfeitas, **ComBat-Seq** (pacote **sva** v3.52+) é aplicado às contagens brutas. O ComBat-Seq foi desenvolvido especificamente para dados de RNA-Seq (contagens negativas binomiais), ao contrário do ComBat original que opera em dados contínuos.

> **Resultado neste experimento:** Coluna `batch` ausente no samplesheet; batch não detectado.
> Contagens originais mantidas sem correção.

Antes do PCA, genes com variância zero (expressos uniformemente em todas as amostras após VST) foram removidos para evitar erros de escalonamento em `prcomp(scale. = TRUE)`:

```r
vst_mat <- vst_mat[apply(vst_mat, 1, var) > 0, , drop = FALSE]
```

---

## 10. Expressão Diferencial — DESeq2

A análise de expressão diferencial foi realizada com **DESeq2** (v1.46+), que usa um modelo de regressão binomial negativo com estimação empírica de Bayes para estabilização dos estimadores de dispersão.

### 10.1 Modelo Estatístico

```r
design = ~ condition   # fórmula unidirecional: BP12 vs. ColWT
```

### 10.2 Estimação e Shrinkage de LFC

O *fold change* foi estimado com *shrinkage* adaptativo pelo método **apeglm** (*approximate posterior estimation for generalized linear models*), com fallback sequencial para **ashr** e para o resultado sem *shrinkage* (método `"normal"`), caso os métodos anteriores falhem.

O *shrinkage* reduz os LFC de genes com baixa expressão média (onde a estimativa bruta é mais ruidosa) em direção a zero, produzindo estimativas mais estáveis e biologicamente interpretáveis para ordenação e visualização.

### 10.3 Normalização

As contagens normalizadas foram obtidas pela transformação **VST** (`vst()`) com `blind = FALSE`, que usa a dispersão estimada a partir do delineamento para estabilizar a variância ao longo dos valores de expressão média. Recomendada para análises exploratórias (PCA, heatmaps, WGCNA) com mais de 30 genes.

### 10.4 Critérios de Significância

| Parâmetro | Valor | Justificativa |
|-----------|-------|---------------|
| FDR (*p*adj) | < 0,05 | Controla taxa de falsos positivos a 5% (método de Benjamini-Hochberg) |
| \|log₂FC\| | ≥ 1,0 | Equivale a ≥ 2× de mudança de expressão |

O ajuste de múltiplos testes é realizado internamente pelo DESeq2 pelo método de **Benjamini-Hochberg (BH)**, que controla a *False Discovery Rate* (FDR).

### 10.5 Visualizações

| Figura | Descrição |
|--------|-----------|
| PCA | PC1 vs PC2 com `plotPCA(intgroup="condition")` |
| Volcano plot | log₂FC × −log₁₀(padj); linhas de corte em \|LFC\| = 1,0 e padj = 0,05 |
| MA plot | log₁₀(baseMean+1) × log₂FC; genes DEG coloridos |
| Heatmap | Top 50 DEGs por padj menor; Z-score por gene; anotação de condição |

---

## 11. Enriquecimento Funcional — GO e KEGG

A análise de enriquecimento funcional foi realizada com **clusterProfiler** (v4.14+), utilizando o banco de dados de anotações de *Arabidopsis thaliana* disponível em `org.At.tair.db` (Bioconductor 3.20+).

### 11.1 Conversão de Identificadores

Os identificadores TAIR10 foram convertidos para Entrez IDs via `bitr()` com `keyType = "TAIR"` e `toType = "ENTREZID"`, pois o `enrichGO()` do clusterProfiler 4.x não mapeia termos GO diretamente a partir de IDs TAIR em `org.At.tair.db`. O KEGG foi executado com os IDs TAIR originais utilizando `keyType = "kegg"`, que é suportado nativamente pelo banco KEGG REST para *Arabidopsis thaliana* (`organism = "ath"`).

### 11.2 Enriquecimento GO (ORA — *Over-Representation Analysis*)

A análise de sobrerepresentação avalia se determinados termos GO aparecem com frequência maior entre os DEGs do que seria esperado por acaso, pelo teste exato de Fisher (hiperteste).

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `OrgDb` | `org.At.tair.db` | Banco de anotações de *A. thaliana* |
| `keyType` | `"ENTREZID"` | Tipo de identificador de gene |
| `ont` | `"BP"`, `"MF"`, `"CC"` | Ontologias testadas |
| `pAdjustMethod` | `"BH"` | Correção de Benjamini-Hochberg |
| `pvalueCutoff` | 0,05 | *p*-valor ajustado máximo |
| `qvalueCutoff` | 0,20 | *q*-valor máximo (FDR local) |
| `readable` | `TRUE` | Converte Entrez IDs para símbolos gênicos |

### 11.3 Enriquecimento KEGG

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `organism` | `"ath"` | Código KEGG para *Arabidopsis thaliana* |
| `keyType` | `"kegg"` | IDs TAIR nativos do KEGG |
| `pAdjustMethod` | `"BH"` | Correção de Benjamini-Hochberg |
| `pvalueCutoff` | 0,05 | *p*-valor ajustado máximo |
| `use_internal_data` | `FALSE` | Busca dados via API REST (dados atualizados) |

O mapa de similaridade semântica entre termos GO foi calculado com `pairwise_termsim()` e visualizado como grafo de enriquecimento (`emapplot()`, top 30 termos GO-BP).

---

## 12. GSEA — *Gene Set Enrichment Analysis*

A análise GSEA (*Gene Set Enrichment Analysis*) foi realizada para GO-BP e KEGG com **gseGO()** e **gseKEGG()** do clusterProfiler, que utiliza o algoritmo baseado em ranking de Subramanian et al. (2005) implementado via **fgsea** (v1.32+).

### 12.1 Lista Ranqueada

Os genes foram ranqueados em ordem decrescente pelo **log₂ Fold Change** obtido do DESeq2, sem filtragem por significância. Isso preserva toda a informação direcional do experimento, incluindo genes com *p*-valor > 0,05 mas com direção biológica consistente.

```
Ranking: genes ordenados por log₂FC (maior → menor)
Estatística de enriquecimento: ES (Kolmogorov-Smirnov modificado)
NES: ES normalizado pelo tamanho do gene set
```

### 12.2 Parâmetros Estatísticos

| Parâmetro | Valor | Análise | Descrição |
|-----------|-------|---------|-----------|
| `keyType` | `"ENTREZID"` | GSEA GO | Identificadores para GO |
| `keyType` | `"kegg"` | GSEA KEGG | IDs TAIR para KEGG |
| `ont` | `"BP"` | GSEA GO | Processo Biológico |
| `pAdjustMethod` | `"BH"` | Ambas | Benjamini-Hochberg |
| `pvalueCutoff` | 0,05 | Ambas | FDR máximo |
| Mínimo de genes | 10 | Ambas | Tamanho mínimo da lista ranqueada |

---

## 13. GOseq — Correção de Viés de Comprimento Gênico

O enriquecimento GO com correção explícita do viés de comprimento gênico foi realizado com **goseq** (v1.58+), seguindo o método de Young et al. (2010). Genes mais longos têm maior probabilidade de serem detectados como DEGs simplesmente por acumularem mais leituras, o que pode inflar artificialmente o enriquecimento de funções associadas a genes longos.

### 13.1 Comprimentos de Genes

Os comprimentos de genes foram calculados como a soma do comprimento não-redundante dos éxons por locus gênico (*union of exons*), utilizando `exonsBy(txdb, by = "gene")` e `reduce()` para remover sobreposições entre éxons alternativos. O objeto TxDb foi construído a partir do GTF com `txdbmaker::makeTxDbFromGFF()` (adaptação para GenomicFeatures ≥ 1.61.1, onde a função foi movida do pacote GenomicFeatures para txdbmaker).

### 13.2 Função de Ponderação (PWF)

A *Probability Weighting Function* (PWF) descreve a relação entre comprimento do gene e probabilidade de detecção como DEG. Quando comprimentos são suficientemente heterogêneos (≥ 6 valores únicos), a PWF é ajustada com `nullp()` pelo método **Wallenius** (distribuição hipergeométrica ponderada). Quando os comprimentos são uniformes, usa-se a aproximação **Hypergeometric** com uma PWF constante (equivalente ao teste de Fisher sem correção de viés).

| Condição | Método | Justificativa |
|---------|--------|---------------|
| ≥ 6 comprimentos únicos | Wallenius | Correção de viés adequada |
| < 6 comprimentos únicos | Hypergeometric | Evita erro de nós do spline no `nullp()` |

### 13.3 Mapeamento Gene–GO

O mapeamento gene→GO foi obtido via `AnnotationDbi::select(org.At.tair.db, keytype = "TAIR", columns = c("GO", "ONTOLOGY"))`, construindo um `gene2cat` personalizado por ontologia (BP, MF, CC). Esta abordagem substitui a dependência do banco interno do goseq (`genome = "tair10", id = "tair"`), que não está disponível em versões recentes do pacote.

O ajuste de múltiplos testes foi realizado pelo método **BH** sobre os *p*-valores da sobrerepresentação:

```r
res_go |> mutate(p.adjust = p.adjust(over_represented_pvalue, method = "BH"))
```

---

## 14. WGCNA — Análise de Co-expressão Gênica

A análise de co-expressão foi realizada com **WGCNA** (*Weighted Gene Co-expression Network Analysis*, v1.72+), que identifica módulos de genes com padrões de expressão correlacionados entre amostras.

### 14.1 Seleção de Genes

Os `min_genes` = 5.000 genes com maior **desvio absoluto mediano** (MAD) das contagens normalizadas por VST foram selecionados para a análise. A MAD é uma estatística robusta a valores extremos, mais adequada que o desvio padrão para seleção de genes informativos.

```r
mad_vals <- apply(datExpr, 2, mad)
datExpr  <- datExpr[, order(mad_vals, decreasing = TRUE)[1:5000]]
```

### 14.2 Validação das Amostras

A função `goodSamplesGenes()` verificou e removeu amostras ou genes com proporção excessiva de valores ausentes.

### 14.3 Topologia de Rede Livre de Escala — Soft Power

O parâmetro de *soft-thresholding* (β, ou *soft power*) foi selecionado automaticamente como o menor valor de β que produz R² ≥ 0,85 no ajuste do modelo de rede livre de escala (*scale-free topology*):

| Critério | Valor |
|---------|-------|
| R² mínimo (*scale-free fit*) | 0,85 |
| Poderes testados | 1–10, 12, 14, 16, 18, 20 |
| Fallback (caso nenhum β ≥ R²) | β = 6 |

### 14.4 Construção dos Módulos

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `TOMType` | `"unsigned"` | Matriz TOM não-dirigida (ignora direção da correlação) |
| `minModuleSize` | 30 | Número mínimo de genes por módulo |
| `mergeCutHeight` | 0,25 | Distância máxima para fusão de módulos similares (correlação de eigengenes ≥ 0,75) |
| `reassignThreshold` | 0 | Limiar para reatribuição de genes entre módulos |
| `pamRespectsDendro` | `FALSE` | PAM sem restrição de dendrograma (mais módulos detectados) |
| `maxBlockSize` | 20.000 | Bloco máximo de genes por processamento (Arabidopsis: ~27k genes) |

O algoritmo emprega a **Topological Overlap Matrix** (TOM) como medida de dissimilaridade entre genes e o agrupamento hierárquico com ligação média (*average linkage*) para detecção de módulos. Cada módulo recebe uma cor como identificador.

### 14.5 Correlação Módulo-Trait

A correlação de Pearson entre os **eigengenes** de cada módulo (primeiro componente principal do módulo, representativo do padrão de expressão do módulo) e a variável de tratamento (0 = controle, 1 = tratamento) foi calculada com `corPvalueStudent()`. O *p*-valor foi calculado pela distribuição *t* de Student com n−2 graus de liberdade.

### 14.6 Hub Genes

Os *hub genes* de cada módulo foram definidos como os 10 genes com maior **conectividade intramodular** (kWithin), calculada a partir da matriz de adjacência ponderada:

```
kWithin_i = Σ_j (adjacência_ij)   para todo j no mesmo módulo
```

---

## 15. rMATS — Splicing Alternativo

O splicing alternativo foi analisado com **rMATS** (v4.1.0+, *replicate Multivariate Analysis of Transcript Splicing*), que utiliza um modelo estatístico hierárquico Bayesiano para detectar eventos de splicing diferenciais.

### 15.1 Tipos de Eventos

rMATS detecta cinco categorias de splicing alternativo:

| Código | Tipo | Descrição |
|--------|------|-----------|
| SE | *Skipped Exon* | Éxon incluído ou pulado |
| A5SS | *Alternative 5' Splice Site* | Sítio alternativo de 5' do doador |
| A3SS | *Alternative 3' Splice Site* | Sítio alternativo de 3' do receptor |
| MXE | *Mutually Exclusive Exons* | Éxons mutuamente exclusivos |
| RI | *Retained Intron* | Retenção de íntron |

### 15.2 Parâmetros de Execução

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `-t` | `paired` | Biblioteca *paired-end* |
| `--readLength` | 100 bp | Comprimento das leituras do experimento |
| `--novelSS` | ativado | Detecta sítios de splicing novos (não presentes na anotação) |

### 15.3 Estatística e Critérios de Significância

rMATS estima a diferença de inclusão (ΔPSI — *Percent Spliced In*) entre condições usando um modelo estatístico com permutações sobre as contagens de junções (*junction counts*). Os resultados são reportados como:

- **PValue:** *p*-valor nominal do teste
- **FDR:** Taxa de Falsos Positivos ajustada (método BH)
- **IncLevelDifference** (ΔPSI): diferença na proporção de inclusão do evento

**Critérios de filtragem (PARSE_RMATS + RMATS_FILTER):**

| Critério | Limiar | Justificativa |
|---------|--------|---------------|
| FDR | < 0,05 | Controle de taxa de falsos positivos |
| \|ΔPSI\| | > 0,10 | Relevância biológica mínima (10% de mudança de inclusão) |
| Cobertura de junction | ≥ 10 leituras totais | Confiabilidade estatística do evento |

> **Resultado neste experimento:** 9 eventos detectados; 0 aprovados no filtro rigoroso.
> A baixa cobertura de junction (≤ 2 leituras na maioria das réplicas) impediu a detecção
> de eventos com FDR < 0,05.

---

## 16. StringTie e GffCompare — Montagem de Transcritos

A montagem *de novo* de transcritos foi realizada com **StringTie** (v2.2.x) em modo guiado por referência, produzindo um GTF de transcritos montados por amostra. Os GTFs individuais foram fundidos com `stringtie --merge` (guiado pelo GTF de referência) para produzir um assembly de transcritos consenso.

A comparação do assembly com a anotação de referência foi realizada com **GffCompare** (v0.12.x), que classifica cada transcrito montado em categorias (*class codes*) com base em sua sobreposição com transcritos conhecidos:

| Class Code | Significado |
|-----------|-------------|
| `=` | Correspondência exata de éxons |
| `c` | Contido em transcrito de referência |
| `k` | Containment reverso |
| `m`, `n` | Retenção de íntron (multi-éxon) |
| `j` | Nova isoforma com pelo menos um éxon compartilhado |
| `e` | Éxon único sobrepoosto a intron de referência |
| `o` | Outro tipo de sobreposição |
| `s` | *Antisense* com éxons sobrepostos |
| `x` | *Antisense* sem éxons sobrepostos |
| `i` | Transcrito intrónico |
| `u` | Intergênico (não anotado) |

Os transcritos com class codes **u** (intergênico), **i** (intrónico), **x** (*antisense*), **s** (*shadow*) e **o** (outra sobreposição) foram considerados candidatos a transcritos novos potencialmente não-codificantes.

---

## 17. Predição de lncRNA

A predição de lncRNAs (*long non-coding RNAs*) nos transcritos novos foi realizada com base nos critérios de codificação do pacote **Biostrings** (v2.72+) do Bioconductor.

### 17.1 Detecção de ORFs

Para cada transcrito novo, todos os possíveis quadros abertos de leitura (ORFs — *Open Reading Frames*) nas 6 fases de tradução (3 diretas + 3 reversas) foram detectados com `matchPattern()` procurando por codons de início (ATG) e terminação (TAA, TAG, TGA).

### 17.2 Critérios de Classificação como lncRNA

| Critério | Limiar | Base |
|---------|--------|------|
| Comprimento do transcrito | ≥ 200 nt | Definição convencional de lncRNA |
| ORF máximo | < 100 aminoácidos (< 300 nt) | Distingue de mRNAs codificantes |

Um transcrito é classificado como candidato a lncRNA quando **ambos** os critérios são satisfeitos simultaneamente: comprimento ≥ 200 nt **E** maior ORF < 100 aminoácidos.

---

## 18. Integração Multi-ômica

Um *integration score* foi calculado para cada gene testado, integrando evidências de múltiplas camadas analíticas:

### 18.1 Fórmula do Score

```
integration_score =
    lfc_score  × 3,0 +    (peso: magnitude do LFC)
    mean_score × 2,0 +    (peso: nível médio de expressão)
    sig_score  × 2,0 +    (peso: significância estatística)
    has_splicing × 1,5 +  (bônus: evento de splicing alternativo)
    in_pathway   × 1,0 +  (bônus: gene em via GO/KEGG enriquecida)
    is_hub       × 0,5     (bônus: hub gene no WGCNA)
```

### 18.2 Normalização dos Componentes

| Componente | Cálculo | Intervalo |
|-----------|---------|-----------|
| `lfc_score` | `min(|log₂FC| / 5, 1)` | [0, 1] |
| `mean_score` | `min(log₁₀(mean_expr + 1) / 5, 1)` | [0, 1] |
| `sig_score` | `min(−log₁₀(padj) / 10, 1)` | [0, 1] |
| `has_splicing` | binário (0 ou 1) | {0, 1} |
| `in_pathway` | binário (GO-BP ou KEGG) | {0, 1} |
| `is_hub` | binário (top kWithin WGCNA) | {0, 1} |

**Intervalo total do score:** 0 a 10 (componentes ponderados e normalizados).

### 18.3 Candidatos Prioritários

Genes com **evidência em ≥ 2 camadas** (DEG + splicing e/ou via e/ou hub) foram classificados como **candidatos prioritários** (*key candidates*). Os pesos foram definidos para priorizar a magnitude do efeito biológico (LFC e expressão) e a significância estatística, com bônus para convergência de múltiplas evidências independentes.

---

## 19. PlantTFDB — Classificação de Fatores de Transcrição

Os genes diferencialmente expressos foram cruzados com o banco **PlantTFDB** (v5.0, *Plant Transcription Factor Database*) para identificação e classificação de fatores de transcrição (TFs) por família.

### 19.1 Fonte de Dados

- **URL:** http://planttfdb.gao-lab.org/download.php
- **Arquivo:** `Ath_TF_list.txt.gz` (*Arabidopsis thaliana* TF list)
- **Identificadores:** IDs TAIR10 (e.g., AT1G01010)

### 19.2 Análise de Enriquecimento de Famílias TF

Para cada família de TFs, a sobrerepresentação entre os DEGs foi testada pelo **teste exato de Fisher** (unilateral), com correção de múltiplos testes pelo método **BH** (FDR < 0,05).

```
Tabela de contingência por família:
       | DEG | não-DEG |
TF     |  a  |    b    |
não-TF |  c  |    d    |
```

---

## 20. GENIE3 — Inferência de Rede Regulatória

A rede de regulação gênica (TF → gene alvo) foi inferida com **GENIE3** (v1.28+, *GEne Network Inference with Ensemble of trees*), que utiliza *Random Forests* para estimar a importância de cada TF na predição da expressão de cada gene alvo.

### 20.1 Configuração

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `treeMethod` | `"RF"` | Random Forest como método de árvore |
| `K` | `"sqrt"` | Número de variáveis testadas por split: √(n_reguladores) |
| `nTrees` | 500 | Número de árvores por gene alvo |
| `nCores` | `task.cpus` | Paralelismo (com fallback para nCores=1) |
| Reguladores | TFs do PlantTFDB presentes na matriz | Restrição biológica de reguladores |
| Alvos | DEGs (padj < 0,05, \|LFC\| > 1,0) | Genes diferencialmente expressos como alvos |

### 20.2 Pós-processamento

Os pesos de importância foram convertidos em uma lista de links com `getLinkList(reportMax = 5000)`, retendo os 5.000 links com maior peso preditivo. Os **hub TFs** foram definidos como os reguladores com maior grau de saída (*out-degree*) ponderado pela soma de pesos.

---

## 21. Machine Learning — Seleção de Biomarcadores

Três algoritmos de aprendizado de máquina supervisionado foram treinados para classificar amostras (BP12 vs. ColWT) e identificar os genes mais preditivos, utilizando o pacote **caret** (v6.0+) como framework unificado.

### 21.1 Feature Selection

As `n_features = 500` features (genes) foram selecionadas como os DEGs com menor FDR (padj), priorizando os genes mais estatisticamente robustos para a classificação.

### 21.2 Validação Cruzada

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `method` | `"repeatedcv"` | Validação cruzada repetida |
| `number` | min(5, n_amostras_por_classe) | k-folds adaptado ao número de amostras |
| `repeats` | 3 | Número de repetições da CV |
| `summaryFunction` | `twoClassSummary` | Métricas: AUC, Sens, Spec |
| `classProbs` | `TRUE` | Estima probabilidades de classe |
| `metric` | `"ROC"` | AUC como critério de seleção de hiperparâmetros |
| Pré-processamento | `c("center", "scale")` | Centralização e escalonamento (Z-score) |

### 21.3 Modelos e Hiperparâmetros

**Random Forest** (`method = "rf"`, pacote randomForest v4.7+):

| Hiperparâmetro | Grid testado |
|---------------|-------------|
| `mtry` | {2, √p, p/3} onde p = número de features |

**SVM com kernel RBF** (`method = "svmRadial"`, pacote kernlab v0.9+):

| Hiperparâmetro | Seleção |
|---------------|---------|
| `C` (custo), `sigma` | Busca automática pelo caret |

**ElasticNet** (`method = "glmnet"`, pacote glmnet v4.1+):

| Hiperparâmetro | Grid testado |
|---------------|-------------|
| `alpha` | {0 (Ridge), 0,5, 1 (LASSO)} |
| `lambda` | 10⁻⁴ a 10⁰ (8 valores em escala log) |

### 21.4 Importância de Variáveis

A importância de cada gene foi calculada com `varImp(scale = TRUE)`, que normaliza as importâncias para o intervalo [0, 100]. Os top 50 genes por modelo foram exportados. O modelo com maior AUC foi utilizado para a visualização dos top 20 biomarcadores.

### 21.5 Curvas ROC

Curvas ROC (*Receiver Operating Characteristic*) foram geradas com o pacote **pROC** (v1.18+) a partir das predições de probabilidade do conjunto de validação da CV, com a classe positiva definida como o grupo controle (ColWT).

> **Nota:** Com n = 6 amostras e k = 3 folds, cada fold de teste contém apenas 1 amostra,
> tornando o AUC = 1,0 trivialmente alcançável por qualquer modelo (overfitting estrutural
> da CV). Os resultados de `feature_importance.tsv` são biologicamente válidos, porém
> as métricas de AUC/Sens/Spec não devem ser usadas como estimativas de generalização
> para este tamanho de amostra.

---

## 22. Rede de Interação Proteína–Proteína — STRINGdb

A rede de interação proteína–proteína (PPI) foi construída com o pacote R **STRINGdb** (v2.8+), que acessa o banco de dados STRING v11.5.

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `species` | 3702 | Taxon ID de *Arabidopsis thaliana* no NCBI |
| `score_threshold` | 400 | *Combined score* mínimo (escala 0–1000) |
| `version` | 11.5 | Versão do banco STRING |

O *combined score* no STRING integra evidências de co-expressão, fusão gênica, co-ocorrência filogenética, experimentos de interação, mineração de texto e banco de vias metabólicas. Pontuação ≥ 400 corresponde ao limiar "medium confidence".

Os **hub genes** da rede PPI foram definidos como os nós no 10º percentil superior de grau (*degree*) da rede. Métricas topológicas adicionais — **betweenness centrality** e **closeness centrality** — foram calculadas com o pacote **igraph** (v2.0+).

---

## 23. Metanálise — Validação Cruzada GEO

A validação cruzada dos DEGs identificados neste experimento foi realizada por metanálise com datasets públicos do **NCBI GEO** (*Gene Expression Omnibus*), obtidos via pacote **GEOquery** (v2.72+).

Para cada dataset GEO com amostras de *Arabidopsis thaliana*:
1. Os dados de expressão bruta foram normalizados (log₂ ou VSN, conforme a plataforma)
2. A análise de expressão diferencial foi realizada com **limma** (v3.62+) via *linear models*
3. Os DEGs do dataset GEO (FDR < 0,05, \|LFC\| > 0,5) foram comparados com os DEGs do presente experimento
4. A sobreposição foi quantificada pelo **coeficiente de Jaccard** e pelo **teste exato de Fisher**

> **Nota:** IDs GEO não foram fornecidos (`geo_accessions = ""`); o módulo gerou
> instruções sobre datasets sugeridos para metanálise (GSE52979, GSE114065, GSE80460).

---

## 24. Relatório e Visualização

O relatório final integrado foi gerado com **Quarto** (versão atual) em formato HTML interativo. O template inclui:

- Sumário executivo com estatísticas principais
- Volcano plot, MA plot e PCA interativos (via plotly)
- Tabelas de DEGs e enriquecimento filtráveis (via DT)
- Módulos WGCNA com heatmap de correlação
- Ranking de candidatos prioritários

O dashboard interativo foi construído com **Shiny** (v1.8+) e **shinydashboard** (v0.7+), permitindo exploração dinâmica dos resultados.

---

## 25. Parâmetros Consolidados

### Parâmetros Globais do Experimento

| Parâmetro | Valor |
|-----------|-------|
| Organismo | *Arabidopsis thaliana* Col-0 |
| Genoma | TAIR10 |
| Contraste | BP12 vs. ColWT |
| Tipo de biblioteca | *Paired-end* |
| Strandedness | 0 (não-stranded) |
| Comprimento das leituras | 100 bp |
| Réplicas biológicas | 3 por condição (n=6 total) |

### Parâmetros de Filtragem e Significância

| Análise | Parâmetro | Valor |
|--------|-----------|-------|
| Pré-filtragem | min.count (filterByExpr) | 10 |
| Pré-filtragem | min.total.count | 15 |
| DESeq2 | FDR (padj) | < 0,05 |
| DESeq2 | \|log₂FC\| | ≥ 1,0 |
| DESeq2 | LFC shrinkage | apeglm → ashr → normal |
| GO ORA | FDR | < 0,05 |
| GO ORA | q-value | < 0,20 |
| KEGG ORA | FDR | < 0,05 |
| GSEA | FDR | < 0,05 |
| GOseq | FDR | < 0,05 |
| rMATS | FDR | < 0,05 |
| rMATS | \|ΔPSI\| | > 0,10 |
| rMATS | Cobertura mínima | ≥ 10 leituras/evento |
| WGCNA | R² mínimo (soft power) | ≥ 0,85 |
| WGCNA | Mínimo de genes/módulo | 30 |
| WGCNA | Merge height | 0,25 |
| PlantTFDB | FDR (família TF) | < 0,05 |
| GENIE3 | Árvores por gene | 500 |
| GENIE3 | Top links exportados | 5.000 |
| STRING | Combined score mínimo | 400 |
| ML — CV | k-folds | min(5, n_classe) |
| ML — CV | Repetições | 3 |
| ML — CV | Métrica | AUC (ROC) |
| ML — ElasticNet | Alpha testado | 0; 0,5; 1 |
| Integração | Evidências mínimas (*key candidates*) | ≥ 2 camadas |

---

## 26. Softwares e Versões

### Ferramentas Bioinformáticas (ambiente `rnaseq-tools`)

| Software | Versão | Referência |
|---------|--------|-----------|
| Nextflow | ≥ 24.04 (DSL2) | Di Tommaso et al., 2017 |
| FastQC | 0.12.1 | Andrews, 2010 |
| MultiQC | 1.21 | Ewels et al., 2016 |
| fastp | 0.23.4 | Chen et al., 2018 |
| STAR | **2.7.10b** | Dobin et al., 2013 |
| Salmon | 1.10.x | Patro et al., 2017 |
| featureCounts (subread) | 2.0.6 | Liao et al., 2014 |
| samtools | 1.18 | Li et al., 2009 |
| rMATS | ≥ 4.1.0 | Shen et al., 2014 |
| StringTie | 2.2.x | Pertea et al., 2015 |
| GffCompare | 0.12.x | Pertea & Pertea, 2020 |
| gffread | 0.12.x | Pertea & Pertea, 2020 |

### Pacotes R/Bioconductor (ambiente `r-analysis`)

| Pacote | Versão | Referência |
|--------|--------|-----------|
| DESeq2 | ≥ 1.46 | Love et al., 2014 |
| edgeR | ≥ 4.4 | Robinson et al., 2010 |
| tximport | ≥ 1.30 | Soneson et al., 2015 |
| sva (ComBat-Seq) | ≥ 3.52 | Zhang et al., 2020 |
| clusterProfiler | ≥ 4.14 | Yu et al., 2012 |
| enrichplot | ≥ 1.24 | Yu, 2022 |
| fgsea | ≥ 1.32 | Korotkevich et al., 2021 |
| org.At.tair.db | Bioconductor | Carlson, 2023 |
| goseq | ≥ 1.58 | Young et al., 2010 |
| GenomicFeatures | ≥ 1.58 | Lawrence et al., 2013 |
| txdbmaker | ≥ 1.2 | Bioconductor, 2024 |
| WGCNA | ≥ 1.72 | Langfelder & Horvath, 2008 |
| GENIE3 | ≥ 1.28 | Huynh-Thu et al., 2010 |
| STRINGdb | ≥ 2.8 | Franceschini et al., 2013 |
| Biostrings | ≥ 2.72 | Pagès et al., 2023 |
| caret | ≥ 6.0 | Kuhn, 2008 |
| randomForest | ≥ 4.7 | Liaw & Wiener, 2002 |
| kernlab | ≥ 0.9 | Karatzoglou et al., 2004 |
| glmnet | ≥ 4.1 | Friedman et al., 2010 |
| pROC | ≥ 1.18 | Robin et al., 2011 |
| igraph | ≥ 2.0 | Csardi & Nepusz, 2006 |
| limma | ≥ 3.62 | Ritchie et al., 2015 |
| GEOquery | ≥ 2.72 | Davis & Meltzer, 2007 |
| ggplot2 | ≥ 3.5 | Wickham, 2016 |
| pheatmap | ≥ 1.0.12 | Kolde, 2019 |
| EnhancedVolcano | ≥ 1.20 | Blighe et al., 2024 |
| Quarto | atual | Allaire et al., 2024 |

---

## Referências

Blighe K et al. (2024). EnhancedVolcano: Publication-ready volcano plots. Bioconductor.

Carlson M (2023). org.At.tair.db: Genome wide annotation for Arabidopsis. Bioconductor.

Chen S et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17):i884–i890.

Csardi G, Nepusz T (2006). The igraph software package for complex network research. *InterJournal*, Complex Systems:1695.

Davis S, Meltzer PS (2007). GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. *Bioinformatics*, 23(14):1846–1847.

Di Tommaso P et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4):316–319.

Dobin A et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1):15–21.

Ewels P et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19):3047–3048.

Franceschini A et al. (2013). STRING v9.1: protein-protein interaction networks, with increased coverage. *Nucleic Acids Research*, 41:D808–D815.

Friedman J et al. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. *Journal of Statistical Software*, 33(1):1–22.

Huynh-Thu VA et al. (2010). Inferring Regulatory Networks from Expression Data Using Tree-Based Methods. *PLoS ONE*, 5(9):e12776.

Karatzoglou A et al. (2004). kernlab — An S4 Package for Kernel Methods in R. *Journal of Statistical Software*, 11(9):1–20.

Korotkevich G et al. (2021). Fast gene set enrichment analysis. *bioRxiv* 2021.

Kuhn M (2008). Building Predictive Models in R Using the caret Package. *Journal of Statistical Software*, 28(5):1–26.

Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9:559.

Lawrence M et al. (2013). Software for Computing and Annotating Genomic Ranges. *PLoS Computational Biology*, 9(8):e1003118.

Li H et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16):2078–2079.

Liao Y et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7):923–930.

Liaw A, Wiener M (2002). Classification and Regression by randomForest. *R News*, 2(3):18–22.

Love MI et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12):550.

Pagès H et al. (2023). Biostrings: Efficient manipulation of biological strings. Bioconductor.

Patro R et al. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*, 14(4):417–419.

Pertea M et al. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. *Nature Biotechnology*, 33(3):290–295.

Pertea G, Pertea M (2020). GFF Utilities: GffRead and GffCompare. *F1000Research*, 9:304.

Ritchie ME et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7):e47.

Robin X et al. (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics*, 12:77.

Robinson MD et al. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1):139–140.

Shen S et al. (2014). rMATS: Robust and Flexible Detection of Differential Alternative Splicing from Replicate RNA-Seq Data. *PNAS*, 111(51):E5593–E5601.

Soneson C et al. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*, 4:1521.

Subramanian A et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43):15545–15550.

Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, New York.

Young MD et al. (2010). Gene ontology analysis for RNA-seq: accounting for selection bias. *Genome Biology*, 11(2):R14.

Yu G et al. (2012). clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS*, 16(5):284–287.

Zhang Y et al. (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics*, 2(3):lqaa078.

---

*Documento gerado em 18-May-2026 | Pipeline v2.0 | Eulalio Santos — UFV*
