# ============================================================
# RNASeq Insight Dashboard – Arabidopsis thaliana TAIR10
# Execute: RESULTS_DIR=results Rscript -e "shiny::runApp('dashboard/app.R', port=3838)"
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(plotly)
  library(DT)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(pheatmap)
  library(scales)
  library(tidyr)
})

RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = "results")

load_tsv_safe <- function(path, ...) {
  full <- file.path(RESULTS_DIR, path)
  if (!file.exists(full)) return(NULL)
  tryCatch(read_tsv(full, show_col_types = FALSE, ...), error = function(e) NULL)
}

# ── UI ────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "RNASeq – A. thaliana"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Visão Geral",          tabName = "overview",    icon = icon("chart-bar")),
      menuItem("Expressão Diferencial",tabName = "de",          icon = icon("dna")),
      menuItem("Enriquecimento",        tabName = "enrichment",  icon = icon("project-diagram")),
      menuItem("Splicing",              tabName = "splicing",    icon = icon("code-branch")),
      menuItem("WGCNA",                 tabName = "wgcna",       icon = icon("network-wired")),
      menuItem("Integração",            tabName = "integration", icon = icon("layer-group")),
      menuItem("Dados Brutos",          tabName = "raw",         icon = icon("table"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f4f4; }
      .box { border-radius: 6px; }
    "))),
    tabItems(

      # ── Visão Geral ──────────────────────────────────────────
      tabItem(tabName = "overview",
        fluidRow(
          valueBoxOutput("vbox_total_de",  width = 3),
          valueBoxOutput("vbox_up",        width = 3),
          valueBoxOutput("vbox_down",      width = 3),
          valueBoxOutput("vbox_splicing",  width = 3)
        ),
        fluidRow(
          box(title = "PCA – Amostras", status = "primary", solidHeader = TRUE,
              width = 6, plotlyOutput("pca_plot", height = "380px")),
          box(title = "DEGs por Regulação", status = "success", solidHeader = TRUE,
              width = 6, plotlyOutput("de_barplot", height = "380px"))
        )
      ),

      # ── Expressão Diferencial ────────────────────────────────
      tabItem(tabName = "de",
        fluidRow(
          box(width = 3, solidHeader = TRUE, status = "warning", title = "Filtros",
            sliderInput("padj_filter", "FDR máximo", min = 0.001, max = 0.2, value = 0.05, step = 0.001),
            sliderInput("lfc_filter",  "|log2FC| mínimo", min = 0, max = 5, value = 1, step = 0.1),
            radioButtons("reg_filter", "Regulação",
                         choices = c("Todos"="all","Up"="up","Down"="down"), selected = "all")
          ),
          box(title = "Volcano Plot", status = "primary", solidHeader = TRUE,
              width = 9, plotlyOutput("volcano_plot", height = "480px"))
        ),
        fluidRow(
          box(title = "Genes Expressos Diferencialmente", width = 12,
              solidHeader = TRUE, status = "info",
              DTOutput("de_table"))
        )
      ),

      # ── Enriquecimento ───────────────────────────────────────
      tabItem(tabName = "enrichment",
        fluidRow(
          box(title = "GO Biological Process", status = "success", solidHeader = TRUE,
              width = 6, plotlyOutput("go_bp_plot", height = "400px")),
          box(title = "KEGG Pathways", status = "warning", solidHeader = TRUE,
              width = 6, plotlyOutput("kegg_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "Resultados GO-BP", width = 6, solidHeader = TRUE, status = "info", DTOutput("go_table")),
          box(title = "Resultados KEGG",  width = 6, solidHeader = TRUE, status = "info", DTOutput("kegg_table"))
        )
      ),

      # ── Splicing ─────────────────────────────────────────────
      tabItem(tabName = "splicing",
        fluidRow(
          box(title = "Eventos de Splicing Alternativo", width = 12,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("splicing_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Tabela de Eventos Significativos", width = 12,
              solidHeader = TRUE, status = "info", DTOutput("splicing_table"))
        )
      ),

      # ── WGCNA ────────────────────────────────────────────────
      tabItem(tabName = "wgcna",
        fluidRow(
          box(title = "Módulos WGCNA – Resumo", status = "primary", solidHeader = TRUE,
              width = 6, plotlyOutput("wgcna_module_plot", height = "380px")),
          box(title = "Hub Genes", status = "success", solidHeader = TRUE,
              width = 6, DTOutput("hub_genes_table"))
        )
      ),

      # ── Integração ───────────────────────────────────────────
      tabItem(tabName = "integration",
        fluidRow(
          box(title = "Top 30 – Integration Score", status = "primary", solidHeader = TRUE,
              width = 8, plotlyOutput("integration_plot", height = "500px")),
          box(title = "Candidatos Chave", status = "success", solidHeader = TRUE,
              width = 4,
              valueBoxOutput("n_candidates", width = 12),
              p("Genes com evidência em ≥ 2 camadas: DE + splicing + vias + hub WGCNA")
          )
        ),
        fluidRow(
          box(title = "Tabela de Candidatos", width = 12,
              solidHeader = TRUE, status = "info", DTOutput("candidates_table"))
        )
      ),

      # ── Dados Brutos ─────────────────────────────────────────
      tabItem(tabName = "raw",
        fluidRow(
          box(title = "Todos os Resultados DESeq2", width = 12,
              solidHeader = TRUE, status = "info",
              DTOutput("raw_deseq2_table"))
        )
      )

    ) # tabItems
  ) # dashboardBody
) # dashboardPage

# ── Server ────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Carregamento lazy dos dados
  de_data      <- reactive({ load_tsv_safe("deseq2/deseq2_results_all.tsv") })
  norm_data    <- reactive({ load_tsv_safe("deseq2/normalized_counts.tsv") })
  meta_data    <- reactive({ load_tsv_safe("counts/sample_metadata.tsv") })
  go_bp_data   <- reactive({ load_tsv_safe("enrichment/go_bp_results.tsv") })
  kegg_data    <- reactive({ load_tsv_safe("enrichment/kegg_results.tsv") })
  splicing_data<- reactive({ load_tsv_safe("splicing/splicing_significant.tsv") })
  wgcna_data   <- reactive({ load_tsv_safe("wgcna/wgcna_modules.tsv") })
  hub_data     <- reactive({ load_tsv_safe("wgcna/wgcna_hub_genes.tsv") })
  ranking_data <- reactive({ load_tsv_safe("integration/gene_ranking.tsv") })
  cand_data    <- reactive({ load_tsv_safe("integration/key_candidates.tsv") })

  de_filtered <- reactive({
    req(de_data())
    df <- de_data() |> filter(!is.na(padj))
    if (input$reg_filter != "all") df <- filter(df, regulation == input$reg_filter)
    filter(df, padj <= input$padj_filter, abs(log2FoldChange) >= input$lfc_filter)
  })

  # Value boxes
  output$vbox_total_de <- renderValueBox({
    n <- if (!is.null(de_data())) sum(de_data()$padj < 0.05, na.rm = TRUE) else 0
    valueBox(n, "DEGs Totais", icon = icon("dna"), color = "blue")
  })
  output$vbox_up <- renderValueBox({
    n <- if (!is.null(de_data())) sum(de_data()$regulation == "up", na.rm = TRUE) else 0
    valueBox(n, "Up-regulated", icon = icon("arrow-up"), color = "red")
  })
  output$vbox_down <- renderValueBox({
    n <- if (!is.null(de_data())) sum(de_data()$regulation == "down", na.rm = TRUE) else 0
    valueBox(n, "Down-regulated", icon = icon("arrow-down"), color = "blue")
  })
  output$vbox_splicing <- renderValueBox({
    n <- if (!is.null(splicing_data())) nrow(splicing_data()) else 0
    valueBox(n, "Eventos Splicing", icon = icon("code-branch"), color = "purple")
  })

  # PCA
  output$pca_plot <- renderPlotly({
    req(norm_data(), meta_data())
    mat    <- norm_data() |> column_to_rownames("gene_id") |> as.matrix()
    pca    <- prcomp(t(mat), scale. = TRUE)
    pct    <- round(summary(pca)$importance[2,] * 100, 1)
    df_pca <- data.frame(
      PC1  = pca$x[,1], PC2 = pca$x[,2],
      sample    = rownames(pca$x),
      condition = meta_data()$condition[match(rownames(pca$x), meta_data()$sample)]
    )
    plot_ly(df_pca, x = ~PC1, y = ~PC2, color = ~condition, text = ~sample,
            type = "scatter", mode = "markers+text", marker = list(size = 12),
            textposition = "top center") |>
      layout(xaxis = list(title = sprintf("PC1 (%.1f%%)", pct[1])),
             yaxis = list(title = sprintf("PC2 (%.1f%%)", pct[2])),
             title = "PCA – Amostras")
  })

  # DEG barplot
  output$de_barplot <- renderPlotly({
    req(de_data())
    df <- de_data() |>
      filter(!is.na(regulation)) |>
      count(regulation) |>
      filter(regulation != "ns")
    plot_ly(df, x = ~regulation, y = ~n, type = "bar",
            color = ~regulation,
            colors = c("up" = "#E41A1C", "down" = "#377EB8")) |>
      layout(title = "DEGs por Regulação", yaxis = list(title = "Contagem"))
  })

  # Volcano
  output$volcano_plot <- renderPlotly({
    req(de_data())
    df <- de_data() |> filter(!is.na(padj), !is.na(log2FoldChange))
    df$color <- case_when(
      df$padj < input$padj_filter & df$log2FoldChange >  input$lfc_filter ~ "Up",
      df$padj < input$padj_filter & df$log2FoldChange < -input$lfc_filter ~ "Down",
      TRUE ~ "ns"
    )
    plot_ly(df, x = ~log2FoldChange, y = ~-log10(padj),
            color = ~color, colors = c("Up"="#E41A1C","Down"="#377EB8","ns"="grey70"),
            text = ~gene_id, type = "scatter", mode = "markers",
            marker = list(size = 4, opacity = 0.7)) |>
      layout(title = "Volcano Plot",
             xaxis = list(title = "log2 Fold Change"),
             yaxis = list(title = "-log10(padj)"),
             shapes = list(
               list(type="line", x0=-input$lfc_filter, x1=-input$lfc_filter, y0=0, y1=15,
                    line=list(dash="dash", color="grey")),
               list(type="line", x0=input$lfc_filter, x1=input$lfc_filter, y0=0, y1=15,
                    line=list(dash="dash", color="grey")),
               list(type="line", x0=-10, x1=10,
                    y0=-log10(input$padj_filter), y1=-log10(input$padj_filter),
                    line=list(dash="dash", color="grey"))
             ))
  })

  # Tabela DE
  output$de_table <- renderDT({
    req(de_filtered())
    df <- de_filtered() |>
      select(gene_id, log2FoldChange, baseMean, padj, regulation) |>
      mutate(across(where(is.numeric), \(x) round(x, 4)))
    datatable(df, extensions = "Buttons",
              options = list(pageLength = 15, dom = "Bfrtip",
                             buttons = c("csv","excel")),
              rownames = FALSE)
  })

  # GO-BP dotplot
  output$go_bp_plot <- renderPlotly({
    req(go_bp_data())
    df <- go_bp_data() |> head(20) |>
      mutate(GeneRatio_n = sapply(GeneRatio, function(r) eval(parse(text = r))))
    plot_ly(df, x = ~GeneRatio_n, y = ~reorder(Description, GeneRatio_n),
            color = ~p.adjust, type = "scatter", mode = "markers",
            marker = list(size = ~Count, sizemode = "area", sizeref = 0.05),
            text = ~paste0(Description, "<br>Count: ", Count)) |>
      layout(yaxis = list(title = ""), xaxis = list(title = "Gene Ratio"),
             title = "GO Biological Process")
  })

  # KEGG dotplot
  output$kegg_plot <- renderPlotly({
    req(kegg_data())
    df <- kegg_data() |> head(20) |>
      mutate(GeneRatio_n = sapply(GeneRatio, function(r) eval(parse(text = r))))
    plot_ly(df, x = ~GeneRatio_n, y = ~reorder(Description, GeneRatio_n),
            color = ~p.adjust, type = "scatter", mode = "markers",
            marker = list(size = ~Count, sizemode = "area", sizeref = 0.05),
            text = ~paste0(Description, "<br>Count: ", Count)) |>
      layout(yaxis = list(title = ""), xaxis = list(title = "Gene Ratio"),
             title = "KEGG Pathways")
  })

  output$go_table   <- renderDT({
    req(go_bp_data())
    datatable(go_bp_data() |> select(ID, Description, GeneRatio, p.adjust, Count),
              options = list(pageLength = 10, dom = "Bfrtip", buttons = c("csv","excel")),
              extensions = "Buttons", rownames = FALSE)
  })
  output$kegg_table <- renderDT({
    req(kegg_data())
    datatable(kegg_data() |> select(ID, Description, GeneRatio, p.adjust, Count),
              options = list(pageLength = 10, dom = "Bfrtip", buttons = c("csv","excel")),
              extensions = "Buttons", rownames = FALSE)
  })

  # Splicing
  output$splicing_plot <- renderPlotly({
    req(splicing_data())
    if (!"event_type" %in% colnames(splicing_data())) return(NULL)
    df <- splicing_data() |> count(event_type)
    plot_ly(df, x = ~event_type, y = ~n, type = "bar",
            color = ~event_type) |>
      layout(title = "Eventos de Splicing por Tipo",
             xaxis = list(title = "Tipo"), yaxis = list(title = "Contagem"))
  })
  output$splicing_table <- renderDT({
    req(splicing_data())
    datatable(splicing_data(),
              options = list(pageLength = 10, scrollX = TRUE, dom = "Bfrtip",
                             buttons = c("csv","excel")),
              extensions = "Buttons", rownames = FALSE)
  })

  # WGCNA
  output$wgcna_module_plot <- renderPlotly({
    req(wgcna_data())
    df <- wgcna_data() |> count(module)
    plot_ly(df, x = ~module, y = ~n, type = "bar",
            marker = list(color = df$module)) |>
      layout(title = "Genes por Módulo WGCNA",
             xaxis = list(title = "Módulo"), yaxis = list(title = "N Genes"))
  })
  output$hub_genes_table <- renderDT({
    req(hub_data())
    datatable(hub_data() |> select(gene_id, module, kWithin) |>
                mutate(kWithin = round(kWithin, 3)),
              options = list(pageLength = 10), rownames = FALSE)
  })

  # Integração
  output$n_candidates <- renderValueBox({
    n <- if (!is.null(cand_data())) nrow(cand_data()) else 0
    valueBox(n, "Candidatos Chave", icon = icon("star"), color = "yellow", width = 12)
  })
  output$integration_plot <- renderPlotly({
    req(ranking_data())
    df <- ranking_data() |> head(30)
    plot_ly(df, x = ~integration_score, y = ~reorder(gene_id, integration_score),
            color = ~regulation,
            colors = c("up" = "#E41A1C", "down" = "#377EB8"),
            type = "bar", orientation = "h",
            text = ~sprintf("LFC: %.2f | padj: %.2e | Camadas: %d",
                             log2FoldChange, padj, evidence_layers)) |>
      layout(yaxis = list(title = ""), xaxis = list(title = "Integration Score"),
             title = "Top 30 Genes – Integration Score")
  })
  output$candidates_table <- renderDT({
    req(cand_data())
    datatable(
      cand_data() |>
        select(gene_id, log2FoldChange, padj, integration_score,
               evidence_layers, has_splicing, in_pathway, is_hub) |>
        mutate(across(c(log2FoldChange, integration_score), \(x) round(x, 3)),
               padj = formatC(padj, format = "e", digits = 2)),
      extensions = "Buttons",
      options = list(pageLength = 20, dom = "Bfrtip", buttons = c("csv","excel")),
      rownames = FALSE
    )
  })

  # Dados brutos
  output$raw_deseq2_table <- renderDT({
    req(de_data())
    datatable(
      de_data() |> mutate(across(where(is.numeric), \(x) round(x, 4))),
      extensions = "Buttons",
      options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
                     buttons = c("csv","excel")),
      rownames = FALSE
    )
  })

}

shinyApp(ui, server)
