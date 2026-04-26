# Módulo: Red de co-ocurrencia microbiana
#
# Inferencia de redes taxon-taxon a partir de correlaciones de Spearman sobre
# datos transformados a CLR (composicionalmente conscientes), con ajuste BH
# del p-valor. Es el approach recomendado en la formación cuando SpiecEasi /
# SparCC no están disponibles, y produce redes interpretables siempre que se
# combine con un threshold de magnitud `|r|` y un threshold de FDR.
#
# Diferencia respecto a `mod_graphs` (red de similaridad sample-sample):
#   - `mod_graphs`     → nodos = muestras, aristas = distancia entre composiciones
#   - este módulo      → nodos = taxones, aristas = co-variación a través de
#                         muestras (qué taxones aparecen juntos)
#
# Implementación:
#   1. Glomerar a un rango taxonómico (default Genus) → tax_glom
#   2. Filtro por prevalencia mínima (default 0.20)
#   3. CLR vía microbiome::transform("clr")
#   4. Matriz de correlación de Spearman + p-valores ajustados BH
#      (psych::corr.test maneja la matriz completa)
#   5. Construir igraph con aristas que cumplen |r| ≥ r_threshold AND
#      p_adj ≤ p_threshold
#   6. Visualizar con ggraph (layout_with_fr, edge color by sign of r,
#      node size by degree)

mod_microbial_network_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Red de co-ocurrencia microbiana", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Grafo donde los nodos son taxones y las aristas representan ",
        "correlación significativa de sus abundancias a través de las ",
        "muestras. Calculado con ", tags$strong("Spearman"),
        " sobre datos transformados a ", tags$strong("CLR"),
        " (composicionalmente correcto), con corrección ",
        tags$strong("BH-FDR"), " del p-valor."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        # ----- Parámetros -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
          bslib::card_body(

            tags$span(class = "andera-eyebrow", "Pre-procesamiento"),
            uiOutput(ns("rank_ui")),
            numericInput(
              ns("min_prevalence"),
              "Prevalencia mínima",
              value = 0.20, min = 0, max = 1, step = 0.05
            ),
            tags$small(class = "andera-form-help",
              "Las redes con miles de taxones son ininterpretables y lentas; ",
              "0.20 = taxón presente en ≥ 20% de las muestras."
            ),

            tags$hr(),
            tags$span(class = "andera-eyebrow", "Aristas"),
            sliderInput(
              ns("r_threshold"),
              HTML("|r| mínimo"),
              min = 0.30, max = 0.95, value = 0.50, step = 0.05
            ),
            numericInput(
              ns("padj_threshold"),
              "p-ajustado (BH-FDR) máximo",
              value = 0.05, min = 0, max = 1, step = 0.01
            ),

            tags$hr(),
            tags$span(class = "andera-eyebrow", "Visualización"),
            uiOutput(ns("color_variable_ui")),
            checkboxInput(ns("show_labels"), "Etiquetar nodos",
                          value = FALSE),

            tags$div(class = "andera-actions",
              actionButton(ns("update_network"), "Construir red",
                           class = "btn btn-primary",
                           icon = icon("play"))
            )
          )
        ),

        # ----- Resultados -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("diagram-3"), " Red"),
          bslib::card_body(

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Resumen"),
              uiOutput(ns("summary"))
            ),

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Grafo"),
              shinycssloaders::withSpinner(
                plotOutput(ns("netPlot"), height = "560px"), type = 5
              )
            ),

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Hubs · top 10 por grado"),
              tableOutput(ns("hubsTable"))
            ),

            tags$div(class = "andera-actions",
              downloadButton(ns("download_net"), "Descargar grafo (.png)",
                             class = "btn-outline-secondary"),
              downloadButton(ns("download_edges"), "Descargar aristas (.csv)",
                             class = "btn-outline-secondary")
            )
          )
        )
      )
    )
  )
}

mod_microbial_network_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    # ----- UI dinámica -----
    output$rank_ui <- renderUI({
      req(physeq())
      tt <- tryCatch(phyloseq::tax_table(physeq(), errorIfNULL = FALSE),
                     error = function(e) NULL)
      ranks <- if (!is.null(tt)) colnames(tt) else character(0)
      choices <- c("(sin glom)" = "__none__", ranks)
      default <- if ("Genus" %in% ranks) "Genus" else "__none__"
      selectInput(
        session$ns("glom_rank"),
        "Glomerar al rango",
        choices  = choices,
        selected = default
      )
    })

    output$color_variable_ui <- renderUI({
      req(physeq())
      tt <- tryCatch(phyloseq::tax_table(physeq(), errorIfNULL = FALSE),
                     error = function(e) NULL)
      ranks <- if (!is.null(tt)) colnames(tt) else character(0)
      selectInput(
        session$ns("color_variable"),
        "Colorear nodos por rango",
        choices  = c("(ninguno)", ranks),
        selected = if ("Phylum" %in% ranks) "Phylum" else "(ninguno)"
      )
    })

    # ----- Estado -----
    net_result <- reactiveVal(NULL)

    observe({
      physeq()
      net_result(NULL)
    })

    # ----- Cómputo -----
    observeEvent(input$update_network, {
      req(physeq(), input$r_threshold, input$padj_threshold)
      ps <- physeq()

      v <- validate_for_count_indices(ps)
      if (!isTRUE(v$ok)) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = "La red de co-ocurrencia requiere conteos enteros (CLR no es bien definido sobre proporciones imputadas con ceros estructurales).",
          type  = "warning"
        )
        return()
      }

      result <- tryCatch(
        withProgress(message = "Red de co-ocurrencia", value = 0, {

          incProgress(0.10, detail = "Pre-procesamiento")
          ps_proc <- ps
          if (input$glom_rank != "__none__") {
            ps_proc <- phyloseq::tax_glom(ps_proc, taxrank = input$glom_rank,
                                          NArm = TRUE)
          }
          if (input$min_prevalence > 0) {
            otu <- as(phyloseq::otu_table(ps_proc), "matrix")
            if (!phyloseq::taxa_are_rows(ps_proc)) otu <- t(otu)
            prev <- rowMeans(otu > 0)
            ps_proc <- phyloseq::prune_taxa(prev >= input$min_prevalence, ps_proc)
          }
          n_taxa <- phyloseq::ntaxa(ps_proc)
          if (n_taxa < 5) {
            stop(sprintf(
              "Tras el pre-procesamiento solo quedan %d taxones; relaja la prevalencia o cambia el rango.",
              n_taxa
            ))
          }

          incProgress(0.20, detail = sprintf("CLR sobre %d taxones", n_taxa))
          ps_clr <- microbiome::transform(ps_proc, "clr")
          otu_clr <- as(phyloseq::otu_table(ps_clr), "matrix")
          if (phyloseq::taxa_are_rows(ps_clr)) otu_clr <- t(otu_clr)
          # Tras el t(): filas = muestras, columnas = taxa.

          incProgress(0.40, detail = "Correlaciones de Spearman + FDR")
          ct <- psych::corr.test(
            otu_clr,
            method  = "spearman",
            adjust  = "none",   # ajustamos a mano para no depender de la
                                  # convención (upper/lower triangle) de psych
            ci      = FALSE
          )
          r_mat <- ct$r
          p_raw <- ct$p
          n <- nrow(r_mat)
          tax_ids <- rownames(r_mat)

          # BH-FDR sobre los pares únicos (upper triangle, excluyendo diagonal)
          ut <- upper.tri(p_raw)
          p_vec <- p_raw[ut]
          p_adj_vec <- stats::p.adjust(p_vec, method = "BH")
          p_adj <- matrix(NA_real_, n, n, dimnames = dimnames(p_raw))
          p_adj[ut] <- p_adj_vec
          p_adj <- pmin(p_adj, t(p_adj), na.rm = TRUE)  # simetrizar

          incProgress(0.20, detail = "Construyendo grafo")
          # Vectorizado: identificar pares (i, j) en upper triangle con
          # |r| ≥ thr AND p_adj ≤ thr. Sin loop.
          mask <- ut &
                  !is.na(r_mat) & !is.na(p_adj) &
                  abs(r_mat) >= input$r_threshold &
                  p_adj      <= input$padj_threshold
          idx <- which(mask, arr.ind = TRUE)

          edges <- if (nrow(idx)) {
            data.frame(
              from  = tax_ids[idx[, "row"]],
              to    = tax_ids[idx[, "col"]],
              r     = r_mat[mask],
              p_adj = p_adj[mask],
              stringsAsFactors = FALSE
            )
          } else {
            data.frame(from = character(0), to = character(0),
                        r = numeric(0), p_adj = numeric(0),
                        stringsAsFactors = FALSE)
          }

          if (nrow(edges) == 0L) {
            stop(sprintf(
              "Ninguna pareja taxon-taxon supera los thresholds (|r|=%.2f, p_adj=%.3f). Relaja los criterios.",
              input$r_threshold, input$padj_threshold
            ))
          }

          # Solo nodos con al menos una arista
          used_taxa <- unique(c(edges$from, edges$to))
          tt <- tryCatch(phyloseq::tax_table(ps_proc, errorIfNULL = FALSE),
                          error = function(e) NULL)
          tax_df <- if (!is.null(tt)) {
            data.frame(name = rownames(tt), as.data.frame(tt@.Data),
                        stringsAsFactors = FALSE, check.names = FALSE)
          } else {
            data.frame(name = used_taxa, stringsAsFactors = FALSE)
          }
          tax_df <- tax_df[tax_df$name %in% used_taxa, , drop = FALSE]

          incProgress(0.10, detail = "Layout y métricas")
          g <- igraph::graph_from_data_frame(
            d = edges, vertices = tax_df, directed = FALSE
          )
          deg <- igraph::degree(g)

          list(
            graph    = g,
            edges    = edges,
            n_taxa   = n_taxa,
            n_nodes  = igraph::vcount(g),
            n_edges  = igraph::ecount(g),
            degrees  = deg,
            meta = list(
              glom_rank      = input$glom_rank,
              min_prevalence = input$min_prevalence,
              r_threshold    = input$r_threshold,
              padj_threshold = input$padj_threshold,
              n_pos          = sum(edges$r > 0),
              n_neg          = sum(edges$r < 0)
            )
          )
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error al construir la red",
            text  = e$message,
            type  = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) net_result(result)
    })

    # ----- Render: resumen -----
    output$summary <- renderUI({
      r <- net_result()
      if (is.null(r)) {
        return(tags$p(class = "andera-muted",
          "Configura los parámetros y pulsa 'Construir red'."))
      }
      m <- r$meta
      tags$dl(class = "andera-meta-dl",
        tags$dt("Pre-procesamiento"),
        tags$dd(sprintf(
          "glom = %s · prevalencia ≥ %.0f%% · %d taxones",
          if (m$glom_rank == "__none__") "ninguno" else m$glom_rank,
          100 * m$min_prevalence, r$n_taxa
        )),
        tags$dt("Thresholds"),
        tags$dd(sprintf("|r| ≥ %.2f · p_adj ≤ %.3f",
                         m$r_threshold, m$padj_threshold)),
        tags$dt("Nodos en la red"),
        tags$dd(sprintf("%d (de %d taxones procesados)", r$n_nodes, r$n_taxa)),
        tags$dt("Aristas"),
        tags$dd(sprintf("%d (positivas %d · negativas %d)",
                         r$n_edges, m$n_pos, m$n_neg))
      )
    })

    # ----- Render: grafo -----
    net_plot <- reactive({
      r <- net_result()
      if (is.null(r)) return(NULL)
      build_network_plot(r$graph, r$edges,
                         color_variable = input$color_variable,
                         show_labels    = isTRUE(input$show_labels))
    })

    output$netPlot <- renderPlot({
      p <- net_plot()
      validate(need(p,
        "Pulsa 'Construir red' para ver la red de co-ocurrencia."))
      p
    })

    # ----- Render: hubs -----
    output$hubsTable <- renderTable({
      r <- net_result()
      validate(need(r, ""))
      deg <- sort(r$degrees, decreasing = TRUE)
      top <- head(deg, 10L)
      tax <- igraph::vertex_attr(r$graph)
      df <- data.frame(
        taxon  = names(top),
        degree = unname(top),
        stringsAsFactors = FALSE
      )
      # Anotar con tax_table si existe
      if (!is.null(tax) && length(tax) > 1L) {
        for (col in setdiff(names(tax), "name")) {
          df[[col]] <- tax[[col]][match(df$taxon, tax$name)]
        }
      }
      df
    }, rownames = FALSE, digits = 0)

    # ----- Descargas -----
    output$download_net <- download_plot(net_plot, "red-microbiana",
                                          width = 10, height = 8)

    output$download_edges <- shiny::downloadHandler(
      filename = function() paste0("red-microbiana-edges-", Sys.Date(), ".csv"),
      content = function(file) {
        r <- net_result()
        shiny::req(r)
        utils::write.csv(r$edges, file, row.names = FALSE)
      }
    )
  })
}

# ===========================================================================
# Helpers fuera del módulo
# ===========================================================================

# Construye el grafo con ggraph. Color de aristas: rojo correlación negativa,
# verde correlación positiva. Tamaño de nodo: degree.
build_network_plot <- function(g, edges, color_variable = "(ninguno)",
                                show_labels = FALSE) {
  use_color <- !is.null(color_variable) && color_variable != "(ninguno)" &&
                nzchar(color_variable)

  # Layout de Fruchterman-Reingold con seed fijo para reproducibilidad
  set.seed(42)
  tg <- tidygraph::as_tbl_graph(g) |>
    tidygraph::mutate(degree = tidygraph::centrality_degree())

  p <- ggraph::ggraph(tg, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(
        edge_alpha = abs(.data$r),
        edge_color = .data$r > 0,
        edge_width = abs(.data$r)
      ),
      show.legend = c(edge_alpha = FALSE, edge_color = TRUE, edge_width = FALSE)
    ) +
    ggraph::scale_edge_color_manual(
      values = c(`TRUE` = "#5B7D8C", `FALSE` = "#C32D28"),
      labels = c(`TRUE` = "positiva (+)", `FALSE` = "negativa (−)"),
      name   = "Correlación"
    ) +
    ggraph::scale_edge_width(range = c(0.3, 1.4)) +
    ggraph::geom_node_point(
      ggplot2::aes(
        size  = .data$degree,
        color = if (use_color) .data[[color_variable]] else I("#4A6D5E")
      )
    ) +
    ggplot2::scale_size_continuous(range = c(2, 7), name = "Grado") +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(legend.position = "right")

  if (use_color) {
    p <- p + ggplot2::labs(color = color_variable)
  }

  if (show_labels) {
    p <- p + ggraph::geom_node_text(
      ggplot2::aes(label = .data$name),
      size = 2.5, repel = TRUE, max.overlaps = 30
    )
  }
  p
}
