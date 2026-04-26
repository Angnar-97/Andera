# Módulo: Mapas de Calor
# plot_heatmap() con pre-filtrado de top-N taxa (un heatmap con miles de
# filas es ilegible y NMDS suele no converger) y selector de método de
# ordenación. phyloseq aplica internamente escala log al color.

mod_heatmaps_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Mapas de calor", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Heatmap de abundancias con ", tags$code("phyloseq::plot_heatmap"),
        ". Las filas se reducen a los ", tags$strong("N taxa más abundantes"),
        " (con miles de filas el plot es ilegible y la ordenación interna ",
        "puede no converger). El color usa escala log internamente."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
          bslib::card_body(
            numericInput(
              ns("top_n"), "Top N taxa por abundancia total",
              value = 50, min = 10, max = 500, step = 10
            ),
            tags$small(class = "andera-form-help",
              "Reducir a 50 taxa hace el heatmap legible y el NMDS estable."
            ),
            radioButtons(
              ns("method"), "Método de ordenación interno",
              choices  = c("PCoA (robusto)" = "PCoA",
                            "NMDS"          = "NMDS"),
              selected = "PCoA",
              inline   = TRUE
            ),
            selectInput(
              ns("distance"), "Distancia para la ordenación",
              choices  = c("bray", "jaccard", "euclidean"),
              selected = "bray"
            ),
            uiOutput(ns("label_rank_ui")),
            tags$small(class = "andera-form-help",
              tags$strong("OTU/ASV"), " etiqueta con identificadores. ",
              tags$strong("Genus / Species"),
              " requieren tax_table con ese rango."
            ),
            tags$div(class = "andera-actions",
              actionButton(ns("update_heatmaps"), "Actualizar",
                           class = "btn btn-primary",
                           icon = icon("arrow-right"))
            )
          )
        ),

        bslib::card(
          bslib::card_header(bsicons::bs_icon("grid-3x3"), " Heatmap"),
          bslib::card_body(
            shinycssloaders::withSpinner(
              plotOutput(ns("heatmap_plot"), height = "640px"), type = 5
            ),
            tags$div(class = "andera-actions",
              downloadButton(ns("download_heatmap"), "Descargar (.png)",
                             class = "btn-outline-secondary")
            )
          )
        )
      )
    )
  )
}

mod_heatmaps_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    output$label_rank_ui <- renderUI({
      req(physeq())
      tt <- tryCatch(phyloseq::tax_table(physeq(), errorIfNULL = FALSE),
                     error = function(e) NULL)
      ranks <- if (!is.null(tt)) colnames(tt) else character(0)
      choices <- c("OTU/ASV" = "OTU", ranks)
      selectInput(session$ns("label_rank"),
                   "Etiquetar taxa por rango",
                   choices = choices, selected = "OTU")
    })

    heatmap_state <- reactiveVal(NULL)

    observe({
      physeq()
      heatmap_state(NULL)
    })

    observeEvent(input$update_heatmaps, {
      req(physeq(), input$top_n, input$method, input$distance,
          input$label_rank)
      ps <- physeq()

      # Validar rango si no es OTU
      if (input$label_rank != "OTU") {
        tt <- tryCatch(phyloseq::tax_table(ps, errorIfNULL = FALSE),
                       error = function(e) NULL)
        if (is.null(tt) || !input$label_rank %in% colnames(tt)) {
          shinyalert::shinyalert(
            title = "Rango no disponible",
            text  = sprintf("El phyloseq no tiene el rango '%s' en tax_table.",
                             input$label_rank),
            type  = "warning"
          )
          return()
        }
      }

      # Pre-filtrar a top N taxa por suma total (taxa_sums)
      sums <- phyloseq::taxa_sums(ps)
      n_keep <- min(as.integer(input$top_n), length(sums))
      top_taxa <- names(sort(sums, decreasing = TRUE))[seq_len(n_keep)]
      ps_top <- phyloseq::prune_taxa(top_taxa, ps)

      heatmap_state(list(
        ps         = ps_top,
        label_rank = input$label_rank,
        method     = input$method,
        distance   = input$distance,
        n_kept     = n_keep,
        n_total    = length(sums)
      ))
    })

    heatmap_plot <- reactive({
      st <- heatmap_state()
      if (is.null(st)) return(NULL)
      # phyloseq::plot_heatmap interpreta `taxa.label` como nombre de columna
      # de tax_table; para etiquetar por OTU/ASV hay que pasar NULL.
      tax_lbl <- if (st$label_rank == "OTU") NULL else st$label_rank
      tryCatch(
        phyloseq::plot_heatmap(
          st$ps,
          method     = st$method,
          distance   = st$distance,
          taxa.label = tax_lbl
        ),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error al generar heatmap",
            text  = paste(
              "Falló la ordenación interna. Prueba con PCoA, otra distance",
              "o un Top N menor. Detalle:", e$message
            ),
            type  = "error"
          )
          NULL
        }
      )
    })

    output$heatmap_plot <- renderPlot({
      p <- heatmap_plot()
      validate(need(p, "Configura los parámetros y pulsa 'Actualizar'."))
      p
    })

    output$download_heatmap <- download_plot(heatmap_plot, "heatmap",
                                              width = 12, height = 9)
  })
}
