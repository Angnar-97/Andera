# Módulo: Filtrado
# Filtra el objeto phyloseq por una variable/valor de metadata con ps_filter.
# Expone un reactive con el phyloseq filtrado (o el original si aún no se filtró).

mod_filtering_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Filtrado", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Restringe el análisis a un subconjunto de muestras mediante ",
        tags$code("microViz::ps_filter"),
        ". El phyloseq filtrado se propaga a las pestañas de análisis."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        bslib::card(
          bslib::card_header(bsicons::bs_icon("funnel"), " Criterios"),
          bslib::card_body(
            numericInput(ns("min_count"), "Número mínimo de conteos",
                         value = 0, min = 0),
            uiOutput(ns("variableui")),
            uiOutput(ns("valueui")),
            tags$div(class = "andera-actions",
              actionButton(ns("filter"), "Aplicar filtro",
                           class = "btn btn-primary",
                           icon = icon("filter")),
              downloadButton(ns("download_filtered"), "Descargar (.rds)",
                             class = "btn-outline-secondary")
            )
          )
        ),

        tagList(
          bslib::card(
            bslib::card_header(bsicons::bs_icon("card-list"), " Resumen del phyloseq filtrado"),
            bslib::card_body(
              shinycssloaders::withSpinner(verbatimTextOutput(ns("filtered_physeq_summary")),
                                            type = 5)
            )
          ),
          bslib::card(
            bslib::card_header(bsicons::bs_icon("diagram-3"), " Árbol filogenético"),
            bslib::card_body(
              shinycssloaders::withSpinner(plotOutput(ns("filtered_phylo_tree"), height = "520px"),
                                            type = 5),
              tags$div(class = "andera-actions",
                downloadButton(ns("download_tree"), "Descargar árbol (.png)",
                               class = "btn-outline-secondary")
              )
            )
          )
        )
      )
    )
  )
}

mod_filtering_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {
    physeq_filtered <- reactiveVal(NULL)

    observe({
      ps <- physeq()
      if (!is.null(ps)) physeq_filtered(ps)
    })

    output$variableui <- renderUI({
      req(physeq_filtered())
      selectInput(
        session$ns("variable"),
        "Variable de metadata",
        choices = colnames(phyloseq::sample_data(physeq_filtered()))
      )
    })

    output$valueui <- renderUI({
      req(input$variable, physeq_filtered())
      selectizeInput(
        session$ns("value"),
        "Valores a conservar",
        choices  = unique(phyloseq::sample_data(physeq_filtered())[, input$variable]),
        multiple = TRUE
      )
    })

    observeEvent(input$filter, {
      req(input$variable, input$value)

      if (!(input$variable %in% colnames(phyloseq::sample_data(physeq_filtered())))) {
        shinyalert::shinyalert(
          title = "Error",
          text  = paste("La variable seleccionada", input$variable,
                        "no existe en los datos de muestra"),
          type  = "error"
        )
        return()
      }

      # Filtrado seguro vía phyloseq::prune_samples — sin eval/parse
      filtered <- tryCatch({
        ps_in <- physeq_filtered()
        sd    <- as(phyloseq::sample_data(ps_in), "data.frame")
        keep  <- as.character(sd[[input$variable]]) %in%
                 as.character(input$value)
        if (!any(keep)) NULL else phyloseq::prune_samples(keep, ps_in)
      }, error = function(e) NULL)

      if (is.null(filtered)) {
        shinyalert::shinyalert(
          title = "Error",
          text  = "El filtro no pudo ser aplicado. Verifica variable y valores.",
          type  = "error"
        )
        return()
      }

      physeq_filtered(filtered)

      showNotification(
        sprintf("Filtro aplicado: %d taxa, %d muestras",
                phyloseq::ntaxa(filtered), phyloseq::nsamples(filtered)),
        type     = "message",
        duration = 4
      )
    })

    tree_plot <- reactive({
      ps <- physeq_filtered()
      if (is.null(ps) || !has_tree(ps)) return(NULL)
      color_col <- if ("SampleType" %in% colnames(phyloseq::sample_data(ps))) "SampleType" else NULL
      phyloseq::plot_tree(ps, color = color_col)
    })

    output$filtered_physeq_summary <- renderPrint({
      ps <- physeq_filtered()
      validate(need(ps, "Carga primero un objeto phyloseq desde 'Carga de Datos'."))
      print(ps)
    })

    output$filtered_phylo_tree <- renderPlot({
      ps <- physeq_filtered()
      validate(need(ps, "Carga primero un objeto phyloseq desde 'Carga de Datos'."))
      validate(need(has_tree(ps),
                    "El objeto phyloseq cargado carece de árbol filogenético."))
      tree_plot()
    })

    output$download_filtered <- download_rds(physeq_filtered, "phyloseq-filtered")
    output$download_tree     <- download_plot(tree_plot, "phylo-tree-filtered",
                                              width = 8, height = 8)

    physeq_filtered
  })
}
