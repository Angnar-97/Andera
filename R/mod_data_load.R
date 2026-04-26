# Módulo: Carga de Datos
# Permite cargar un objeto phyloseq desde archivo o desde los datasets de
# ejemplo del paquete phyloseq. Expone un reactive con el phyloseq "cargado"
# (es decir, aceptado por el usuario al pulsar 'Actualizar').

mod_data_load_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Carga de Datos", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Sube un objeto phyloseq en formato ", tags$code(".rds"),
        " o elige uno de los datasets de ejemplo del paquete phyloseq."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Fuente"),
          bslib::card_body(
            radioButtons(
              ns("data_source"), "Origen del phyloseq",
              choices  = c("Archivo", "Conjunto de datos de ejemplo"),
              selected = "Conjunto de datos de ejemplo"
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Archivo'", ns("data_source")),
              fileInput(ns("file1"), "Selecciona un archivo .rds")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Conjunto de datos de ejemplo'", ns("data_source")),
              selectInput(
                ns("dataset"), "Dataset de ejemplo",
                choices  = c("Global Patterns", "Enterotype", "Soil"),
                selected = "Global Patterns"
              )
            ),
            tags$div(class = "andera-actions",
              actionButton(ns("update"), "Actualizar",
                           class = "btn btn-primary",
                           icon = icon("arrow-right")),
              downloadButton(ns("downloadData"), "Descargar (.rds)",
                             class = "btn-outline-secondary")
            )
          )
        ),

        tagList(
          bslib::card(
            bslib::card_header(bsicons::bs_icon("card-list"), " Resumen del phyloseq"),
            bslib::card_body(
              shinycssloaders::withSpinner(verbatimTextOutput(ns("physeq_summary")),
                                            type = 5)
            )
          ),
          bslib::card(
            bslib::card_header(bsicons::bs_icon("diagram-3"), " Árbol filogenético"),
            bslib::card_body(
              shinycssloaders::withSpinner(plotOutput(ns("phylo_tree"), height = "520px"),
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

mod_data_load_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    pending_ps <- reactive({
      if (input$data_source == "Archivo") {
        req(input$file1)
        # Validar que el RDS contiene un phyloseq
        ps_candidate <- tryCatch(
          readRDS(input$file1$datapath),
          error = function(e) {
            shinyalert::shinyalert(
              title = "RDS ilegible",
              text  = paste(
                "No se pudo leer el archivo. Verifica que es un .rds válido",
                "exportado con saveRDS()."
              ),
              type  = "error"
            )
            NULL
          }
        )
        if (is.null(ps_candidate)) return(NULL)
        if (!inherits(ps_candidate, "phyloseq")) {
          shinyalert::shinyalert(
            title = "Objeto no compatible",
            text  = sprintf(
              "El archivo contiene un objeto de clase '%s'; se esperaba un phyloseq. Exporta tu objeto con saveRDS(ps, 'archivo.rds') antes de subirlo.",
              paste(class(ps_candidate), collapse = "/")
            ),
            type  = "warning"
          )
          return(NULL)
        }
        ps_candidate
      } else {
        load_example_dataset(input$dataset)
      }
    })

    loaded_ps <- reactiveVal(NULL)
    observeEvent(input$update, loaded_ps(pending_ps()))

    tree_plot <- reactive({
      ps <- loaded_ps()
      if (is.null(ps) || !has_tree(ps)) return(NULL)
      color_col <- if ("SampleType" %in% colnames(phyloseq::sample_data(ps))) "SampleType" else NULL
      phyloseq::plot_tree(ps, color = color_col)
    })

    output$physeq_summary <- renderPrint({
      ps <- loaded_ps()
      validate(need(ps, "Pulsa 'Actualizar' para cargar el objeto phyloseq."))
      print(ps)
    })

    output$phylo_tree <- renderPlot({
      ps <- loaded_ps()
      validate(need(ps, "Pulsa 'Actualizar' para cargar el objeto phyloseq."))
      validate(need(has_tree(ps),
                    "El objeto phyloseq cargado carece de árbol filogenético."))
      tree_plot()
    })

    output$downloadData  <- download_rds(loaded_ps, "phyloseq")
    output$download_tree <- download_plot(tree_plot, "phylo-tree",
                                          width = 8, height = 8)

    loaded_ps
  })
}
