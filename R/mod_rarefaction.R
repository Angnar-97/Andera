# Módulo: Curvas de rarefacción
#
# Visualiza la riqueza observada (taxa con conteo > 0) como función de la
# profundidad de secuenciación, para evaluar si las muestras alcanzan
# saturación. Imprescindible antes de comparar diversidad alfa: si las curvas
# no se aplanan, las diferencias entre muestras pueden ser artefactos del
# esfuerzo de muestreo (depth bias).
#
# Implementación:
#   - vegan::rarecurve(otu, step, tidy = TRUE) → data.frame largo
#   - ggplot con una línea por muestra, coloreada opcionalmente por metadata
#   - Línea de referencia vertical (umbral de profundidad)
#   - Resumen de profundidades de secuenciación por muestra

mod_rarefaction_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Curvas de rarefacción", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Riqueza observada como función de la profundidad de secuenciación. ",
        "Si las curvas no alcanzan asíntota (saturación), las comparaciones ",
        "de diversidad alfa entre muestras pueden ser artefactos del esfuerzo ",
        "de muestreo. Pre-flight check obligatorio antes de Diversidad."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        # ----- Parámetros -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
          bslib::card_body(
            uiOutput(ns("color_variable_ui")),
            numericInput(
              ns("step_size"), "Paso (lecturas por punto)",
              value = 100, min = 10, step = 50
            ),
            tags$small(class = "andera-form-help",
              "Resolución del muestreo. Valores bajos producen curvas más ",
              "suaves pero el cómputo es más lento."
            ),
            uiOutput(ns("threshold_ui")),
            tags$div(class = "andera-actions",
              actionButton(ns("update_rarefaction"), "Calcular curvas",
                           class = "btn btn-primary",
                           icon = icon("play"))
            )
          )
        ),

        # ----- Resultados -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("graph-up"),
                             " Curvas de rarefacción"),
          bslib::card_body(
            shinycssloaders::withSpinner(
              plotOutput(ns("rarePlot"), height = "440px"), type = 5
            ),
            uiOutput(ns("depthSummary")),
            tags$div(class = "andera-actions",
              downloadButton(ns("download_rarefaction"), "Descargar (.png)",
                             class = "btn-outline-secondary")
            )
          )
        )
      )
    )
  )
}

mod_rarefaction_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    # ---- UI dinámica ----
    output$color_variable_ui <- renderUI({
      req(physeq())
      selectInput(
        session$ns("color_variable"),
        "Colorear por variable",
        choices  = c("(ninguna)", colnames(phyloseq::sample_data(physeq()))),
        selected = "(ninguna)"
      )
    })

    output$threshold_ui <- renderUI({
      req(physeq())
      depths <- sample_depths(physeq())
      max_d  <- max(depths)
      min_d  <- min(depths)
      sliderInput(
        session$ns("threshold"),
        "Umbral de referencia (línea vertical)",
        min   = 0,
        max   = max_d,
        value = min_d,
        step  = max(1, round(max_d / 100))
      )
    })

    rare_data <- reactiveVal(NULL)

    observe({
      physeq()
      rare_data(NULL)
    })

    observeEvent(input$update_rarefaction, {
      req(physeq(), input$step_size)
      ps <- physeq()

      # ---- Validación: rarefacción solo tiene sentido con conteos enteros
      v <- validate_for_count_indices(ps)
      if (!isTRUE(v$ok)) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = paste(
            "Las curvas de rarefacción requieren conteos enteros (conteos de",
            "lecturas, no proporciones)."
          ),
          type  = "warning"
        )
        return()
      }

      result <- tryCatch(
        withProgress(message = "Rarefacción", value = 0, {
          incProgress(0.20, detail = "Preparando matriz")
          otu <- as(phyloseq::otu_table(ps), "matrix")
          if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)
          # Tras este t(), filas = muestras, columnas = taxa.

          incProgress(0.55,
                      detail = sprintf("vegan::rarecurve (step = %d)",
                                       input$step_size))
          rc <- vegan::rarecurve(
            otu,
            step  = max(1L, as.integer(input$step_size)),
            label = FALSE,
            tidy  = TRUE
          )
          # Columnas devueltas: Site, Sample, Species

          incProgress(0.25, detail = "Empaquetando")
          list(
            df          = rc,
            depths      = rowSums(otu),
            sample_ids  = rownames(otu)
          )
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en rarefacción",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) rare_data(result)
    })

    # ---- Construcción del ggplot reactivo ----
    rare_plot <- reactive({
      r <- rare_data()
      if (is.null(r)) return(NULL)

      df <- r$df
      ps <- physeq()
      sd <- as(phyloseq::sample_data(ps), "data.frame")

      color_var <- input$color_variable
      use_color <- !is.null(color_var) && nzchar(color_var) &&
                    color_var != "(ninguna)" && color_var %in% colnames(sd)

      if (use_color) {
        df$.colorvar <- sd[[color_var]][match(df$Site, rownames(sd))]
      } else {
        df$.colorvar <- "Todas"
      }

      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x     = .data$Sample,
          y     = .data$Species,
          group = .data$Site,
          color = .data$.colorvar
        )
      ) +
        ggplot2::geom_line(alpha = 0.7, linewidth = 0.6) +
        ggplot2::labs(
          x     = "Lecturas (sub-muestreadas)",
          y     = "Riqueza observada (taxa)",
          color = if (use_color) color_var else NULL
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.position  = if (use_color) "right" else "none"
        )

      if (!is.null(input$threshold) && input$threshold > 0) {
        p <- p + ggplot2::geom_vline(
          xintercept = input$threshold,
          linetype   = "dashed",
          color      = "#5D263E",
          alpha      = 0.5
        ) + ggplot2::annotate(
          "text", x = input$threshold,
          y = max(df$Species), label = "umbral",
          vjust = -0.5, hjust = -0.1,
          family = "sans", size = 3, color = "#5D263E"
        )
      }
      p
    })

    output$rarePlot <- renderPlot({
      p <- rare_plot()
      validate(need(p, "Pulsa 'Calcular curvas' para ver la rarefacción."))
      p
    })

    # ---- Resumen de profundidades ----
    output$depthSummary <- renderUI({
      r <- rare_data()
      if (is.null(r)) return(NULL)
      depths    <- r$depths
      threshold <- input$threshold
      n_below   <- if (!is.null(threshold) && threshold > 0) {
        sum(depths < threshold)
      } else {
        NA_integer_
      }

      tags$div(class = "andera-result-section",
        tags$span(class = "andera-eyebrow",
                  "Profundidad de secuenciación (lecturas / muestra)"),
        tags$dl(class = "andera-meta-dl",
          tags$dt("Muestras"),
          tags$dd(format(length(depths), big.mark = ",")),
          tags$dt("Mínimo"),
          tags$dd(format(min(depths), big.mark = ",")),
          tags$dt("Mediana"),
          tags$dd(format(stats::median(depths), big.mark = ",")),
          tags$dt("Media"),
          tags$dd(format(round(mean(depths)), big.mark = ",")),
          tags$dt("Máximo"),
          tags$dd(format(max(depths), big.mark = ",")),
          if (!is.na(n_below)) tagList(
            tags$dt("Muestras bajo el umbral"),
            tags$dd(sprintf(
              "%d / %d (%.0f%%)",
              n_below, length(depths),
              100 * n_below / length(depths)
            ))
          )
        )
      )
    })

    output$download_rarefaction <- download_plot(
      rare_plot, "rarefaccion",
      width = 10, height = 6
    )
  })
}

# ---------------------------------------------------------------------------
# Helper: profundidades de secuenciación por muestra (suma de lecturas).
# ---------------------------------------------------------------------------
sample_depths <- function(ps) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) colSums(otu) else rowSums(otu)
}
