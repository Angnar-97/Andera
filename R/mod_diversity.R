# Módulo: Diversidad Alfa
#
# Estimaciones de diversidad alfa con phyloseq::plot_richness. Los índices
# basados en conteos (Chao1, ACE, Observed) requieren conteos enteros y se
# validan antes de ejecutar el cálculo. Shannon, Simpson, InvSimpson y Fisher
# son robustos a proporciones.

# Índices que dependen estrictamente de singletons / doubletons o conteos.
ALPHA_COUNT_INDICES <- c("Observed", "Chao1", "ACE")

mod_diversity_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Diversidad alfa", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Estimaciones de diversidad ", tags$em("dentro"),
        " de cada muestra (", tags$em("alpha diversity"),
        "): riqueza observada, Chao1, ACE, Shannon, Simpson, InvSimpson y Fisher, ",
        "calculadas con ", tags$code("phyloseq::plot_richness"), "."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
          bslib::card_body(
            checkboxGroupInput(
              ns("diversity"), "Índices de diversidad",
              choices  = c("Observed", "Chao1", "ACE",
                           "Shannon", "Simpson", "InvSimpson", "Fisher"),
              selected = c("Shannon", "Simpson"),
              inline   = FALSE
            ),
            tags$small(class = "andera-form-help",
              tags$strong("Observed / Chao1 / ACE"),
              " requieren conteos enteros (singletons, doubletons). ",
              tags$strong("Shannon / Simpson / InvSimpson / Fisher"),
              " funcionan también con proporciones."
            ),
            uiOutput(ns("alpha_variableui")),
            tags$div(class = "andera-actions",
              actionButton(ns("update_diversity"), "Actualizar",
                           class = "btn btn-primary",
                           icon  = icon("arrow-right"))
            )
          )
        ),

        bslib::card(
          bslib::card_header(bsicons::bs_icon("bar-chart-line"),
                             " Riqueza / diversidad"),
          bslib::card_body(
            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Boxplot"),
              shinycssloaders::withSpinner(
                plotOutput(ns("diversityPlot"), height = "520px"),
                type = 5
              )
            ),
            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow",
                        "Test estadístico · Welch t-test (2 grupos) o Kruskal-Wallis (≥3) · BH-FDR"),
              tableOutput(ns("diversityStats"))
            ),
            tags$div(class = "andera-actions",
              downloadButton(ns("download_diversity"), "Descargar (.png)",
                             class = "btn-outline-secondary")
            )
          )
        )
      )
    )
  )
}

mod_diversity_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    output$alpha_variableui <- renderUI({
      req(physeq())
      selectInput(
        session$ns("variable_alpha"),
        "Variable de metadata",
        choices = colnames(phyloseq::sample_data(physeq()))
      )
    })

    diversity_result <- reactiveVal(NULL)

    observe({
      physeq()
      diversity_result(NULL)
    })

    observeEvent(input$update_diversity, {
      req(physeq(), input$diversity, input$variable_alpha)
      ps <- physeq()

      # ---- Validaciones ----
      uses_count_indices <- any(input$diversity %in% ALPHA_COUNT_INDICES)
      if (uses_count_indices) {
        v <- validate_for_count_indices(ps)
        if (!isTRUE(v$ok)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = v$reason, type = "warning"
          )
          return()
        }
      }
      v <- validate_for_grouping(ps, input$variable_alpha,
                                  min_n_per_group = 3L)
      group_ok <- isTRUE(v$ok)
      # No se rechaza si la variable no cumple — el plot se hace igual,
      # pero no calculamos test estadístico.

      result <- tryCatch({
        rich_df <- phyloseq::estimate_richness(
          ps, measures = input$diversity
        )
        rich_df$.sample <- rownames(rich_df)
        sd <- as(phyloseq::sample_data(ps), "data.frame")
        rich_df$.group <- sd[[input$variable_alpha]][
          match(rich_df$.sample, rownames(sd))
        ]

        # Long-format para ggplot
        long_df <- stats::reshape(
          rich_df,
          varying    = setdiff(colnames(rich_df), c(".sample", ".group")),
          v.names    = "value",
          times      = setdiff(colnames(rich_df), c(".sample", ".group")),
          timevar    = "index",
          direction  = "long"
        )
        long_df$index <- factor(long_df$index, levels = input$diversity)

        # Test estadístico por índice si hay grupos válidos
        tests <- NULL
        if (group_ok) {
          n_levels <- length(unique(stats::na.omit(long_df$.group)))
          test_name <- if (n_levels == 2L) "Welch t-test" else "Kruskal-Wallis"
          tests <- do.call(rbind, lapply(input$diversity, function(idx) {
            sub <- long_df[long_df$index == idx & !is.na(long_df$value), ]
            if (length(unique(sub$.group)) < 2L) {
              return(data.frame(index = idx, p = NA_real_, test = test_name,
                                  stringsAsFactors = FALSE))
            }
            p <- tryCatch({
              if (n_levels == 2L) {
                stats::t.test(value ~ .group, data = sub, var.equal = FALSE)$p.value
              } else {
                stats::kruskal.test(value ~ .group, data = sub)$p.value
              }
            }, error = function(e) NA_real_)
            data.frame(index = idx, p = p, test = test_name,
                        stringsAsFactors = FALSE)
          }))
          tests$p_adj_BH <- stats::p.adjust(tests$p, method = "BH")
        }

        list(
          long_df    = long_df,
          tests      = tests,
          variable   = input$variable_alpha,
          n_levels   = if (group_ok)
            length(unique(stats::na.omit(long_df$.group)))
          else NA_integer_
        )
      }, error = function(e) {
        shinyalert::shinyalert(
          title = "Error al calcular la diversidad",
          text  = e$message, type = "error"
        )
        NULL
      })
      if (!is.null(result)) diversity_result(result)
    })

    diversity_plot <- reactive({
      r <- diversity_result()
      if (is.null(r)) return(NULL)
      df <- r$long_df

      p <- ggplot2::ggplot(
        df, ggplot2::aes(x = .data$.group, y = .data$value,
                          fill = .data$.group)
      ) +
        ggplot2::geom_boxplot(alpha = 0.65, outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.20, alpha = 0.7, size = 1.6,
                              ggplot2::aes(color = .data$.group)) +
        ggplot2::facet_wrap(~ index, scales = "free_y") +
        ggplot2::labs(x = r$variable, y = "Diversidad alfa",
                       fill = r$variable, color = r$variable) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1),
          legend.position  = "none",
          strip.text       = ggplot2::element_text(face = "bold",
                                                    color = "#3B2A4A")
        )

      # Anotar p-values cuando hay test
      if (!is.null(r$tests)) {
        ann <- r$tests
        ann$label <- ifelse(
          is.na(ann$p_adj_BH), "—",
          sprintf("%s\np_adj = %s", ann$test,
                   ifelse(ann$p_adj_BH < 0.001, "<0.001",
                           formatC(ann$p_adj_BH, digits = 3, format = "f")))
        )
        p <- p + ggplot2::geom_text(
          data = ann,
          ggplot2::aes(label = .data$label),
          x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
          size = 3.2, color = "#3B2A4A",
          inherit.aes = FALSE
        )
      }
      p
    })

    output$diversityPlot <- renderPlot({
      p <- diversity_plot()
      validate(need(p, "Pulsa 'Actualizar' para calcular la diversidad."))
      p
    })

    output$diversityStats <- renderTable({
      r <- diversity_result()
      validate(need(r, ""))
      validate(need(r$tests,
        "Configura una variable con ≥2 niveles y ≥3 muestras/grupo para ver tests estadísticos."))
      r$tests
    }, rownames = FALSE, digits = 4, na = "—")

    output$download_diversity <- download_plot(
      diversity_plot, "diversidad-alfa",
      width = 10, height = 6
    )
  })
}
