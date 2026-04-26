# Módulo: Vistas avanzadas (paquete microbiome)
#
# Tres visualizaciones del paquete microbiome (Lahti & Shetty) que
# complementan las del paquete phyloseq base:
#
#   - Core microbiome — heatmap prevalencia × abundancia (microbiome::plot_core).
#     Útil para identificar el "core" estable de la comunidad y elegir
#     thresholds de filtrado de forma justificada.
#
#   - Landscape — proyección 2D (PCoA o NMDS) con densidad de muestras
#     superpuesta como contornos (microbiome::plot_landscape). Más informativo
#     que un PCoA pelado cuando hay muchas muestras o solapamiento de grupos.
#
#   - Prevalencia por taxón — barras horizontales de prevalencia,
#     organizadas por un rango taxonómico superior (microbiome::plot_taxa_prevalence).
#     Diagnóstico rápido sobre qué taxones son ubicuos vs raros.
#
# Cada vista tiene su propio panel dentro de un navset_card_tab para no
# forzar al usuario a recorrer 3 pestañas en el navbar principal.

mod_microbiome_views_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Vistas avanzadas (microbiome)", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Tres visualizaciones canónicas del paquete ",
        tags$a(href = "https://microbiome.github.io",
                target = "_blank", rel = "noopener", "microbiome"),
        " que complementan las gráficas de phyloseq base. Selecciona la ",
        "pestaña que quieras explorar; cada vista tiene parámetros propios."
      ),

      bslib::navset_card_tab(
        id = ns("view_tabs"),

        # =================================================================
        # Core microbiome
        # =================================================================
        bslib::nav_panel(
          "Core microbiome",
          bslib::layout_columns(
            col_widths = c(4, 8),
            bslib::card(
              bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
              bslib::card_body(
                tags$p(class = "andera-form-help",
                  "Heatmap prevalencia × abundancia: ¿qué fracción de tus ",
                  "muestras tienen un taxón con al menos X% de abundancia? ",
                  "Identifica el ", tags$em("core"), " estable de la comunidad."
                ),
                uiOutput(ns("core_rank_ui")),
                sliderInput(ns("core_prev_min"),
                            "Prevalencia mínima a mostrar",
                            min = 0.05, max = 1, value = 0.50, step = 0.05),
                tags$div(class = "andera-actions",
                  actionButton(ns("update_core"), "Calcular core",
                               class = "btn btn-primary",
                               icon = icon("play"))
                )
              )
            ),
            bslib::card(
              bslib::card_header(bsicons::bs_icon("grid-3x3"),
                                 " Heatmap del core"),
              bslib::card_body(
                shinycssloaders::withSpinner(
                  plotOutput(ns("corePlot"), height = "480px"), type = 5
                ),
                tags$div(class = "andera-actions",
                  downloadButton(ns("download_core"), "Descargar (.png)",
                                 class = "btn-outline-secondary")
                )
              )
            )
          )
        ),

        # =================================================================
        # Landscape
        # =================================================================
        bslib::nav_panel(
          "Landscape",
          bslib::layout_columns(
            col_widths = c(4, 8),
            bslib::card(
              bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
              bslib::card_body(
                tags$p(class = "andera-form-help",
                  "Ordenación 2D con contornos de densidad — útil para ver ",
                  "regiones del espacio donde se acumulan muestras."
                ),
                radioButtons(ns("landscape_method"), "Método",
                              choices = c("PCoA" = "PCoA", "NMDS" = "NMDS"),
                              selected = "PCoA", inline = TRUE),
                selectInput(ns("landscape_distance"), "Distancia",
                             choices = ANDERA_DISTANCES, selected = "bray"),
                uiOutput(ns("landscape_color_ui")),
                tags$div(class = "andera-actions",
                  actionButton(ns("update_landscape"), "Generar landscape",
                               class = "btn btn-primary",
                               icon = icon("play"))
                )
              )
            ),
            bslib::card(
              bslib::card_header(bsicons::bs_icon("graph-up-arrow"),
                                 " Landscape"),
              bslib::card_body(
                shinycssloaders::withSpinner(
                  plotOutput(ns("landscapePlot"), height = "480px"), type = 5
                ),
                tags$div(class = "andera-actions",
                  downloadButton(ns("download_landscape"), "Descargar (.png)",
                                 class = "btn-outline-secondary")
                )
              )
            )
          )
        ),

        # =================================================================
        # Taxa prevalence
        # =================================================================
        bslib::nav_panel(
          "Prevalencia por taxón",
          bslib::layout_columns(
            col_widths = c(4, 8),
            bslib::card(
              bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
              bslib::card_body(
                tags$p(class = "andera-form-help",
                  "Prevalencia (fracción de muestras donde aparece) por taxón, ",
                  "facetado por un rango taxonómico superior. Complementa al ",
                  "core: aquí se ve el detalle por taxón individual."
                ),
                uiOutput(ns("prev_level_ui")),
                numericInput(ns("prev_detection"),
                             "Detección mínima (abundancia relativa)",
                             value = 0.001, min = 0, max = 1, step = 0.001),
                tags$small(class = "andera-form-help",
                  "Un taxón cuenta como 'presente' en una muestra si su ",
                  "abundancia relativa supera este umbral."
                ),
                tags$div(class = "andera-actions",
                  actionButton(ns("update_prevalence"), "Calcular prevalencia",
                               class = "btn btn-primary",
                               icon = icon("play"))
                )
              )
            ),
            bslib::card(
              bslib::card_header(bsicons::bs_icon("bar-chart-line"),
                                 " Prevalencia por taxón"),
              bslib::card_body(
                shinycssloaders::withSpinner(
                  plotOutput(ns("prevalencePlot"), height = "560px"), type = 5
                ),
                tags$div(class = "andera-actions",
                  downloadButton(ns("download_prevalence"),
                                 "Descargar (.png)",
                                 class = "btn-outline-secondary")
                )
              )
            )
          )
        )
      )
    )
  )
}

mod_microbiome_views_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    # =====================================================================
    # Core microbiome
    # =====================================================================

    output$core_rank_ui <- renderUI({
      req(physeq())
      tt <- tryCatch(phyloseq::tax_table(physeq(), errorIfNULL = FALSE),
                     error = function(e) NULL)
      ranks <- if (!is.null(tt)) colnames(tt) else character(0)
      choices <- c("(sin glom)" = "__none__", ranks)
      default <- if ("Genus" %in% ranks) "Genus" else "__none__"
      selectInput(session$ns("core_rank"),
                   "Glomerar al rango",
                   choices = choices, selected = default)
    })

    core_plot_obj <- reactiveVal(NULL)
    observe({ physeq(); core_plot_obj(NULL) })

    observeEvent(input$update_core, {
      req(physeq(), input$core_prev_min)
      ps <- physeq()

      v <- validate_for_count_indices(ps)
      if (!isTRUE(v$ok)) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = "plot_core requiere conteos enteros para calcular prevalencia y detección.",
          type  = "warning"
        )
        return()
      }

      p <- tryCatch(
        withProgress(message = "Core microbiome", value = 0, {
          incProgress(0.30, detail = "Pre-procesamiento")
          ps_proc <- ps
          if (input$core_rank != "__none__") {
            ps_proc <- phyloseq::tax_glom(ps_proc, taxrank = input$core_rank,
                                          NArm = TRUE)
          }
          ps_rel <- microbiome::transform(ps_proc, "compositional")

          incProgress(0.50, detail = "plot_core")
          # microbiome::plot_core construye un heatmap; le damos
          # gradientes razonables.
          prev <- seq(input$core_prev_min, 1, length.out = 10)
          det  <- 10^seq(log10(0.001), log10(0.2), length.out = 10)
          microbiome::plot_core(
            ps_rel,
            prevalences      = prev,
            detections       = det,
            plot.type        = "heatmap",
            colours          = grDevices::colorRampPalette(
              c("#F1ECDF", "#C4A962", "#C32D28", "#5D263E"))(50)
          ) + ggplot2::theme_minimal(base_size = 11)
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en core microbiome",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(p)) core_plot_obj(p)
    })

    output$corePlot <- renderPlot({
      p <- core_plot_obj()
      validate(need(p, "Pulsa 'Calcular core' para ver el heatmap."))
      p
    })

    output$download_core <- download_plot(core_plot_obj, "core-microbiome",
                                           width = 10, height = 7)

    # =====================================================================
    # Landscape
    # =====================================================================

    output$landscape_color_ui <- renderUI({
      req(physeq())
      vars <- colnames(phyloseq::sample_data(physeq()))
      selectInput(session$ns("landscape_color"),
                   "Colorear por variable",
                   choices = c("(ninguna)", vars))
    })

    landscape_plot_obj <- reactiveVal(NULL)
    observe({ physeq(); landscape_plot_obj(NULL) })

    observeEvent(input$update_landscape, {
      req(physeq(), input$landscape_method, input$landscape_distance)
      ps <- physeq()
      method   <- input$landscape_method
      dist_m   <- input$landscape_distance
      color_v  <- input$landscape_color
      use_color <- !is.null(color_v) && color_v != "(ninguna)" && nzchar(color_v)

      if (distance_requires_tree(dist_m)) {
        v <- validate_for_unifrac(ps)
        if (!isTRUE(v$ok)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = v$reason, type = "warning"
          )
          return()
        }
      }

      p <- tryCatch(
        withProgress(message = "Landscape", value = 0, {
          incProgress(0.40, detail = sprintf("%s + %s", method, dist_m))
          # microbiome::plot_landscape acepta un phyloseq o una matriz
          # de coordenadas. Para mantener control sobre la distancia,
          # ordenamos primero y luego construimos el plot landscape sobre
          # las coordenadas resultantes.
          d <- compute_distance(ps, dist_m)
          ord <- phyloseq::ordinate(ps, method = method, distance = d)

          incProgress(0.40, detail = "Generando contornos de densidad")
          coords <- if (method == "PCoA") {
            ord$vectors[, 1:2, drop = FALSE]
          } else {
            ord$points[, 1:2, drop = FALSE]
          }
          coords_df <- data.frame(
            x = coords[, 1], y = coords[, 2],
            sample = rownames(coords)
          )
          if (use_color) {
            sd <- as(phyloseq::sample_data(ps), "data.frame")
            coords_df[[color_v]] <- sd[[color_v]][match(coords_df$sample,
                                                          rownames(sd))]
          }

          incProgress(0.20, detail = "ggplot2")
          p <- ggplot2::ggplot(coords_df, ggplot2::aes(x = .data$x,
                                                          y = .data$y)) +
            ggplot2::stat_density_2d(
              ggplot2::aes(fill = ggplot2::after_stat(.data$level)),
              geom = "polygon", alpha = 0.30, color = NA
            ) +
            ggplot2::scale_fill_gradient(low = "#E5DDC9", high = "#5D263E",
                                          guide = "none")
          if (use_color) {
            p <- p + ggplot2::geom_point(
              ggplot2::aes(color = .data[[color_v]]), size = 3, alpha = 0.85
            )
          } else {
            p <- p + ggplot2::geom_point(color = "#4A6D5E",
                                          size = 3, alpha = 0.85)
          }
          p +
            ggplot2::labs(
              x = sprintf("%s 1", method),
              y = sprintf("%s 2", method),
              color = if (use_color) color_v else NULL
            ) +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en landscape",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(p)) landscape_plot_obj(p)
    })

    output$landscapePlot <- renderPlot({
      p <- landscape_plot_obj()
      validate(need(p, "Pulsa 'Generar landscape' para ver el plot."))
      p
    })

    output$download_landscape <- download_plot(
      landscape_plot_obj, "landscape", width = 9, height = 7
    )

    # =====================================================================
    # Taxa prevalence
    # =====================================================================

    output$prev_level_ui <- renderUI({
      req(physeq())
      tt <- tryCatch(phyloseq::tax_table(physeq(), errorIfNULL = FALSE),
                     error = function(e) NULL)
      ranks <- if (!is.null(tt)) colnames(tt) else character(0)
      if (!length(ranks)) {
        return(tags$p(class = "andera-muted",
          "El phyloseq no tiene tax_table; esta vista no es aplicable."))
      }
      default <- if ("Phylum" %in% ranks) "Phylum" else ranks[1]
      selectInput(session$ns("prev_level"),
                   "Facetar por rango",
                   choices = ranks, selected = default)
    })

    prevalence_plot_obj <- reactiveVal(NULL)
    observe({ physeq(); prevalence_plot_obj(NULL) })

    observeEvent(input$update_prevalence, {
      req(physeq(), input$prev_level, input$prev_detection)
      ps <- physeq()

      v <- validate_for_count_indices(ps)
      if (!isTRUE(v$ok)) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = "Esta vista requiere conteos enteros para calcular prevalencia.",
          type  = "warning"
        )
        return()
      }

      p <- tryCatch(
        withProgress(message = "Prevalencia por taxón", value = 0, {
          incProgress(0.40, detail = "Transformando a abundancias relativas")
          ps_rel <- microbiome::transform(ps, "compositional")

          incProgress(0.50,
                      detail = sprintf("plot_taxa_prevalence (%s)",
                                       input$prev_level))
          microbiome::plot_taxa_prevalence(
            ps_rel,
            level     = input$prev_level,
            detection = input$prev_detection
          ) +
            ggplot2::scale_color_manual(values = andera_palette) +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::theme(legend.position = "right",
                            strip.text       = ggplot2::element_text(
                              face = "bold", color = "#5D263E"))
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en prevalencia por taxón",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(p)) prevalence_plot_obj(p)
    })

    output$prevalencePlot <- renderPlot({
      p <- prevalence_plot_obj()
      validate(need(p,
        "Pulsa 'Calcular prevalencia' para ver la distribución por taxón."))
      p
    })

    output$download_prevalence <- download_plot(
      prevalence_plot_obj, "prevalencia-por-taxon",
      width = 10, height = 8
    )
  })
}
