# Módulo: Ordenación
#
# Reducción de dimensionalidad de matrices de distancia entre muestras.
# Tres métodos disponibles, con diagnóstico específico para cada uno:
#
#   - PCoA  (Principal Coordinates Analysis): proyección métrica de la matriz
#           de distancias. Reporta el porcentaje de varianza explicada por
#           los primeros dos ejes.
#   - NMDS  (Non-metric MDS): minimiza el `stress` entre rangos de distancia
#           y posición proyectada. Robusto a no linealidad. Reporta `stress`
#           con código de calidad (verde < 0.10 · amarillo 0.10–0.20 ·
#           rojo > 0.20, Kruskal 1964).
#   - CAP   (Constrained Analysis of Principal Coordinates): PCoA constrained
#           a una covariable de diseño. Reporta R² constrained y un test de
#           permutación (anova.cca, vegan).
#
# Todos aceptan cualquier distancia del catálogo `ANDERA_DISTANCES`.

mod_ordination_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Ordenación", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Proyección 2D de las muestras a partir de una matriz de distancia ",
        "ecológica. ", tags$strong("PCoA"), " preserva métrica; ",
        tags$strong("NMDS"), " preserva rangos (más robusto a no linealidad); ",
        tags$strong("CAP"), " aísla la varianza atribuible a una covariable ",
        "de diseño."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        # ----- Parámetros -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
          bslib::card_body(
            radioButtons(
              ns("method"), "Método",
              choices = c(
                "PCoA"  = "PCoA",
                "NMDS"  = "NMDS",
                "CAP"   = "CAP"
              ),
              selected = "PCoA",
              inline   = TRUE
            ),
            selectInput(
              ns("distance"), "Distancia ecológica",
              choices  = ANDERA_DISTANCES,
              selected = "bray"
            ),
            tags$small(class = "andera-form-help",
              tags$strong("UniFrac"), " requiere árbol filogenético. ",
              tags$strong("Aitchison"),
              " es la opción composicionalmente correcta sin árbol."
            ),
            uiOutput(ns("color_variable_ui")),
            uiOutput(ns("constraint_variable_ui")),  # solo CAP
            tags$div(class = "andera-actions",
              actionButton(ns("update_ordination"), "Actualizar",
                           class = "btn btn-primary",
                           icon = icon("arrow-right"))
            )
          )
        ),

        # ----- Resultados -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("braces"), " Ordenación"),
          bslib::card_body(

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Plot"),
              shinycssloaders::withSpinner(
                plotOutput(ns("ordPlot"), height = "480px"), type = 5
              )
            ),

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Diagnóstico"),
              uiOutput(ns("diagnostic"))
            ),

            r_code_panel(ns),

            tags$div(class = "andera-actions",
              downloadButton(ns("download_ord"), "Descargar (.png)",
                             class = "btn-outline-secondary")
            )
          )
        )
      ),

      # ===== Comparador de distancias =====
      tags$br(),
      bslib::card(
        bslib::card_header(bsicons::bs_icon("layout-three-columns"),
                           " Comparador de distancias"),
        bslib::card_body(
          tags$p(class = "andera-form-help",
            "PCoA con distintas distancias en una sola vista (Bray–Curtis, ",
            "Jaccard, Aitchison y, si hay árbol, weighted UniFrac). Si la ",
            "separación entre grupos se mantiene en todas las distancias, la ",
            "conclusión es robusta. Si depende mucho de la métrica, conviene ",
            "reportar varias en el manuscrito."
          ),
          bslib::layout_columns(
            col_widths = c(4, 8),
            tagList(
              uiOutput(ns("compare_color_ui")),
              tags$div(class = "andera-actions",
                actionButton(ns("compare_distances"),
                             "Generar comparación",
                             class = "btn btn-primary btn-sm",
                             icon = icon("layer-group"))
              )
            ),
            tagList(
              shinycssloaders::withSpinner(
                plotOutput(ns("comparePlot"), height = "560px"), type = 5
              ),
              tags$div(class = "andera-actions",
                downloadButton(ns("download_compare"), "Descargar (.png)",
                               class = "btn-outline-secondary")
              )
            )
          )
        )
      )
    )
  )
}

mod_ordination_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    # ----- UI dinámica -----
    output$color_variable_ui <- renderUI({
      req(physeq())
      selectInput(
        session$ns("color_variable"),
        "Colorear por variable",
        choices = c("(ninguna)", colnames(phyloseq::sample_data(physeq())))
      )
    })

    output$constraint_variable_ui <- renderUI({
      req(physeq(), input$method)
      if (input$method != "CAP") return(NULL)
      vars <- colnames(phyloseq::sample_data(physeq()))
      selectInput(
        session$ns("constraint_variable"),
        "Variable de diseño (constraint)",
        choices = vars
      )
    })

    # ----- Estado -----
    ord_result <- reactiveVal(NULL)

    observe({
      physeq()
      ord_result(NULL)
    })

    # ----- Cómputo -----
    observeEvent(input$update_ordination, {
      req(physeq(), input$distance, input$method)
      ps          <- physeq()
      method      <- input$method
      dist_method <- input$distance

      if (distance_requires_tree(dist_method)) {
        v <- validate_for_unifrac(ps)
        if (!isTRUE(v$ok)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = v$reason, type = "warning"
          )
          return()
        }
      }

      if (method == "CAP") {
        cv <- input$constraint_variable
        if (is.null(cv) || !nzchar(cv)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = "CAP requiere una variable de diseño.",
            type  = "warning"
          )
          return()
        }
        v <- validate_for_grouping(ps, cv, min_n_per_group = 3L)
        if (!isTRUE(v$ok)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = v$reason, type = "warning"
          )
          return()
        }
      }

      result <- tryCatch(
        withProgress(
          message = sprintf("%s + %s", method, dist_method),
          value   = 0, {

            incProgress(0.30, detail = "Calculando matriz de distancias")
            d <- compute_distance(ps, dist_method)

            incProgress(0.50, detail = sprintf("Ordenación %s", method))
            ord <- if (method == "CAP") {
              # phyloseq::ordinate envuelve vegan::capscale para CAP.
              # Necesita el formula constraint.
              cv  <- input$constraint_variable
              frm <- stats::reformulate(cv)
              phyloseq::ordinate(
                ps, method = "CAP", distance = d, formula = frm
              )
            } else {
              phyloseq::ordinate(ps, method = method, distance = d)
            }

            incProgress(0.20, detail = "Empaquetando")
            color_var <- input$color_variable
            if (is.null(color_var) || color_var == "(ninguna)") {
              color_var <- NULL
            }

            list(
              ps                  = ps,
              ord                 = ord,
              method              = method,
              distance            = dist_method,
              color               = color_var,
              constraint_variable = if (method == "CAP")
                input$constraint_variable else NULL
            )
          }
        ),
        error = function(e) {
          shinyalert::shinyalert(
            title = sprintf("Error al calcular %s", method),
            text  = e$message,
            type  = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) ord_result(result)
    })

    # ----- Plot reactivo -----
    ord_plot <- reactive({
      r <- ord_result()
      if (is.null(r)) return(NULL)

      p <- phyloseq::plot_ordination(
        r$ps, r$ord, type = "samples", color = r$color
      ) +
        ggplot2::geom_point(size = 3, alpha = 0.85)

      # Stat ellipses si hay color discreto
      if (!is.null(r$color)) {
        sd <- as(phyloseq::sample_data(r$ps), "data.frame")
        if (is.character(sd[[r$color]]) || is.factor(sd[[r$color]])) {
          p <- p + ggplot2::stat_ellipse(type = "norm", linetype = "dashed",
                                          alpha = 0.4)
        }
      }

      p +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    })

    output$ordPlot <- renderPlot({
      p <- ord_plot()
      validate(need(p, "Pulsa 'Actualizar' para calcular la ordenación."))
      p
    })

    # ----- Diagnóstico (cambia según método) -----
    output$diagnostic <- renderUI({
      r <- ord_result()
      if (is.null(r)) {
        return(tags$p(class = "andera-muted",
          "Ejecuta la ordenación para ver el diagnóstico."))
      }

      if (r$method == "PCoA") {
        eig <- ord_eigenvalues(r$ord, "PCoA")
        if (is.null(eig)) {
          return(tags$p(class = "andera-muted", "Sin diagnóstico disponible."))
        }
        var_pct <- eig / sum(eig) * 100
        tags$dl(class = "andera-meta-dl",
          tags$dt("Método"),       tags$dd("PCoA"),
          tags$dt("Distancia"),    tags$dd(r$distance),
          tags$dt("Eje 1"),        tags$dd(sprintf("%.1f%% varianza", var_pct[1])),
          tags$dt("Eje 2"),        tags$dd(sprintf("%.1f%% varianza", var_pct[2])),
          tags$dt("Eje 1 + 2"),    tags$dd(sprintf("%.1f%% acumulado", sum(var_pct[1:2]))),
          tags$dt("Total ejes"),   tags$dd(length(eig))
        )

      } else if (r$method == "NMDS") {
        stress <- r$ord$stress
        cat <- nmds_stress_category(stress)
        alert <- tags$div(
          class = paste("andera-alert", cat$css),
          tags$strong(cat$icon, " ", cat$label),
          tags$p(cat$message)
        )
        tagList(
          tags$dl(class = "andera-meta-dl",
            tags$dt("Método"),     tags$dd("NMDS (Kruskal)"),
            tags$dt("Distancia"),  tags$dd(r$distance),
            tags$dt("Stress"),     tags$dd(sprintf("%.4f", stress)),
            tags$dt("Convergencia"), tags$dd(if (isTRUE(r$ord$converged)) "sí" else "no")
          ),
          alert
        )

      } else if (r$method == "CAP") {
        # CAP devuelve un objeto cca/capscale (envuelto por phyloseq).
        d_obj <- compute_distance(r$ps, r$distance)
        sd    <- as(phyloseq::sample_data(r$ps), "data.frame")
        cv    <- r$constraint_variable

        # Test de permutación con anova.cca, sólo si la formula tiene términos
        anova_res <- tryCatch(
          vegan::anova.cca(r$ord, permutations = 999),
          error = function(e) NULL
        )

        constrained_var <- r$ord$CCA$tot.chi
        total_var       <- r$ord$tot.chi
        r2 <- if (!is.null(constrained_var) && !is.null(total_var) && total_var > 0)
          constrained_var / total_var else NA_real_

        children <- list(
          tags$dl(class = "andera-meta-dl",
            tags$dt("Método"),               tags$dd("CAP"),
            tags$dt("Distancia"),            tags$dd(r$distance),
            tags$dt("Variable de diseño"),    tags$dd(cv),
            tags$dt("R² constrained"),       tags$dd(if (is.na(r2)) "—" else sprintf("%.3f (%.1f%%)", r2, 100 * r2))
          )
        )
        if (!is.null(anova_res)) {
          F_val <- anova_res$F[1]
          p_val <- anova_res$`Pr(>F)`[1]
          children <- c(children, list(
            tags$dl(class = "andera-meta-dl",
              tags$dt("F (anova.cca)"),  tags$dd(sprintf("%.3f", F_val)),
              tags$dt("p-valor"),       tags$dd(format_p(p_val)),
              tags$dt("Permutaciones"), tags$dd("999")
            )
          ))
        }
        do.call(tagList, children)
      }
    })

    output$download_ord <- download_plot(ord_plot, "ordinacion",
                                          width = 8, height = 6)

    # ----- Reproducir en R -----
    r_code_text <- reactive({
      r <- ord_result()
      if (is.null(r)) return(NULL)
      r_code_ordination(r)
    })
    r_code_handlers(output, r_code_text, "ordinacion")

    # =====================================================================
    # Comparador de distancias
    # =====================================================================
    output$compare_color_ui <- renderUI({
      req(physeq())
      selectInput(
        session$ns("compare_color"),
        "Colorear por variable",
        choices = c("(ninguna)", colnames(phyloseq::sample_data(physeq())))
      )
    })

    compare_plot_obj <- reactiveVal(NULL)
    observe({
      physeq()
      compare_plot_obj(NULL)
    })

    observeEvent(input$compare_distances, {
      req(physeq())
      ps <- physeq()

      distances <- c("bray", "jaccard")
      if (has_tree(ps)) distances <- c(distances, "wunifrac")
      distances <- c(distances, "aitchison")
      distances <- distances[seq_len(min(4L, length(distances)))]

      color_var  <- input$compare_color
      use_color  <- !is.null(color_var) && color_var != "(ninguna)" && nzchar(color_var)

      p <- tryCatch(
        withProgress(message = "Comparador de distancias", value = 0, {
          sd <- as(phyloseq::sample_data(ps), "data.frame")
          step <- 0.85 / length(distances)

          results <- lapply(seq_along(distances), function(i) {
            dm <- distances[i]
            incProgress(step,
                        detail = sprintf("PCoA · %s (%d/%d)",
                                         dm, i, length(distances)))
            d <- compute_distance(ps, dm)
            ord <- phyloseq::ordinate(ps, method = "PCoA", distance = d)
            coords <- ord$vectors[, 1:2, drop = FALSE]
            eig <- ord$values$Eigenvalues
            eig_pos <- eig[eig > 0]
            var_pct <- (eig[1:2] / sum(eig_pos)) * 100
            data.frame(
              sample   = rownames(coords),
              PC1      = coords[, 1],
              PC2      = coords[, 2],
              distance = sprintf("%s · (%.0f%% + %.0f%%)",
                                  dm, var_pct[1], var_pct[2]),
              stringsAsFactors = FALSE
            )
          })
          df <- do.call(rbind, results)
          # Mantener orden de las distancias
          df$distance <- factor(df$distance, levels = unique(df$distance))

          if (use_color) {
            df$.colorvar <- sd[[color_var]][match(df$sample, rownames(sd))]
          }

          incProgress(0.15, detail = "Renderizando")
          p <- if (use_color) {
            ggplot2::ggplot(df,
              ggplot2::aes(.data$PC1, .data$PC2, color = .data$.colorvar))
          } else {
            ggplot2::ggplot(df, ggplot2::aes(.data$PC1, .data$PC2))
          }
          p +
            ggplot2::geom_point(size = 2.5, alpha = 0.85) +
            ggplot2::facet_wrap(~ distance, scales = "free", ncol = 2) +
            ggplot2::labs(color = if (use_color) color_var else NULL) +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(
              panel.grid.minor = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(face = "bold",
                                                  color = "#5D263E"),
              strip.background = ggplot2::element_blank()
            )
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en comparador de distancias",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(p)) compare_plot_obj(p)
    })

    output$comparePlot <- renderPlot({
      p <- compare_plot_obj()
      validate(need(p,
        "Pulsa 'Generar comparación' para superponer 2-4 distancias."))
      p
    })

    output$download_compare <- download_plot(
      compare_plot_obj, "comparador-distancias",
      width = 12, height = 9
    )
  })
}

# ===========================================================================
# Helpers locales
# ===========================================================================

# Extrae los autovalores positivos de una ordenación PCoA.
ord_eigenvalues <- function(ord, method) {
  if (method != "PCoA") return(NULL)
  eig <- tryCatch(ord$values$Eigenvalues,    error = function(e) NULL)
  if (is.null(eig)) eig <- tryCatch(ord$eig,  error = function(e) NULL)
  eig <- eig[!is.na(eig) & eig > 0]
  if (!length(eig)) return(NULL)
  eig
}

# Clasifica el stress NMDS según convención de Kruskal/Clarke.
nmds_stress_category <- function(stress) {
  if (is.na(stress)) {
    return(list(css = "andera-alert-warning",
                icon = "?",
                label = "Stress no disponible",
                message = "No se pudo extraer el valor de stress."))
  }
  if (stress < 0.10) {
    list(css = "andera-alert-success",
         icon = "✓",
         label = sprintf("Stress = %.3f · ajuste excelente", stress),
         message = "El ordenamiento de rangos refleja con alta fidelidad las distancias originales (Kruskal: stress < 0.10).")
  } else if (stress < 0.20) {
    list(css = "andera-alert-info",
         icon = "·",
         label = sprintf("Stress = %.3f · ajuste aceptable", stress),
         message = "Lectura general posible; cuidado con interpretaciones detalladas (Kruskal: 0.10–0.20).")
  } else {
    list(css = "andera-alert-warning",
         icon = "⚠",
         label = sprintf("Stress = %.3f · ajuste deficiente", stress),
         message = "Stress > 0.20: la representación 2D distorsiona apreciablemente las distancias. Considera más dimensiones (k = 3), otra distancia, o interpretar con cautela.")
  }
}
