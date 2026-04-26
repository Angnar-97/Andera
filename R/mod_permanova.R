# Módulo: PERMANOVA
#
# Test multivariante de varianza por permutaciones (vegan::adonis2) sobre la
# matriz de distancia elegida, con dos extensiones canónicas:
#
#   1. Diagnóstico de homogeneidad de dispersiones (vegan::betadisper +
#      permutest), siguiendo Anderson 2001 — sin él, un PERMANOVA significativo
#      podría reflejar varianzas heterogéneas y no centroides distintos.
#
#   2. Comparaciones pareadas post-hoc cuando la variable de agrupación
#      tiene >2 niveles. Para cada par (g_i, g_j) se subsetea la matriz de
#      distancias y se vuelve a correr adonis2; los p-valores se ajustan
#      con BH-FDR (Benjamini-Hochberg).
#
# Salida: tabla adonis2 + tabla permutest(betadisper) + alerta interpretativa
# + tabla pairwise (cuando aplica) + bloque de reproducibilidad + script R.

mod_permanova_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("PERMANOVA", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Test multivariante por permutaciones con ",
        tags$code("vegan::adonis2"),
        " + diagnóstico de dispersiones con ", tags$code("betadisper"),
        ". Cuando la variable tiene >2 niveles, se añaden comparaciones ",
        "pareadas post-hoc con corrección BH-FDR."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        # ----- Parámetros -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
          bslib::card_body(
            selectInput(
              ns("distance_permanova"), "Distancia ecológica",
              choices  = ANDERA_DISTANCES,
              selected = "bray"
            ),
            tags$small(class = "andera-form-help",
              tags$strong("UniFrac"), " requiere árbol filogenético. ",
              tags$strong("Aitchison"),
              " es la opción composicionalmente correcta cuando no hay árbol."
            ),
            uiOutput(ns("grouping_variable_ui")),
            numericInput(
              ns("permutations"), "Número de permutaciones",
              value = 999, min = 99, step = 100
            ),
            tags$small(class = "andera-form-help",
              "Mínimo recomendado: 999. Para p ≤ 0.001 fiable, 9999."
            ),
            checkboxInput(
              ns("do_pairwise"),
              "Pairwise PERMANOVA si hay >2 niveles (BH-FDR)",
              value = TRUE
            ),
            tags$div(class = "andera-actions",
              actionButton(ns("update_permanova"), "Ejecutar PERMANOVA",
                           class = "btn btn-primary",
                           icon = icon("play"))
            )
          )
        ),

        # ----- Resultados -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("table"), " Resultados"),
          bslib::card_body(

            # PERMANOVA global
            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow",
                        "PERMANOVA · diferencia de centroides"),
              shinycssloaders::withSpinner(
                tableOutput(ns("permanovaResults")), type = 5
              )
            ),

            # betadisper
            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow",
                        "Diagnóstico de dispersión · betadisper"),
              tableOutput(ns("dispersionResults")),
              uiOutput(ns("dispersionAlert"))
            ),

            # Pairwise (condicional)
            uiOutput(ns("pairwiseSection")),

            # Reproducibilidad
            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Reproducibilidad"),
              uiOutput(ns("metaInfo"))
            ),

            r_code_panel(ns),

            tags$div(class = "andera-actions",
              downloadButton(ns("download_permanova"),
                             "Descargar todo (.csv)",
                             class = "btn-outline-secondary")
            )
          )
        )
      )
    )
  )
}

mod_permanova_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    # Selector dinámico de variable
    output$grouping_variable_ui <- renderUI({
      req(physeq())
      selectInput(
        session$ns("grouping_variable"),
        "Variable de agrupación",
        choices = colnames(phyloseq::sample_data(physeq()))
      )
    })

    # Resultado consolidado: list(adonis, disp, disp_test, pairwise, meta)
    permanova_result <- reactiveVal(NULL)

    # Reset al cambiar el phyloseq activo
    observe({
      physeq()
      permanova_result(NULL)
    })

    observeEvent(input$update_permanova, {
      req(physeq(), input$distance_permanova, input$grouping_variable,
          input$permutations)
      ps          <- physeq()
      dist_method <- input$distance_permanova
      var         <- input$grouping_variable
      n_perm      <- input$permutations
      do_pw       <- isTRUE(input$do_pairwise)

      # ---- Validación ----
      checks <- list(
        validate_for_grouping(ps, var, min_n_per_group = 3L)
      )
      if (distance_requires_tree(dist_method)) {
        checks <- c(checks, list(validate_for_unifrac(ps)))
      }
      v <- first_invalid(checks)
      if (!isTRUE(v$ok)) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = v$reason,
          type  = "warning"
        )
        return()
      }

      # ---- Cómputo ----
      result <- tryCatch(
        withProgress(message = "PERMANOVA + betadisper", value = 0, {
          incProgress(0.15, detail = "Calculando matriz de distancias")
          d <- compute_distance(ps, dist_method)

          incProgress(0.05, detail = "Preparando modelo")
          sdat   <- as(phyloseq::sample_data(ps), "data.frame")
          frm    <- stats::reformulate(var, response = "d")
          groups <- factor(sdat[[var]])

          incProgress(0.30,
                      detail = sprintf("adonis2 (%s permutaciones)",
                                       format(n_perm, big.mark = ",")))
          adonis <- vegan::adonis2(frm, data = sdat, permutations = n_perm)

          incProgress(0.20, detail = "betadisper + permutest")
          # betadisper requiere un factor sin NA. Filtrar muestras válidas.
          keep   <- !is.na(groups)
          d_sub  <- if (all(keep)) d else stats::as.dist(as.matrix(d)[keep, keep])
          g_sub  <- droplevels(groups[keep])
          sdat_kept <- sdat[keep, , drop = FALSE]
          bd     <- vegan::betadisper(d_sub, g_sub)
          bd_t   <- vegan::permutest(bd, permutations = n_perm)

          # ---- Pairwise PERMANOVA (si aplica) ----
          pairwise <- NULL
          n_levels <- nlevels(g_sub)
          if (do_pw && n_levels > 2L) {
            n_pairs <- choose(n_levels, 2L)
            incProgress(0.20,
                        detail = sprintf("Pairwise PERMANOVA (%d pares · BH-FDR)",
                                         n_pairs))
            pairs <- utils::combn(levels(g_sub), 2L, simplify = FALSE)

            pw_rows <- lapply(pairs, function(pair) {
              keep_pair <- as.character(g_sub) %in% pair
              d_p <- stats::as.dist(as.matrix(d_sub)[keep_pair, keep_pair])
              sdat_p <- sdat_kept[keep_pair, , drop = FALSE]
              sdat_p[[var]] <- factor(sdat_p[[var]], levels = pair)
              a <- tryCatch(
                vegan::adonis2(
                  stats::reformulate(var, response = "d_p"),
                  data = sdat_p, permutations = n_perm
                ),
                error = function(e) NULL
              )
              if (is.null(a)) {
                data.frame(
                  group_a = pair[1], group_b = pair[2],
                  n_a = sum(sdat_p[[var]] == pair[1]),
                  n_b = sum(sdat_p[[var]] == pair[2]),
                  R2 = NA_real_, F = NA_real_, p = NA_real_,
                  stringsAsFactors = FALSE
                )
              } else {
                data.frame(
                  group_a = pair[1], group_b = pair[2],
                  n_a = sum(sdat_p[[var]] == pair[1]),
                  n_b = sum(sdat_p[[var]] == pair[2]),
                  R2  = a$R2[1],
                  F   = a$F[1],
                  p   = a$`Pr(>F)`[1],
                  stringsAsFactors = FALSE
                )
              }
            })
            pairwise <- do.call(rbind, pw_rows)
            pairwise$p_adj_BH <- stats::p.adjust(pairwise$p, method = "BH")
            pairwise <- pairwise[order(pairwise$p_adj_BH, na.last = TRUE), ,
                                  drop = FALSE]
            rownames(pairwise) <- NULL
          }

          incProgress(0.10, detail = "Empaquetando")
          list(
            adonis    = adonis,
            disp      = bd,
            disp_test = bd_t,
            pairwise  = pairwise,
            meta      = list(
              distance      = dist_method,
              variable      = var,
              formula       = paste0("d ~ ", var),
              permutations  = n_perm,
              n_per_group   = as.list(table(g_sub)),
              n_levels      = n_levels,
              R2            = adonis$R2[1],
              adonis_F      = adonis$F[1],
              adonis_p      = adonis$`Pr(>F)`[1],
              disp_F        = bd_t$tab[1, "F"],
              disp_p        = bd_t$tab[1, "Pr(>F)"],
              has_pairwise  = !is.null(pairwise),
              n_sig_pw      = if (!is.null(pairwise))
                sum(pairwise$p_adj_BH < 0.05, na.rm = TRUE) else NA_integer_,
              vegan_v       = as.character(utils::packageVersion("vegan"))
            )
          )
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error al ejecutar PERMANOVA",
            text  = e$message,
            type  = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) permanova_result(result)
    })

    # ---- Render: tabla adonis2 ----
    output$permanovaResults <- renderTable({
      res <- permanova_result()
      validate(need(res, "Pulsa 'Ejecutar PERMANOVA' para obtener resultados."))
      as.data.frame(res$adonis)
    }, rownames = TRUE, digits = 4, na = "—")

    # ---- Render: tabla betadisper ----
    output$dispersionResults <- renderTable({
      res <- permanova_result()
      validate(need(res, ""))
      as.data.frame(res$disp_test$tab)
    }, rownames = TRUE, digits = 4, na = "—")

    # ---- Render: alerta interpretativa ----
    output$dispersionAlert <- renderUI({
      res <- permanova_result()
      if (is.null(res)) return(NULL)

      ap <- res$meta$adonis_p
      dp <- res$meta$disp_p
      sig <- function(x) !is.na(x) && x < 0.05

      if (!sig(ap)) {
        tags$div(class = "andera-alert andera-alert-info",
          tags$strong("PERMANOVA no significativo."),
          " Sin evidencia de diferencias en composición entre grupos al nivel α = 0.05."
        )
      } else if (sig(ap) && sig(dp)) {
        tags$div(class = "andera-alert andera-alert-warning",
          tags$strong("⚠ Atención · dispersión heterogénea."),
          tags$p(
            "Tanto PERMANOVA (p = ", format_p(ap), ") como ",
            "betadisper (p = ", format_p(dp), ") son significativos. ",
            "El efecto observado puede deberse a diferencias de centroide ",
            tags$em("o"), " a diferencias de dispersión multivariante entre ",
            "grupos. ", tags$strong("No"), " se puede atribuir el resultado ",
            "únicamente a las medias grupales (Anderson 2001). Considera ",
            "reportar ambos resultados o usar un test alternativo (e.g. ",
            "ANOSIM, MRPP)."
          )
        )
      } else {
        tags$div(class = "andera-alert andera-alert-success",
          tags$strong("✓ Dispersiones homogéneas."),
          tags$p(
            "PERMANOVA significativo (p = ", format_p(ap), ") y betadisper ",
            "no significativo (p = ", format_p(dp), "): el efecto es ",
            "atribuible a diferencias de centroide entre grupos."
          )
        )
      }
    })

    # ---- Render: sección pairwise (condicional) ----
    output$pairwiseSection <- renderUI({
      res <- permanova_result()
      if (is.null(res)) return(NULL)
      if (is.null(res$pairwise)) {
        # No hay pairwise — muestra hint cuando hay solo 2 niveles
        if (!is.null(res$meta$n_levels) && res$meta$n_levels == 2L) {
          return(tags$div(class = "andera-result-section",
            tags$span(class = "andera-eyebrow",
                      "Comparaciones pareadas"),
            tags$p(class = "andera-muted",
              "Solo 2 niveles: el PERMANOVA global ya es la única comparación posible."
            )
          ))
        }
        return(NULL)
      }
      tags$div(class = "andera-result-section",
        tags$span(class = "andera-eyebrow",
                  sprintf("Comparaciones pareadas · post-hoc (%d significativas)",
                           res$meta$n_sig_pw)),
        tableOutput(session$ns("pairwiseResults")),
        tags$small(class = "andera-form-help",
          "Pairwise PERMANOVA con ajuste BH-FDR (Benjamini-Hochberg). ",
          tags$code("p_adj_BH"), " es el p-valor corregido por múltiples ",
          "comparaciones; usa este valor para decidir significancia."
        )
      )
    })

    output$pairwiseResults <- renderTable({
      res <- permanova_result()
      validate(need(res, ""))
      validate(need(res$pairwise, ""))
      res$pairwise
    }, rownames = FALSE, digits = 4, na = "—")

    # ---- Render: bloque de reproducibilidad ----
    output$metaInfo <- renderUI({
      res <- permanova_result()
      if (is.null(res)) {
        return(tags$p(class = "andera-muted",
          "Ejecuta PERMANOVA para ver los parámetros usados."))
      }
      m <- res$meta
      tags$dl(class = "andera-meta-dl",
        tags$dt("Distancia"),
        tags$dd(m$distance),

        tags$dt("Fórmula"),
        tags$dd(tags$code(m$formula)),

        tags$dt("Permutaciones"),
        tags$dd(format(m$permutations, big.mark = ",")),

        tags$dt("R²"),
        tags$dd(sprintf("%.4f", m$R2)),

        tags$dt("p-valor (adonis2)"),
        tags$dd(format_p(m$adonis_p)),

        tags$dt("p-valor (betadisper)"),
        tags$dd(format_p(m$disp_p)),

        tags$dt("n por grupo"),
        tags$dd(paste(
          names(m$n_per_group),
          unlist(m$n_per_group),
          sep = " = ", collapse = ", "
        )),

        if (isTRUE(m$has_pairwise)) tagList(
          tags$dt("Pares pareados"),
          tags$dd(sprintf("%d pares · %d significativos (p_adj < 0.05)",
                           choose(m$n_levels, 2L), m$n_sig_pw))
        ),

        tags$dt("vegan"),
        tags$dd(paste0("v", m$vegan_v))
      )
    })

    # ---- Descarga: CSV con metadatos + todas las tablas ----
    output$download_permanova <- shiny::downloadHandler(
      filename = function() paste0("permanova-", Sys.Date(), ".csv"),
      content  = function(file) {
        res <- permanova_result()
        shiny::req(res)
        m <- res$meta

        con <- file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)

        meta_lines <- c(
          "# PERMANOVA results - Andera",
          sprintf("# Date: %s", Sys.Date()),
          sprintf("# Distance: %s", m$distance),
          sprintf("# Formula: %s", m$formula),
          sprintf("# Permutations: %s", format(m$permutations, big.mark = ",")),
          sprintf("# n per group: %s", paste(
            names(m$n_per_group), unlist(m$n_per_group),
            sep = " = ", collapse = ", "
          )),
          sprintf("# R-squared: %.6f", m$R2),
          sprintf("# adonis2 p-value:    %.6f", m$adonis_p),
          sprintf("# betadisper p-value: %.6f", m$disp_p),
          sprintf("# Pairwise PERMANOVA: %s",
                   if (isTRUE(m$has_pairwise))
                     sprintf("%d pairs (BH-FDR)", choose(m$n_levels, 2L))
                   else "no (only 2 levels)"),
          sprintf("# vegan version: %s", m$vegan_v),
          ""
        )
        writeLines(meta_lines, con)

        writeLines("# adonis2 (PERMANOVA)", con)
        utils::write.table(
          as.data.frame(res$adonis),
          file = con, sep = ",",
          row.names = TRUE, col.names = NA, quote = TRUE
        )
        writeLines("", con)

        writeLines("# permutest(betadisper) - homogeneity of dispersions", con)
        utils::write.table(
          as.data.frame(res$disp_test$tab),
          file = con, sep = ",",
          row.names = TRUE, col.names = NA, quote = TRUE
        )

        if (!is.null(res$pairwise)) {
          writeLines("", con)
          writeLines("# Pairwise PERMANOVA (BH-FDR)", con)
          utils::write.table(
            res$pairwise, file = con, sep = ",",
            row.names = FALSE, col.names = TRUE, quote = TRUE
          )
        }
      }
    )

    # ----- Reproducir en R -----
    r_code_text <- reactive({
      r <- permanova_result()
      if (is.null(r)) return(NULL)
      r_code_permanova(r$meta)
    })
    r_code_handlers(output, r_code_text, "permanova")
  })
}

# ---------------------------------------------------------------------------
# Helper local: formateo "humano" de p-valores
# ---------------------------------------------------------------------------
format_p <- function(p) {
  if (is.na(p))  return("—")
  if (p < 0.001) return("< 0.001")
  formatC(p, digits = 3, format = "f")
}
