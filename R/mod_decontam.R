# Módulo: Decontaminación (decontam)
#
# QC pre-análisis para detectar taxones contaminantes — secuencias procedentes
# del kit de extracción, agua del laboratorio, polvo del aire — distinguibles
# de los taxones reales de las muestras gracias a controles negativos o a la
# concentración de DNA medida por muestra (Davis 2018, Microbiome 6:226).
#
# Dos métodos disponibles:
#   - "prevalence": usa controles negativos (columna logical en sample_data
#     que marca cada control como TRUE). Un taxón contaminante aparecerá con
#     mayor prevalencia en los controles que en las muestras reales.
#   - "frequency": usa concentración de DNA por muestra (columna numérica).
#     Un taxón contaminante muestra una correlación negativa abundancia-DNA
#     (más visible cuando hay menos DNA real).
#   - "combined": combinación de Fisher de ambos.
#
# Output: lista de taxones marcados como contaminantes con su `p` y
# decisión (TRUE/FALSE). El módulo NO altera el phyloseq activo: el usuario
# revisa la lista y filtra manualmente si lo desea.

mod_decontam_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Detección de contaminantes (decontam)",
              class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "QC para identificar taxones procedentes de contaminación de ",
        "laboratorio (kit de extracción, reactivos, ambiente). Usa el ",
        "método de Davis 2018 implementado en ",
        tags$a(href = "https://benjjneb.github.io/decontam/",
                target = "_blank", rel = "noopener", "decontam"),
        ". Requiere haber etiquetado los controles negativos o disponer de ",
        "concentraciones de DNA por muestra en el ", tags$code("sample_data"),
        "."
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
                "Prevalencia (controles negativos)" = "prevalence",
                "Frecuencia (concentración DNA)"     = "frequency",
                "Combinado"                          = "combined"
              ),
              selected = "prevalence",
              inline   = FALSE
            ),

            uiOutput(ns("control_var_ui")),  # solo prevalence / combined
            uiOutput(ns("dna_var_ui")),      # solo frequency / combined

            sliderInput(
              ns("threshold"),
              "Umbral de probabilidad (P)",
              min = 0.01, max = 0.50, value = 0.10, step = 0.01
            ),
            tags$small(class = "andera-form-help",
              "P es la probabilidad bajo la hipótesis nula de que el taxón ",
              "no es contaminante. Davis 2018 recomienda 0.10 como umbral ",
              "general; 0.50 es muy permisivo."
            ),

            tags$div(class = "andera-actions",
              actionButton(ns("update_decontam"), "Detectar contaminantes",
                           class = "btn btn-primary",
                           icon = icon("play"))
            )
          )
        ),

        # ----- Resultados -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("shield-check"),
                             " Resultados"),
          bslib::card_body(

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Resumen"),
              uiOutput(ns("summary"))
            ),

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow",
                        "Distribución de P bajo H₀"),
              shinycssloaders::withSpinner(
                plotOutput(ns("pPlot"), height = "300px"), type = 5
              )
            ),

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Taxones contaminantes"),
              shinycssloaders::withSpinner(
                DT::DTOutput(ns("contaminantTable")), type = 5
              )
            ),

            tags$div(class = "andera-actions",
              downloadButton(ns("download_contam"),
                             "Descargar lista (.csv)",
                             class = "btn-outline-secondary")
            )
          )
        )
      )
    )
  )
}

mod_decontam_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    # ----- UI dinámica -----
    output$control_var_ui <- renderUI({
      req(physeq(), input$method)
      if (input$method == "frequency") return(NULL)
      sd <- as(phyloseq::sample_data(physeq()), "data.frame")
      # Detectar columnas logical o factor con 2 niveles plausibles
      logical_cols <- names(sd)[vapply(sd, is.logical, logical(1))]
      tagList(
        selectInput(
          session$ns("control_var"),
          "Variable que marca los controles negativos (logical)",
          choices = c("(elige una)", logical_cols),
          selected = if (length(logical_cols)) logical_cols[1] else "(elige una)"
        ),
        if (!length(logical_cols)) tags$small(class = "andera-form-help",
          "Ninguna columna lógica encontrada. Necesitas una columna ",
          tags$code("logical"), " donde TRUE marque cada muestra como ",
          "control negativo."
        )
      )
    })

    output$dna_var_ui <- renderUI({
      req(physeq(), input$method)
      if (input$method == "prevalence") return(NULL)
      sd <- as(phyloseq::sample_data(physeq()), "data.frame")
      numeric_cols <- names(sd)[vapply(sd, is.numeric, logical(1))]
      tagList(
        selectInput(
          session$ns("dna_var"),
          "Variable de concentración de DNA (numérica)",
          choices = c("(elige una)", numeric_cols),
          selected = if (length(numeric_cols)) numeric_cols[1] else "(elige una)"
        ),
        if (!length(numeric_cols)) tags$small(class = "andera-form-help",
          "Ninguna columna numérica encontrada en sample_data. ",
          "Necesitas una columna numérica con la concentración medida en cada muestra."
        )
      )
    })

    # ----- Estado -----
    decontam_result <- reactiveVal(NULL)

    observe({
      physeq()
      decontam_result(NULL)
    })

    # ----- Cómputo -----
    observeEvent(input$update_decontam, {
      req(physeq(), input$method, input$threshold)
      ps <- physeq()

      v <- validate_for_count_indices(ps)
      if (!isTRUE(v$ok)) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = "decontam requiere conteos enteros. El phyloseq cargado contiene proporciones.",
          type  = "warning"
        )
        return()
      }

      # ---- Validación de la columna seleccionada según el método ----
      sd <- as(phyloseq::sample_data(ps), "data.frame")
      args <- list(
        seqtab    = ps,
        method    = input$method,
        threshold = input$threshold
      )

      if (input$method %in% c("prevalence", "combined")) {
        cv <- input$control_var
        if (is.null(cv) || cv == "(elige una)" || !nzchar(cv) ||
            !cv %in% names(sd)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = "Selecciona la variable que marca los controles negativos.",
            type  = "warning"
          )
          return()
        }
        if (!is.logical(sd[[cv]])) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = sprintf("La columna '%s' no es lógica (TRUE/FALSE).", cv),
            type  = "warning"
          )
          return()
        }
        if (!any(sd[[cv]], na.rm = TRUE)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = sprintf("La columna '%s' no contiene ningún TRUE: no hay controles marcados.", cv),
            type  = "warning"
          )
          return()
        }
        args$neg <- cv
      }

      if (input$method %in% c("frequency", "combined")) {
        dv <- input$dna_var
        if (is.null(dv) || dv == "(elige una)" || !nzchar(dv) ||
            !dv %in% names(sd)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = "Selecciona la variable de concentración de DNA.",
            type  = "warning"
          )
          return()
        }
        if (!is.numeric(sd[[dv]])) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = sprintf("La columna '%s' no es numérica.", dv),
            type  = "warning"
          )
          return()
        }
        args$conc <- dv
      }

      result <- tryCatch(
        withProgress(message = "decontam", value = 0, {
          incProgress(0.30, detail = "Calculando isContaminant")
          ic <- do.call(decontam::isContaminant, args)

          incProgress(0.30, detail = "Anotando taxonomía")
          tt <- tryCatch(phyloseq::tax_table(ps, errorIfNULL = FALSE),
                         error = function(e) NULL)
          if (!is.null(tt)) {
            tax_df <- as.data.frame(tt@.Data, stringsAsFactors = FALSE)
            tax_df$taxon <- rownames(tax_df)
            ic$taxon <- rownames(ic)
            res_df <- merge(ic, tax_df, by = "taxon", all.x = TRUE)
          } else {
            ic$taxon <- rownames(ic)
            res_df <- ic
          }
          # Reordenar: contaminantes primero, luego por p
          res_df <- res_df[order(-res_df$contaminant, res_df$p), ]

          incProgress(0.40, detail = "Empaquetando")
          list(
            df = res_df,
            meta = list(
              method      = input$method,
              threshold   = input$threshold,
              control_var = if (input$method != "frequency") input$control_var else NA_character_,
              dna_var     = if (input$method != "prevalence") input$dna_var     else NA_character_,
              n_total     = nrow(res_df),
              n_contam    = sum(res_df$contaminant, na.rm = TRUE),
              decontam_v  = as.character(utils::packageVersion("decontam"))
            )
          )
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en decontam",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) decontam_result(result)
    })

    # ----- Render: resumen -----
    output$summary <- renderUI({
      r <- decontam_result()
      if (is.null(r)) {
        return(tags$p(class = "andera-muted",
          "Configura los parámetros y pulsa 'Detectar contaminantes'."))
      }
      m <- r$meta
      tags$dl(class = "andera-meta-dl",
        tags$dt("Método"),
        tags$dd(m$method),
        if (!is.na(m$control_var)) tagList(
          tags$dt("Variable de controles"), tags$dd(m$control_var)
        ),
        if (!is.na(m$dna_var)) tagList(
          tags$dt("Variable DNA"), tags$dd(m$dna_var)
        ),
        tags$dt("Umbral P"),
        tags$dd(sprintf("%.2f", m$threshold)),
        tags$dt("Taxones evaluados"),
        tags$dd(format(m$n_total, big.mark = ",")),
        tags$dt("Contaminantes detectados"),
        tags$dd(tags$strong(sprintf("%d (%.1f%%)",
                                      m$n_contam,
                                      100 * m$n_contam / max(m$n_total, 1L)))),
        tags$dt("decontam"),
        tags$dd(paste0("v", m$decontam_v))
      )
    })

    # ----- Render: histograma de P -----
    output$pPlot <- renderPlot({
      r <- decontam_result()
      validate(need(r, "Sin resultados aún."))
      df <- r$df
      df <- df[!is.na(df$p), ]
      ggplot2::ggplot(df, ggplot2::aes(x = .data$p,
                                          fill = .data$contaminant)) +
        ggplot2::geom_histogram(bins = 30, color = "white") +
        ggplot2::geom_vline(xintercept = r$meta$threshold,
                             linetype = "dashed", color = "#5D263E") +
        ggplot2::scale_fill_manual(
          values = c(`TRUE` = "#C32D28", `FALSE` = "#6D8A7A"),
          labels = c(`TRUE` = "contaminante", `FALSE` = "no"),
          name   = NULL
        ) +
        ggplot2::labs(x = "P (decontam)", y = "Taxones") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    })

    # ----- Render: tabla -----
    output$contaminantTable <- DT::renderDT({
      r <- decontam_result()
      validate(need(r, ""))
      tbl <- r$df
      keep_cols <- c("taxon", "p", "contaminant",
                     intersect(c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
                               colnames(tbl)))
      tbl <- tbl[, keep_cols]
      tbl$p <- signif(tbl$p, 4)
      DT::datatable(
        tbl,
        rownames = FALSE,
        options = list(pageLength = 15, scrollX = TRUE, dom = "lfrtip")
      )
    })

    # ----- Descarga -----
    output$download_contam <- shiny::downloadHandler(
      filename = function() paste0("decontam-", Sys.Date(), ".csv"),
      content = function(file) {
        r <- decontam_result()
        shiny::req(r)
        utils::write.csv(r$df, file, row.names = FALSE)
      }
    )
  })
}
