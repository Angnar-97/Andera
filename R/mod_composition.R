# Módulo: Composición taxonómica (stacked bars)
#
# Visualización canónica del microbioma: barras apiladas con la abundancia
# relativa de los N taxones más frecuentes a un rango taxonómico dado, con el
# resto agregado como "Other". Soporta dos vistas:
#   - Por muestra: una barra por muestra
#   - Por grupo: media de abundancia relativa por grupo de metadata
#
# Implementación:
#   - tax_glom(ps, taxrank) colapsa al rango elegido
#   - transform_sample_counts → proporciones (TSS)
#   - Top N por abundancia media; el resto colapsado a "Other"
#   - psmelt + ggplot stacked bar

mod_composition_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Composición taxonómica", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Barras apiladas de abundancia relativa por rango taxonómico. ",
        "Vista canónica para identificar taxones dominantes y diferencias ",
        "compositivas entre muestras o grupos. Internamente: ",
        tags$code("tax_glom"), " + ", tags$code("transform_sample_counts"),
        " + agregación de los menos abundantes en “Other”."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        # ----- Parámetros -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
          bslib::card_body(
            uiOutput(ns("rank_ui")),
            numericInput(
              ns("top_n"), "Top N taxones (resto = \"Other\")",
              value = 10, min = 3, max = 30, step = 1
            ),
            radioButtons(
              ns("aggregation"), "Agregación",
              choices = c(
                "Por muestra"    = "sample",
                "Por grupo (media)" = "group"
              ),
              selected = "sample",
              inline   = FALSE
            ),
            uiOutput(ns("group_variable_ui")),
            tags$div(class = "andera-actions",
              actionButton(ns("update_composition"), "Actualizar",
                           class = "btn btn-primary",
                           icon = icon("arrow-right"))
            )
          )
        ),

        # ----- Plot + tabla -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("bar-chart"),
                             " Barras apiladas"),
          bslib::card_body(
            shinycssloaders::withSpinner(
              plotOutput(ns("compositionPlot"), height = "560px"), type = 5
            ),
            r_code_panel(ns),
            tags$div(class = "andera-actions",
              downloadButton(ns("download_composition"),
                             "Descargar gráfico (.png)",
                             class = "btn-outline-secondary"),
              downloadButton(ns("download_composition_csv"),
                             "Descargar tabla (.csv)",
                             class = "btn-outline-secondary")
            )
          )
        )
      )
    )
  )
}

mod_composition_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    # ----- UI dinámico: rangos disponibles en el phyloseq -----
    output$rank_ui <- renderUI({
      req(physeq())
      tt <- tryCatch(phyloseq::tax_table(physeq(), errorIfNULL = FALSE),
                     error = function(e) NULL)
      ranks <- if (!is.null(tt)) colnames(tt) else character(0)

      if (!length(ranks)) {
        return(tags$p(class = "andera-muted",
          "El phyloseq cargado no tiene tax_table; la composición no puede calcularse."
        ))
      }

      default <- if ("Phylum" %in% ranks) "Phylum" else ranks[1]
      selectInput(
        session$ns("rank"), "Rango taxonómico",
        choices = ranks, selected = default
      )
    })

    output$group_variable_ui <- renderUI({
      req(physeq())
      vars <- colnames(phyloseq::sample_data(physeq()))
      selectInput(
        session$ns("group_variable"),
        "Variable de agrupación",
        choices = c("(ninguna)", vars),
        selected = "(ninguna)"
      )
    })

    # ----- Resultado -----
    composition_data <- reactiveVal(NULL)

    observe({
      physeq()
      composition_data(NULL)
    })

    observeEvent(input$update_composition, {
      req(physeq(), input$rank, input$top_n)
      ps <- physeq()

      if (input$aggregation == "group") {
        gv <- input$group_variable
        if (is.null(gv) || gv == "(ninguna)" || !nzchar(gv)) {
          shinyalert::shinyalert(
            title = "Configuración inválida",
            text  = "Para agregar por grupo necesitas seleccionar una variable de metadata.",
            type  = "warning"
          )
          return()
        }
      }

      result <- tryCatch(
        withProgress(message = "Composición taxonómica", value = 0, {
          incProgress(0.30, detail = sprintf("Glomerando a %s", input$rank))
          ps_glom <- phyloseq::tax_glom(ps, taxrank = input$rank, NArm = TRUE)

          incProgress(0.20, detail = "Transformando a proporciones")
          ps_rel <- phyloseq::transform_sample_counts(
            ps_glom, function(x) {
              s <- sum(x)
              if (s == 0) x else x / s
            }
          )

          incProgress(0.30, detail = "Agregando \"Other\"")
          # Top N calculado sobre la matriz de proporciones (taxa × samples),
          # NO sobre el data.frame melted: psmelt expande cada taxón a una fila
          # por muestra, lo que double-cuenta muestras con muchos ceros y
          # sesga el ranking. Usamos rowMeans de la otu_table directamente.
          rank_col <- input$rank
          otu_rel <- as(phyloseq::otu_table(ps_rel), "matrix")
          if (!phyloseq::taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)
          mean_per_otu <- rowMeans(otu_rel, na.rm = TRUE)
          tt_glom <- as.data.frame(
            phyloseq::tax_table(ps_glom)@.Data, stringsAsFactors = FALSE
          )
          # Tras tax_glom, cada OTU representa un único taxón al rango pedido.
          # Mapeamos cada OTU a su nombre de rango y agregamos por nombre
          # (puede haber varios OTUs con el mismo nombre si el rango es alto).
          taxon_per_otu <- tt_glom[[rank_col]]
          mean_per_taxon <- tapply(mean_per_otu, taxon_per_otu, sum,
                                     na.rm = TRUE)
          mean_per_taxon <- sort(mean_per_taxon, decreasing = TRUE,
                                  na.last = TRUE)
          top_n_taxa <- head(names(mean_per_taxon), input$top_n)

          df <- phyloseq::psmelt(ps_rel)
          # `df` tiene columnas: OTU, Sample, Abundance, + sample_data + tax_table

          df$.fill <- ifelse(df[[rank_col]] %in% top_n_taxa,
                             as.character(df[[rank_col]]),
                             "Other")
          # Mantén el orden de top_n_taxa con "Other" al final
          df$.fill <- factor(df$.fill, levels = c(top_n_taxa, "Other"))

          incProgress(0.20, detail = "Empaquetando")
          list(
            df              = df,
            rank            = rank_col,
            top_n           = input$top_n,
            aggregation     = input$aggregation,
            group_variable  = if (input$aggregation == "group") input$group_variable else NULL,
            n_taxa_total    = phyloseq::ntaxa(ps_glom),
            n_samples       = phyloseq::nsamples(ps_glom)
          )
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en composición",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) composition_data(result)
    })

    # ----- ggplot reactivo -----
    composition_plot <- reactive({
      r <- composition_data()
      if (is.null(r)) return(NULL)
      df <- r$df

      if (r$aggregation == "group") {
        gv <- r$group_variable
        # Agregar por (grupo, taxon-fill)
        agg <- aggregate(
          df$Abundance,
          by  = list(group = df[[gv]], taxon = df$.fill),
          FUN = function(x) mean(x, na.rm = TRUE)
        )
        colnames(agg) <- c("group", "taxon", "abundance")
        # Re-aplicar el orden factorial con Other al final
        agg$taxon <- factor(agg$taxon,
                            levels = levels(df$.fill))

        p <- ggplot2::ggplot(
          agg,
          ggplot2::aes(x = .data$group, y = .data$abundance, fill = .data$taxon)
        ) +
          ggplot2::geom_col(position = "fill") +
          ggplot2::labs(
            x = gv,
            y = "Abundancia relativa media",
            fill = r$rank
          )
      } else {
        p <- ggplot2::ggplot(
          df,
          ggplot2::aes(
            x    = .data$Sample,
            y    = .data$Abundance,
            fill = .data$.fill
          )
        ) +
          ggplot2::geom_col(position = "fill") +
          ggplot2::labs(
            x    = "Muestra",
            y    = "Abundancia relativa",
            fill = r$rank
          ) +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
          )
      }

      n_levels <- length(levels(df$.fill))
      fill_palette <- composition_palette(n_levels)
      p +
        ggplot2::scale_fill_manual(values = fill_palette) +
        ggplot2::scale_y_continuous(labels = scales::percent_format()) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor   = ggplot2::element_blank(),
          legend.position    = "right",
          legend.text        = ggplot2::element_text(size = 9)
        )
    })

    output$compositionPlot <- renderPlot({
      p <- composition_plot()
      validate(need(p,
        "Pulsa 'Actualizar' para calcular la composición taxonómica."))
      p
    })

    # ----- Descargas -----
    output$download_composition <- download_plot(
      composition_plot, "composicion-taxonomica",
      width = 12, height = 7
    )

    # ----- Reproducir en R -----
    r_code_text <- reactive({
      r <- composition_data()
      if (is.null(r)) return(NULL)
      r_code_composition(r)
    })
    r_code_handlers(output, r_code_text, "composicion-taxonomica")

    output$download_composition_csv <- shiny::downloadHandler(
      filename = function() {
        paste0("composicion-taxonomica-", Sys.Date(), ".csv")
      },
      content = function(file) {
        r <- composition_data()
        shiny::req(r)
        if (r$aggregation == "group") {
          gv <- r$group_variable
          tbl <- aggregate(
            r$df$Abundance,
            by  = list(group = r$df[[gv]], taxon = r$df$.fill),
            FUN = function(x) mean(x, na.rm = TRUE)
          )
          colnames(tbl) <- c(gv, "taxon", "rel_abundance")
        } else {
          tbl <- data.frame(
            sample        = r$df$Sample,
            taxon         = r$df$.fill,
            rel_abundance = r$df$Abundance
          )
        }
        utils::write.csv(tbl, file, row.names = FALSE)
      }
    )
  })
}

# ---------------------------------------------------------------------------
# Helper: paleta para el stacked bar
# Usa la paleta editorial de Andera para los top-N y un gris neutro para
# "Other" — visualmente "Other" pesa menos que los taxones nombrados.
# ---------------------------------------------------------------------------
composition_palette <- function(n_levels) {
  if (n_levels <= 0) return(character(0))
  base <- c(
    "#4A6D5E", "#C32D28", "#C4A962", "#5B7D8C",
    "#8FA872", "#6B4A7A", "#8C3B3B", "#2E5D5C",
    "#A8907F", "#7E897F", "#D97757", "#5D8A75",
    "#9B7A4A", "#6B8E7F", "#A88B5F", "#3E5C4F",
    "#C9AF53", "#7D5A6F", "#5A7D6E", "#B89968"
  )
  # Si hay más niveles que colores en `base`, repetimos con
  # opacidad implícita vía colorRampPalette para mantener variedad.
  cols <- if (n_levels - 1L <= length(base)) {
    base[seq_len(n_levels - 1L)]
  } else {
    grDevices::colorRampPalette(base)(n_levels - 1L)
  }
  c(cols, "#9C9C9C")  # último = "Other" en gris neutro
}
