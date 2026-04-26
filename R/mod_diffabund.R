# Módulo: Abundancia diferencial
#
# Identificación de taxones con abundancia diferencial entre dos grupos
# usando DESeq2 (NB GLM, Love 2014) y ALDEx2 (Monte Carlo Dirichlet + CLR,
# Fernandes 2014). Reporta ambos métodos en paralelo y la intersección,
# siguiendo la recomendación de Nearing 2022 (Nat Commun 13:342) de no
# confiar en un único método dado que producen resultados discordantes.
#
# Entrada: phyloseq con conteos enteros + variable de agrupación con ≥2 niveles.
# Pre-procesamiento configurable: glomerado a un rango taxonómico (default
# Genus, reduce drásticamente el coste computacional sin perder señal
# biológica) + filtrado por prevalencia mínima.

mod_diffabund_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("Abundancia diferencial", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Taxones con abundancia diferencial entre dos grupos. Calcula en ",
        "paralelo ", tags$b("DESeq2"), " (modelo NB con shrinkage) y ",
        tags$b("ALDEx2"), " (composicional, basado en CLR). La intersección ",
        "de ambos métodos es el conjunto más conservador y reproducible — ",
        "métodos de DA discrepan con frecuencia, así que reportar uno solo es ",
        "engañoso (Nearing 2022)."
      ),

      bslib::layout_columns(
        col_widths = c(4, 8),

        # ----- Parámetros -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("sliders"), " Parámetros"),
          bslib::card_body(

            tags$span(class = "andera-eyebrow", "Comparación"),
            uiOutput(ns("grouping_variable_ui")),
            uiOutput(ns("group_a_ui")),
            uiOutput(ns("group_b_ui")),

            tags$hr(),
            tags$span(class = "andera-eyebrow", "Pre-procesamiento"),
            uiOutput(ns("rank_ui")),
            numericInput(
              ns("min_prevalence"),
              "Prevalencia mínima (fracción de muestras con conteo > 0)",
              value = 0.10, min = 0, max = 1, step = 0.05
            ),
            tags$small(class = "andera-form-help",
              "Filtrar taxones presentes en menos del X % de muestras antes ",
              "de testear reduce ruido y testeo múltiple. 0 desactiva el filtro."
            ),

            tags$hr(),
            tags$span(class = "andera-eyebrow", "Test"),
            checkboxGroupInput(
              ns("methods"), "Métodos",
              choices  = c("DESeq2", "ALDEx2"),
              selected = c("DESeq2", "ALDEx2"),
              inline   = TRUE
            ),
            numericInput(
              ns("padj_threshold"),
              "Umbral de p ajustado (BH-FDR)",
              value = 0.05, min = 0, max = 1, step = 0.01
            ),
            numericInput(
              ns("log2fc_threshold"),
              "Umbral |log₂FC| (DESeq2 only, 0 = sin filtro)",
              value = 1, min = 0, max = 10, step = 0.5
            ),

            tags$div(class = "andera-actions",
              actionButton(ns("update_da"), "Ejecutar test",
                           class = "btn btn-primary",
                           icon = icon("play"))
            )
          )
        ),

        # ----- Resultados -----
        bslib::card(
          bslib::card_header(bsicons::bs_icon("table"), " Resultados"),
          bslib::card_body(

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Resumen"),
              uiOutput(ns("summary"))
            ),

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow",
                        "Volcano plot · log₂FC vs −log₁₀(p ajustado)"),
              shinycssloaders::withSpinner(
                plotOutput(ns("volcanoPlot"), height = "440px"), type = 5
              )
            ),

            tags$div(class = "andera-result-section",
              tags$span(class = "andera-eyebrow", "Tabla de resultados"),
              shinycssloaders::withSpinner(
                DT::DTOutput(ns("daTable")), type = 5
              )
            ),

            r_code_panel(ns),

            tags$div(class = "andera-actions",
              downloadButton(ns("download_volcano"),
                             "Descargar volcano (.png)",
                             class = "btn-outline-secondary"),
              downloadButton(ns("download_da_csv"),
                             "Descargar tabla (.csv)",
                             class = "btn-outline-secondary")
            )
          )
        )
      )
    )
  )
}

mod_diffabund_server <- function(id, physeq) {
  moduleServer(id, function(input, output, session) {

    # ----- UI dinámica -----
    output$grouping_variable_ui <- renderUI({
      req(physeq())
      vars <- colnames(phyloseq::sample_data(physeq()))
      selectInput(
        session$ns("grouping_variable"),
        "Variable de agrupación",
        choices = vars
      )
    })

    group_levels <- reactive({
      req(physeq(), input$grouping_variable)
      sd <- as(phyloseq::sample_data(physeq()), "data.frame")
      sort(unique(as.character(sd[[input$grouping_variable]])))
    })

    output$group_a_ui <- renderUI({
      req(group_levels())
      selectInput(
        session$ns("group_a"), "Grupo A",
        choices  = group_levels(),
        selected = group_levels()[1]
      )
    })

    output$group_b_ui <- renderUI({
      req(group_levels())
      lvls <- group_levels()
      selectInput(
        session$ns("group_b"), "Grupo B",
        choices  = lvls,
        selected = if (length(lvls) >= 2) lvls[2] else lvls[1]
      )
    })

    output$rank_ui <- renderUI({
      req(physeq())
      tt <- tryCatch(phyloseq::tax_table(physeq(), errorIfNULL = FALSE),
                     error = function(e) NULL)
      ranks <- if (!is.null(tt)) colnames(tt) else character(0)
      choices <- c("(sin glom)" = "__none__", ranks)
      default <- if ("Genus" %in% ranks) "Genus" else "__none__"
      selectInput(
        session$ns("glom_rank"),
        "Glomerar al rango taxonómico",
        choices  = choices,
        selected = default
      )
    })

    # ----- Estado -----
    da_result <- reactiveVal(NULL)

    observe({
      physeq()
      da_result(NULL)
    })

    # ----- Cómputo -----
    observeEvent(input$update_da, {
      req(physeq(), input$grouping_variable, input$group_a, input$group_b,
          input$methods)
      ps     <- physeq()
      var    <- input$grouping_variable
      gA     <- input$group_a
      gB     <- input$group_b
      meths  <- input$methods
      padj_t <- input$padj_threshold
      lfc_t  <- input$log2fc_threshold

      # ---- Validaciones ----
      if (gA == gB) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = "Los grupos A y B deben ser distintos.", type = "warning"
        )
        return()
      }
      v <- validate_for_count_indices(ps)
      if (!isTRUE(v$ok)) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = paste(
            "DESeq2 y ALDEx2 requieren conteos enteros (modelo NB / Dirichlet).",
            "El phyloseq cargado contiene proporciones — no aplica."
          ),
          type  = "warning"
        )
        return()
      }
      v <- validate_for_grouping(ps, var, min_n_per_group = 3L)
      if (!isTRUE(v$ok)) {
        shinyalert::shinyalert(
          title = "Configuración inválida",
          text  = v$reason, type = "warning"
        )
        return()
      }

      # ---- Cómputo ----
      result <- tryCatch(
        withProgress(message = "Abundancia diferencial", value = 0, {

          incProgress(0.05, detail = "Filtrando muestras de los 2 grupos")
          sd <- as(phyloseq::sample_data(ps), "data.frame")
          keep_samples <- sd[[var]] %in% c(gA, gB) & !is.na(sd[[var]])
          ps_2 <- phyloseq::prune_samples(keep_samples, ps)

          # Garantizar que hay ≥3 por cada grupo retenido
          tab <- table(as(phyloseq::sample_data(ps_2), "data.frame")[[var]])
          if (any(tab[c(gA, gB)] < 3)) {
            stop(sprintf(
              "Tras seleccionar A=%s vs B=%s, algún grupo tiene <3 muestras (n_A=%d, n_B=%d).",
              gA, gB, tab[gA], tab[gB]
            ))
          }

          # Glomerar si se pidió
          if (input$glom_rank != "__none__") {
            incProgress(0.05,
                        detail = sprintf("Glomerando a %s", input$glom_rank))
            ps_2 <- phyloseq::tax_glom(ps_2, taxrank = input$glom_rank,
                                       NArm = TRUE)
          }

          # Filtro de prevalencia
          if (input$min_prevalence > 0) {
            incProgress(0.05, detail = "Filtrado por prevalencia")
            otu <- as(phyloseq::otu_table(ps_2), "matrix")
            if (!phyloseq::taxa_are_rows(ps_2)) otu <- t(otu)
            prev <- rowMeans(otu > 0)
            keep_taxa <- prev >= input$min_prevalence
            ps_2 <- phyloseq::prune_taxa(keep_taxa, ps_2)
          }

          n_taxa_tested <- phyloseq::ntaxa(ps_2)
          if (n_taxa_tested < 5) {
            stop(sprintf(
              "Tras el pre-procesamiento solo quedan %d taxa para testear; relaja la prevalencia o usa otro rango.",
              n_taxa_tested
            ))
          }

          # ----- DESeq2 -----
          deseq_tbl <- NULL
          if ("DESeq2" %in% meths) {
            incProgress(0.30, detail = "DESeq2 (NB GLM)")
            deseq_tbl <- run_deseq2(ps_2, var, gA, gB)
          }

          # ----- ALDEx2 -----
          aldex_tbl <- NULL
          if ("ALDEx2" %in% meths) {
            incProgress(0.30, detail = "ALDEx2 (Monte Carlo CLR)")
            aldex_tbl <- run_aldex2(ps_2, var, gA, gB)
          }

          # ----- Merge ambos resultados -----
          incProgress(0.20, detail = "Combinando resultados")
          combined <- merge_da_results(deseq_tbl, aldex_tbl, padj_t, lfc_t)

          # ----- Anotar tax_table del phyloseq glomerado -----
          tt <- tryCatch(phyloseq::tax_table(ps_2, errorIfNULL = FALSE),
                          error = function(e) NULL)
          if (!is.null(tt)) {
            tax_df <- as.data.frame(tt@.Data, stringsAsFactors = FALSE)
            tax_df$taxon <- rownames(tax_df)
            combined <- merge(combined, tax_df, by = "taxon", all.x = TRUE)
          }

          incProgress(0.05, detail = "Empaquetando")
          n_sig_d <- if (!is.null(deseq_tbl))
            sum(combined$sig_DESeq2, na.rm = TRUE) else NA_integer_
          n_sig_a <- if (!is.null(aldex_tbl))
            sum(combined$sig_ALDEx2, na.rm = TRUE) else NA_integer_
          n_sig_both <- if (!is.null(deseq_tbl) && !is.null(aldex_tbl))
            sum(combined$sig_DESeq2 & combined$sig_ALDEx2, na.rm = TRUE)
          else NA_integer_

          list(
            table = combined,
            meta  = list(
              variable          = var,
              group_a           = gA,
              group_b           = gB,
              n_a               = unname(tab[gA]),
              n_b               = unname(tab[gB]),
              glom_rank         = input$glom_rank,
              min_prevalence    = input$min_prevalence,
              n_taxa_tested     = n_taxa_tested,
              methods           = meths,
              padj_threshold    = padj_t,
              log2fc_threshold  = lfc_t,
              n_sig_deseq2      = n_sig_d,
              n_sig_aldex2      = n_sig_a,
              n_sig_both        = n_sig_both,
              deseq2_v          = if ("DESeq2" %in% meths)
                as.character(utils::packageVersion("DESeq2")) else NA_character_,
              aldex2_v          = if ("ALDEx2" %in% meths)
                as.character(utils::packageVersion("ALDEx2")) else NA_character_
            )
          )
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en abundancia diferencial",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) da_result(result)
    })

    # ----- Renders -----
    output$summary <- renderUI({
      r <- da_result()
      if (is.null(r)) {
        return(tags$p(class = "andera-muted",
          "Configura los parámetros y pulsa 'Ejecutar test'."))
      }
      m <- r$meta
      tags$dl(class = "andera-meta-dl",
        tags$dt("Comparación"),
        tags$dd(sprintf("%s = %s vs %s", m$variable, m$group_a, m$group_b)),

        tags$dt("n por grupo"),
        tags$dd(sprintf("A = %d, B = %d", m$n_a, m$n_b)),

        tags$dt("Pre-procesamiento"),
        tags$dd(sprintf(
          "glom = %s · prevalencia ≥ %.0f%% · %d taxa testeados",
          if (m$glom_rank == "__none__") "ninguno" else m$glom_rank,
          100 * m$min_prevalence,
          m$n_taxa_tested
        )),

        tags$dt("Umbral de significación"),
        tags$dd(sprintf("p_adj ≤ %.3f%s",
          m$padj_threshold,
          if (m$log2fc_threshold > 0)
            sprintf(" · |log₂FC| ≥ %.1f", m$log2fc_threshold)
          else ""
        )),

        if ("DESeq2" %in% m$methods) tagList(
          tags$dt("DESeq2 significativos"),
          tags$dd(sprintf("%d taxa (v%s)", m$n_sig_deseq2, m$deseq2_v))
        ),
        if ("ALDEx2" %in% m$methods) tagList(
          tags$dt("ALDEx2 significativos"),
          tags$dd(sprintf("%d taxa (v%s)", m$n_sig_aldex2, m$aldex2_v))
        ),
        if (length(m$methods) == 2L) tagList(
          tags$dt("Intersección"),
          tags$dd(tags$strong(sprintf("%d taxa", m$n_sig_both)))
        )
      )
    })

    volcano_plot <- reactive({
      r <- da_result()
      if (is.null(r)) return(NULL)
      build_volcano(r$table, r$meta)
    })

    output$volcanoPlot <- renderPlot({
      p <- volcano_plot()
      validate(need(p,
        "Pulsa 'Ejecutar test' para ver el volcano plot."))
      p
    })

    output$daTable <- DT::renderDT({
      r <- da_result()
      validate(need(r, "Pulsa 'Ejecutar test' para ver la tabla."))
      tbl <- r$table
      # Ordenar por padj_DESeq2 si existe, si no por padj_ALDEx2
      sort_col <- if ("padj_DESeq2" %in% colnames(tbl)) "padj_DESeq2" else "padj_ALDEx2"
      tbl <- tbl[order(tbl[[sort_col]], na.last = TRUE), ]

      # Redondear numéricas
      num_cols <- vapply(tbl, is.numeric, logical(1))
      for (col in names(num_cols)[num_cols]) {
        tbl[[col]] <- signif(tbl[[col]], 4)
      }
      DT::datatable(
        tbl,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          pageLength = 15,
          scrollX    = TRUE,
          dom        = "lfrtip"
        )
      )
    })

    output$download_volcano <- download_plot(
      volcano_plot, "volcano-da",
      width = 9, height = 6
    )

    # ----- Reproducir en R -----
    r_code_text <- reactive({
      r <- da_result()
      if (is.null(r)) return(NULL)
      r_code_diffabund(r$meta)
    })
    r_code_handlers(output, r_code_text, "abundancia-diferencial")

    output$download_da_csv <- shiny::downloadHandler(
      filename = function() paste0("abundancia-diferencial-", Sys.Date(), ".csv"),
      content = function(file) {
        r <- da_result()
        shiny::req(r)
        m <- r$meta

        con <- file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)

        meta_lines <- c(
          "# Differential abundance results - Andera",
          sprintf("# Date: %s", Sys.Date()),
          sprintf("# Variable: %s", m$variable),
          sprintf("# Group A vs B: %s vs %s (n = %d vs %d)",
                  m$group_a, m$group_b, m$n_a, m$n_b),
          sprintf("# Glom rank: %s",
                  if (m$glom_rank == "__none__") "none" else m$glom_rank),
          sprintf("# Min prevalence: %.2f", m$min_prevalence),
          sprintf("# Taxa tested: %d", m$n_taxa_tested),
          sprintf("# Methods: %s", paste(m$methods, collapse = ", ")),
          sprintf("# padj threshold: %.3f", m$padj_threshold),
          sprintf("# log2FC threshold: %.2f", m$log2fc_threshold),
          if (!is.na(m$n_sig_deseq2)) sprintf("# DESeq2 significant: %d", m$n_sig_deseq2),
          if (!is.na(m$n_sig_aldex2)) sprintf("# ALDEx2 significant: %d", m$n_sig_aldex2),
          if (!is.na(m$n_sig_both))   sprintf("# Intersection: %d",     m$n_sig_both),
          if (!is.na(m$deseq2_v)) sprintf("# DESeq2 version: %s", m$deseq2_v),
          if (!is.na(m$aldex2_v)) sprintf("# ALDEx2 version: %s", m$aldex2_v),
          ""
        )
        writeLines(meta_lines, con)
        utils::write.table(
          r$table, file = con, sep = ",",
          row.names = FALSE, col.names = TRUE, quote = TRUE
        )
      }
    )
  })
}

# ===========================================================================
# Helpers fuera del módulo
# ===========================================================================

# Geometric mean sobre conteos POSITIVOS (Anders 2010 / phyloseq vignette).
# La media se hace sobre x[x > 0], no sobre length(x), para evitar el sesgo a
# la baja en filas sparse. Si todos son cero o NA, devuelve NA (la fila
# debería filtrarse aguas arriba, no entrar al modelo NB).
gm_mean_pos <- function(x, na.rm = TRUE) {
  pos <- x[is.finite(x) & x > 0]
  if (!length(pos)) return(NA_real_)
  exp(mean(log(pos), na.rm = na.rm))
}

# Ejecuta DESeq2 sobre un phyloseq con la fórmula simple ~ var.
# Devuelve data.frame con: taxon, baseMean, log2FoldChange, lfcSE, stat,
# pvalue, padj. El contraste es B vs A (B en el numerador).
run_deseq2 <- function(ps, var, group_a, group_b) {
  sd <- as(phyloseq::sample_data(ps), "data.frame")
  sd[[var]] <- factor(sd[[var]], levels = c(group_a, group_b))
  phyloseq::sample_data(ps) <- phyloseq::sample_data(sd)

  dds <- phyloseq::phyloseq_to_deseq2(ps, stats::reformulate(var))
  geo <- apply(DESeq2::counts(dds), 1, gm_mean_pos)
  # Filtrar taxa con geometric mean inválido (todos los conteos a 0 → NA)
  keep_taxa <- is.finite(geo) & geo > 0
  if (!all(keep_taxa)) dds <- dds[keep_taxa, , drop = FALSE]
  geo <- geo[keep_taxa]
  dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geo)
  dds <- DESeq2::DESeq(dds, fitType = "local", quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c(var, group_b, group_a))

  out <- as.data.frame(res)
  out$taxon <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c("taxon", "baseMean", "log2FoldChange", "lfcSE",
                 "stat", "pvalue", "padj")]
  colnames(out) <- c("taxon", "baseMean_DESeq2", "log2FC_DESeq2",
                     "lfcSE_DESeq2", "stat_DESeq2",
                     "pvalue_DESeq2", "padj_DESeq2")
  out
}

# Ejecuta ALDEx2. Devuelve data.frame con: taxon, effect_ALDEx2,
# diff_btw_ALDEx2, padj_ALDEx2 (usa we.eBH, Welch t-test BH-ajustada).
run_aldex2 <- function(ps, var, group_a, group_b) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  storage.mode(otu) <- "integer"
  sd <- as(phyloseq::sample_data(ps), "data.frame")
  # ALDEx2 usa el primer nivel encontrado como referencia.
  # Forzamos el factor a c(group_a, group_b) para que `effect` y `diff.btw`
  # se interpreten como B vs A (mismo sentido que DESeq2 contrast = B/A).
  conds <- as.character(
    factor(sd[[var]], levels = c(group_a, group_b))
  )
  ax <- ALDEx2::aldex(
    reads      = otu,
    conditions = conds,
    test       = "t",
    effect     = TRUE,
    verbose    = FALSE,
    denom      = "all"
  )

  out <- data.frame(
    taxon            = rownames(ax),
    effect_ALDEx2    = ax$effect,
    diff_btw_ALDEx2  = ax$diff.btw,
    pvalue_ALDEx2    = ax$we.ep,
    padj_ALDEx2      = ax$we.eBH,
    stringsAsFactors = FALSE
  )
  out
}

# Combina las dos tablas (any o ambas) en una con un único taxon por fila
# y banderas de significación según los thresholds dados.
merge_da_results <- function(deseq_tbl, aldex_tbl, padj_t, lfc_t) {
  combined <- if (!is.null(deseq_tbl) && !is.null(aldex_tbl)) {
    merge(deseq_tbl, aldex_tbl, by = "taxon", all = TRUE)
  } else if (!is.null(deseq_tbl)) {
    deseq_tbl
  } else {
    aldex_tbl
  }

  if ("padj_DESeq2" %in% colnames(combined)) {
    sig_d <- !is.na(combined$padj_DESeq2) & combined$padj_DESeq2 <= padj_t
    if (lfc_t > 0 && "log2FC_DESeq2" %in% colnames(combined)) {
      sig_d <- sig_d & abs(combined$log2FC_DESeq2) >= lfc_t
    }
    combined$sig_DESeq2 <- sig_d
  }
  if ("padj_ALDEx2" %in% colnames(combined)) {
    combined$sig_ALDEx2 <- !is.na(combined$padj_ALDEx2) &
                           combined$padj_ALDEx2 <= padj_t
  }
  combined
}

# Volcano plot. Si hay DESeq2 → x = log2FC_DESeq2. Si solo ALDEx2 →
# x = effect_ALDEx2 (efecto de Cohen-like). Y siempre = -log10(p_adj).
build_volcano <- function(tbl, meta) {
  has_d <- "padj_DESeq2"  %in% colnames(tbl)
  has_a <- "padj_ALDEx2"  %in% colnames(tbl)

  if (has_d) {
    df <- data.frame(
      x        = tbl$log2FC_DESeq2,
      y        = -log10(tbl$padj_DESeq2),
      sig_d    = if ("sig_DESeq2" %in% colnames(tbl)) tbl$sig_DESeq2 else FALSE,
      sig_a    = if (has_a && "sig_ALDEx2" %in% colnames(tbl)) tbl$sig_ALDEx2 else FALSE,
      taxon    = tbl$taxon
    )
    x_lab <- "log₂FC (DESeq2: B / A)"
  } else {
    df <- data.frame(
      x        = tbl$effect_ALDEx2,
      y        = -log10(tbl$padj_ALDEx2),
      sig_d    = FALSE,
      sig_a    = if ("sig_ALDEx2" %in% colnames(tbl)) tbl$sig_ALDEx2 else FALSE,
      taxon    = tbl$taxon
    )
    x_lab <- "Effect size (ALDEx2)"
  }

  df <- df[is.finite(df$x) & is.finite(df$y), ]

  df$category <- with(df, ifelse(sig_d & sig_a, "Ambos",
                          ifelse(sig_d,         "DESeq2",
                          ifelse(sig_a,         "ALDEx2", "n.s."))))
  df$category <- factor(df$category,
                         levels = c("Ambos", "DESeq2", "ALDEx2", "n.s."))

  cat_palette <- c(
    "Ambos"  = "#5D263E",  # burgundy
    "DESeq2" = "#C32D28",  # terracota
    "ALDEx2" = "#5B7D8C",  # nord blue
    "n.s."   = "#B0B0B0"
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                          color = .data$category)) +
    ggplot2::geom_point(alpha = 0.75, size = 2) +
    ggplot2::geom_hline(yintercept = -log10(meta$padj_threshold),
                         linetype = "dashed", color = "#666666",
                         alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                         color = "#666666", alpha = 0.4) +
    ggplot2::scale_color_manual(values = cat_palette, drop = FALSE) +
    ggplot2::labs(
      x     = x_lab,
      y     = "−log₁₀ (p ajustado)",
      color = "Significación",
      title = sprintf("Volcano: %s = %s vs %s",
                       meta$variable, meta$group_a, meta$group_b)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right"
    )

  if (meta$log2fc_threshold > 0 && has_d) {
    p <- p + ggplot2::geom_vline(
      xintercept = c(-meta$log2fc_threshold, meta$log2fc_threshold),
      linetype = "dotted", color = "#666666", alpha = 0.4
    )
  }
  p
}
