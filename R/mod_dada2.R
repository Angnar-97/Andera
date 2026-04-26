# Módulo: DADA2 — FASTQ a phyloseq (small-batch)
#
# Pipeline interactivo para construir un objeto phyloseq directamente desde
# FASTQs paired-end de secuenciación de amplicones (16S rRNA, ITS, 18S).
#
# Pre-condiciones que el módulo asume y comunica al usuario:
#   1. Lecturas paired-end Illumina con nombres `*_R1.fastq.gz` /
#      `*_R2.fastq.gz` (o `_1` / `_2`).
#   2. Primers ya recortados (cutadapt previo al pipeline).
#   3. Tamaño manejable: cap explícito de 20 muestras × 50 MB por fichero
#      para mantener interactividad. Para datasets reales, descargar el
#      script R generado y correrlo con Rscript.
#
# Pipeline (Callahan 2016, Nature Methods 13:581):
#   plotQualityProfile → filterAndTrim → learnErrors → dada → mergePairs →
#   makeSequenceTable → removeBimeraDenovo → phyloseq(otu_table)
#
# Salida: phyloseq sin tax_table (la asignación taxonómica con SILVA exige
# descargar referencias multi-GB, fuera del scope de un dashboard).
# El objeto se descarga como .rds y se carga después en la pestaña de
# "Carga" usando la opción de archivo.

# Constantes
DADA2_MAX_SAMPLES   <- 20L
DADA2_MAX_FILE_MB   <- 50L
DADA2_MAX_TOTAL_MB  <- 500L

mod_dada2_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "andera-page",
      tags$h2("DADA2 — FASTQ a phyloseq", class = "andera-page-title"),
      tags$p(class = "andera-page-lead",
        "Pipeline ", tags$a(href = "https://benjjneb.github.io/dada2/",
                              target = "_blank", rel = "noopener", "DADA2"),
        " para inferir ", tags$strong("ASVs"),
        " (Amplicon Sequence Variants) a partir de FASTQs paired-end. ",
        "Pensado para ", tags$strong("datasets pequeños"),
        " (≤ 20 muestras, ≤ 50 MB/fichero) — para producción, descarga el ",
        "script R generado y córrelo con Rscript en local."
      ),

      # ----- Pre-condiciones -----
      tags$div(class = "andera-alert andera-alert-info",
        tags$strong("Antes de subir:"),
        tags$ul(
          tags$li("Lecturas ", tags$b("paired-end Illumina"),
                  " con nombres ", tags$code("*_R1.fastq.gz"), " y ",
                  tags$code("*_R2.fastq.gz"),
                  " (o ", tags$code("_1"), "/", tags$code("_2"), ")."),
          tags$li("Primers ", tags$b("ya recortados"),
                  " — DADA2 no los elimina (usa cutadapt previo)."),
          tags$li(sprintf(
            "Cap por sesión interactiva: %d muestras, %d MB/fichero, %d MB totales.",
            DADA2_MAX_SAMPLES, DADA2_MAX_FILE_MB, DADA2_MAX_TOTAL_MB)),
          tags$li("La asignación taxonómica (SILVA, GreenGenes) ",
                  tags$b("no"), " se incluye — exige referencias multi-GB; ",
                  "el phyloseq generado tiene OTU table + secuencias.")
        )
      ),

      # =====================================================================
      # 1. Upload
      # =====================================================================
      bslib::card(
        bslib::card_header(bsicons::bs_icon("upload"),
                            " 1 · Subir FASTQs paired-end"),
        bslib::card_body(
          fileInput(
            ns("fastq_files"),
            "Selecciona todos los FASTQ (R1 y R2 mezclados; el módulo los empareja)",
            multiple   = TRUE,
            accept     = c(".fastq", ".fastq.gz", ".fq", ".fq.gz")
          ),
          uiOutput(ns("upload_summary"))
        )
      ),

      # =====================================================================
      # 2. Quality profile
      # =====================================================================
      bslib::card(
        bslib::card_header(bsicons::bs_icon("activity"),
                            " 2 · Quality profile (R1 + R2)"),
        bslib::card_body(
          tags$p(class = "andera-form-help",
            "Inspecciona el perfil de calidad para decidir ",
            tags$code("truncLen"),
            " — corta donde la calidad cae por debajo de Q20."
          ),
          actionButton(ns("compute_qc"), "Calcular quality profiles",
                       class = "btn btn-primary btn-sm",
                       icon = icon("play")),
          tags$div(class = "andera-result-section",
            shinycssloaders::withSpinner(
              plotOutput(ns("qcPlotF"), height = "260px"), type = 5
            ),
            shinycssloaders::withSpinner(
              plotOutput(ns("qcPlotR"), height = "260px"), type = 5
            )
          )
        )
      ),

      # =====================================================================
      # 3. Filter + ASV inference
      # =====================================================================
      bslib::card(
        bslib::card_header(bsicons::bs_icon("filter"),
                            " 3 · Filter, infer ASVs, merge & remove chimeras"),
        bslib::card_body(
          bslib::layout_columns(
            col_widths = c(6, 6),
            tagList(
              numericInput(ns("trunclen_f"),
                           "truncLen R1 (lecturas forward)",
                           value = 240, min = 100, max = 300, step = 10),
              numericInput(ns("trunclen_r"),
                           "truncLen R2 (lecturas reverse)",
                           value = 200, min = 100, max = 300, step = 10),
              tags$small(class = "andera-form-help",
                "Recorta cada lectura a esta longitud antes de inferir ASVs. ",
                "Ajusta según el quality profile: cuando Q baja de 20-25, ",
                "trunca por encima."
              )
            ),
            tagList(
              numericInput(ns("max_ee_f"),
                           "maxEE R1 (errores esperados máx.)",
                           value = 2, min = 0.5, max = 10, step = 0.5),
              numericInput(ns("max_ee_r"),
                           "maxEE R2",
                           value = 2, min = 0.5, max = 10, step = 0.5),
              tags$small(class = "andera-form-help",
                tags$code("maxEE = 2"),
                " es el default Callahan 2016. Aumentar si las lecturas son ",
                "de alta calidad (HiSeq); bajar si calidades dudosas (NovaSeq)."
              )
            )
          ),
          tags$div(class = "andera-actions",
            actionButton(ns("run_pipeline"),
                          "Ejecutar pipeline DADA2",
                          class = "btn btn-primary",
                          icon = icon("play"))
          )
        )
      ),

      # =====================================================================
      # 4. Results
      # =====================================================================
      bslib::card(
        bslib::card_header(bsicons::bs_icon("table"), " 4 · Resultados"),
        bslib::card_body(
          tags$div(class = "andera-result-section",
            tags$span(class = "andera-eyebrow", "Resumen del pipeline"),
            uiOutput(ns("pipeline_summary"))
          ),
          tags$div(class = "andera-result-section",
            tags$span(class = "andera-eyebrow",
                      "Tracking de lecturas por etapa"),
            tableOutput(ns("trackTable"))
          ),
          tags$div(class = "andera-actions",
            downloadButton(ns("download_phyloseq"),
                            "Descargar phyloseq (.rds)",
                            class = "btn-outline-secondary"),
            downloadButton(ns("download_track"),
                            "Tracking (.csv)",
                            class = "btn-outline-secondary"),
            downloadButton(ns("download_r_pipeline"),
                            "Script R reproducible (.R)",
                            class = "btn-outline-secondary")
          )
        )
      )
    )
  )
}

mod_dada2_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # ----- Estado -----
    qc_data       <- reactiveVal(NULL)
    pipeline_data <- reactiveVal(NULL)

    # ----- Parsear ficheros subidos: emparejar R1/R2 -----
    paired_files <- reactive({
      f <- input$fastq_files
      if (is.null(f) || nrow(f) == 0L) return(NULL)
      pair_fastq_files(f)
    })

    # ----- Reset al subir nuevos ficheros -----
    observeEvent(input$fastq_files, {
      qc_data(NULL)
      pipeline_data(NULL)
    })

    # ----- Resumen del upload con validaciones -----
    output$upload_summary <- renderUI({
      pf <- paired_files()
      if (is.null(pf)) {
        return(tags$p(class = "andera-muted",
          "Esperando ficheros."))
      }
      issues <- pf$issues
      ok     <- length(issues) == 0L

      content <- list(
        tags$dl(class = "andera-meta-dl",
          tags$dt("Ficheros subidos"),
          tags$dd(sprintf("%d", nrow(input$fastq_files))),
          tags$dt("Pares R1/R2 detectados"),
          tags$dd(sprintf("%d", length(pf$samples))),
          tags$dt("Tamaño total"),
          tags$dd(sprintf("%.1f MB", pf$total_mb))
        )
      )

      if (!ok) {
        content <- c(content, list(tags$div(
          class = "andera-alert andera-alert-warning",
          tags$strong("Avisos de validación:"),
          tags$ul(lapply(issues, tags$li))
        )))
      } else {
        content <- c(content, list(tags$div(
          class = "andera-alert andera-alert-success",
          tags$strong(sprintf("✓ %d muestras emparejadas correctamente",
                                length(pf$samples)))
        )))
      }
      do.call(tagList, content)
    })

    # ----- 2. Quality profile -----
    observeEvent(input$compute_qc, {
      pf <- paired_files()
      req(pf, length(pf$issues) == 0L)

      result <- tryCatch(
        withProgress(message = "Quality profile", value = 0, {
          incProgress(0.40, detail = "plotQualityProfile R1")
          # Limitar a las primeras 4 muestras para mantener interactividad
          n_show <- min(4L, length(pf$samples))
          plot_F <- dada2::plotQualityProfile(pf$forward[seq_len(n_show)])
          incProgress(0.40, detail = "plotQualityProfile R2")
          plot_R <- dada2::plotQualityProfile(pf$reverse[seq_len(n_show)])
          incProgress(0.20, detail = "Empaquetando")
          list(forward = plot_F, reverse = plot_R, n_shown = n_show)
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en quality profile",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) qc_data(result)
    })

    output$qcPlotF <- renderPlot({
      r <- qc_data()
      validate(need(r, "Pulsa 'Calcular quality profiles'."))
      r$forward
    })
    output$qcPlotR <- renderPlot({
      r <- qc_data()
      validate(need(r, ""))
      r$reverse
    })

    # ----- 3. Run full pipeline -----
    observeEvent(input$run_pipeline, {
      pf <- paired_files()
      req(pf, length(pf$issues) == 0L)
      req(input$trunclen_f, input$trunclen_r,
          input$max_ee_f, input$max_ee_r)

      result <- tryCatch(
        withProgress(message = "Pipeline DADA2", value = 0, {

          # Crear directorio temporal para los filtrados
          incProgress(0.05, detail = "Preparando directorios")
          tdir <- tempfile("dada2_filt_")
          dir.create(tdir, recursive = TRUE)
          on.exit(unlink(tdir, recursive = TRUE), add = TRUE)
          fnFs.filt <- file.path(tdir, paste0(pf$samples, "_F_filt.fastq.gz"))
          fnRs.filt <- file.path(tdir, paste0(pf$samples, "_R_filt.fastq.gz"))

          # 1. filterAndTrim
          incProgress(0.15,
                      detail = sprintf("filterAndTrim (truncLen=%d/%d, maxEE=%g/%g)",
                                       input$trunclen_f, input$trunclen_r,
                                       input$max_ee_f,   input$max_ee_r))
          filt_out <- dada2::filterAndTrim(
            pf$forward, fnFs.filt,
            pf$reverse, fnRs.filt,
            truncLen     = c(input$trunclen_f, input$trunclen_r),
            maxEE        = c(input$max_ee_f,   input$max_ee_r),
            truncQ       = 2, rm.phix = TRUE,
            compress     = TRUE, multithread = FALSE,
            verbose      = FALSE
          )

          # Comprobar que sobreviven muestras
          surviving <- filt_out[, "reads.out"] > 0L
          if (!any(surviving)) {
            stop("Tras filterAndTrim no quedan lecturas en ninguna muestra. ",
                 "Relaja truncLen o maxEE.")
          }
          fnFs.filt <- fnFs.filt[surviving]
          fnRs.filt <- fnRs.filt[surviving]
          samples_kept <- pf$samples[surviving]

          # 2. learnErrors
          incProgress(0.20, detail = "learnErrors R1")
          errF <- dada2::learnErrors(fnFs.filt, multithread = FALSE,
                                       verbose = FALSE)
          incProgress(0.10, detail = "learnErrors R2")
          errR <- dada2::learnErrors(fnRs.filt, multithread = FALSE,
                                       verbose = FALSE)

          # 3. dada
          incProgress(0.15, detail = "dada R1")
          dadaFs <- dada2::dada(fnFs.filt, err = errF,
                                  multithread = FALSE, verbose = 0)
          incProgress(0.10, detail = "dada R2")
          dadaRs <- dada2::dada(fnRs.filt, err = errR,
                                  multithread = FALSE, verbose = 0)

          # 4. mergePairs
          incProgress(0.10, detail = "mergePairs")
          mergers <- dada2::mergePairs(dadaFs, fnFs.filt,
                                          dadaRs, fnRs.filt, verbose = FALSE)

          # 5. makeSequenceTable
          incProgress(0.05, detail = "makeSequenceTable")
          seqtab <- dada2::makeSequenceTable(mergers)

          # 6. removeBimeraDenovo
          incProgress(0.05, detail = "removeBimeraDenovo")
          seqtab.nochim <- dada2::removeBimeraDenovo(
            seqtab, method = "consensus", multithread = FALSE, verbose = FALSE
          )

          # Track de lecturas por etapa
          getN <- function(x) sum(dada2::getUniques(x))
          track <- cbind(
            input    = filt_out[surviving, "reads.in"],
            filtered = filt_out[surviving, "reads.out"],
            denoisedF = vapply(if (is.list(dadaFs)) dadaFs else list(dadaFs), getN, numeric(1)),
            denoisedR = vapply(if (is.list(dadaRs)) dadaRs else list(dadaRs), getN, numeric(1)),
            merged   = vapply(if (is.list(mergers)) mergers else list(mergers), getN, numeric(1)),
            nonchim  = rowSums(seqtab.nochim)
          )
          rownames(track) <- samples_kept

          # 7. phyloseq object
          incProgress(0.05, detail = "Construyendo phyloseq")
          ps <- phyloseq::phyloseq(
            phyloseq::otu_table(seqtab.nochim, taxa_are_rows = FALSE),
            phyloseq::sample_data(data.frame(
              sample_id = samples_kept,
              row.names = samples_kept,
              stringsAsFactors = FALSE
            ))
          )

          list(
            phyloseq      = ps,
            seqtab.nochim = seqtab.nochim,
            track         = track,
            params = list(
              n_samples = length(samples_kept),
              n_asv     = ncol(seqtab.nochim),
              trunclen  = c(input$trunclen_f, input$trunclen_r),
              maxEE     = c(input$max_ee_f, input$max_ee_r),
              dada2_v   = as.character(utils::packageVersion("dada2"))
            )
          )
        }),
        error = function(e) {
          shinyalert::shinyalert(
            title = "Error en el pipeline DADA2",
            text  = e$message, type = "error"
          )
          NULL
        }
      )
      if (!is.null(result)) pipeline_data(result)
    })

    # ----- Render: resumen del pipeline -----
    output$pipeline_summary <- renderUI({
      r <- pipeline_data()
      if (is.null(r)) {
        return(tags$p(class = "andera-muted",
          "Ejecuta el pipeline para ver los resultados."))
      }
      p <- r$params
      tags$dl(class = "andera-meta-dl",
        tags$dt("Muestras retenidas"),
        tags$dd(sprintf("%d", p$n_samples)),
        tags$dt("ASVs únicos"),
        tags$dd(format(p$n_asv, big.mark = ",")),
        tags$dt("truncLen"),
        tags$dd(sprintf("R1 = %d, R2 = %d", p$trunclen[1], p$trunclen[2])),
        tags$dt("maxEE"),
        tags$dd(sprintf("R1 = %g, R2 = %g", p$maxEE[1], p$maxEE[2])),
        tags$dt("dada2"),
        tags$dd(paste0("v", p$dada2_v))
      )
    })

    # ----- Tabla de tracking -----
    output$trackTable <- renderTable({
      r <- pipeline_data()
      validate(need(r, ""))
      df <- as.data.frame(r$track)
      df$sample <- rownames(df)
      df[, c("sample", colnames(r$track))]
    }, rownames = FALSE, digits = 0)

    # ----- Descargas -----
    output$download_phyloseq <- shiny::downloadHandler(
      filename = function() paste0("dada2-phyloseq-", Sys.Date(), ".rds"),
      content = function(file) {
        r <- pipeline_data()
        shiny::req(r)
        saveRDS(r$phyloseq, file)
      }
    )

    output$download_track <- shiny::downloadHandler(
      filename = function() paste0("dada2-tracking-", Sys.Date(), ".csv"),
      content = function(file) {
        r <- pipeline_data()
        shiny::req(r)
        df <- as.data.frame(r$track)
        df$sample <- rownames(df)
        utils::write.csv(df[, c("sample", colnames(r$track))],
                         file, row.names = FALSE)
      }
    )

    output$download_r_pipeline <- shiny::downloadHandler(
      filename = function() paste0("dada2-pipeline-", Sys.Date(), ".R"),
      content = function(file) {
        r <- pipeline_data()
        if (is.null(r)) {
          writeLines("# Ejecuta el pipeline para generar el script reproducible.",
                     file)
          return()
        }
        writeLines(r_code_dada2_pipeline(r$params), file)
      }
    )
  })
}

# ===========================================================================
# Helpers
# ===========================================================================

# Empareja R1/R2 a partir de los ficheros subidos. Devuelve una lista con:
#   - samples  : character vector con nombres de muestra
#   - forward  : datapaths a los ficheros R1 (ordenados por samples)
#   - reverse  : datapaths a los ficheros R2 (idem)
#   - total_mb : tamaño total en MB
#   - issues   : character vector de problemas detectados (vacío si OK)
pair_fastq_files <- function(file_df) {
  total_mb <- sum(file_df$size) / 1024^2

  # Regex permisiva: cubre ambas convenciones
  #   short:  sample_R1.fastq.gz, sample_1.fq.gz
  #   illumina: Sample_S1_L001_R1_001.fastq.gz (sufijo lane/index variable)
  fwd_pat <- "_R?1(_\\d+)?\\.(fastq|fq)(\\.gz)?$"
  rev_pat <- "_R?2(_\\d+)?\\.(fastq|fq)(\\.gz)?$"

  is_forward <- grepl(fwd_pat, file_df$name, ignore.case = TRUE)
  is_reverse <- grepl(rev_pat, file_df$name, ignore.case = TRUE)

  forward_files <- file_df[is_forward, , drop = FALSE]
  reverse_files <- file_df[is_reverse, , drop = FALSE]

  # Extraer stem de la muestra (todo lo que precede al _R?[12](_NNN)?.fastq...)
  forward_files$sample <- sub(fwd_pat, "", forward_files$name,
                                ignore.case = TRUE)
  reverse_files$sample <- sub(rev_pat, "", reverse_files$name,
                                ignore.case = TRUE)

  common <- intersect(forward_files$sample, reverse_files$sample)
  common <- sort(common)

  # Reordenar por sample
  fnFs <- forward_files$datapath[match(common, forward_files$sample)]
  fnRs <- reverse_files$datapath[match(common, reverse_files$sample)]

  issues <- character(0)
  unmatched_F <- setdiff(forward_files$sample, common)
  unmatched_R <- setdiff(reverse_files$sample, common)
  if (length(unmatched_F))
    issues <- c(issues, sprintf("R1 sin pareja R2: %s",
                                  paste(unmatched_F, collapse = ", ")))
  if (length(unmatched_R))
    issues <- c(issues, sprintf("R2 sin pareja R1: %s",
                                  paste(unmatched_R, collapse = ", ")))
  if (length(common) > DADA2_MAX_SAMPLES)
    issues <- c(issues, sprintf("Más de %d muestras (recibidas %d) — el cap interactivo es %d.",
                                  DADA2_MAX_SAMPLES, length(common),
                                  DADA2_MAX_SAMPLES))
  big_files <- which(file_df$size > DADA2_MAX_FILE_MB * 1024^2)
  if (length(big_files))
    issues <- c(issues,
      sprintf("Ficheros >%d MB: %s",
              DADA2_MAX_FILE_MB,
              paste(file_df$name[big_files], collapse = ", ")))
  if (total_mb > DADA2_MAX_TOTAL_MB)
    issues <- c(issues,
      sprintf("Total >%d MB (%.1f) — recorta el batch.",
              DADA2_MAX_TOTAL_MB, total_mb))

  list(
    samples  = common,
    forward  = fnFs,
    reverse  = fnRs,
    total_mb = total_mb,
    issues   = issues
  )
}

# Genera el script R reproducible con los parámetros usados.
r_code_dada2_pipeline <- function(params) {
  c(
    "# Andera · Pipeline DADA2 reproducible",
    sprintf("# Generado: %s", Sys.Date()),
    sprintf("# dada2 version: %s", params$dada2_v),
    "",
    "library(dada2)",
    "library(phyloseq)",
    "",
    "# Ajusta estas rutas a tus FASTQs paired-end:",
    "fastq_dir <- '<RUTA_A_TUS_FASTQS>'",
    "fnFs <- sort(list.files(fastq_dir,",
    "                pattern = '_R?1(_\\\\d+)?\\\\.(fastq|fq)(\\\\.gz)?$',",
    "                full.names = TRUE))",
    "fnRs <- sort(list.files(fastq_dir,",
    "                pattern = '_R?2(_\\\\d+)?\\\\.(fastq|fq)(\\\\.gz)?$',",
    "                full.names = TRUE))",
    "samples <- sub('_R?1(_\\\\d+)?\\\\.(fastq|fq)(\\\\.gz)?$', '',",
    "                 basename(fnFs))",
    "",
    "# Directorio para los ficheros filtrados",
    "filt_dir <- file.path(fastq_dir, 'filtered')",
    "dir.create(filt_dir, showWarnings = FALSE)",
    "fnFs.filt <- file.path(filt_dir, paste0(samples, '_F_filt.fastq.gz'))",
    "fnRs.filt <- file.path(filt_dir, paste0(samples, '_R_filt.fastq.gz'))",
    "",
    "# 1. filter & trim",
    "filt_out <- dada2::filterAndTrim(",
    "  fnFs, fnFs.filt, fnRs, fnRs.filt,",
    sprintf("  truncLen = c(%d, %d),", params$trunclen[1], params$trunclen[2]),
    sprintf("  maxEE    = c(%g, %g),", params$maxEE[1],    params$maxEE[2]),
    "  truncQ   = 2, rm.phix = TRUE,",
    "  compress = TRUE, multithread = TRUE",
    ")",
    "",
    "# 2. error model",
    "errF <- dada2::learnErrors(fnFs.filt, multithread = TRUE)",
    "errR <- dada2::learnErrors(fnRs.filt, multithread = TRUE)",
    "",
    "# 3. ASV inference",
    "dadaFs <- dada2::dada(fnFs.filt, err = errF, multithread = TRUE)",
    "dadaRs <- dada2::dada(fnRs.filt, err = errR, multithread = TRUE)",
    "",
    "# 4. merge pairs",
    "mergers <- dada2::mergePairs(dadaFs, fnFs.filt, dadaRs, fnRs.filt)",
    "",
    "# 5. sequence table + chimera removal",
    "seqtab        <- dada2::makeSequenceTable(mergers)",
    "seqtab.nochim <- dada2::removeBimeraDenovo(",
    "  seqtab, method = 'consensus', multithread = TRUE",
    ")",
    "",
    "# 6. asignación taxonómica (descarga SILVA / Greengenes2 antes)",
    "# taxa <- dada2::assignTaxonomy(seqtab.nochim, '<silva_nr99_v138.fa.gz>',",
    "#                                multithread = TRUE)",
    "",
    "# 7. construir phyloseq",
    "ps <- phyloseq::phyloseq(",
    "  phyloseq::otu_table(seqtab.nochim, taxa_are_rows = FALSE)",
    "  # phyloseq::tax_table(taxa),                 # si hay taxonomía",
    "  # phyloseq::sample_data(metadata)            # si hay metadata",
    ")",
    "",
    "saveRDS(ps, 'phyloseq.rds')"
  )
}
