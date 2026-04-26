# ----- Interfaz de la Aplicación ----
# Top navigation bar plano con bslib::page_navbar.
# Las pestañas se agrupan visualmente por secciones (DATOS · ANÁLISIS · SOBRE)
# usando `nav_item` con la clase `andera-nav-section` como divisor no
# clickeable. Esto da la jerarquía sin esconder pestañas tras submenús.

ui <- tagList(
  # ----- Recursos externos (favicon, CSS, meta) -----
  tags$head(
    tags$link(rel = "icon",       type = "image/png", href = "celta.png"),
    tags$link(rel = "stylesheet", type = "text/css",  href = "andera.css"),

    # ----- SEO: meta básicas -----
    tags$meta(
      name    = "description",
      content = paste(
        "Andera · Explora microbiomas online sin código. Dashboard Shiny",
        "gratuito y de código abierto para el análisis exploratorio de",
        "microbiomas con phyloseq: carga de datos de secuenciación 16S rRNA,",
        "ITS y metagenómica, filtrado por metadata, diversidad alfa",
        "(Shannon, Simpson, Chao1, ACE, Fisher), diversidad beta",
        "(Bray-Curtis, Jaccard, UniFrac), ordenación PCoA,",
        "PERMANOVA con vegan::adonis2, mapas de calor y redes de",
        "co-ocurrencia entre muestras y taxa."
      )
    ),
    tags$meta(
      name    = "keywords",
      content = paste(
        "microbioma, análisis microbioma, metagenómica, metagenómica 16S,",
        "secuenciación 16S rRNA, ITS, 18S, amplicon sequencing,",
        "phyloseq, microViz, vegan, microbiome R, Bioconductor, Shiny, R,",
        "diversidad alfa, índice Shannon, índice Simpson, Chao1, ACE, Fisher,",
        "diversidad beta, Bray-Curtis, Jaccard, UniFrac, UniFrac ponderada,",
        "PCoA, NMDS, PERMANOVA, adonis2, ecología microbiana,",
        "comunidades microbianas, bioinformática, DADA2, QIIME 2"
      )
    ),
    tags$meta(name = "author",   content = "Alejandro Navas González"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$meta(name = "robots",   content = "index, follow"),

    # ----- Open Graph -----
    tags$meta(property = "og:title",
              content  = "Andera · Explora microbiomas online sin código"),
    tags$meta(property = "og:description",
              content  = paste(
                "Dashboard Shiny gratuito para el análisis exploratorio de",
                "microbiomas con phyloseq: diversidad alfa y beta,",
                "PCoA, PERMANOVA, mapas de calor y redes de co-ocurrencia."
              )),
    tags$meta(property = "og:image",  content = "celta.png"),
    tags$meta(property = "og:type",   content = "website"),
    tags$meta(property = "og:locale", content = "es_ES"),

    # ----- Twitter Card -----
    tags$meta(name = "twitter:card",        content = "summary"),
    tags$meta(name = "twitter:title",
              content = "Andera · Explora microbiomas online sin código"),
    tags$meta(name = "twitter:description",
              content = "Análisis interactivo de microbiomas con phyloseq, PCoA, PERMANOVA y redes — sin escribir una línea de R.")
  ),

  bslib::page_navbar(
    id    = "tabs",
    title = tags$span(
      class = "andera-brand",
      tags$img(src = "celta.png", class = "andera-logo", alt = "Andera — árbol celta"),
      tags$span(class = "andera-brand-text", "Andera")
    ),
    window_title = "Andera · Explora microbiomas online — phyloseq, PCoA, PERMANOVA",
    theme        = andera_theme(),
    fillable     = FALSE,
    padding      = 0,
    navbar_options = bslib::navbar_options(
      position    = "fixed-top",
      collapsible = TRUE
    ),

    # ----- Header global (banner del dataset activo) -----
    header = tagList(
      tags$div(class = "andera-navbar-spacer"),
      mod_dataset_banner_ui("banner")
    ),

    # ----- Footer editorial -----
    footer = tags$footer(
      class = "andera-footer",
      tags$div(class = "container-xl",
        tags$div(class = "andera-footer-grid",
          tags$div(class = "andera-footer-col",
            tags$span(class = "andera-eyebrow", "Andera"),
            tags$p(class = "andera-footer-lead",
              "Análisis exploratorio e interactivo de microbiomas con phyloseq, microViz y vegan.")
          ),
          tags$div(class = "andera-footer-col",
            tags$span(class = "andera-eyebrow", "Autoría"),
            tags$a(href = "https://github.com/Angnar-97",
                   target = "_blank", rel = "noopener",
                   "Alejandro Navas González"),
            tags$br(),
            tags$a(href = "mailto:angnar@telaris.es", "angnar@telaris.es")
          ),
          tags$div(class = "andera-footer-col",
            tags$span(class = "andera-eyebrow", "Código"),
            tags$a(href = "https://github.com/Angnar-97",
                   target = "_blank", rel = "noopener", "GitHub ↗"),
            tags$br(),
            tags$span(class = "andera-muted", "Proyecto independiente © 2026")
          )
        )
      )
    ),

    # =====================================================================
    # Pestañas (plano · agrupado por secciones via nav_item divisores)
    # =====================================================================

    bslib::nav_panel(
      title = "Inicio", value = "home",
      mod_home_ui("home")
    ),

    # ----- DATOS -----
    bslib::nav_item(tags$span(class = "andera-nav-section", "Datos")),
    bslib::nav_panel(
      title = "Carga", value = "data_load",
      mod_data_load_ui("data_load")
    ),
    bslib::nav_panel(
      title = "DADA2 (FASTQ)", value = "dada2",
      mod_dada2_ui("dada2")
    ),
    bslib::nav_panel(
      title = "Filtrado", value = "filtering",
      mod_filtering_ui("filtering")
    ),
    bslib::nav_panel(
      title = "Decontam", value = "decontam",
      mod_decontam_ui("decontam")
    ),

    # ----- ANÁLISIS -----
    bslib::nav_item(tags$span(class = "andera-nav-section", "Análisis")),
    bslib::nav_panel(
      title = "Rarefacción", value = "rarefaction",
      mod_rarefaction_ui("rarefaction")
    ),
    bslib::nav_panel(
      title = "Composición", value = "composition",
      mod_composition_ui("composition")
    ),
    bslib::nav_panel(
      title = "Diversidad", value = "diversity",
      mod_diversity_ui("diversity")
    ),
    bslib::nav_panel(
      title = "Heatmaps", value = "heatmaps",
      mod_heatmaps_ui("heatmaps")
    ),
    bslib::nav_panel(
      title = "Ordenación", value = "ordination",
      mod_ordination_ui("ordination")
    ),
    bslib::nav_panel(
      title = "PERMANOVA", value = "permanova",
      mod_permanova_ui("permanova")
    ),
    bslib::nav_panel(
      title = "Abundancia diferencial", value = "diffabund",
      mod_diffabund_ui("diffabund")
    ),
    bslib::nav_panel(
      title = "Vistas avanzadas", value = "microbiome_views",
      mod_microbiome_views_ui("microbiome_views")
    ),
    bslib::nav_panel(
      title = "Red muestras", value = "grafos",
      mod_graphs_ui("grafos")
    ),
    bslib::nav_panel(
      title = "Red taxa", value = "microbial_net",
      mod_microbial_network_ui("microbial_net")
    ),

    bslib::nav_spacer(),

    bslib::nav_panel(
      title = "Contacto", value = "contact",
      mod_contact_ui("contact")
    )
  )
)
