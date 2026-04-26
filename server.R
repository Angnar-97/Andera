# ----- Servidor de la Aplicación ----
# Orquesta los módulos de R/ y comparte el phyloseq activo entre pestañas.

server <- function(input, output, session) {
  # ----- Portada + Contacto -----
  mod_home_server   ("home", session)
  mod_contact_server("contact")

  # ----- Carga + filtrado -----
  physeq          <- mod_data_load_server("data_load")
  physeq_filtered <- mod_filtering_server("filtering", physeq)

  # ----- DADA2 stub -----
  mod_dada2_server("dada2")

  # ----- QC: detección de contaminantes -----
  mod_decontam_server("decontam", physeq)

  # ----- Banner global -----
  mod_dataset_banner_server("banner", physeq, physeq_filtered)

  # Phyloseq "activo": el filtrado si existe, si no el original.
  active_physeq <- reactive({
    filt <- physeq_filtered()
    if (!is.null(filt)) filt else physeq()
  })

  # ----- Módulos de análisis -----
  mod_rarefaction_server("rarefaction", active_physeq)
  mod_composition_server("composition", active_physeq)
  mod_diversity_server  ("diversity",   active_physeq)
  mod_heatmaps_server   ("heatmaps",    active_physeq)
  mod_ordination_server ("ordination",  active_physeq)
  mod_permanova_server  ("permanova",   active_physeq)
  mod_diffabund_server  ("diffabund",   active_physeq)
  mod_graphs_server     ("grafos",      active_physeq)
  mod_microbial_network_server("microbial_net", active_physeq)
  mod_microbiome_views_server ("microbiome_views", active_physeq)
}
