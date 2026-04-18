# ----- Application Interface ----
ui <- dashboardPage(
  
  # ----- Dashboard Header -----   
  dashboardHeader(title = "Andera - Exploración de Microbiomas"),
  
  # ----- Dashboard Panel ----- 
  dashboardSidebar(
    # Código para la barra lateral
    sidebarMenu(
      menuItem("Carga de Datos", tabName = "data_load", icon = icon("database")),
      menuItem("Filtrado", tabName = "filtering", icon = icon("yin-yang")),
      menuItem("Diversidad Alfa y Beta", tabName = "diversity", icon = icon("bank")),
      menuItem("Mapas de Calor", tabName = "heatmaps", icon = icon("snowflake")),
      menuItem("PCoA", tabName = "pcoa", icon = icon("hurricane")),
      menuItem("PERMANOVA", tabName = "permanova", icon = icon("flask-vial")),
      menuItem("Grafos", tabName = "grafos", icon = icon("spider", lib = "font-awesome"))
    )
  ),
  
  # ----- Cuerpo del Dashboard ----- 
  dashboardBody(
    
    # Se añade el tema que se ha creado en fresh para la Dashboard. 
    includeCSS("themes/theme_cerulean.css"),
    
    
    tabItems(
      # Tab de Carga de Datos
      tabItem(tabName = "data_load",
              fluidRow(
                box(
                  title = "Carga de Datos",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  radioButtons(
                    "data_source", 
                    "Seleccione la Fuente de Datos",
                    choices = c("Archivo", "Conjunto de datos de ejemplo"),
                    selected = "Conjunto de datos de ejemplo"
                    ),
                  conditionalPanel(
                    condition = "input.data_source == 'Archivo'",
                    fileInput("file1", "Elija un archivo Phyloseq")
                    ),
                  conditionalPanel(
                    condition = "input.data_source == 'Conjunto de datos de ejemplo'", 
                    selectInput("dataset", "Seleccione un objeto phyloseq de muestra",
                                choices = c("Global Patterns", "Enterotype", "Soil"), 
                                selected = "Global Patterns")),
                  actionButton("update", "Actualizar"),
                  downloadButton("downloadData", "Descargar objeto Phyloseq")
                ),
                box(
                  title = "Resumen del Objeto Phyloseq",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  verbatimTextOutput("physeq_summary")
                ),
                box(
                  title = "Árbol Filogenético",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  plotOutput("phylo_tree")
                )
              )
      ),
      # Tab de Filtrado
      tabItem(tabName = "filtering",
              fluidRow(
                box(
                  title = "Filtrado",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  numericInput("min_count", "Número mínimo de conteos", value = 0, min = 0),
                  uiOutput("variableui"),
                  uiOutput("valueui"),
                  #numericInput("min_prevalence", "Prevalencia mínima", value = 0, min = 0, max = 1, step = 0.01),
                  actionButton("filter", "Filtrar")
                ),
                box(
                  title = "Resumen del Objeto Phyloseq",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  verbatimTextOutput("filtered_physeq_summary")
                ),
                box(
                  title = "Árbol Filogenético",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  plotOutput("filtered_phylo_tree")
                )
              )
      ),
      # Tab de Diversidad Alfa y Beta
      tabItem(tabName = "diversity",
              fluidRow(
                box(
                  title = "Estudio de la Riqueza",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  selectInput("diversity", "Seleccione el tipo de diversidad", 
                              choices = c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), 
                              multiple = TRUE, selected = c('Shannon', 'Simpson')),
                  uiOutput("alpha_variableui"),
                  actionButton("update_diversity", "Actualizar"),
                  #selectInput("distance", "Seleccione el tipo de distancia", choices = c("bray", "jaccard", "unifrac", "wunifrac")),
                ),
                box(
                  title = "Diversidad Alfa & Beta",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  plotOutput("diversityPlot")
                )
              )
      ),
      # Tab de Mapas de Calor
      tabItem(tabName = "heatmaps",
              fluidRow(
                box(
                  title = "Mapas de Calor",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  actionButton("update_heatmaps", "Actualizar"),
                ),
                box(
                  title = "Mapa de Calor de Abundancia de OTUs",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  plotOutput("heatmap_otus")
                ),
                box(
                  title = "Mapa de Calor de Abundancia de Especies",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  plotOutput("heatmap_species")
                )
              )
      ),
      # Tab de PCoA
      tabItem(tabName = "pcoa",
              fluidRow(
                box(
                  title = "PCoA",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  selectInput("distance_pcoa", "Seleccione el tipo de distancia", choices = c("bray", "jaccard", "unifrac", "wunifrac")),
                  plotOutput("pcoaPlot")
                )
              )
      ),
      # Tab de PERMANOVA
      tabItem(tabName = "permanova",
              fluidRow(
                box(
                  title = "PERMANOVA",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  selectInput("distance_permanova", "Seleccione el tipo de distancia", choices = c("bray", "jaccard", "unifrac", "wunifrac")),
                  tableOutput("permanovaResults")
                )
              )
      ),
      # Tab de Grafos
      tabItem(tabName = "grafos",
              fluidRow(
                box(
                  title = "Grafos",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  actionButton("update_graph", "Actualizar"),
                )        
      )
    )
  )
))