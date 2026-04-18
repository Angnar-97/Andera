# ----- Interfaz de la Aplicación ----

ui <- dashboardPage(
    
    
    # ----- Encabezado de la Dashboard -----  
    dashboardHeader(title =  h3(strong("The Doge Whisperer"), align = "center"),
                    # Aquí se añade el logo de la aplicación. 
                    tags$li(a(img(src = "logocelta3.png",
                                  title = "Whisperer", height = "40px"), 
                              style = "padding-top:5px; padding-bottom:5px;"),
                            class = "dropdown"),
                    titleWidth = 250
    ),
    
    
    # ----- Panel del Dashboard -----  
    dashboardSidebar(
        
        width = 200, 
        
        # Crear un menú
        sidebarMenu(id = "mySideMenu",
                    
                    # TabName es el identificador de cada elemento del menú.
                    # Con selected especificas que ese ítem del menú sea el que aparezca al iniciar la app. 
                    menuItem("Estudio", tabName = "EstudioID", icon = icon("cloud", lib = "glyphicon"), selected = TRUE,
                             menuSubItem("Marco de Datos", tabName = "dataTableID"),
                             menuSubItem("Resumen Estadístico", tabName = "summaryID")),
                    
                    menuItem("Análisis", tabName = "AnalysisID", icon = icon("cd", lib = "glyphicon"),
                             menuSubItem("Visualización", tabName = "graphicsID"),
                             menuSubItem("Detección de Valores Atípicos", tabName = "outlierID"),
                             menuSubItem("Correlación", tabName = "corrID")),
                    
                    menuItem("Predicciones", tabName = "PredictionsID", icon = icon("grain", lib = "glyphicon"),
                             menuSubItem("Modelo Seleccionado", tabName = "modelID"),
                             menuSubItem("Predicción", tabName = "predictionID"),
                             menuSubItem("Evaluación", tabName = "EvaluateID")),
                    
                    menuItem("Aprendizaje Profundo", tabName = "DeepID", icon = icon("cog", lib = "glyphicon"),
                             menuSubItem("LSTM Seleccionada", tabName = "kerasMID"),
                             menuSubItem("Predicción", tabName = "predictionKerasID"),
                             menuSubItem("Evaluación", tabName = "EvaluateKerasID"))
                    
        ),
        
        # Barra horizontal para separar secciones en el panel. 
        hr(),
        hr(),
        
        # Selección de las criptomonedas. 
        checkboxGroupButtons("cryptocoin", h5(strong("Selecciona una o varias Criptodivisas"), style = "color : #2E3440;"), 
                             choices = list("Bitcoin"   = "BTC/USDT",
                                            "Ethereum"  = "ETH/USDT",
                                            "Polkadot"  = "DOT/USDT",
                                            "Cardano"   = "ADA/USDT",
                                            "Ripple"    = "XRP/USDT",
                                            "Litecoin"  = "LTC/USDT",
                                            "Solana"    = "SOL/USDT",
                                            "Binance"   = "BNB/USDT",
                                            "Chainlink" = "LINK/USDT",
                                            "Uniswap"   = "UNI/USDT"
                             ), 
                             justified = FALSE, direction = "horizontal", individual = TRUE, width = "450px",
                             selected = c("SOL/USDT", "ADA/USDT", "DOT/USDT")),
        
        hr(),
        hr()
        
    ),
    
    # ----- Cuerpo del Dashboard ----- 
    dashboardBody(
        
        # Se añade el tema que se ha creado en fresh para la Dashboard. 
        use_theme(mytheme),
        
        # Pestañas: se pone tabItems para ligar cada elemento a los ítems del menú creado. 
        tabItems(
            
            # -----  Primera pestaña -----
            tabItem(tabName = "dataTableID",
                    
                    # Está en desuso, pero aquí se ponen las alertas shiny que más adelante se añaden en el servidor. 
                    useShinyalert(),
                    
                    fluidRow(
                        
                        # Entrada para el límite inferior de la fecha. 
                        column(4, textInput("FromDateTime", "Desde:", value = "2021-06-23 00:00:00")), 
                        
                        # Entrada para el límite superior de la fecha.
                        column(4, textInput("TillDateTime", "Hasta:", value = "2021-06-23 03:00:00")),
                        
                    ),
                    
                    # Para dar un salto de línea. 
                    br(),
                    
                    fluidRow(
                        
                        
                        # Botón para procesar los datos y ejecutar la estandarización. 
                        column(1, offset = 6, actionButton("scaleID", "Estandarizar",
                                                           icon = icon("qrcode", lib = "glyphicon"))),
                        
                        # Botón para reescalar los datos. 
                        column(1, offset = 1, actionButton("rescaleID", "Reescalar",
                                                           icon = icon("qrcode", lib = "glyphicon"))),
                        
                        # Botón para el guardado del marco de datos. 
                        column(1, offset = 1, downloadButton("downloadDF", "Descargar",
                                                             icon = icon("download-alt", lib = "glyphicon")))
                        
                    ),
                    
                    br(),
                    
                    fluidRow(
                        
                        # Tabla que muestre el marco de datos. 
                        column(12, DT::dataTableOutput("dataframeID"))
                        
                    )
            ), 
            
            # ----- Segunda pestaña -----
            
            tabItem(tabName = "summaryID",
                    
                    tabBox(id = "tabset1", height = "1000px", width = 12, 
                           
                           tabPanel("Estadísticos",
                                    
                                    # Tabla con los principales estadísticos del marco de datos.
                                    fluidRow(column(12, DT::dataTableOutput("SummaryF")))
                           ),
                           
                           tabPanel("Resumen Skim",
                                    
                                    # Resumen estadístico de los datos.
                                    fluidRow(verbatimTextOutput("SkimF"), width = 12)
                           )
                    )
            ),
            
            
            # ----- Tercera pestaña -----
            
            tabItem(tabName = "graphicsID",
                    
                    tabBox(
                        
                        id = "tabset2", height = "1000px", width = 12, 
                        
                        tabPanel("Serie Temporal",
                                 
                                 fluidRow(
                                     
                                     # Serie temporal de las criptodivisas seleccionadas. 
                                     column(width = 12, plotOutput(outputId = "closeLineID"))
                                     
                                 ),
                                 
                                 fluidRow(
                                     
                                     # Histograma con el número de transacciones en base al tiempo. 
                                     column(width = 12, plotOutput(outputId = "transactionsHistogramID"))
                                 )
                                 
                        ),
                        
                        tabPanel("Volumen",
                                 
                                 fluidRow(
                                     
                                     # Distribución del volumen de criptodivisa.
                                     column(width = 12, plotOutput(outputId = "volumeDensityID")),
                                     
                                     
                                 )
                        )
                        
                    )
                    
            ),
            
            # ----- Cuarta pestaña -----
            
            tabItem(tabName = "outlierID",
                    
                    
                    tabBox(id = "tabset3", height = "1000px", width = 12, 
                           
                           tabPanel("Tabla de Valores Atípicos",
                                    
                                    
                                    fluidRow(
                                        
                                        # Lista desplegable para seleccionar la variable en la que estudiar los valores atípicos. 
                                        column(3, selectInput("selectVarOutlier", h5(strong("Selecciona una variable:"), style = "color : #2E3440;"), 
                                                              choices = varquant, selected = "tradecount")), 
                                        
                                        # Botón para el guardado de la tabla de los valores atípicos.
                                        column(1, offset = 7, downloadButton("downloadOutlierDF", "Descargar",
                                                                             icon = icon("download-alt", lib = "glyphicon")))
                                        
                                        
                                    ),
                                    
                                    # Tabla de los valores atípicos. 
                                    fluidRow(column(12, DT::dataTableOutput("outlierDFID")))
                           ),
                           
                           tabPanel("Diagramas de Cajas",
                                    
                                    # Gráfica con los diagramas de bigotes de las criptotransacciones.  
                                    fluidRow(column(12, plotOutput(outputId = "boxplotsOutlierID")))
                           ),
                           
                           tabPanel("Histogramas",
                                    
                                    # Gráfica con histogramas de cada criptomoneda. Compara la media del conjunto con la de cada moneda.
                                    fluidRow(column(12, plotOutput(outputId = "histogramOutlierID")))
                           ) 
                    )
            ),  
            
            # ----- Quinta pestaña -----
            
            tabItem(tabName = "corrID",
                    
                    tabBox(id = "tabset4", height = "1000px", width = 12, 
                           
                           tabPanel("Matriz de Correlaciones",
                                    
                                    fluidRow(
                                        
                                        # Botón para el guardado de la matriz de correlaciones. 
                                        column(1, offset = 10, downloadButton("downloadCorrDF", "Descargar",
                                                                              icon = icon("download-alt", lib = "glyphicon")))
                                    ),
                                    
                                    br(),
                                    
                                    # Marco de datos con las matrices de correlaciones. 
                                    fluidRow(column(12, DT::dataTableOutput(outputId = "corrDFID")))
                           ),
                           
                           tabPanel("Diagrama de Correlación",
                                    
                                    # Diagramas de las matrices de correlaciones de cada criptomoneda. 
                                    fluidRow(column(12, plotOutput(outputId = "corrPlotID")))
                           ),
                           
                           tabPanel("Mapa de Calor",
                                    
                                    # Mapa de calor de las correlaciones del marco general con .  
                                    fluidRow(column(12, d3heatmapOutput(outputId = "heatPlotID")))
                           )
                    )
            ),
            
            # ----- Sexta pestaña -----
            
            tabItem(tabName = "modelID",
                    
                    
                    fluidRow(
                        
                        # Lista desplegable para seleccionar un modelo. 
                        column(5, selectInput("selectionCaret", h5(strong("Selecciona un Modelo de Regresión Lineal:"), 
                                                                   style = "color : #2E3440;"), 
                                              choices = models, selected = "models/modeloSolana.rds")), 
                        
                        # Botón para ejecutar la predicción. 
                        column(2, offset = 4, actionButton("predictButtonCaret", "Predecir",
                                                           icon = icon("screenshot", lib = "glyphicon")))
                        
                    ),
                    
                    # Rseumen del modelo que incluye sus principales parámetros. 
                    fluidRow(column(12, verbatimTextOutput("modelSummary")))
                    
            ),
            
            # ----- Séptima pestaña -----
            
            tabItem(tabName = "predictionID",
                    
                    fluidRow(
                        
                        # Botón para el guardado de las predicciones obtenidas por los modelos de regresión. 
                        column(1, offset = 10, downloadButton("downloadCaretDF", "Descargar",
                                                              icon = icon("download-alt", lib = "glyphicon")))
                    ),
                    
                    br(),
                    
                    # Marco de datos con las predicciones de los modelos de regresión, los valores reales y la fecha. 
                    fluidRow(column(12, DT::dataTableOutput(outputId = "predictionsCaretDF")))
                    
                    
            ),
            
            # ----- Octava pestaña -----
            tabItem(tabName = "EvaluateID",
                    
                    fluidRow(
                        
                        # Caja que recoge el valor del coeficiente de determinación. 
                        valueBoxOutput(width = 3, outputId = "R2CaretID"),
                        
                        # Caja que recoge el valor del coeficiente de determinación con factor de corrección. 
                        valueBoxOutput(width = 3, outputId = "BR2CaretID"),
                        
                        # Caja que recoge el valor del error absoluto medio. 
                        valueBoxOutput(width = 3, outputId = "MAECaretID"), 
                        
                        # Caja que recoge el valor del error cuadrático medio.  
                        valueBoxOutput(width = 3, outputId = "MSECaretID") 
                        
                    ),
                    
                    br(),
                    
                    # Gráfica que enseña la serie temporal de valores reales frente a las predicciones. 
                    fluidRow(column(12, plotOutput(outputId = "predictionsPlotCaretID")))
                    
                    
            ),
            
            # ----- Novena pestaña -----
            # Aquí se añadirán los LSTM
            tabItem(tabName = "kerasMID",
                    
                    
                    fluidRow(
                        
                        # Lista desplegable para seleccionar un modelo de Aprendizaje Profundo. 
                        column(5, selectInput("selectionKeras", h5(strong("Selecciona un Modelo de Aprendizaje Profundo:"),
                                                                   style = "color : #2E3440;"), 
                                              choices = modelsKeras, selected = "models/LSTMSolana.h5")), 
                        
                        # Botón para ejecutar la predicción. 
                        column(2, offset = 4, actionButton("predictButtonKeras", "Predecir",
                                                           icon = icon("screenshot", lib = "glyphicon")))
                        
                    ),
                    
                    # Se muestra un resumen del modelo LSTM obtenido por Keras. 
                    fluidRow(column(12, verbatimTextOutput("modelKerasSummary")))
                    
            ),
            
            # ----- Décima pestaña -----
            tabItem(tabName = "predictionKerasID",
                    
                    fluidRow(
                        
                        # Botón para el guardado de las predicciones obtenidas por los modelos de Keras.  
                        column(1, offset = 10, downloadButton("downloadKerasDF", "Descargar",
                                                              icon = icon("download-alt", lib = "glyphicon")))
                    ),
                    
                    br(),
                    
                    # Marco de datos con las predicciones obtenidas por los modelos de Keras, los valores reales y la fecha. 
                    fluidRow(column(12, DT::dataTableOutput(outputId = "predictionsKerasDF")))
                    
            ),
            
            # ----- Undécima pestaña -----
            tabItem(tabName = "EvaluateKerasID",
                    
                    fluidRow(
                        
                        # Caja que recoge el valor del coeficiente de determinación.
                        valueBoxOutput(width = 3, outputId = "R2KerasID"), 
                        
                        # Caja que recoge el valor del coeficiente de determinación con factor de corrección.
                        valueBoxOutput(width = 3, outputId = "BR2KerasID"),
                        
                        # Caja que recoge el valor del error absoluto medio.
                        valueBoxOutput(width = 3, outputId = "MAEKerasID"), 
                        
                        # Caja que recoge el valor del error cuadrático medio.
                        valueBoxOutput(width = 3, outputId = "MSEKerasID") 
                        
                    ),
                    
                    br(),
                    
                    # Serie temporal de los valores reales enfrentados a las predicciones dadas por las LSTM. 
                    fluidRow(column(12, plotOutput(outputId = "predictionsPlotKerasID")))
                    
            )
            
            
            
        )
        
        
    )
    
)