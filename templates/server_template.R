# ----- Servidor de la Aplicación ----
server <- function(input, output, session) {
    
    
    # ----- Variables Reactivas -----
    # Aplicación de filtros
    # Se crea una variable reactiva que responda a los inputs de entrada. 
    cryptoReactive <- reactive({
        
        newCrypto <- crypto %>% 
            filter(symbol %in% input$cryptocoin) %>%
            filter(date > input$FromDateTime & date < input$TillDateTime) 
        
        # Condición para evitar que el usuario elija intervalos de tiempo superiores a las seis horas. 
        if (abs(as.numeric(as.POSIXct(input$FromDateTime) - as.POSIXct(input$TillDateTime))) > 6L |
            abs(as.numeric(as.Date(input$FromDateTime) - as.Date(input$TillDateTime))) >= 1L){
            
            # Alerta de shiny que le indica al usuario la restricción del intervalo temporal. 
            shinyalert("Advertencia", "¡Se han de seleccionar intervalos inferiores a seis horas!", type = "warning")
            # Se vacía el marco de datos para evitar que el usuario colapse la aplicación. 
            newCrypto <- newCrypto[0,]
            return(newCrypto)
            
        } else {
            
            return(newCrypto)
            
        }
        
    })
    
    
    # Se crea una variable de valores reactivos (esto sirve para los observeEvents). 
    cryptoDF <- reactiveValues(data = NULL)
    
    # Se guarda en la variable de valores reactivos el marco resultante de los filtros. 
    observe({cryptoDF$data <- cryptoReactive()})
    
    # Valor reactivo que sirva de bandera para la estandarización. Se inicializa como falso. 
    dataScaled <- reactiveVal(FALSE)
    
    # Evento de estandarizado. 
    observeEvent(input$scaleID, {
        
        if(dataScaled() == FALSE){
            
            # Estandarización del marco de datos. 
            cryptoDF$data <- cryptoDF$data %>% dplyr::mutate_at(c("open", "high", "low", "close",
                                                                  "tradecount", "volume", "USDT"), scale)
            
            # Cambio de valor en la variable bandera. 
            dataScaled(TRUE)
        }
        
    })
    
    # Evento para regresar al estado normal tras estandarizar. 
    observeEvent(input$rescaleID, {
        
        if(dataScaled() == TRUE){
            
            # Reescalado del marco de datos. 
            cryptoDF$data <- cryptoDF$data %>% dplyr::mutate_at(c("open", "high", "low", "close",
                                                                  "tradecount", "volume", "USDT"), unscale)
            
            # Cambio de valor en la variable bandera. 
            dataScaled(FALSE)
            
        }
        
    })
    
    # Variable reactiva para manipular estéticamente la tabla de salida. 
    cryptotable <- reactive({
        
        # Guardar en una variable el marco de datos reactivo.
        # Esto sirve para poder cambiar el nombre de las columnas.
        # Si se utiliza el marco reactivo directamente para el renderizado da error. 
        cryptotable <- cryptoDF$data
        
        # Poner los nombres de cabecera que se quiera para la tabla. 
        colnames(cryptotable) <- c("Fecha", "Cripto", "Apertura",
                                   "Máximo", "Mínimo", "Cierre",
                                   "Volumen", "USDT", "Transacciones")
        
        # Se devuelve el objeto. 
        return(cryptotable)
        
    })
    
    
    # ----- Renderización de la tabla sobre criptomonedas ----- 
    output$dataframeID <- DT::renderDataTable({
        
        # Pasar a formato de data.table. Este paso no es necesario. 
        cryptotable() %>% data.table()
        
        
    }, rownames = FALSE, # Eliminar los números de las filas.
    extensions = 'Buttons', # Extensión para añadir botones en la propia DT. 
    # filter = list(position = 'top'),
    options = list(initComplete = JS("function(settings, json) {", 
                                     "$(this.api().table().header()).css({
                                        'background-color': '#3B4A52',
                                        'margin':'100px',
                                        'text-align': 'center',
                                        'font-family': 'Cambria',
                                        'color': '#B4C7D8'});",
                                     "}"), # Cabecera para la tabla con formato CSS. 
                   # container = sketch(),
                   searching = TRUE, # Para mostrar el buscador de la tabla.
                   searchHighlight = TRUE,
                   displayStart = 1,
                   autoWidth = TRUE,
                   drop = FALSE, # Con drop = FALSE se evita tener una tabla vacía.
                   orderClasses = TRUE, # Se permite ordenar las clases.
                   lengthMenu = c(15, 25, 35, 45, 55, 65, 75, 85), # Así se permite editar el número de sujetos de la tabla. 
                   pageLength = 25, # Número de filas mostradas por defecto. 
                   dom = "tplf", # Sirve para configurar el modo de visualización. 
                   # Con t sólo deja visible la tabla. 
                   # Con i muestras la info de la tabla en la parte inferior. 
                   # Con p muestra además la paginación.
                   # Con l se permite el input de control de entrada de filas por tabla.
                   # Con f se da entrada al filtro.
                   # Todas las letras anteriores se pueden combinar en el número y orden que se quiera para la configuración.
                   scrollY = "500px", # Permite crear una barra de scroll.
                   # Opciones de guardado para la tabla. Al tener un marco de datos grande esta opción no sirve. 
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print') 
    )
    )
    
    
    
    # Acción de guardado del marco de datos.
    output$downloadDF <- downloadHandler(
        
        filename = function(){
            
            # Nombre del archivo. 
            paste0("CryptoDataTable", as.numeric(as.Date(input$FromDateTime)), as.numeric(as.Date(input$TillDateTime)), ".xlsx")
            
        },
        
        content = function(file) {
            
            # Guardar el archivo como xlsx. Se pueden añadir más formatos.  
            write.xlsx(cryptoDF$data, file, row.names = FALSE)
            
        }
    )
    
    # ----- Resumen Estadístico -----
    output$SummaryF <- DT::renderDataTable({
        
        statsDF <- cryptoDF$data
        
        # Obtención del resumen estadístico y pivotaje de la tabla para mejorar la presentación. 
        statsDF <- statsDF %>% dplyr::select(-c(date, symbol)) %>% describe() %>%
            rownames_to_column(var = "variable") %>% dplyr::select(-c(vars, n)) %>%
            pivot_longer(-variable) %>% pivot_wider(names_from = variable, values_from = value)
        
        # Cambio de los nombres de las variables a castellano. 
        colnames(statsDF) <- c(" ", "Apertura", "Máximo", "Mínimo", "Cierre",
                               "Volumen", "USDT", "Transacciones")
        
        # Se redondean los resultados para mejorar la claridad de la presentación. 
        statsDF %>% mutate_if(is.numeric, round, 3)
        
    }, rownames = FALSE,
    options = list(scrollX = TRUE, searching = FALSE, dom = "ti",
                   initComplete = JS("function(settings, json) {", 
                                     "$(this.api().table().header()).css({
                                        'background-color': '#3B4A52',
                                        'margin':'100px',
                                        'text-align': 'center',
                                        'font-family': 'Cambria',
                                        'color': '#B4C7D8'});",
                                     "}"),
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))
    
    
    )
    
    output$SkimF <-  metaRender(renderPrint, ({
        
        # Resumen skim del marco de datos sin histogramas (en ocasiones hay glitches y los muestra igual). 
        cryptoDF$data  %>% skim_without_charts()
        
    }))
    
    
    # Salida de la serie temporal de la ventana seleccionada.  
    output$closeLineID <- renderPlot({
        
        ggplot(data = cryptoDF$data,
               aes(x = date, y = close)) +
            geom_line(size = 0.6) +
            labs(title = "Valores de Cierre de las Criptotransacciones",
                 subtitle = paste0("Fecha: ", input$FromDateTime, " / ", input$TillDateTime),
                 y = "", x = " ") +
            stat_smooth(aes(col = symbol), method = 'gam') +
            scale_color_poke(pokemon = "Totodile") +
            facet_wrap(symbol~.,  scales = "free") +
            theme(legend.position = "none")
        
    })
    
    # Salida del histograma del número de criptotransacciones de la ventana seleccionada.  
    output$transactionsHistogramID <- renderPlot({
        
        ggplot(data = cryptoDF$data, 
               aes(x = date, y = tradecount, fill = symbol)) +
            geom_bar(stat = "identity") +
            labs(title = "Número de Criptotransacciones",
                 subtitle = paste0("Fecha: ", input$FromDateTime, " / ", input$TillDateTime),
                 y = "", x = " ") +
            scale_fill_poke(pokemon = "Totodile") +
            facet_wrap(symbol~.,  scales = "free", ncol = 3) +
            theme(legend.position = "none")
        
    })
    
    # Salida de la distribución del volumen de las criptomonedas de la ventana seleccionada. 
    output$volumeDensityID <- renderPlot({
        
        ggplot(data = cryptoDF$data) +
            geom_area(aes(x = date, y = volume, fill = symbol), alpha = 0.5) +
            labs(title = "Distribución del Volumen de las Criptotransacciones",
                 subtitle = paste0("Fecha: ", input$FromDateTime, " / ", input$TillDateTime),
                 y = "", x = " ", fill = " ") + 
            scale_y_log10() +
            scale_fill_poke(pokemon = "Totodile") +
            theme(legend.position = "top")
        
    }, height = 700)
    
    
    # ----- Detección de Valores Atípicos ----- 
    # Variable reactiva que detecta los valores atípicos para una variable dada.
    outlierDetection <- reactive({
        
        # La sintaxis de boxplot exige que sea un dataframe, no otro objeto. 
        detectionDF <- cryptoDF$data %>% as.data.frame()
        
        # Obtención de los valores de los índices donde se hallan valores atípicos. 
        index <- boxplot(detectionDF[, input$selectVarOutlier], plot = FALSE)$out
        
        return(index)
    })
    
    # Variable reactiva que contiene el subconjunto de los valores atípicos de la variable seleccionada.
    outlierDFR <- reactive({
        
        # La función which también solicita que se le pase un dataframe. 
        outliersDF <- cryptoDF$data %>% as.data.frame()
        
        # Obtención del subconjunto donde se encuentran valores atípicos.  
        outliersDF <- outliersDF[which(outliersDF[, input$selectVarOutlier] %in% outlierDetection()),]
        
        # Cambio del nombre de las variables a castellano para mejorar la presentación. 
        colnames(outliersDF) <- c("Fecha", "Cripto", "Apertura",
                                  "Máximo", "Mínimo", "Cierre",
                                  "Volumen", "USDT", "Transacciones")
        
        return(outliersDF)
        
        
    })
    
    # Salida de la variable reactiva que contiene el marco de datos con los valores atípicos. 
    output$outlierDFID <- DT::renderDataTable({
        
        outliersDF <- outlierDFR()
        
        outliersDF
        
    }, rownames = FALSE, options = list(initComplete = JS("function(settings, json) {", 
                                                          "$(this.api().table().header()).css({
                                        'background-color': '#3B4A52',
                                        'margin':'100px',
                                        'text-align': 'center',
                                        'font-family': 'Cambria',
                                        'color': '#B4C7D8'});",
                                                          "}"),
                                        searching = TRUE,  searchHighlight = TRUE,
                                        displayStart = 1, autoWidth = FALSE, drop = FALSE,
                                        orderClasses = TRUE, pageLength = 20, dom = "tplf", 
                                        lengthMenu = c(5, 10, 15, 20),
                                        scrollY = TRUE, scrollX = TRUE)
    )
    
    
    
    # Acción de guardado del marco de datos de los valores atípicos. 
    output$downloadOutlierDF <- downloadHandler(
        
        filename = function(){
            
            # Nombre del archivo. 
            paste0("ValoresAtípicos", input$selectVarOutlier, as.numeric(as.Date(input$FromDateTime)), as.numeric(as.Date(input$TillDateTime)), ".xlsx")
            
        },
        
        content = function(file) {
            
            # Guardar el archivo como xlsx. 
            write.xlsx(outlierDFR(), file, row.names = FALSE)
            
        }
    )
    
    # Salida de los diagramas de bigotes de la variable elegida.  
    output$boxplotsOutlierID <- renderPlot({
        
        # Test T Paramétrico. 
        p1 <- ggbetweenstats(
            data = cryptoDF$data, outlier.tagging = TRUE,
            # Se pone la doble exclamación para que entienda que la cadena de caracteres de la entrada es una variable (paquete rlang). 
            x = symbol, y = !!input$selectVarOutlier,
            xlab = "Criptodivisa", ylab = "Transacciones",
            plot.type = "box",
            type = "p", conf.level = 0.99,
            title = "Test Paramétrico",
            package = "palettetown",
            palette = "charizard"
        )
        
        # Test U de Mann-Whitney.
        p2 <- ggbetweenstats(
            data = cryptoDF$data, outlier.tagging = TRUE,
            x = symbol, y = !!input$selectVarOutlier,
            xlab = "Criptodivisa", ylab = "Transacciones",
            plot.type = "violin",
            type = "np", conf.level = 0.99,
            title = "Test No Paramétrico",
            package = "palettetown",
            palette = "bulbasaur"
        )
        
        # Test T Robusto. 
        p3 <- ggbetweenstats(
            data = cryptoDF$data, outlier.tagging = TRUE,
            x = symbol, y = !!input$selectVarOutlier, 
            xlab = "Criptodivisa", ylab = "Transacciones",
            plot.type = "boxviolin",
            type = "r", conf.level = 0.99,
            title = "Test Robusto",
            tr = 0.005, k = 3,
            package = "palettetown",
            palette = "porygon"
        )
        
        # Test T Paramétrico con factor de Bayes
        p4 <- ggbetweenstats(
            data = cryptoDF$data, outlier.tagging = TRUE,
            x = symbol, y = !!input$selectVarOutlier,
            xlab = "Criptodivisa", ylab = "Transacciones",
            type = "bayes",
            plot.type = "box",
            title = "Test Bayesiano",
            package = "palettetown",
            palette = "surskit"
        )
        
        # Combinar los cuatro gráficos en un solo marco.
        combine_plots(list(p1, p2, p3, p4),
                      plotgrid.args = list(nrow = 2),
                      annotation.args = list(
                          title = "Comparación del Número de Transacciones entre Criptodivisas",
                          caption = "Fuente: Binance")
        )
        
    }, height = 900)
    
    # Salida del histograma de la variable seleccionada. Compara la media del conjunto con la de cada criptomoneda. 
    output$histogramOutlierID <- renderPlot({
        
        grouped_gghistostats(
            data = cryptoDF$data,
            x = !!input$selectVarOutlier, xlab = "Número de Transacciones", 
            grouping.var = symbol, normal.curve = FALSE,
            type = "robust",
            test.value = mean(cryptoDF$data$tradecount), 
            centrality.line.args = list(color = "#D55E00", linetype = "dashed"),
            ggtheme = ggthemes::theme_pander(),
            annotation.args = list(
                title = "Distribución del Número de Transacciones entre Criptodivisas",
                caption = "Fuente: Binance",
                tag_levels = "I"
            ), plotgrid.args = list(ncol = 2)
        )
    }, height = 700)
    
    # ----- Estudio de la correlación entre las variables -----
    # Variable reactiva con la matriz de correlaciones. 
    corrDFR <- reactive({
        
        # Se extrae un marco de datos con información estadística de la matriz de correlaciones. 
        corrDF <- grouped_ggcorrmat(
            data = cryptoDF$data, 
            grouping.var = symbol,
            type = "robust",
            # digits = 4, No se ve redondeado en la aplicación Shiny. Es un glitch rarillo... 
            output = "dataframe")
        
        return(corrDF)
        
    })
    
    # Salida de la tabla que recoge la matriz de correlaciones. 
    output$corrDFID <- DT::renderDataTable({
        
        corrDF <- corrDFR()
        
        # Se redondean todos los valores para mejorar la claridad de la exposición. 
        corrDF %>% mutate_if(is.numeric, round, 4)
        
    }, rownames = FALSE, options = list(initComplete = JS("function(settings, json) {", 
                                                          "$(this.api().table().header()).css({
                                        'background-color': '#3B4A52',
                                        'margin':'100px',
                                        'text-align': 'center',
                                        'font-family': 'Cambria',
                                        'color': '#B4C7D8'});",
                                                          "}"),
                                        searching = TRUE,  searchHighlight = TRUE,
                                        displayStart = 1, autoWidth = FALSE, drop = FALSE,
                                        orderClasses = TRUE, pageLength = 10, dom = "tplf", 
                                        lengthMenu = c(5, 10, 15, 20, 25),
                                        scrollY = TRUE, scrollX = TRUE)
    )
    
    
    # Guardado de la matriz de correlaciones. 
    output$downloadCorrDF <- downloadHandler(
        
        filename = function(){
            
            # Nombre del archivo. 
            paste0("MatrizCorrelaciones", as.numeric(as.Date(input$FromDateTime)), as.numeric(as.Date(input$TillDateTime)), ".xlsx")
            
        },
        
        content = function(file) {
            
            corrDF <- corrDFR() %>% as.data.frame()
            
            corrDF <- corrDF %>% mutate_if(is.numeric, round, 4)
            
            # Guardar el archivo como xlsx. 
            write.xlsx(corrDF, file, row.names = FALSE)
            
        }
    )
    
    
    # Salida de los diagramas de correlaciones. 
    output$corrPlotID <- renderPlot({
        
        # Se aplica una separación en base a la criptodivisa. 
        grouped_ggcorrmat(
            data = cryptoDF$data,
            type = "robust", 
            grouping.var = symbol,
            plotgrid.args = list(ncol = 2),
            colors = c("#B4C7D8", "white", "#7A578E"),
            annotation.args = list(
                tag_levels = "A",
                title = " ", # Diagrama de Correlación
                caption = "Fuente: Binance"
            ), 
            
            ggcorrplot.args = list(outline.color = "black", hc.order = TRUE)
            
        )
    }, height = 850)
    
    
    # Creación de un marco de datos que añada las variables de bitcoin y bitcoin con salto temporal de menos treinta minutos. 
    kryptoBTC <- reactive({
        
        bitcoinLag30 <- lag(bitcoin, 30)
        
        bitcoinLag30$date <- bitcoinLag30$date-(30*60)
        
        bitcoinLag30 <- bitcoinLag30 %>% dplyr::select(date, close) %>% rename(closeBTC_lag30 = close)
        
        bitcoinNew <- merge(bitcoin, bitcoinLag30, by = 'date')
        
        kryptoBTC <- merge(cryptoDF$data, bitcoinNew %>% dplyr::select(-unix), by = 'date', suffixes = c("","_BTC")) 
        
        return(kryptoBTC)
        
    })
    
    # Salida de los mapas de calor obtenidos mediante D3Heatmap. 
    output$heatPlotID <- renderD3heatmap({
        
        kryptoBTCRenamed <- kryptoBTC()
        
        # Se castellaniza el nombre de las variables para mejorar la calidad de la presentación. 
        colnames(kryptoBTCRenamed) <- c("Fecha", "Cripto", "Apertura",
                                        "Máximo", "Mínimo", "Cierre",
                                        "Volumen", "USDT", "Transacciones",
                                        "Cripto_BTC", "Apertura_BTC",
                                        "Máximo_BTC", "Mínimo_BTC", "Cierre_BTC",
                                        "Volumen_BTC", "USDT_BTC", 
                                        "Transacciones_BTC", "CierreLag30_BTC")
        
        d3heatmap::d3heatmap(cor(kryptoBTCRenamed %>% select_if(is.numeric)), dendrogram = 'none')
        
    })
    
    
    # ----- Aplicación de los Modelos de Regresión -----
    
    
    # Variables reactivas con los valores de la estandarización. 
    # Servirán al propósito de reescalar las predicciones de los LSTM. 
    # Variable reactiva que toma el atributo centro de la estandarización. 
    centerSt <- reactive({
        
        if(dataScaled() == FALSE){
            
            scaledDF <- cryptoDF$data %>% select_if(is.numeric) %>% scale()
            
        } else {
            
            scaledDF <- cryptoDF$data %>% dplyr::mutate_at(c("open", "high", "low", "close",
                                                             "tradecount", "volume", "USDT"), unscale)
            
            scaledDF <- scaledDF %>% select_if(is.numeric) %>% scale()
            
        }
        
        return(attr(scaledDF, 'scaled:center')[4])
        
    })
    
    # Variable reactiva que extrae el atributo escala de la estandarización.
    scaleSt <- reactive({
        
        if(dataScaled() == FALSE){
            
            scaledDF <- cryptoDF$data %>% select_if(is.numeric) %>% scale()
            
        } else {
            
            scaledDF <- cryptoDF$data %>% dplyr::mutate_at(c("open", "high", "low", "close",
                                                             "tradecount", "volume", "USDT"), unscale)
            
            scaledDF <- scaledDF %>% select_if(is.numeric) %>% scale()
            
        }
        
        return(attr(scaledDF, 'scaled:scale')[4])
        
    })
    
    # Función para reescalar el marco de datos.  
    unscale2 <- function(x){
        
        return(x *  scaleSt() + centerSt())
        
    }
    
    
    # Variable reactiva sobre la que se efectúan las predicciones en base al modelo seleccionado. 
    # Se crea esta nueva variable para jugar con la estandarización sin que afecte a otros componentes de la aplicación. 
    predictionDF <- reactive({
        
        if(dataScaled() == TRUE){
            
            predictionDF <- kryptoBTC()
            
            # Siempre va a estar sin estandarizar de base, haga lo que haga el usuario con el apartado de análisis.  
            predictionDF <- predictionDF %>% dplyr::mutate_at(c("open", "high", "low", "close",
                                                                "tradecount", "volume", "USDT"),
                                                              unscale2)
        } else {
            
            predictionDF <- kryptoBTC()
            
        }
        
        return(predictionDF)
        
    })
    
    
    # Carga del modelo de regresión seleccionado.
    modelCaret <- reactive({
        
        # Se utiliza la función de carga de modelos con extensión rds. 
        model <- readRDS(input$selectionCaret)
        
        return(model)
        
    })
    
    # Resumen del modelo de regresión cargado. 
    output$modelSummary <- renderPrint({
        
        summary(modelCaret())
        
    })
    
    # Variables reactivas con las métricas de los modelos.
    R2Caret <- reactiveVal(0)
    
    BR2Caret <- reactiveVal(0)
    
    MAECaret <- reactiveVal(0)
    
    MSECaret <- reactiveVal(0)
    
    
    # Salidas de las métricas en las cajas de valores de la interfaz. 
    output$R2CaretID <- renderValueBox({
        
        valueBox(subtitle = "R2", value = R2Caret(), color = "light-blue")
        
    })
    
    output$BR2CaretID <- renderValueBox({
        
        valueBox(subtitle = "BR2", value = BR2Caret(), color = "blue")
        
    })
    
    output$MAECaretID <- renderValueBox({
        
        valueBox(subtitle = "MAE", value = MAECaret(), color = "navy")
        
    })
    
    output$MSECaretID <- renderValueBox({
        
        valueBox(subtitle = "MSE", value = MSECaret(), color = "purple")
        
    })
    
    # Variable reactiva con el marco de datos de las predicciones. 
    predictionsCaretDF <- reactive({
        
        kryptoCaret <- predictionDF() 
        
        # Se seleccionan las variables de la fecha, la criptodivisa y el valor de cierre.
        # Además, se crea una variable de predicciones inicializada en un vector de ceros. 
        kryptoCaret <- kryptoCaret %>% dplyr::select(date, symbol, close) %>% add_column(predictions = 0)
        
        return(kryptoCaret)
        
    })
    
    # Se crea una variable de valores reactivos para el evento de predicción. 
    predictionsCaret <- reactiveValues(data = NULL)
    
    # Se guarda en la variable de valores reactivos el marco resultante. 
    observe({predictionsCaret$data <- predictionsCaretDF()})
    
    
    # Salida del marco de datos con las predicciones. 
    output$predictionsCaretDF <- DT::renderDataTable({
        
        predictions <- predictionsCaret$data
        
        colnames(predictions) <- c("Fecha", "Criptodivisa", "Valor Real del Cierre", "Predicción de Cierre")
        
        predictions %>% mutate_if(is.numeric, round, 4)
        
    }, rownames = FALSE, options = list(initComplete = JS("function(settings, json) {", 
                                                          "$(this.api().table().header()).css({
                                        'background-color': '#3B4A52',
                                        'margin':'100px',
                                        'text-align': 'center',
                                        'font-family': 'Cambria',
                                        'color': '#B4C7D8'});",
                                                          "}"),
                                        searching = TRUE,  searchHighlight = TRUE,
                                        displayStart = 1, autoWidth = FALSE, drop = FALSE,
                                        orderClasses = TRUE, pageLength = 20, dom = "tplf", 
                                        lengthMenu = c(5, 10, 15, 20, 25, 30),
                                        scrollY = TRUE, scrollX = TRUE)
    )
    
    
    # Guardado del marco de datos con las predicciones obtenidas del modelo de regresión. 
    output$downloadCaretDF <- downloadHandler(
        
        filename = function(){
            
            # Nombre del archivo. 
            paste0("PrediccionesModeloLineal", as.numeric(as.Date(input$FromDateTime)), as.numeric(as.Date(input$TillDateTime)), ".xlsx")
            
        },
        
        content = function(file) {
            
            # Guardar el archivo como xlsx. 
            write.xlsx(predictionsCaret$data, file, row.names = FALSE)
            
        }
    )
    
    # Gráfica que recoge la serie temporal real frente a la predicción realizada. 
    output$predictionsPlotCaretID <- renderPlot({
        
        ggplot(data = predictionsCaret$data,
               aes(x = date, y = close, col = symbol)) +
            geom_line(size = 0.7) +
            geom_line(aes(y = predictions), size = 0.7, linetype = "dashed", color = "#750025") +
            labs(title = "",
                 subtitle = "",
                 y = "", x = "",
                 caption = paste0("Fecha: ", input$FromDateTime, " / ", input$TillDateTime)) +
            scale_colour_brewer(palette = "Dark2") +
            facet_wrap(symbol~.,  scales = "free") +
            theme(legend.position = "none")
        
    })
    
    
    # Acción de efectuar la predicción. 
    observeEvent(input$predictButtonCaret, {
        
        # Alerta para la recomendación de que solamente se apliquen los modelos 
        # para efectuar las predicciones de una criptodivisa. 
        if(length(unique(as.character(predictionsKeras$data$symbol))) != 1){
            
            shinyalert("Aviso", "¡Es mejor elegir solamente una criptodivisa!", type = "info")
            
        }
        
        # Se realizan las predicciones en en virtud del modelo elegido y se guardan en el marco de datos. 
        predictionsCaret$data$predictions <- predict(modelCaret(), predictionDF())
        
        # Se recalculan las métricas que se verán en los resultados.  
        R2Caret(round(caret::R2(predictionsCaret$data$predictions, predictionsCaretDF()$close), digits = 2))
        BR2Caret(round(hydroGOF::br2(predictionsCaret$data$predictions, predictionsCaretDF()$close), digits = 2))
        MAECaret(round(hydroGOF::mae(predictionsCaret$data$predictions, predictionsCaretDF()$close), digits = 2))
        MSECaret(round(hydroGOF::mse(predictionsCaret$data$predictions, predictionsCaretDF()$close), digits = 2))
        
        
    })
    
    
    # ----- Aplicación de los modelos de Aprendizaje Profundo -----
    
    # Cargar el modelo seleccionado que ha sido creado con Keras
    modelKeras <- reactive({
        
        # Se utiliza la función de carga de keras de modelos con extensión h5. 
        modelKeras <- load_model_hdf5(input$selectionKeras)
        
    }) 
    
    # Resumen del modelo cargado. 
    output$modelKerasSummary <- renderPrint({
        
        summary(modelKeras())
        
    })
    
    # Variables reactivas con las métricas de los modelos.
    R2Keras <- reactiveVal(0)
    
    BR2Keras <- reactiveVal(0)
    
    MAEKeras <- reactiveVal(0)
    
    MSEKeras <- reactiveVal(0)
    
    # Salidas de las métricas de Keras en las cajas de valores de la interfaz. 
    output$R2KerasID <- renderValueBox({
        
        valueBox(subtitle = "R2", value = R2Keras(), color = "light-blue")
        
    })
    
    output$BR2KerasID <- renderValueBox({
        
        valueBox(subtitle = "BR2", value = BR2Keras(), color = "blue")
        
    })
    
    output$MAEKerasID <- renderValueBox({
        
        valueBox(subtitle = "MAE", value = MAEKeras(), color = "navy")
        
    })
    
    output$MSEKerasID <- renderValueBox({
        
        valueBox(subtitle = "MSE", value = MSEKeras(), color = "purple")
        
    })
    
    # Variable reactiva con el marco de datos de las predicciones para los modelos de Keras. 
    predictionsKerasDF <- reactive({
        
        kryptoKeras <- kryptoBTC() 
        
        # Se seleccionan las variables de la fecha, la criptodivisa y el valor de cierre.
        # También se crea una variable de predicciones inicializada en un vector de ceros.
        kryptoKeras <- kryptoKeras %>% dplyr::select(date, symbol, close) %>% add_column(predictions = 0)
        
        return(kryptoKeras)
        
    })
    
    # Se crea una variable de valores reactivos para el evento de predicción. 
    predictionsKeras <- reactiveValues(data = NULL)
    
    # Se guarda en la variable de valores reactivos el marco resultante. 
    observe({predictionsKeras$data <- predictionsKerasDF()})
    
    # Salida de la tabla de resultados de las predicciones realizadas. 
    output$predictionsKerasDF <- DT::renderDataTable({
        
        predictions <- predictionsKeras$data
        
        colnames(predictions) <- c("Fecha", "Criptodivisa", "Valor Real del Cierre", "Predicción de Cierre")
        
        predictions %>% mutate_if(is.numeric, round, 4)
        
    }, rownames = FALSE, options = list(initComplete = JS("function(settings, json) {", 
                                                          "$(this.api().table().header()).css({
                                        'background-color': '#3B4A52',
                                        'margin':'100px',
                                        'text-align': 'center',
                                        'font-family': 'Cambria',
                                        'color': '#B4C7D8'});",
                                                          "}"),
                                        searching = TRUE,  searchHighlight = TRUE,
                                        displayStart = 1, autoWidth = FALSE, drop = FALSE,
                                        orderClasses = TRUE, pageLength = 20, dom = "tplf", 
                                        lengthMenu = c(5, 10, 15, 20, 25, 30),
                                        scrollY = TRUE, scrollX = TRUE)
    )
    
    # Guardado del marco de datos con las predicciones obtenidas del modelo LSTM. 
    output$downloadCaretDF <- downloadHandler(
        
        filename = function(){
            
            # Nombre del archivo. 
            paste0("PrediccionesLSTM", as.numeric(as.Date(input$FromDateTime)), as.numeric(as.Date(input$TillDateTime)), ".xlsx")
            
        },
        
        content = function(file) {
            
            # Guardar el archivo como xlsx. 
            write.xlsx(predictionsKeras$data, file, row.names = FALSE)
            
        }
    )
    
    # Gráfica con la serie temporal real frente a la predicción obtenida por los LSTM. 
    output$predictionsPlotKerasID <- renderPlot({
        
        ggplot(data = predictionsKeras$data,
               aes(x = date, y = close, col = symbol)) +
            geom_line(size = 0.7) +
            geom_line(aes(y = predictions), size = 0.7, linetype = "dashed", color = "#750025") +
            labs(title = "",
                 subtitle = "",
                 y = "", x = " ",
                 caption = paste0("Fecha: ", input$FromDateTime, " / ", input$TillDateTime)) +
            scale_colour_brewer(palette = "Dark2") +
            facet_wrap(symbol~.,  scales = "free") +
            theme(legend.position = "none")
        
    })
    
    # Variables reactivas que sirven de entrada para los modelos. 
    # Se ha de transformar el marco de entrada en un array. 
    # Antes de la creación del array se de ha de ejecutar la estandarización. 
    arrayX1 <- reactive({
        
        test  <- predictionDF()
        
        test <- test %>% mutate_if(is.numeric, scale)
        
        arrayX <- array(data = lag(cbind(test$close,
                                         test$tradecount,
                                         test$volume),
                                   datalags)[-(1:datalags), ],
                        dim = c(nrow(test), datalags, 3))
        
        return(arrayX)
        
    })
    
    # Hay modelos que toman las variables alto y bajo en lugar del volumen y las transacciones. 
    # Esto exige la creación de otra variable reactiva que devuelva una estructura de array apropiada. 
    arrayX2 <- reactive({
        
        test  <- predictionDF()
        
        test <- test %>% mutate_if(is.numeric, scale)
        
        arrayX <- array(data = lag(cbind(test$close,
                                         test$high,
                                         test$low),
                                   datalags)[-(1:datalags), ],
                        dim = c(nrow(test), datalags, 3))
        
        return(arrayX)
        
    })
    
    
    # Acción de efectuar la predicción. 
    observeEvent(input$predictButtonKeras, {
        
        # Advertencia emergente para que el usuario seleccione sólo una criptodivisa. 
        if(length(unique(as.character(predictionsKeras$data$symbol))) != 1){
            
            shinyalert("Aviso", "¡Es mejor elegir solamente una criptodivisa!", type = "info")
            
        }
        
        # Condición para la elección del array de entrada en virtud del modelo seleccionado. 
        if(input$selectionKeras %in% c("models/LSTMSolana.h5", "models/LSTMCardano.h5")){
            
            predictionsKeras$data$predictions <- modelKeras() %>% predict(arrayX1(), batch_size = 1) %>% .[,1]
            predictionsKeras$data$predictions <- unscale2(predictionsKeras$data$predictions)
            
        } else{
            
            predictionsKeras$data$predictions <- modelKeras() %>% predict(arrayX2(), batch_size = 1) %>% .[,1]
            predictionsKeras$data$predictions <- unscale2(predictionsKeras$data$predictions)
            
        }
        
        # Cambio de las métricas en función de las predicciones. 
        R2Keras(round(caret::R2(predictionsKeras$data$predictions, predictionsKerasDF()$close), digits = 2))
        BR2Keras(round(hydroGOF::br2(predictionsKeras$data$predictions, predictionsKerasDF()$close), digits = 2))
        MAEKeras(round(hydroGOF::mae(predictionsKeras$data$predictions, predictionsKerasDF()$close), digits = 2))
        MSEKeras(round(hydroGOF::mse(predictionsKeras$data$predictions, predictionsKerasDF()$close), digits = 2))
        
        
    })
    
}
