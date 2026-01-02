  # Librerías requeridas por la aplicación web interactiva. 
  library(shiny)
  library(data.table)
  library(xgboost)
  library(ranger)
  library(tools)
  
  
  # Carga de los datos anteriores. 
  load("models_with_genes.RData")
  
  # Se genera una longitud para cada gen, lo cual es necesario para cada gen.
  gene_length <- rep(1000, length(selected_genes))
  names(gene_length) <- selected_genes
  
  # Tamaño máximo permitido para los archivos subidos.
  options(shiny.maxRequestSize = 200*1024^2)
  
  # Interfaz de usuario.
  ui <- fluidPage(
    titlePanel("Clasificación jerárquica de muestras normales y tumorales"),
    
    sidebarLayout(
      sidebarPanel(
        fileInput("file", "Subir muestras", 
                  accept = c(".csv", ".txt", ".tsv", ".gz")),
        downloadButton("downloadExample", "Descargar archivo de ejemplo"),
        actionButton("predict", "Clasificar")),
      
      mainPanel(
        tabsetPanel(
          # Panel de la clasificación de las muestas. 
          tabPanel("Clasificación",
                   h3("Resultados de clasificación"),
                   tableOutput("resultsTable"),
                   downloadButton("downloadResults", "Descargar resultados")),
          
          # Panel de ayuda al usuario. 
          tabPanel("Requisitos / Ayudas al usuario",
                   h3("Formato del archivo"),
                   p("La aplicación acepta archivos que se encuentren en formato CSV, TXT, TSV o GZ."),
                   p("Cada fila del documento debe corresponder a una muestra, y cada columna a un gen."),
                   p("Los valores deben de ser numéricos, exceptuando la columna de los identificadores."),
                   
                   h3("Identificador de muestra"),
                   p("El archivo puede contener una columna denominada 'Sample' o 'ID' que contenga los identificadores de muestra."),
                   p("En caso de no existir esta columna, la aplicación asigna nombre de manera automática a cada muestra."),
                   
                   h3("Genes requeridos"),
                   p("El archivo debe contener todas las columnas correspondientes a los genes utilizados por el modelo. Si falta alguno, la clasificación no podrá realizarse."),
                   p("El nombre de los genes debe de aparecer escrito tal y como se muestra en el siguiente listado."),
                   verbatimTextOutput("geneList"),
                   p("Los nombres de los genes deben seguir la nomenclatura oficial del", 
                     a("HUGO Gene Nomenclature Comittee (HGNC)",
                       href = "https://www.genenames.org/",
                       target = "_blank")),
                   
                   h3("Número de muestras"),
                   p("No hay un límite estricto de muestras, pero archivos muy grandes pueden realizar el proceso."),
                   
                   h3("Escala de los datos"),
                   p("La aplicación acepta datos crudos o datos normalizados por TPM."),
                   p("En caso de que los datos no estén en TPM, la aplicación los normaliza automáticamente."),
                   
                   h3("Valores esperados"),
                   p("No debe haber valores negativos."),
                   p("Valores extremadamente altos pueden indicar un error en el archivo."),
                   
                   h3("Archivo de ejemplo"),
                   p("Se puede descargar un archivo de ejemplo desde el panel lateral para ver el formato esperado de fichero.")),
          
          # Panel con especificaciones tecnicas. 
          tabPanel("Especificaciones técnicas",
                   h3("Modelos utilizados"),
                   p("Random Forest entrenado para la clasificación de muestras normales y tumorales."),
                   p("XGBoost entrenado para la clasificación de muestras tumorales respecto al tipo de cáncer del que se trata."),
                   
                   h3("Métricas de evaluación del Random Forest"),
                   p("Precisión muestras normales: 0.9841"),
                   p("Precisión muestras tumorales: 0.9968"),
                   p("Precisión global: 0.9905"),
                   p("Recall muestras normales: 0.9764"),
                   p("Recall global: 0.9871"),
                   p("F1-score muestras normales: 0.9802"),
                   p("F1-score global: 0.9888"),
                   p("Kappa: 0.9775"),
                   p("Sensibilidad: 0.9978"),
                   p("Especificidad: 0.9764"),
                   p("Balanced Accuracy: 0.9871"),
                   
                   h3("Métricas de evaluación del XGBoost"),
                   p("Precisión: 0.9779"),
                   p("Recall: 0.9788"),
                   p("Especificidad: 0.9957"),
                   p("F1 macro: 0.9782"),
                   p("Accuracy: 0.9785"),
                   p("Kappa: 0.9742"),
                   p("Balanced Accuracy: 0.9873"),
                   shiny::h4("Valores AUC-ROC para cada tipo de cáncer"),
                   p("BRCA: 0.9999"),
                   p("KIRC: 1"),
                   p("LUAD: 0.9975"),
                   p("LUSC: 0.9973"),
                   p("PRAD: 1"),
                   p("THCA: 1"),
                   p("Global: 0.9991"),
                   
                   h3("Origen de los datos"),
                   p("Los datos provienen del estudio 'Alternatively processed and compiled RNA-Sequencing and clinical data for thousands of samples from The Cancer Genome Atlas', que se puede encontrar en NCBI - Gene Expression Omnibus mediante el identificador ",
                     a("GSE62944", 
                       href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944",
                       target = "_blank"),"."),
                   p("Este estudio, como su nombre indica, trabaja con algunas muestras obtenidas a su vez del The Cancer Genome Atlas, disponible a través de la web del National Cancer Institute - NIH."))
        )
      )
    )
  )
  
  # Función para comprobar si los datos se encuentran en formato TPM.
  is_tpm <- function(data, tolerance = 0.20) {
    rs_med <- median(rowSums(data, na.rm = TRUE), na.rm = TRUE)
    pct_diff <- abs(rs_med - 1e6) / 1e6
    return(pct_diff < tolerance)
  }
  
  # Función que aplica la normlización TPM. 
  tpm_normalize <- function(counts, gene_length) {
    rpk <- counts / (gene_length / 1000)
    tpm <- t(t(rpk) / colSums(rpk)) * 1e6
    return(tpm)
  }
  
  # Función para comprobar que los datos se encuentran en formato TPM.
  ensure_tpm <- function(data, gene_length, tolerance = 0.20) {
    if (is_tpm(data, tolerance)) {
      return(data)
    } else {
      return(tpm_normalize(data, gene_length))
    }
  }
  
  # Servidor. 
  server <- function(input, output) {
    
    # Mostrar la lista de los genes requeridos en la ayuda al usuario. 
    output$geneList <- renderPrint({
      selected_genes
    })
    
    # Descarga del archivo de ejemplo. 
    output$downloadExample <- downloadHandler(filename = function() { 
      "test_samples.csv" 
    },
    content = function(file) {
      file.copy("test_samples.csv", file)
    })
    
    # Evento de clasificación. 
    observeEvent(input$predict, {
      req(input$file)
      
      # Diferentes formas de lectura del arcgivo según de su extensión. 
      ext <- tools::file_ext(input$file$name)
      if (ext == "csv") {
        samples <- fread(input$file$datapath)
      } else if (ext %in% c("txt", "tsv", "gz")) {
        samples <- fread(input$file$datapath, sep = "\t")
      } else {
        showNotification("Formato de archivo no soportado. Usar CSV, TXT o TSV.", 
                         type = "error")
        return(NULL)
      }
      
      # Reconocimiento de identificadores de muestra, si es que los hay. 
      if ("Sample" %in% colnames(samples)) {
        sample_names <- samples$Sample
        samples[, Sample := NULL]
      } else if ("ID" %in% colnames(samples)) {
        sample_names <- samples$ID
        samples[, ID := NULL]
      } else {
        sample_names <- paste0("Sample ", seq_len(nrow(samples)))
      }
      
      # Comprobación de que los valores son numéricos. 
      if (!all(sapply(samples, is.numeric))) {
        showNotification("El archivo contiene valores no numéricos.", 
                         type = "error")
        return(NULL)
      }
      
      # Comprobación de que todos los genes requeridos se encuentran en los datos. 
      missing_genes <- setdiff(selected_genes, colnames(samples))
      if (length(missing_genes) > 0) {
        showNotification(
          paste("El archivo no contiene todos los genes necesarios. Faltan: ", 
                paste(missing_genes, collapse = ", ")), type = "error", 
          duration = NULL)
        return(NULL)
      }
      
      # Normalización TPM de los datos (en caso que sea necesario). 
      samples_mat <- as.matrix(samples)
      samples_mat <- ensure_tpm(samples_mat, gene_length)
      samples <- as.data.table(samples_mat)
      
      # Transformación log2(x + 1).
      samples <- log2(samples + 1)
      
      # Comprobación de si el archivo contiene valores negativos o demasiado
      # elevados.
      if (any(samples < 0, na.rm = TRUE)) {
        showNotification("El archivo contiene valores negativos.")
      }
      
      if (any(samples > 50, na.rm = TRUE)) {
        showNotification("El archivo contiene valores extremadamente altos.")
      }
      
      # Filtrado de los genes para mantener solo los genes necesitados por los
      # modelos. 
      filtered_samples <- samples[, ..selected_genes]
      filtered_samples <- as.matrix(filtered_samples)
      
      # Predicciones del modelo Random Forest para clasificar según el tipo de 
      # muestra. 
      rf_pred <- predict(rf_final, data.frame(filtered_samples))$predictions
      rf_pred <- as.matrix(rf_pred)
      
      if (ncol(rf_pred) == 1) {
        rf_pred <- cbind(1 - rf_pred, rf_pred)
      }
      
      rf_prob <- rf_pred[, 2]
      rf_classif <- ifelse(rf_prob > 0.5, "Tumoral", "Normal")
      
      prob_rf <- ifelse(rf_classif == "Tumoral", rf_pred[, 2], rf_pred[, 1])
      
      # Inicialización de los resultados del modelo XGBoost. 
      xgb_classif <- rep("-", nrow(filtered_samples))
      cancer_prob <- rep(NA, nrow(filtered_samples))
      
      # Obtención de las muestras tumorales predichas por el Random Forest. 
      tumor_idx <- which(rf_classif == "Tumoral")
      
      # En caso de que haya alguna, se realiza la predicción sobre estas muestras. 
      if (length(tumor_idx) > 0) {
        prob <- predict(xgb_final, newdata = filtered_samples[tumor_idx, ])
        prob <- matrix(prob, ncol = length(levels_cancer), byrow = TRUE)
        colnames(prob) <- levels_cancer
        
        pred_classif <- apply(prob, 1, function(x) colnames(prob)[which.max(x)])
        xgb_classif[tumor_idx] <- pred_classif
        
        cancer_prob[tumor_idx] <- apply(prob, 1, max)
      }
      
      # Generación de un dataframe con los identificadores de muestra, la 
      # clasificación según el tipo de muestra y tipo de cáncer y la probabilidad
      # de la clase obtenida en cada caso. 
      results <- data.frame(Muestra = sample_names, 
                            "Tipo de Muestra" = rf_classif, 
                            "Prob. de Muestra" = round(prob_rf, 4),
                            "Tipo de Cáncer" = xgb_classif,
                            "Prob. de Cáncer" = round(cancer_prob, 4),
                            stringsAsFactors = FALSE, check.names = FALSE)
      
      # Renderizado del dataframe de resultados para la obtención de la tabla 
      # en la pantalla de la aplicación. 
      output$resultsTable <- renderTable({ results }, rownames = FALSE)
      
      # Descarga de los resultados en formato CSV. 
      output$downloadResults <- downloadHandler(filename = function() {
        paste0("resultados_clasificacion", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(results, file, row.names = FALSE)
      })
    })
  }
  
  # Lanzamiento de la aplicación. 
  shinyApp(ui = ui, server = server)