# TFM-Raul-Garcia-Gonzalo

Este repositorio contiene el desarrollo completo del Trabajo de Fin de Máster (TFM) de Raúl García Gonzalo, centrado en la implementación de un pipeline reproducible para la detección y tipificación de cáncer a partir de datos de expresión génica obtenidos mediante RNA-Seq. El objetivo principal del trabajo es evaluar diferentes modelos de aprendizaje automático para clasificar muestras biológicas en función de si se trata de muestras sanas (normales) o tumorales y, en caso de ser tumorales, predecir el tipo específico de cáncer.

El repositorio incluye tanto el código fuente como los modelos entrenados y una aplicación web desarrollada con *Shiny*, la cual permite aplicar los modelos a nuevas muestras. El documento principal del trabajo es el fichero **TFM_RaúlGarcíaGonzalo.Rmd**, cuyo contenido se corresponde con el fichero **TFM_RaúlGarcíaGonzalo.pdf**. En dicho documento se encuentra todo el código utilizado y comentado paso a paso. 

Es importante destacar que el código incluido en el documento **TFM_RaúlGarcíaGonzalo.Rmd**  para la generación de la aplicación web interactiva no ejecuta directamente la aplicación *Shiny*. La aplicación debe ejecutarse mediante el fichero **app.R**, que constituye el punto de entrada real de la interfaz web. Por tanto, para utilizar la aplicación localmente es necesario lanzar directamente este fichero y no el RMarkdown.

Para que la ejecución sea correcta, el fichero **app.R** debe encontrarse en el mismo directorio de trabajo que los ficheros **models_with_genes.RData** y **test_samples.csv**, ya que estos contienen los modelos entrenados y ejemplos de muestras de entrada que son cargados automáticamente por la aplicación. Si estos ficheros no están disponibles en el mismo directorio, la aplicación no podrá inicializarse correctamente.

De forma análoga, si se desea ejecutar el análisis completo desde el fichero **TFM_RaúlGarcíaGonzalo.Rmd**, este debe situarse en el mismo directorio que los ficheros de datos originales **GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz** y **GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz**, ya que el *RMarkdown* asume que dichos ficheros están disponibles localmente para realizar los pasos de carga y preprocesamiento de datos.

El proyecto ha sido desarrollado en el entorno R versión 4.4.1, sobre un sistema operativo Windows 11 de 64 bits. La aplicación y el análisis se han implementado utilizando los siguientes paquetes y versiones, obtenidos mediante la función *sessionInfo()* en *RStudio*: 

  •	**rsconnect** (1.7.0) 

  •	**shiny** (1.10.0) 

  •	**iml** (0.11.4)

  •	**lime** (0.5.3)

  •	**ranger** (0.17.0) 

  •	**pROC** (1.18.5)

  •	**MLmetrics** (1.1.3)

  •	**lightgbm** (4.6.0)

  •	**kknn** (1.4.1)

  •	**class** (7.3-23)

  •	**reticulate** (1.42.0)

  •	**tensorflow** (2.16.0)

  •	**keras3** (1.4.0)

  •	**xgboost** (1.7.11.1)

  •	**e1071** (1.7-16)

  •	**randomForest** (4.7-1.2)

  •	**glmnet** (4.1-8) 

  •	**Matrix** (1.7-3)

  •	**tidyr** (1.3.1)

  •	**MASS** (7.3-65)

  •	**tibble** (3.2.1)

  •	**knitr** (1.50)

  •	**caret** (7.0-1)

  •	**lattice** (0.22-6)

  •	**ggplot2** (3.5.2)

  •	**kableExtra** (1.4.0)

  •	**dplyr** (1.1.4)

  •	**curl** (6.2.2)

  •	**data.table** (1.17.2)

  •	**limma** (3.62.2)

Finalmente, el repositorio incluye una versión desplegada de la aplicación en la plataforma *ShinyApps.io*, cuyo enlace se encuentra especificado en el documento del trabajo. Esta versión permite utilizar la aplicación directamente desde el navegador sin necesidad de configurar un entorno local.
