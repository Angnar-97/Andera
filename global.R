# ----- Project Description ----
#
# Shiny application for the study of microbiomes
# Date: 18/04/2026
# Author: Alejandro Navas González
# 

# ----- General Zone ----

# Loading of the packages for shiny development
library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(shinyalert)
library(shinyWidgets)
library(shinymeta)
library(htmlwidgets)

# Loading of the package used to edit the interface theme
library(fresh)

# Load the packages to be used for the tables in the application. 
library(tidyverse)
library(data.table)
library(DT)

# Loading of the packages for the statistical study
library(skimr)
library(ggcorrplot)
library(ggstatsplot)
library(PMCMRplus)
library(psych)

# Loading of bioconductor packages
library(phyloseq)
library(vegan)
library(microbiome)
library(ComplexHeatmap)
# library(microViz)

# Loading of packages to edit the graphics palettes
library(palettetown)
library(ggsci)
library(viridis)
library(RColorBrewer)


# Establish a common theme for the graphics.
theme_set(
  theme_light() + 
    theme(
      text = element_text(size = 16, family = "serif")
      )
  )
