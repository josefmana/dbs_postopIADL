# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", "tidyverse", # for data wrangling
  "brms", "bayestestR", # for prediction calculations and summaries 
  "ggplot2", "patchwork" # for plotting
)

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# read model and scaling values
m <- readRDS("shiny/stat_model.rds")$m
scl <- readRDS("shiny/stat_model.rds")$scl
