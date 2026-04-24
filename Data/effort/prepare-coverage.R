library(readxl)
library(dplyr)
source(file.path("support_scripts", "data_helpers.R"))

#### EFFORT ####
dir <- './Data'
effort <- read_effort_coverage(dir) %>% 
  rbind(., data.frame(
    survey = paste(thresholds, 3, sep='_'),
    coverage = 'Good'
  ))
