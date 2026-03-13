library(readxl)
#### EFFORT ####
get_effort <- function(dir){
  atlas1 <- read_excel(file.path(dir, 'effort', 'dækning atlas1+2.xlsx'), sheet = 'kvadrater atlas1')
  atlas2 <- read_excel(file.path(dir, 'effort', 'dækning atlas1+2.xlsx'), sheet = 'kvadrater atlas2')
  
  atlas2 <- atlas2 %>% 
    select(-UNDS_KODE_) %>% 
    rename(UNDS_KODE_ = UNDS_KOD_1)
  
  # Load design (Note: 'design' object is overwritten here in original code)
  # load('HmscOutputs/2026-02-10_07-03-53_All_All_Atlas3_MinOccs5/Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata')
  # design <- read.csv('Data/data/1_preprocessing/design/studyDesign.csv')
  # xycoords <- design[design$atlas == 2, c('site', 'lat', 'lon')]
  
  key <- read_excel(file.path(dir, 'effort', 'dækning atlas1+2.xlsx'), sheet = 'dækning') %>% 
    mutate(coverage = case_when(
      UNDSKODE == 0 ~ 'Not understood',
      UNDSKODE == 1 ~ 'Good',
      UNDSKODE == 2 ~ 'Average',
      UNDSKODE == 3 ~ 'Bad'
    ))
  
  effort <- rbind(
    atlas1 %>% mutate(atlas = 1),
    atlas2 %>% mutate(atlas = 2)
  ) %>% 
    left_join(key %>% rename(UNDS_KODE_ = UNDSKODE), by = "UNDS_KODE_") %>% 
    rename(site = KVADRATKOD) #%>% 
    #left_join(xycoords, by = 'site') 
  
  effort$coverage <- factor(effort$coverage, levels = c('Good','Average','Bad','Not understood'))
  return(effort)
}

dir <- './Data'
effort <- get_effort(dir) %>% 
  mutate(survey = paste(.$site, .$atlas, sep='_')) %>% 
  select(survey, coverage) %>% 
  rbind(., data.frame(
    survey = paste(thresholds, 3, sep='_'),
    coverage = 'Good'
  ))
