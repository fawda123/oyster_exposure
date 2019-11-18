library(tidyverse)
library(readxl)
library(here)

# consistent treatment terminology names
trts <- c("7.7A0.2", "7.7A0.5", "7.7C", "8.0A0.2", "8.0C")

# weight data processing --------------------------------------------------

# none of these files include treatments, just individual ids and week

# initial weights week 0
iniwts <- read.csv(here('data/raw', 'Oyster initial weights_Rfile.csv'), stringsAsFactors = F) %>% 
  rename(
    id = Individual.ID
  ) %>% 
  mutate(
    Week = 0
  ) 

# weights after week 0
trtwts <- read.csv(here('data/raw', 'Weight Data_all weeks_Rfile.csv'), stringsAsFactors = F) %>% 
  rename(
    id = Individual.ID
  ) %>% 
  select(-Date, -Comments2, -Unable.to.confirm.ID, -Respiration, -Feeding.calcification)

# length data processing --------------------------------------------------

alllen <- read_excel(here('data/raw', 'preliminary_oyster_data.xlsx'), 
                     col_types = c('text', 'numeric', 'text', 'text', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'), 
                     na = 'NA') %>% 
  rename(
    id = individual_id
  ) %>% 
  mutate(
    lendif = 
  )
