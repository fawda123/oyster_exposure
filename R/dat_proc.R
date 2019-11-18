library(tidyverse)
library(readxl)
library(here)
library(janitor)

# consistent treatment terminology names
trts <- tibble(
  shrtlab = c("7.7A0.2", "7.7A0.5", "7.7C", "8.0A0.2", "8.0C"),
  lngslab = c("7.7 Fluctuating 0.2A", "7.7 Fluctuating 0.5A", "7.7 Constant", "8.0 Fluctuating 0.2A", "8.0 Constant")
  )

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


# respiration data --------------------------------------------------------

# respoiration, in umol_hr_g
allres <- read.csv(here::here("data/raw", "Respiration_oyster_alldata.csv"), header = T, stringsAsFactors = F) %>% 
  clean_names %>% 
  na.omit %>% 
  select(week, treatment, jar, id = oyster_id, resp = respiration_rate_umol_hr_g, species) %>% 
  mutate(
    ph = gsub('(^[0-9])*\\s.*$', '\\1', treatment), 
    amp = gsub('^[0-9]\\.[0-9]\\s', '', treatment)
  )

save(allres, file = here::here('data', 'allres.RData'), compress = 'xz')
