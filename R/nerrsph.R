# pull ten years of ph data for select NERR sites
# uses data from ignore folder (not online) of swmp_rats repo

library(tidyverse)
library(here)

fls <- list.files('C:/proj/swmp_rats/ignore/proc1/', full.names = T)
fls <-  c('pdbgs', 'pdbby', 'sosch', 'sosva', 'sfbcc', 'sfbgc', 'elkvm', 'elknm', 'tjrbr', 'tjros', 'kacss', 'kach3') %>% 
  paste0(., 'wq') %>% 
  paste0(., collapse = '|') %>% 
  grep(., fls, value = T) 

for(fl in fls){
  
  cat(fl, '\n')
  
  load(fl)
  nm <- basename(fl) %>% 
    gsub('\\.RData$', '', .)
  
  dat <- nm %>% 
    get %>% 
    select(DateTimeStamp = datetimestamp, pH = ph) %>% 
    mutate(StationCode = gsub('wq$', '', nm)) %>% 
    filter(lubridate::year(DateTimeStamp) >= 2010)
  
  tosv <- gsub('wq$', 'ph.csv', nm) %>% 
    paste0('data/raw/', .)
  
  write.csv(dat, here(tosv), row.names = F)
  
}
