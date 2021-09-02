library(SWMPr)
library(tidyverse)
library(lubridate)
library(here)

pdbby <- import_local(here('~/Desktop/116853.zip'), 'pdbbywq', trace = T)
nocec <- import_local(here('~/Desktop/116853.zip'), 'nocecwq', trace = T)
nartb <- import_local(here('~/Desktop/116853.zip'), 'nartbwq', trace = T)
elkvm <- import_local(here('~/Desktop/116853.zip'), 'elkvmwq', trace = T)

pdbbyph <- qaqc(pdbby, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'pdbby') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)
nocecph <- qaqc(nocec, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'nocec') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)
nartbph <- qaqc(nartb, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'nartb') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)
elkvmph <- qaqc(elkvm, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'elkvm') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)

write.csv(pdbbyph, here('data/raw/pdbbyph.csv'), row.names = F)
write.csv(nocecph, here('data/raw/nocecph.csv'), row.names = F)
write.csv(nartbph, here('data/raw/nartbph.csv'), row.names = F)
write.csv(elkvmph, here('data/raw/elkvmph.csv'), row.names = F)

nerrsph <- bind_rows(pdbbyph, nocecph, nartbph, elkvmph)

p <- ggplot(nerrsph, aes(x = DateTimeStamp, y = pH)) + 
  geom_line() + 
  facet_wrap(~StationCode, ncol = 1) + 
  theme_minimal() +
  labs(
    x = NULL
  )

png(here('figs/nerrsph.png'), height = 8, width = 12, res = 400, units = 'in', family = 'serif') 
print(p)
dev.off()
