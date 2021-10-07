library(SWMPr)
library(tidyverse)
library(lubridate)
library(here)

kach3 <- import_local(here('data/raw/810300.zip'), 'kach3wq', trace = T)
kacss <- import_local(here('data/raw/657751.zip'), 'kacsswq', trace = T)
sfbcc <- import_local(here('data/raw/810300.zip'), 'sfbccwq', trace = T)
sfbgc <- import_local(here('data/raw/810300.zip'), 'sfbgcwq', trace = T)
tjrbr <- import_local(here('data/raw/810300.zip'), 'tjrbrwq', trace = T)

kach3ph <- qaqc(kach3, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'kach3') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)
kacssph <- qaqc(kacss, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'kacss') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)
sfbccph <- qaqc(sfbcc, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'sfbcc') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)
sfbgcph <- qaqc(sfbgc, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'sfbgc') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)
tjrbrph <- qaqc(tjrbr, qaqc_keep = c('0', '1', '2', '3', '4', '5')) %>% 
  select(DateTimeStamp = datetimestamp, pH = ph) %>% 
  mutate(StationCode = 'tjrbr') %>% 
  filter(lubridate::year(DateTimeStamp) >= 2010)

write.csv(kach3ph, here('data/raw/kach3ph.csv'), row.names = F)
write.csv(kacssph, here('data/raw/kacssph.csv'), row.names = F)
write.csv(sfbccph, here('data/raw/sfbccph.csv'), row.names = F)
write.csv(sfbgcph, here('data/raw/sfbgcph.csv'), row.names = F)
write.csv(tjrbrph, here('data/raw/tjrbrph.csv'), row.names = F)

nerrsph <- bind_rows(kach3ph, kacssph, sfbccph, sfbgcph, tjrbrph)

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
