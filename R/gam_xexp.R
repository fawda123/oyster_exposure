library(knitr)
library(tidyverse)
library(readxl)
library(here)
library(lubridate)
library(janitor)
library(WtRegDO)

opts_chunk$set(warning = FALSE, message = FALSE)

# buoy, farm WDNR
obsdat1 <- read.csv(here('data/raw/2018 WDNR pH sensor data - 2018 WDNR pH sensor data.csv.csv')) %>% 
  select(Date, TimePST, BuoypHAG, MidpHP1) %>% 
  gather('loc', 'pH', -Date, -TimePST) %>% 
  mutate(
    loc = case_when(
      grepl('Buoy', loc) ~ 'WDNR Buoy', 
      grepl('Mid', loc) ~ 'WDNR Farm center'
    ), 
    loc = factor(loc, levels = c('WDNR Buoy', 'WDNR Farm edge', 'WDNR Farm center'))
  ) %>% 
  unite('datetime', Date, TimePST, sep = ' ') %>% 
  mutate(
    datetime = mdy_hms(datetime, tz = 'Pacific/Pitcairn')
  ) %>% 
  filter(as.Date(datetime) <= ymd('2018-07-07') & as.Date(datetime) >= ymd('2018-05-25')) %>% 
  filter(minute(datetime) == 0)

subdts <- ymd(c('2018-05-25', '2018-07-06'))

# SB1 1m
obsdat2 <- read.csv(here('data/raw/SB1_1m_PMEL_2017-2018.csv')) %>% 
  group_by(year, month, day, hour) %>% 
  sample_n(1) %>% # sometimes there's more than one measurement per hour
  ungroup %>% 
  unite('date', month, day, year, sep = '-') %>% 
  unite('time', hour, min, sec, sep = ':') %>% 
  unite('datetime', date, time, sep = ' ') %>% 
  mutate(
    datetime = mdy_hms(datetime, tz = 'UTC'), 
    datetime = with_tz(datetime, tz = 'Pacific/Pitcairn'),
    datetime = floor_date(datetime, unit = 'hour'),
    loc = 'SB1 1m'
  ) %>% 
  rename(pH = pH_tot) %>% 
  filter(as.Date(datetime) <= subdts[2] & as.Date(datetime) >= subdts[1])

# NB1 3m
obsdat3 <- read.csv(here('data/raw/NB1_3m_SeaFET001_pH_tot.csv')) %>% 
  group_by(year, month, day, hour) %>% 
  sample_n(1) %>% # sometimes there's more than one measurement per hour
  ungroup %>% 
  unite('date', month, day, year, sep = '-') %>% 
  unite('time', hour, min, sec, sep = ':') %>% 
  unite('datetime', date, time, sep = ' ') %>% 
  mutate(
    datetime = mdy_hms(datetime, tz = 'UTC'), 
    datetime = with_tz(datetime, tz = 'Pacific/Pitcairn'),
    datetime = floor_date(datetime, unit = 'hour'),
    loc = 'NB1 3m'
  ) %>% 
  rename(pH = pH_tot) %>% 
  filter(as.Date(datetime) <= subdts[2] & as.Date(datetime) >= subdts[1])

# combine
obsdat <- bind_rows(obsdat1, obsdat2, obsdat3)

# tidal data from ADCP
tiddat <- read.csv(here('data/raw/SL1_Lander_depth_UTC.csv')) %>%
  unite('datetime', year, month, day, hour, min, sec) %>% 
  # clean_names() %>%
  rename( 
    Tide = depth
  ) %>% 
  mutate(
    datetime = ymd_hms(datetime), 
    datetime = with_tz(datetime, tzone = 'Pacific/Pitcairn')
  ) %>% 
  filter(minute(datetime) == 0 & second(datetime) == 0) %>% 
  filter(datetime >= min(obsdat$datetime) & datetime <= max(obsdat$datetime))

# combine ph and tidal data
alldat <- obsdat %>% 
  full_join(tiddat, by = 'datetime') %>% 
  rename(DateTimeStamp = datetime) 

lat <- 47.87167
long <- -122.60491
tz <- 'Pacific/Pitcairn'


tmp <- alldat %>% 
  filter(loc %in% 'WDNR Farm center')

library(mgcv)

tomod <- alldat %>% 
  filter(loc %in% 'WDNR Farm center') %>% 
  mutate(
    date = as.Date(DateTimeStamp, tz = attr(DateTimeStamp, 'tzone')),
    sec = second(DateTimeStamp),
    min = minute(DateTimeStamp),
    hr = hour(DateTimeStamp), 
    mindate = date - min(date),
    secs = hr * 3600 + min * 60 + sec,
    dec_hr = secs / (60 * 60),
    dec_day = as.numeric(mindate) + (secs / 86400)
  ) %>% 
  select(DateTimeStamp, pH, Tide, dec_hr, dec_day)

modtmp <- gam(pH ~ s(dec_day, k = 40) + te(dec_hr, Tide, bs = c("cc", "tp")), data = tomod, 
              knots = list(doy = c(0, 24), Tide = 40))

modtmp <- gam(pH ~ te(dec_day, dec_hr, Tide, bs = c("tp", "cc", "tp")), data = tomod, 
              knots = list(dec_day = 40, doy = c(0, 24), Tide = 40))

toplo <- tomod %>% 
  mutate(
    pH_prd = predict(modtmp)
  )

ggplot(toplo, aes(x = DateTimeStamp)) + 
  geom_point(aes(y = pH)) + 
  geom_line(aes(y = pH_prd), col = 'red')

# create prediction matrix

