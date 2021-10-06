library(tidyverse)
library(lubridate)
library(janitor)

trts <- tibble(
  shrtlab = c("7.7C", "7.7A0.2", "7.7A0.5", "8.0C", "8.0A0.2"),
  lngslab = c("7.7 Constant", "7.7 Fluctuating 0.2A", "7.7 Fluctuating 0.5A", "8.0 Constant", "8.0 Fluctuating 0.2A")
)


# dissolution data
dissdat <- read.csv('data/raw/SEM scoring datasheet_MRVERSION.csv') %>% 
  clean_names() %>% 
  filter(week != 0) %>% 
  select(
    jar, 
    id = individual_id, 
    species, 
    week, 
    trt = treatment,
    rep, 
    scr = dissolution_score1
  ) %>% 
  group_by(jar, id, species, week, trt) %>% 
  summarise(
    val = mean(scr, na.rm = T),
    .groups = 'drop'
  ) %>% 
  mutate(
    jar = factor(jar),
    trt = factor(trt, levels = trts$shrtlab)
  ) %>% 
  group_by(species, week, trt) %>% 
  summarise(
    val = mean(val, na.rm = T), 
    .groups = 'drop'
    )


# function for creating sine wave
sinefunc <- function(time_in, alpha = 0, beta = 1, freq = 1, phi = 0){
  
  # timestep per hour
  time_step <- unique(diff(time_in))
  
  # set phi as difference in hours from start of time_in
  phi <- min(time_in) + phi * 3600
  phi<- as.numeric(difftime(phi, min(time_in)))
  phi <- phi / time_step
  
  # get input values to cos func
  in_vals <- seq(0, length(time_in), length = length(time_in))
  in_vals <- in_vals / time_step
  in_vals <- 2 * pi * in_vals * 1 / freq
  
  # wave
  y <- alpha + beta * sin(in_vals + phi)
  
  return(y)
  
}

trtdat <- '2017-07-01 00:00:00' %>% 
  ymd_hms(tz = Sys.timezone()) %>% 
  seq.POSIXt(., to = . + weeks(6), by = 'hours') %>% 
  tibble(
    datetime = .
  ) %>% 
  mutate(
    `8.0C` = sinefunc(datetime, freq = 12, beta = 0, alpha = 8),
    `8.0A0.2` = sinefunc(datetime, freq = 12, beta = 0.2, alpha = 8),
    `7.7C` = sinefunc(datetime, freq = 24, beta = 0, alpha = 7.7),
    `7.7A0.2` = sinefunc(datetime, freq = 24, beta = 0.2, alpha = 7.7),
    `7.7A0.5` = sinefunc(datetime, freq = 24, beta = 0.5, alpha = 7.7)
  ) %>% 
  gather('trt', 'ph', -datetime) %>% 
  mutate(
    trt = factor(trt, levels = trts$shrtlab)
  )

ggplot(trtdat, aes(x = datetime, y = ph)) + 
  geom_line() + 
  facet_wrap(~trt, ncol = 1)

lim <- 7.7

expdat <- trtdat %>% 
  group_by(trt) %>% 
  mutate(
    `7.0` = cumsum(ph <= 7), 
    `7.5` = cumsum(ph <= 7.5),
    `8.0` = cumsum(ph <= 8),
    `8.5` = cumsum(ph <= 85)
  )

toplo <- expdat %>% 
  select(-ph) %>% 
  gather('thr', 'hours', -datetime, -trt) %>% 
  mutate(
    trt = factor(trt, levels = trts$shrtlab)
  )

ggplot(toplo, aes(x = datetime, y = hours, col = trt)) + 
  geom_line() +
  facet_grid(trt~thr)

tocmp <- toplo %>% 
  group_by(trt, thr) %>% 
  mutate(
    week = week(datetime) - min(week(datetime)), 
    dups = duplicated(week, fromLast = T)
  ) %>% 
  filter(week %in% c(2, 4, 6) & !dups) %>% 
  select(trt, thr, hours, week) %>% 
  full_join(dissdat, by = c('trt', 'week')) %>% 
  filter()

ggplot(tocmp, aes(x = hours, y = val, color = trt, group = species)) + 
  geom_point() + 
  facet_grid(thr~species)
  

