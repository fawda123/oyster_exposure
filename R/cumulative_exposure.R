library(tidyverse)
library(lubridate)

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

tm <- '2017-07-01 00:00:00' %>% 
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
  )

plot(`7.0C` ~ datetime, tm, type = 'l')
