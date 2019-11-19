library(tidyverse)
library(readxl)
library(here)
library(janitor)
library(lmerTest)
library(psycho)

# consistent treatment terminology names
trts <- tibble(
  shrtlab = c("7.7C", "7.7A0.2", "7.7A0.5", "8.0C", "8.0A0.2"),
  lngslab = c( "7.7 Constant", "7.7 Fluctuating 0.2A", "7.7 Fluctuating 0.5A", "8.0 Constant", "8.0 Fluctuating 0.2A")
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
  ) 


# respiration data --------------------------------------------------------

# respoiration, in umol_hr_g
allres <- read.csv(here::here("data/raw", "Respiration_oyster_alldata.csv"), header = T, stringsAsFactors = F) %>% 
  clean_names %>% 
  na.omit %>% 
  select(week, trt = treatment, jar, id = oyster_id, resp = respiration_rate_umol_hr_g, species) %>% 
  mutate(
    ph = gsub('(^[0-9])*\\s.*$', '\\1', trt), 
    ph = factor(ph),
    amp = gsub('^[0-9]\\.[0-9]\\s', '', trt), 
    amp = factor(amp),
    trt = factor(trt, levels = trts$lngslab, labels = trts$shrtlab), 
    jar = factor(jar)
  )

save(allres, file = here::here('data', 'allres.RData'), compress = 'xz')

# respiration models ------------------------------------------------------

# fit separate mixed models by species, week, jar as random effect
resmods <- allres %>%
  group_by(week, species) %>% 
  nest %>% 
  mutate(
    mixmod = map(data, function(x){
      
      lmerTest::lmer(resp ~ trt + (1|jar), data = x)
      
    }),
    anomod = map(mixmod, anova), 
    summod = map(mixmod, analyze),
    plomod = pmap(list(week, species, mixmod), function(week, species, mixmod){
      
      # estimates
      mnsval <- get_means(mixmod, 'trt')
      cnsval <- get_contrasts(mixmod, 'trt') %>% 
        mutate(
          sig = ifelse(p < 0.05, 'sig', 'notsig'),
          sig = factor(sig,levels = c('notsig', 'sig'), labels = c(' not significant', 'significant'))
        )
      
      # sample size
      n <- mixmod@frame %>% nrow
      
      # labels
      subttl <- paste0('Respiration (umol/hr/g), week ', week, ', ', species, ' oyster')
      captns <- paste0('Significance where CI does not include zero, alpha = 0.05, total n = ', n)
      
      # mean estimate plots

      p1 <- ggplot(mnsval, aes(x = trt, y = Mean)) + 
        geom_point(size = 2) + 
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher), colour = 'black', size = 1) + 
        labs(x = NULL, y = 'Estimated means (+/- 95% CI)', title = 'Treatment estimates', subtitle = subttl) + 
        theme_ipsum() + 
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
        ) + 
        coord_flip()
      
      # contrast plots
      p2 <- ggplot(cnsval, aes(x = Contrast, y = Difference, colour = sig)) + 
        geom_point(aes(colour = sig), size = 2) + 
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher), colour = 'black', size = 1) + 
        labs(x = NULL, y = 'Estimated differences (+/- 95% CI)', title = 'Treatment differences', subtitle = subttl, caption = captns) +
        theme_ipsum() + 
        theme(
          legend.title = element_blank(), 
          legend.position = 'bottom', 
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
        ) +
        geom_hline(yintercept = 0, linetype = 'dotted', size = 1) + 
        scale_colour_manual(drop = F, values = c('black', 'tomato1')) + 
        coord_flip()
      
      p1 + p2 + plot_layout(ncol  = 2)   
      
      
    })
  )

save(resmods, file = here::here('data', 'resmods.RData'), compress = 'xz')


