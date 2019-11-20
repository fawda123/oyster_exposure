library(tidyverse)
library(readxl)
library(here)
library(janitor)
library(lmerTest)
library(psycho)
library(hrbrthemes)
library(patchwork)

# consistent treatment terminology names
trts <- tibble(
  shrtlab = c("7.7C", "7.7A0.2", "7.7A0.5", "8.0C", "8.0A0.2"),
  lngslab = c( "7.7 Constant", "7.7 Fluctuating 0.2A", "7.7 Fluctuating 0.5A", "8.0 Constant", "8.0 Fluctuating 0.2A")
  )

# consistent dissolution terminology names
diss <- tibble(
  shrtlab = c("dissolution_score", "folds", "irreg", "total"),
  lngslab = c("Dissolution score (0-3)", "Number of folds", "Number of irregularities", "Total observations")
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

data(allres)

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
      captns <- paste0('Significance is where CI does not include zero, alpha = 0.05, total n = ', n)
      
      # mean estimate plots

      p1 <- ggplot(mnsval, aes(x = trt, y = Mean)) + 
        geom_point(size = 3) + 
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
        geom_point(aes(colour = sig), size = 3) + 
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher, colour = sig), size = 1) + 
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

# respiration models with week ----------------------------------------------

data(allres)

# treatment colors
cls <- RColorBrewer::brewer.pal(6, 'BrBG')
cls <- c(cls[1:3], cls[c(6, 4)])

# fit separate mixed models by species, jar as random effect
reswksmods <- allres %>%
  gather('var', 'val', resp) %>% 
  group_by(species) %>% 
  nest %>% 
  mutate(
    mixmod = map(data, function(x){
      
      tomod <- x %>% 
        mutate(
          week = factor(week), 
          jar = fct_drop(jar)
        )
      
      lmerTest::lmer(val ~ trt * week + (1|jar), data = tomod)
      
    }),
    anomod = map(mixmod, anova), 
    summod = map(mixmod, analyze),
    plomod = pmap(list(species, mixmod), function(species, mixmod){
      
      
      # estimates
      mnsval <- get_means(mixmod, 'trt * week')
      
      # sample size
      n <- mixmod@frame %>% nrow
      
      # labels
      subttl <- paste0('Respiration (umol/hr/g), ', species, ' oyster')
      captns <- paste0('Significance is where CI does not include zero, alpha = 0.05, total n = ', n)
      
      # mean esimate plots
      p1 <- ggplot(mnsval, aes(x = week, y = Mean, group = trt, colour = trt, fill = trt)) + 
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher, colour = trt), width = 0, position = position_dodge(0.3), size = 1) + 
        geom_line(position = position_dodge(0.3)) + 
        geom_point(size = 3, position = position_dodge(0.3), pch = 21, colour = 'black') + 
        labs(x = 'Exposure week', y = 'Estimated means (+/- 95% CI)', title = 'Treatment estimates', subtitle = subttl, caption = captns) + 
        scale_colour_manual(values = cls) + 
        scale_fill_manual(values = cls) + 
        theme_ipsum() + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.title = element_blank()
        )
      
      p1
      
      
    })
  )

save(reswksmods, file = here::here('data', 'reswksmods.RData'), compress = 'xz')

# dissolution data --------------------------------------------------------

# there are multiple pictures per individual - dissolution averaged across pictures by jar, treatment, species, week, id
alldis <- read.csv(here::here("data/raw", "SEM scoring datasheet_MRVERSION.csv"), header = T, stringsAsFactors = F) %>% 
  clean_names %>% 
  select(week, trt = treatment, jar, id = individual_id, species, dissolution_score, folds = folds_in_shell, irreg = irregular_growth, total, picture) %>% 
  gather('var', 'val', dissolution_score, folds, irreg, total) %>% 
  mutate(
    trt = factor(trt, levels = trts$shrtlab, labels = trts$shrtlab), 
    var = factor(var, levels = diss$shrtlab, labels = diss$lngslab),
    jar = factor(jar)
  ) %>% 
  group_by(week, trt, jar, id, species, var) %>% 
  summarise(
    val = mean(val, na.rm = T)
  ) %>% 
  ungroup %>% 
  filter(week != 0)

save(alldis, file = here::here('data', 'alldis.RData'), compress = 'xz')

# dissolution models ------------------------------------------------------

data(alldis)

# fit separate mixed models by species, week, jar as random effect
dismods <- alldis %>%
  group_by(week, species, var) %>% 
  nest %>% 
  mutate(
    mixmod = map(data, function(x){
      
      lmerTest::lmer(val ~ trt + (1|jar), data = x)
      
    }),
    anomod = map(mixmod, anova), 
    summod = map(mixmod, analyze),
    plomod = pmap(list(week, species, var, mixmod), function(week, species, var, mixmod){
      
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
      subttl <- paste0(var, ', week ', week, ', ', species, ' oyster')
      captns <- paste0('Significance is where CI does not include zero, alpha = 0.05, total n = ', n)
      
      # mean esimate plots
      p1 <- ggplot(mnsval, aes(x = trt, y = Mean)) + 
        geom_point(size = 3) + 
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
        geom_point(aes(colour = sig), size = 3) + 
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher, colour = sig), size = 1) + 
        labs(x = NULL, y = 'Estimated differences (+/- 95% CI)', title = 'Treatment differences', subtitle = subttl,
             caption = captns) +
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

save(dismods, file = here::here('data', 'dismods.RData'), compress = 'xz')


# dissolution models with week --------------------------------------------

data(alldis)

# treatment colors
cls <- RColorBrewer::brewer.pal(6, 'BrBG')
cls <- c(cls[1:3], cls[c(6, 4)])

# fit separate mixed models by species, week, jar as random effect
diswksmods <- alldis %>%
  group_by(species, var) %>% 
  nest %>% 
  mutate(
    mixmod = map(data, function(x){

      tomod <- x %>% 
        mutate(
          week = factor(week), 
          jar = fct_drop(jar)
        )
      
      lmerTest::lmer(val ~ trt * week + (1|jar), data = tomod)
      
    }),
    anomod = map(mixmod, anova), 
    summod = map(mixmod, analyze),
    plomod = pmap(list(species, var, mixmod), function(species, var, mixmod){

      # estimates
      mnsval <- get_means(mixmod, 'trt * week')

      # sample size
      n <- mixmod@frame %>% nrow
      
      # labels
      subttl <- paste0(var, ', ', species, ' oyster')
      captns <- paste0('Significance is where CI does not include zero, alpha = 0.05, total n = ', n)
      
      # mean esimate plots
      p1 <- ggplot(mnsval, aes(x = week, y = Mean, group = trt, colour = trt, fill = trt)) + 
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher, colour = trt), width = 0, position = position_dodge(0.3), size = 1) + 
        geom_line(position = position_dodge(0.3)) + 
        geom_point(size = 3, position = position_dodge(0.3), pch = 21, colour = 'black') + 
        labs(x = 'Exposure week', y = 'Estimated means (+/- 95% CI)', title = 'Treatment estimates', subtitle = subttl, caption = captns) + 
        scale_colour_manual(values = cls) + 
        scale_fill_manual(values = cls) + 
        theme_ipsum() + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.title = element_blank()
        )
     
    })
  )

save(diswksmods, file = here::here('data', 'diswksmods.RData'), compress = 'xz')


