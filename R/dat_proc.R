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

# consistent response variable terminology names
rsps <- tibble(
  shrtlab = c("dissolution_score", "folds", "irreg", "length_cm", "resp", "shell_weight", "tissue_weight", "total", "whole_organism_weight", "width_cm"),
  lngslab = c("Dissolution score (0-3)", "Number of folds", "Number of irregularities", "Length (cm)", "Respiration (umol/hr/g)", "Shell weight (g)", "Tissue weight (g)", "Total observations",  "Whole weight (g)", "Width (cm)")
)

# data processing to tidy -------------------------------------------------
 
##
# respiration, in umol_hr_g

allres <- read.csv(here::here("data/raw", "Respiration_oyster_alldata.csv"), header = T, stringsAsFactors = F) %>% 
  clean_names %>% 
  na.omit %>% 
  select(week, trt = treatment, jar, id = oyster_id, resp = respiration_rate_umol_hr_g, species) %>% 
  mutate(
    trt = factor(trt, levels = trts$lngslab, labels = trts$shrtlab), 
    trt = as.character(trt)
  ) %>% 
  gather('var', 'val', resp)

##
# length data
# some week have initial/final measurements, initial measurements span all four weeks

alllen <- read.csv(here::here('data/raw', 'oysterdata-combined-no blank ids.csv'), stringsAsFactors = F) %>% 
  clean_names() %>% 
  rename(
    id = individual_id,
    trt = treatment
  ) %>% 
  filter(initial_final == 'initial') %>% 
  filter(week != 0) %>%
  select(week, trt, jar, id, species, length_cm, width_cm) %>% 
  gather('var', 'val', length_cm, width_cm)

##
# weight data

# treatment by jar info (each jar is unique to a week) 
trtinfo <- read.csv(here::here('data/raw', 'oysterdata-combined-no blank ids.csv'), stringsAsFactors = F) %>% 
  rename(trt = treatment) %>% 
  select(trt, jar) %>% 
  unique

allwts <- read.csv(here::here('data/raw', 'Weight Data_all weeks_Rfile.csv'), stringsAsFactors = F) %>% 
  clean_names() %>% 
  filter(week != 0) %>% 
  filter(species != '') %>% 
  mutate(species = gsub('\\s*$', '', species)) %>% 
  rename(
    id = individual_id
  ) %>%
  left_join(trtinfo, by = c('jar')) %>% 
  select(week, trt, jar, id, species, whole_organism_weight, shell_weight, tissue_weight) %>% 
  gather('var', 'val', whole_organism_weight, shell_weight, tissue_weight) %>% 
  na.omit

##
# dissolution data
# there are multiple pictures per individual - dissolution averaged across pictures by jar, treatment, species, week, id
alldis <- read.csv(here::here("data/raw", "SEM scoring datasheet_MRVERSION.csv"), header = T, stringsAsFactors = F) %>% 
  clean_names %>% 
  select(week, trt = treatment, jar, id = individual_id, species, dissolution_score, folds = folds_in_shell, irreg = irregular_growth, total, picture) %>% 
  gather('var', 'val', dissolution_score, folds, irreg, total) %>% 
  group_by(week, trt, jar, id, species, var) %>% 
  summarise(
    val = mean(val, na.rm = T)
  ) %>% 
  ungroup %>% 
  filter(week != 0)

# ##
# # ingestion data
# alling <- read.csv(here::here('data/raw', 'Ingestionrate_Rfile.csv'), stringsAsFactors = F)

##
# combine all exposure data
allexp <- bind_rows(allres, alllen, alldis, allwts) %>% 
  mutate(
    trt = factor(trt, levels = trts$shrtlab, labels = trts$shrtlab), 
    jar = factor(jar),
    var = factor(var, levels = rsps$shrtlab, labels = rsps$lngslab)
  )

save(allexp, file = here::here('data', 'allexp.RData'), compress = 'xz')

# models by week ----------------------------------------------------------

data(allexp)

# fit separate mixed models by species, week, jar as random effect
bymods <- allexp %>%
  group_by(week, species, var) %>% 
  nest %>% 
  mutate(
    mixmod = pmap(list(var, data), function(var, data){
      
      tomod <- data %>%
        mutate(
          trt = fct_drop(trt),
          jar = fct_drop(jar)
        )
        
      # no replicate jars by treatment for whole weight
      if(var == 'Whole weight (g)')
        out <- glm(val ~ trt, data = tomod)
      
      if(var != 'Whole weight (g)')
        out <- lmerTest::lmer(val ~ trt + (1|jar), data = tomod)
      
      return(out)
      
    }),
    anomod = map(mixmod, function(x){
      
      # if(inherits(x, 'glm'))
      #   out <- summary(x)$coefficients
      # if(inherits(x, 'lmerModLmerTest'))
        out <- anova(x)
      
      return(out)
      
    }), 
    # summod = map(mixmod, analyze),
    plomod = pmap(list(week, species, var, mixmod), function(week, species, var, mixmod){
      
      # estimates
      mnsval <- get_means(mixmod, 'trt')
      cnsval <- get_contrasts(mixmod, 'trt') %>% 
        mutate(
          sig = ifelse(p < 0.05, 'sig', 'notsig'),
          sig = factor(sig,levels = c('notsig', 'sig'), labels = c(' not significant', 'significant'))
        )

      # sample size
      if(inherits(mixmod, 'glm'))
        n <- mixmod$model %>% nrow
      if(inherits(mixmod, 'lmerModLmerTest'))
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

save(bymods, file = here::here('data', 'bymods.RData'), compress = 'xz')

# models with week --------------------------------------------------------

data(allexp)

# treatment colors
cls <- RColorBrewer::brewer.pal(6, 'BrBG')
cls <- c(cls[1:3], cls[c(6, 4)])

# fit separate mixed models by species, week, jar as random effect
wtmods <- allexp %>%
  group_by(species, var) %>% 
  nest %>% 
  mutate(
    mixmod = pmap(list(var, data), function(var, data){

      tomod <- data %>% 
        mutate(
          week = factor(week), 
          jar = fct_drop(jar), 
          trt = fct_drop(trt)
        )
      
      # no replicate jars by treatment for whole weight
      if(var == 'Whole weight (g)')
        out <- glm(val ~ trt * week, data = tomod)
      
      if(var != 'Whole weight (g)')
        out <- lmerTest::lmer(val ~ trt * week + (1|jar), data = tomod)
      
      return(out)
      
    }),
    mixmodnoint = pmap(list(var, data), function(var, data){
      
      tomod <- data %>% 
        mutate(
          week = factor(week), 
          jar = fct_drop(jar), 
          trt = fct_drop(trt)
        )
      
      # no replicate jars by treatment for whole weight
      if(var == 'Whole weight (g)')
        out <- glm(val ~ trt + week, data = tomod)
      
      if(var != 'Whole weight (g)')
        out <- lmerTest::lmer(val ~ trt + week + (1|jar), data = tomod)
      
      return(out)
      
    }),
    anomod = map(mixmod, function(x){
      
      # if(inherits(x, 'glm'))
      #   out <- summary(x)$coefficients
      # if(inherits(x, 'lmerModLmerTest'))
        out <- anova(x)
      
      return(out)
      
    }), 
    anomodnoint = map(mixmodnoint, function(x){
      
      # if(inherits(x, 'glm'))
      #   out <- summary(x)$coefficients
      # if(inherits(x, 'lmerModLmerTest'))
      out <- anova(x)
      
      return(out)
      
    }), 
    # summod = map(mixmod, analyze),
    plomod = pmap(list(species, var, mixmod), function(species, var, mixmod){

      # estimates
      mnsval <- get_means(mixmod, 'trt * week')

      # sample size
      if(inherits(mixmod, 'glm'))
        n <- mixmod$model %>% nrow
      if(inherits(mixmod, 'lmerModLmerTest'))
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
      
      p1
     
    })
  )

save(wtmods, file = here::here('data', 'wtmods.RData'), compress = 'xz')


