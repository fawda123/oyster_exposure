library(tidyverse)
library(readxl)
library(here)
library(janitor)
library(lmerTest)
library(psycho)
library(hrbrthemes)
library(patchwork)
library(modelbased)

# consistent treatment terminology names
trts <- tibble(
  shrtlab = c("7.7C", "7.7A0.2", "7.7A0.5", "8.0C", "8.0A0.2"),
  lngslab = c( "7.7 Constant", "7.7 Fluctuating 0.2A", "7.7 Fluctuating 0.5A", "8.0 Constant", "8.0 Fluctuating 0.2A")
  )

# consistent response variable terminology names
rsps <- tibble(
  shrtlab = c("delt_size", "final_size", "rate_size", "delt_length", "delt_width", "dissolution_score", "final_length", "final_width", "folds", "irreg", "rate_length", "rate_width", "resp", "shell_weight", "tissue_weight", "total", "whole_organism_weight"),
  lngslab = c("Delta size (cm2)", "Final size (cm2)", "Size rate (cm2/d)", "Delta length (cm)", "Delta width (cm)", "Dissolution score (0-3)", "Final length (cm)", "Final width (cm)", "Number of folds", "Number of irregularities", "Length rate (cm/d)", "Width rate (cm/d)", "Respiration (umol/hr/g)", "Shell weight (g)", "Tissue weight (g)", "Total observations",  "Whole weight (g)")
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

# initial column measured a few days prior to week zero
alllen <- read_excel(here::here('data/raw/Length_kelp_4_3_2020.xlsx')) %>%
  clean_names() %>%
  select(week, trt = treatment, jar, id = individual_id, species, final_length = length_cm_final, initial_length = length_cm_initial, final_width = width_cm_final, initial_width = width_cm_initial) %>%
  mutate(
    final_length = case_when(
      week == 0 & is.na(final_length) ~ initial_length,
      T ~ final_length
    ),
    final_width = case_when(
      week == 0 & is.na(final_width) ~ initial_width,
      T ~ final_width
    ),
    delt_length = final_length - initial_length,
    delt_width = final_width - initial_width,
    rate_length = case_when(
      week == 2 ~ delt_length / 14,
      week == 4 ~ delt_length / 28
    ),
    rate_width = case_when(
      week == 2 ~ delt_width / 14,
      week == 4 ~ delt_width / 28
    )
  ) %>%
  gather('var', 'val', -week, -trt, -jar, -id, -species) %>%
  filter(week != 6) %>%
  filter(!(week %in% 0 & var %in% c('delt_length', 'delt_width', 'rate_width', 'rate_length'))) %>%
  filter(!var %in% c('initial_length', 'initial_width'))

##
# weight data

allwts <- read.csv(here::here('data/raw/Weight_with dead but not multiple_Kelp_4_3.csv'), stringsAsFactors = F) %>% 
  clean_names() %>% 
  filter(week != 0) %>% 
  filter(species != '') %>% 
  mutate(species = gsub('\\s*$', '', species)) %>% 
  rename(
    id = individual_id, 
    jar = i_jar,
    trt = treatment
  ) %>%
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

# body size, cm2
allsiz <- read_csv(here::here('data/raw/Oyster Data with Surface Area from Photoshopped Images.csv')) %>% 
    clean_names() %>%
    select(week, trt = treatment, jar, id = individual_id, species, final_size = surface_area_cm2_final, initial_size = surface_area_cm2_initial) %>%
    mutate(
      final_size = case_when(
        week == 0 & is.na(final_size) ~ initial_size,
        T ~ final_size
      ),
    delt_size = final_size - initial_size,
    rate_size = case_when(
      week == 2 ~ delt_size / 14,
      week == 4 ~ delt_size / 28
    )
  ) %>%
  gather('var', 'val', -week, -trt, -jar, -id, -species) %>%
  filter(week != 6) %>%
  filter(!(week %in% 0 & var %in% c('delt_size', 'rate_size'))) %>%
  filter(!var %in% c('initial_size'))

##
# combine all exposure data
allexp <- bind_rows(allres, alldis, allsiz, alllen, allwts) %>% 
  mutate(
    trt = factor(trt, levels = trts$shrtlab, labels = trts$shrtlab), 
    jar = factor(jar),
    var = factor(var, levels = rsps$shrtlab, labels = rsps$lngslab)
  ) %>% 
  arrange(species, var, week)

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
      mnsval <- estimate_means(mixmod, 'trt')
      cnsval <- estimate_contrasts(mixmod, 'trt') %>% 
        mutate(
          sig = ifelse(p < 0.05, 'sig', 'notsig'),
          sig = factor(sig,levels = c('notsig', 'sig'), labels = c(' not significant', 'significant'))
        ) %>% 
        unite('Contrast', Level1, Level2, sep = '-')

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
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high), colour = 'black', size = 1) + 
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
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high, colour = sig), size = 1) + 
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
      mnsval <- estimate_means(mixmod)

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
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high, colour = trt), width = 0, position = position_dodge(0.3), size = 1) + 
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

# data tracking individuals -----------------------------------------------

# consistent treatment terminology names
trts <- tibble(
  shrtlab = c("7.7C", "7.7A0.2", "7.7A0.5", "8.0C", "8.0A0.2"),
  lngslab = c( "7.7 Constant", "7.7 Fluctuating 0.2A", "7.7 Fluctuating 0.5A", "8.0 Constant", "8.0 Fluctuating 0.2A")
)

# consistent response variable names for weight
rsps <- tibble(
  shrtlab = c("final_size", "shell_weight", "tissue_weight"),
  lngslab = c("Final size (cm2)", "Shell weight (g)", "Tissue weight (g)")
)

# body size, cm2
allsiz <- read_csv(here::here('data/raw/Oyster Data with Surface Area from Photoshopped Images.csv')) %>% 
  clean_names() %>%
  select(week, trt = treatment, jar, id = individual_id, species, final_size = surface_area_cm2_final, initial_size = surface_area_cm2_initial) %>%
  mutate(
    final_size = case_when(
      week == 0 & is.na(final_size) ~ initial_size,
      T ~ final_size
    )
  ) %>%
  gather('var', 'val', -week, -trt, -jar, -id, -species) %>%
  filter(!var %in% c('initial_size'))

# weight data
allwts <- read.csv(here::here('data/raw/Weight_with dead but not multiple_Kelp_4_3.csv'), stringsAsFactors = F) %>% 
  clean_names() %>% 
  # filter(week != 0) %>% 
  filter(species != '') %>% 
  mutate(species = gsub('\\s*$', '', species)) %>% 
  rename(
    id = individual_id, 
    jar = i_jar,
    trt = treatment
  ) %>%
  select(week, trt, jar, id, species, shell_weight, tissue_weight) %>% 
  gather('var', 'val', shell_weight, tissue_weight) %>% 
  na.omit 

allind <- bind_rows(allsiz, allwts) %>% 
  mutate(
    trt = factor(trt, levels = trts$shrtlab, labels = trts$shrtlab), 
    var = factor(var, levels = rsps$shrtlab, labels = rsps$lngslab), 
    jar = factor(jar)
  )

save(allind, file = here::here('data', 'allind.RData'), compress = 'xz')

# models by week tracking individual differences --------------------------

trts <- tibble(
  shrtlab = c("7.7C", "7.7A0.2", "7.7A0.5", "8.0C", "8.0A0.2"),
  lngslab = c( "7.7 Constant", "7.7 Fluctuating 0.2A", "7.7 Fluctuating 0.5A", "8.0 Constant", "8.0 Fluctuating 0.2A")
)

# consistent response variable terminology names
rsps <- tibble(
  shrtlab = c('size', 'weight'),
  lngslab = c("Final size (cm2)", "Whole weight (g)")
)

# body size, cm2
szdat<- read_csv(here::here('data/raw/Oyster Data with Surface Area from Photoshopped Images.csv')) %>% 
  clean_names() %>%
  select(week, trt = treatment, jar, id = individual_id, species, final_size = surface_area_cm2_final, initial_size = surface_area_cm2_initial) %>%
  mutate(
    delt_size = final_size - initial_size, 
    init_size = 0,
  ) %>% 
  filter(week != 0) %>% 
  select(week, trt, jar, id, species, delt_size, init_size) %>% 
  gather('var', 'val', delt_size, init_size) %>% 
  separate(var, c('period', 'var')) %>% 
  group_by(week, trt, jar, id, species, period, var) %>% 
  summarise(val = mean(val, na.rm = T)) %>% 
  ungroup %>% 
  filter(!week == 0) %>% 
  mutate(
    week = case_when(
      period == 'init' ~ 0, 
      T ~ week
    )
  ) %>% 
  filter(week != 0)

wtdat <- read.csv(here::here('data/raw/Weight_with dead but not multiple_Kelp_4_3.csv'), stringsAsFactors = F) %>% 
  clean_names() %>% 
  # filter(week != 0) %>% 
  filter(species != '') %>% 
  mutate(
    species = gsub('\\s*$', '', species),
    whole_organism_weight = case_when(
      is.na(whole_organism_weight) &  !is.na(shell_weight) & !is.na(tissue_weight) ~ shell_weight + tissue_weight, 
      T ~ whole_organism_weight
    )
  ) %>% 
  filter(!is.na(initial_wetweight) & !is.na(whole_organism_weight)) %>% 
  select(week, trt = treatment, jar = i_jar, id = individual_id, species, initial_weight = initial_wetweight, 
         final_weight = whole_organism_weight) %>% 
  mutate(
    delt_weight = final_weight - initial_weight, 
    init_weight = 0,
  ) %>% 
  filter(week != 0) %>% 
  select(week, trt, jar, id, species, delt_weight, init_weight) %>% 
  gather('var', 'val', delt_weight, init_weight) %>% 
  separate(var, c('period', 'var')) %>% 
  group_by(week, trt, jar, id, species, period, var) %>% 
  summarise(val = mean(val, na.rm = T)) %>% 
  ungroup %>% 
  filter(!week == 0) %>% 
  mutate(
    week = case_when(
      period == 'init' ~ 0L, 
      T ~ week
    )
  ) %>% 
  filter(week != 0)

tomod <- bind_rows(wtdat, szdat) %>% 
  select(-period) %>% 
  mutate(
    trt = factor(trt, levels = trts$shrtlab), 
    var = factor(var, levels = rsps$shrtlab, rsps$lngslab), 
    jar = factor(jar)
  )

# fit separate mixed models by species, week, jar as random effect
indbymods <- tomod %>%
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
      mnsval <- estimate_means(mixmod, 'trt') %>% 
        data.frame
      cnsval <- estimate_contrasts(mixmod, 'trt') %>%
        mutate(
          sig = ifelse(p < 0.05, 'sig', 'notsig'),
          sig = factor(sig,levels = c('notsig', 'sig'), labels = c(' not significant', 'significant'))
        ) %>% 
        unite('Contrast', Level1, Level2, sep = '-')
      
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
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high), colour = 'black', size = 1) +
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
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high, colour = sig), size = 1) +
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

save(indbymods, file = here('data/indbymods.RData'), compress = 'xz')

#  models by week tracking individual differences, pos/neg ----------------

dirtomod <- tomod %>% 
  mutate(
    direction = sign(val)
  ) %>% 
  filter(direction != 0)

# fit separate mixed models by species, week, jar as random effect, separate by direction
dirindbymods <- dirtomod %>%
  group_by(week, species, var, direction) %>% 
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
      mnsval <- estimate_means(mixmod, 'trt') %>% 
        data.frame
      cnsval <- estimate_contrasts(mixmod, 'trt') %>%
        mutate(
          sig = ifelse(p < 0.05, 'sig', 'notsig'),
          sig = factor(sig,levels = c('notsig', 'sig'), labels = c(' not significant', 'significant'))
        ) %>% 
        unite('Contrast', Level1, Level2, sep = '-')
      
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
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high), colour = 'black', size = 1) +
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
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high, colour = sig), size = 1) +
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

save(dirindbymods, file = here('data/dirindbymods.RData'), compress = 'xz')
