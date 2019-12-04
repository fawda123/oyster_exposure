
output_func <- function(yvar){
  
  out <- paste0("
## ", yvar, " {.tabset .tabset-pills}

### Raw data

```{r, fig.height= 8, fig.width = 5, fig.align = 'center', out.width = '65%'}

yvar <- '", yvar, "'

toplo <- allexp %>% 
  filter(var %in% yvar)

p <- ggplot(toplo, aes(x = factor(week), y = val, colour = species, fill = species)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.1) + 
  geom_point(aes(group = species), position = position_jitterdodge(jitter.width = 0.1), size = 2, alpha = 0.6) + 
  theme_ipsum_rc() + 
  facet_grid(trt~ .) + 
  scale_colour_ipsum() + 
  scale_fill_ipsum() + 
  theme(legend.title = element_blank()) + 
  labs(x = 'Exposure week', y = yvar)

p
```

### Models by week

#### Olympia oyster

```{r, fig.height = 5, fig.width = 10, out.width = '90%', results = 'hide'}

bymods %>% 
  filter(species == 'Olympia') %>% 
  filter(var == yvar) %>% 
  pull(plomod)
```

#### Pacific oyster

```{r, fig.height = 5, fig.width = 10, out.width = '90%', results = 'hide'}

bymods %>% 
  filter(species == 'Pacific') %>% 
  filter(var == yvar) %>% 
  pull(plomod)
```

### Models with week

#### Olympia oyster

```{r}

dat <- wtmods %>% 
  filter(species == 'Olympia') %>% 
  filter(var == yvar)

dat %>% 
  filter(species == 'Olympia') %>% 
  pull(anomodnoint) %>% 
  .[[1]] %>% 
  knitr::kable(format = 'html', digits = 2, caption = paste0(yvar, ' vs week + treatment')) %>%
  kable_styling(full_width = T, font_size = 14) %>% 
  HTML

dat %>% 
  filter(species == 'Olympia') %>% 
  pull(anomod) %>% 
  .[[1]] %>% 
  knitr::kable(format = 'html', digits = 2, caption = paste0(yvar, ' vs week x treatment')) %>%
  kable_styling(full_width = T, font_size = 14) %>% 
  HTML
```

```{r, fig.height = 5, fig.width = 6, out.width = '70%', results = 'hide', fig.align = 'center'}

dat %>% 
  pull(plomod)
```

#### Pacific oyster

```{r}

dat <- wtmods %>% 
  filter(species == 'Pacific') %>% 
  filter(var == yvar)

dat %>% 
  filter(species == 'Pacific') %>% 
  pull(anomodnoint) %>% 
  .[[1]] %>% 
  knitr::kable(format = 'html', digits = 2, caption = paste0(yvar, ' vs week + treatment')) %>%
  kable_styling(full_width = T, font_size = 14) %>% 
  HTML

dat %>% 
  filter(species == 'Pacific') %>% 
  pull(anomod) %>% 
  .[[1]] %>% 
  knitr::kable(format = 'html', digits = 2, caption = paste0(yvar, ' vs week x treatment')) %>%
  kable_styling(full_width = T, font_size = 14) %>% 
  HTML
```

```{r, fig.height = 5, fig.width = 6, out.width = '70%', results = 'hide', fig.align = 'center'}

dat %>% 
  pull(plomod)
```
")
  
  return(out)
  
}

