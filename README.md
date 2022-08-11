# HiContacts

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/js2264/HiContacts/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/js2264/HiContacts/actions/workflows/check-bioc.yml)
[![pkgdown](https://github.com/js2264/HiContacts/workflows/pkgdown/badge.svg)](https://github.com/js2264/HiContacts/actions)
<!-- badges: end -->

HiContacts provides tools to import `(m)cool` matrices in R and work with them there. 

Rather than creating redundant classes, it relies on pre-existing Bioconductor objects, namely `InteractionSet`, `GenomicInterations` and `ContactMatrix` (`Lun, Perry & Ing-Simmons, F1000Research 2016`).

## Installation

```r
remotes::install_github('js2264/HiContacts')
```

## Import a .(m)cool file as `contacts`

`contacts()` can be used to import a Hi-C matrix. 

```r
cool_file <- system.file("extdata", 'CH112.cool', package = 'HiContacts')
range <- 'I:20000-80000' # range of interest
contacts <- contacts(cool_file, focus = range)
contacts
```

`cool2gi` works with `.mcool` files as well: in this case, the user needs to specify the resolution at which count values are recovered. 

```r
mcool_file <- system.file("extdata", 'CH112.mcool', package = 'HiContacts')
range <- 'I:20000-80000' # range of interest
lsCoolResolutions(mcool_file)
contacts <- contacts(mcool_file, focus = range) # This throws an error!
contacts <- contacts(mcool_file, focus = range, res = 1000)
contacts
```

## Plotting matrices 

### Plot matrix heatmaps

> Diagonal style

```r
plotMatrix(contacts, use.assay = 'raw')
plotMatrix(contacts, use.assay = 'balanced', limits = c(-4, -1))
```

> Plot matrix of correlation 

```r
contacts(mcool_file, focus = 'II', res = 1000) %>% 
    autocorrelate() %>% 
    plotMatrix(scale = 'linear', limits = c(-1, 1))
```

> Plot signal over expected 

```r
contacts(mcool_file, focus = 'II', res = 1000) %>% 
    detrend() %>% 
    plotMatrix(use.assay = 'expected')
contacts(mcool_file, focus = 'II', res = 1000) %>% 
    detrend() %>% 
    plotMatrix(use.assay = 'detrended', scale = 'linear', limits = c(-2, 2))
```

### Plot loops

```r
loops <- system.file("extdata", 'CH112_1kb.bedpe', package = 'HiContacts') %>% 
    rtracklayer::import() %>% 
    InteractionSet::makeGInteractionsFromGRangesPairs()
contacts(mcool_file, focus = 'IV', res = 1000) %>% 
    detrend() %>% 
    plotMatrix(use.assay = 'detrended', loops = loops, scale = 'linear', limits = c(-2, 2))
```

### Plot borders

```r
borders <- rtracklayer::import('AT507_borders.bed')
contacts(mcool_file, focus = 'IV', res = 1000) %>% 
    detrend() %>% 
    plotMatrix(use.assay = 'detrended', loops = loops, borders = borders, scale = 'linear', limits = c(-2, 2))
```

### Many plots at once

```r
p <- cowplot::plot_grid(
    contacts(mcool_file, focus = 'IV', res = 1000) %>% plotMatrix(use.assay = 'raw'), 
    contacts(mcool_file, focus = 'IV', res = 1000) %>% plotMatrix(use.assay = 'balanced'), 
    contacts(mcool_file, focus = 'IV', res = 1000) %>% detrend() %>% plotMatrix(use.assay = 'expected', scale = 'linear'), 
    contacts(mcool_file, focus = 'IV', res = 1000) %>% detrend() %>% plotMatrix(use.assay = 'detrended', scale = 'linear', limits = c(-3, 3)), 
    contacts(mcool_file, focus = 'IV', res = 1000) %>% detrend() %>% plotMatrix(use.assay = 'detrended', scale = 'linear', limits = c(-3, 3), loops = loops), 
    contacts(mcool_file, focus = 'IV', res = 8000) %>% autocorrelate() %>% plotMatrix(scale = 'linear', limits = c(-1, 1))
)
ggsave('tmp.pdf', width = 20, height = 10)
```


## Analysis in R

### Auto-correlation 

```{r}

```

### Ratio of maps

```{r}
divide() %>% smoothen()
```

### Virtual 4C

```r
contacts <- contacts(mcool_file, res = 1000)
v4C <- virtual4C(contacts, viewpoint = GRanges('V:150000-170000'))
gg4C(v4C, aes = aes(x = center, y = score, col = chr)) + 
    facet_wrap(~chr, scales = 'free_x')
```

### Cis-trans ratios

```r
contacts <- contacts(mcool_file, res = 1000)
cis_trans(contacts)
```

### Compute P(s)

```r
pairs <- system.file("extdata", 'CH112.pairs', package = 'HiContacts') %>% 
    pairs2gi()
ps <- getPs(pairs)
ggPs(ps, aes(x = binned_distance, y = norm_p), xlim = c(5000, 5e5), ylim = c(1e-10, 1e-4))
```

### Compute scalograms

```{r}
pairs <- pairs2gi('tmp2.gz')
df <- getScalogram(pairs, ylim = c(1e4, 1e7))
p <- ggScalogram(df, aes(x = binned_pos, ymin = y0, ymax = dist_quantile, group = prob, fill = prob), ylim = c(1e4, 1e7))
```

### APA : Aggregated Plot Analysis

```{r}

```
