# HiContacts

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/js2264/HiContacts/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/js2264/HiContacts/actions/workflows/check-bioc.yml)
[![pkgdown](https://github.com/js2264/HiContacts/workflows/pkgdown/badge.svg)](https://github.com/js2264/HiContacts/actions)
<!-- badges: end -->

HiContacts provides tools to import `(m)cool` matrices in R and work with them there. 

It creates a new `contacts` class of objects, built on pre-existing Bioconductor objects, namely `InteractionSet`, `GenomicInterations` and `ContactMatrix` (`Lun, Perry & Ing-Simmons, F1000Research 2016`), and provides **analytical** and **visualization** tools to investigate contact maps. 

## Installation

```r
remotes::install_github('js2264/HiContacts')
```

## Import a .(m)cool file as `contacts`

```r
mcool_file <- HiContactsData::HiContactsData('yeast_wt', format = 'mcool')
range <- 'I:20000-80000' # range of interest
lsCoolResolutions(mcool_file)
contacts <- contacts(mcool_file, focus = range, res = 1000)
contacts
```

## Plotting matrices 

```r
plotMatrix(contacts, use.assay = 'raw')
plotMatrix(contacts, use.assay = 'balanced', limits = c(-4, -1))
```

## P(s)

```r
contacts <- contacts(
    mcool_file, 
    pairs = HiContactsData::HiContactsData('yeast_wt', format = 'pairs')
)
ps <- getPs(contacts)
ggPs(ps, aes(x = binned_distance, y = norm_p))
```

## Virtual 4C

```r
contacts <- contacts(mcool_file, res = 1000)
v4C <- virtual4C(contacts, viewpoint = GRanges('V:150000-170000'))
gg4C(v4C, aes(x = center, y = score, col = chr)) + 
    facet_wrap(~chr, scales = 'free_x')
```

## Cis-trans ratios

```r
contacts <- contacts(mcool_file, res = 1000)
cis_trans(contacts)
```
