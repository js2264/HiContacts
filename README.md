# coolerr

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![pkgdown](https://github.com/js2264/coolerr/workflows/pkgdown/badge.svg)](https://github.com/js2264/coolerr/actions)
<!-- badges: end -->

coolerr provides tools to import `cool` matrices in R and work with them there. 

Rather than creating redundant classes, it relies on pre-existing Bioconductor objects, namely `InteractionSet`, `GenomicInterations` and `ContactMatrix` (`Lun, Perry & Ing-Simmons, F1000Research 2016`).

## Installation

```r
remotes::install_github('js2264/coolerr')
```
## Import `cool` file as `GenomicInterations`

`cool2gi` can be used to import contact counts as `GenomicInterations`. 

```r
file <- 'path/to/file.cool'
range <- 'chrI:2000-8000'
cool2seqinfo(file)
cool2gi(file, coords = range)
```

`cool2gi` works with `.mcool` files as well: in this case, the user needs to specify the resolution at which count values are recovered. 

```r
file <- 'path/to/file.mcool'
range <- 'chrI:2000-8000'
listCoolResolutions(file)
cool2gi(file, coords = range, res = 1000)
```

