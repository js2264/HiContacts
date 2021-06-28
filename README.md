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

## Plot `gis` heatmap

> Diagonal style

```r
file <- 'path/to/file.mcool'
range <- 'chr13:100000000-120000000'
res <- 40000
gis <- cool2gi(file, coords = range, res = res)
p <- plotMatrix(gis, limits = c(-3, -1), dpi = 1000)
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 1000)
```

> Plot signal over expected 

```r
file <- 'path/to/file.mcool'
range <- 'chr13:50000000-120000000'
res <- 40000
gis <- cool2gi(file, coords = range, res = res)
p <- plotOverExpected(gis, dpi = 1000)
ggplot2::ggsave('plot2.png', width = 10, height = 10, dpi = 1000)
```

> Horizontal style

```r
p <- plotTriangularMatrix(gis, limits = c(-3, -1), truncate_tip = 0.2)
ggplot2::ggsave('plot3.png', width = 10, height = 10, dpi = 1000)
```

> Horizontal style with a list of multiple `GenomicInterations`

```r
gis1 <- cool2gi(file, coords = range, res = 40000)
gis2 <- cool2gi(file, coords = range, res = 80000)
gis3 <- cool2gi(file, coords = range, res = 160000)
p <- plotMatrixList(
    ls = list('res40kb' = gis1, 'res80kb' = gis2, 'res160kb' = gis3)
)
ggplot2::ggsave('plot4.png', width = 10, height = 10, dpi = 300)
```
