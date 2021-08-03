# coolerr

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/js2264/coolerr/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/js2264/coolerr/actions/workflows/check-bioc.yml)
[![pkgdown](https://github.com/js2264/coolerr/workflows/pkgdown/badge.svg)](https://github.com/js2264/coolerr/actions)
<!-- badges: end -->

coolerr provides tools to import `(m)cool` matrices in R and work with them there. 

Rather than creating redundant classes, it relies on pre-existing Bioconductor objects, namely `InteractionSet`, `GenomicInterations` and `ContactMatrix` (`Lun, Perry & Ing-Simmons, F1000Research 2016`).

## Installation

```r
remotes::install_github('js2264/coolerr')
```

## Import a .cool file as GenomicInterations

`cool2gi` can be used to import contact counts as `GenomicInterations`. 

```r
file <- 'path/to/file.cool'
range <- 'ranges_of_interest' # e.g. range <- 'chrI:2000-8000'
gis <- cool2gi(file, coords = range)
gis
```

`cool2gi` works with `.mcool` files as well: in this case, the user needs to specify the resolution at which count values are recovered. 

```r
file <- 'path/to/file.mcool'
range <- 'ranges_of_interest' # e.g. range <- 'chrI:2000-8000'
lsCoolResolutions(file)
gis <- cool2gi(file, coords = range, res = 1000)
gis
```

## Plot matrix heatmaps

> Diagonal style

```r
gis <- cool2gi(file, coords = range, res = res)
p <- plotMatrix(gis, limits = c(-3, -1), dpi = 400)
```

> Two-entry diagonal style

```r
gis <- cool2gi(file, coords = range1, coords2 = range2, res = res)
p <- plotMatrix(gis, dpi = 400, symmetrical = FALSE)
```

> Horizontal style

```r
gis <- cool2gi(file, coords = range, res = res)
p <- plotTriangularMatrix(gis, dist_max = 0.4)
```

> Horizontal style with a list of multiple GenomicInterations

```r
gis1 <- cool2gi(file, coords = range, res = 40000)
gis2 <- cool2gi(file, coords = range, res = 80000)
gis3 <- cool2gi(file, coords = range, res = 160000)
p <- plotMatrixList(
    ls = list('res40kb' = gis1, 'res80kb' = gis2, 'res160kb' = gis3)
)
```

## Special matrices

> Plot matrix of correlation 

```r
gis <- cool2gi(file, coords = range, res = res)
p <- plotCorrelatedMatrix(gis, dpi = 500)
```

> Plot signal over expected 

```r
gis <- cool2gi(file, coords = range, res = res)
p <- plotOverExpected(gis, dpi = 500)
```

## Plot aggregated matrices

> At borders (e.g. TAD boundaries)

```r
file <- 'path/to/file.mcool'
# - Get TAD borders
coords <- rtracklayer::import('path/to/tads.bed') %>% 
    GenomicRanges::resize(fix = 'start', width = 1) %>% 
    GenomicRanges::resize(fix = 'center', width = 2000000)
# - Trim extended borders to only keep full-size GRanges
seqlevels(coords) <- seqlevels(cool2seqinfo(file, res = res))
seqinfo(coords) <- cool2seqinfo(file, res = res)
coords <- trim(coords)
coords <- coords[width(coords) == 2000000]
# - Plot aggregated signal over borders
p <- plotAggregatedMatrix(
    file, 
    coords[1:10], 
    res = res, 
    BPPARAM = BiocParallel::MulticoreParam(workers = 2, tasks = 200, progressbar = TRUE)
)
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 300)
```

> On pairs of coordinates (e.g. structural loops)

```r
file <- 'path/to/file.mcool'
# - Get loop coordinates from `bedpe` file as an `InteractionSet` object
loops <- rtracklayer::import('path/to/loops.bedpe') %>% ## interactions as `Pairs` object
    GenomicRanges::resize(width = 30000, fix = 'center')
# - Only retain loops for which extended anchors are both within the genome
loops <- loops[InteractionSet::anchors(loops)[[1]] %within% GRanges(cool2seqinfo(file)) & InteractionSet::anchors(loops)[[2]] %within% GRanges(cool2seqinfo(file))]
# - Plot aggregated signal over sets of loop coordinates
p <- plotAggregatedMatrix(file, coords = loops, symmetrical = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 16, tasks = 200, progressbar = TRUE))
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 300)
```

## Display features on matrices

> Plot TADs on the matrix

```r
gis <- cool2gi(file, coords = range, res = res)
tads <- rtracklayer::import('path/to/tads.bed')
p <- plotMatrix(gis, dpi = 500) %>% 
    addTads(tads, coords = range)
```

> Plot loops on the matrix

```r
gis <- cool2gi(file, coords = range, res = res)
loops <- rtracklayer::import('path/to/loops.bedpe')
p <- plotMatrix(gis, dpi = 500) %>% 
    addTads(tads, coords = range) %>% 
    addLoops(loops, coords = range)
```

> Plot genebodies and specific coordinates underneath a matrix

```r
## -- Generate heatmap
file <- 'path/to/file.mcool'
range <- 'chr13:110000000-115000000'
res <- 20000
gis <- cool2gi(file, coords = range, res = res)
p <- plotMatrix(gis, limits = c(-3, -1), dpi = 500)

## -- Import gene annotations
library(plyranges)
mm10_genes <- AnnotationHub::query(AnnotationHub::AnnotationHub(), c('Mus_musculus.GRCm39.104.gtf'))[[1]] %>% 
    filter(gene_biotype == 'protein_coding', source == 'ensembl_havana', type == 'gene') %>% 
    select(gene_id, gene_name) %>% 
    mutate(ID = gene_name)
seqlevelsStyle(mm10_genes) <- 'UCSC'

## -- Import CTCF binding sites
# mm10_CTCF <- rtracklayer::import('http://hgdownload.soe.ucsc.edu/gbdb/mm10/encode3/ccre/encodeCcreCombined.bb') %>% 
mm10_CTCF <- rtracklayer::import('~/encodeCcreCombined.bb') %>% 
    filter(grepl('CTCF-bound', ccre), encodeLabel == 'CTCF-only') 

## -- Import profiles
comps <- rtracklayer::import('~/Projects/20210602_MCCs_HiC-deuts/compartments/AT409_100kb.cis.bw', as = 'Rle')
insul <- rtracklayer::import('~/Projects/20210602_MCCs_HiC-deuts/tads/AT409_25kb_insulation-scores.bw', as = 'Rle')

## -- Combine all
p_withTracks <- addTracks(p, range, 
    annotations = list(genes = mm10_genes, CTCF = mm10_CTCF), 
    profiles = list(eigen = comps, insulation = insul)
)
ggplot2::ggsave('plot.png', plot = p_withTracks, width = 10, height = 10, dpi = 500)
```

## Complex examples 

```r
tads <- rtracklayer::import('path/to/tads.bed')
loops <- rtracklayer::import('path/to/loops.bedpe')
p <- cowplot::plot_grid(
    plotMatrix(cool2gi(file, coords = 'chr13', res = 160000), dpi = 500, limits = c(-3, 0)), 
    plotCorrelatedMatrix(cool2gi(file, coords = 'chr13', res = 160000), dpi = 500), 
    plotMatrixList(
        ls = list(
            'chr13' = cool2gi(file, coords = 'chr13', res = 40000), 
            '50-115Mb' = cool2gi(file, coords = 'chr13:50000000-115000000', res = 40000), 
            '100-115Mb' = cool2gi(file, coords = 'chr13:100000000-115000000', res = 40000),
            '111-115Mb' = cool2gi(file, coords = 'chr13:110000000-115000000', res = 40000)
        ), 
        limits = c(-3, 0)
    ) %>% 
        addTads(tads, coords = 'chr13:110000000-115000000') %>% 
        addLoops(loops, coords = 'chr13:110000000-115000000'), 
    plotOverExpected(cool2gi(file, coords = 'chr13:110000000-115000000', res = 10000), limits = c(-1, 1))
)
ggsave('HiC.png', width = 12, height = 12, dpi = 500)
```
