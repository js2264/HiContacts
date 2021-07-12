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
gis <- cool2gi(file, coords = range)
gis
```

`cool2gi` works with `.mcool` files as well: in this case, the user needs to specify the resolution at which count values are recovered. 

```r
file <- 'path/to/file.mcool'
range <- 'chrI:2000-8000'
listCoolResolutions(file)
gis <- cool2gi(file, coords = range, res = 1000)
gis
```

## Plot matrix heatmaps

> Diagonal style

```r
file <- 'path/to/file.mcool'
range <- 'chr13:100000000-120000000'
res <- 40000
gis <- cool2gi(file, coords = range, res = res)
p <- plotMatrix(gis, limits = c(-3, -1), dpi = 1000)
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 1000)
```

> Two-entry diagonal style

```r
file <- 'path/to/file.mcool'
range1 <- 'chr13:101000000-107000000'
range2 <- 'chr13:104000000-110000000'
res <- 40000
gis <- cool2gi(file, coords = range1, coords2 = range2, res = res)
p <- plotMatrix(gis, limits = c(-3, -1), dpi = 1000)
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 1000)
```

> Horizontal style

```r
file <- 'path/to/file.mcool'
range <- 'chr13:100000000-120000000'
res <- 40000
gis <- cool2gi(file, coords = range, res = res)
p <- plotTriangularMatrix(gis, truncate_tip = 0.2)
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 1000)
```

> Horizontal style with a list of multiple `GenomicInterations`

```r
gis1 <- cool2gi(file, coords = range, res = 40000)
gis2 <- cool2gi(file, coords = range, res = 80000)
gis3 <- cool2gi(file, coords = range, res = 160000)
p <- plotMatrixList(
    ls = list('res40kb' = gis1, 'res80kb' = gis2, 'res160kb' = gis3)
)
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 300)
```

> Plot signal over expected 

```r
file <- 'path/to/file.mcool'
range <- 'chr13:50000000-120000000'
res <- 40000
gis <- cool2gi(file, coords = range, res = res)
p <- plotOverExpected(gis, dpi = 1000)
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 1000)
```

> Plot genebodies and specific coordinates underneath a matrix

```r
## -- Generate heatmap
file <- 'path/to/file.mcool'
range <- 'chr13:110000000-115000000'
res <- 20000
gis <- cool2gi(file, coords = range, res = res)
p <- plotMatrix(gis, limits = c(-3, -1), dpi = 1000)

## -- Import gene annotations
library(plyranges)
mm10_genes <- AnnotationHub::query(AnnotationHub::AnnotationHub(), c('Mus_musculus.GRCm39.104.gtf'))[[1]] %>% 
    filter(gene_biotype == 'protein_coding', source == 'ensembl_havana', type == 'gene') %>% 
    select(gene_id, gene_name) %>% 
    mutate(ID = gene_name)
seqlevelsStyle(mm10_genes) <- 'UCSC'

## -- Import CTCF binding sites
mm10_CTCF <- rtracklayer::import('http://hgdownload.soe.ucsc.edu/gbdb/mm10/encode3/ccre/encodeCcreCombined.bb') %>% 
    filter(grepl('CTCF-bound', ccre), encodeLabel == 'CTCF-only') 

## -- Import profiles
comps <- rtracklayer::import('~/Documents/PostDoc_Koszul/__Bioinfo/Projects/20210602_MCCs_HiC-deuts/compartments/AT409.cis.bw', as = 'Rle')
insul <- rtracklayer::import('~/Documents/PostDoc_Koszul/__Bioinfo/Projects/20210602_MCCs_HiC-deuts/tads/AT409_insulation-scores.bw', as = 'Rle')

## -- Combine all
p_withTracks <- addTracks(p, range, 
    annotations = list(genes = mm10_genes, CTCF = mm10_CTCF), 
    profiles = list(eigen = comps, insulation = insul)
)
ggplot2::ggsave('plot.png', plot = p_withTracks, width = 10, height = 10, dpi = 1000)
```

## Plot aggregated matrices

> On diagonal

```r
file <- 'path/to/file.mcool'
res <- 20000
coords <- rtracklayer::import('path/to/insulation-score.mcool') %>% 
    GenomicRanges::resize(fix = 'start', width = 1) %>% 
    GenomicRanges::resize(fix = 'center', width = 2000000)
p <- plotAggregatedMatrix(file, res, coords, BPPARAM = BiocParallel::MulticoreParam(workers = 16, tasks = 200, progressbar = TRUE))
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 300)
```

> On pairs of coordinates

```r
file <- 'path/to/file.mcool'
res <- 1000
loops <- rtracklayer::import('path/to/loops.bedpe') %>% ## interactions as `Pairs` object
loops <- rtracklayer::import('~/Downloads/pairs.bedpe') %>% ## interactions as `Pairs` object
    InteractionSet::makeGInteractionsFromGRangesPairs() %>% ## interactions as `InteractionSet` object
    GenomicRanges::resize(width = 30000, fix = 'center')
seqlevelsStyle(loops) <- 'NCBI'
loops <- loops[InteractionSet::anchors(loops)[[1]] %within% GRanges(cool2seqinfo(file, res = 1000)) & InteractionSet::anchors(loops)[[2]] %within% GRanges(cool2seqinfo(file, res = 1000))]
p <- plotAggregatedMatrix(file, res, loops, BPPARAM = BiocParallel::MulticoreParam(workers = 16, tasks = 200, progressbar = TRUE))
ggplot2::ggsave('plot.png', width = 10, height = 10, dpi = 300)
```
