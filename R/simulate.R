simulateBorders <- function(
    x = contacts_yeast(), 
    n = 10, 
    binwidth = 30,
    factor = 100,
    SEED = 2021
) {
    set.seed(SEED)
    ints <- interactions(x)
    sc <- scores(x, use.scores)
    pts <- GRanges(
         seqnames = 'II', 
         IRanges(start(sample(regions(x), n)), width = resolution(x)*binwidth)
    )
    borders <- resize(pts, width = 1, fix = 'start')
    xx <- lapply(seq_along(pts), function(K) {
        gr <- pts[K]
        sc[anchors(x)[['first']] %over% gr & anchors(x)[['second']] %over% gr] <- sc[anchors(x)[['first']] %over% gr & anchors(x)[['second']] %over% gr] * factor
        sc
    }) 
    scores(x, 'sim') <- rowMeans(do.call(cbind, xx), na.rm = TRUE)
    topologicalFeatures(x, 'borders') <- borders
    x
}
