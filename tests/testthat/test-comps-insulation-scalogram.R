test_that("compartments works", {
    library(BSgenome.Scerevisiae.UCSC.sacCer3)
    genome <- BSgenome.Scerevisiae.UCSC.sacCer3
    GenomeInfoDb::seqlevelsStyle(genome) <- "NCBI"
    full_contacts_yeast <- HiCExperiment::contacts_yeast(full = TRUE)
    comps <- getCompartments(full_contacts_yeast, genome = genome)
    comps_VI <- getCompartments(full_contacts_yeast, chromosomes = "VI")
    expect_no_error(getCompartments(full_contacts_yeast, 
        genome = genome, 
        chromosomes = "VI"
    ))
    expect_s4_class(
        comps_VI, 
        'HiCExperiment'
    )
    expect_s4_class(
        HiCExperiment::topologicalFeatures(comps_VI, 'compartments'), 
        'GRanges'
    )
    expect_s4_class(
        S4Vectors::metadata(comps_VI)$eigens, 
        'GRanges'
    )
    expect_s3_class(
        plotSaddle(comps_VI), 'gg'
    )
})

test_that("compartments can be phased with an RleList track", {
    library(BSgenome.Scerevisiae.UCSC.sacCer3)
    genome <- BSgenome.Scerevisiae.UCSC.sacCer3
    GenomeInfoDb::seqlevelsStyle(genome) <- "NCBI"
    full_contacts_yeast <- HiCExperiment::contacts_yeast(full = TRUE)
    chr_ids <- as.character(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(full_contacts_yeast)))
    gc_cov <- lapply(Biostrings::getSeq(genome, chr_ids), function(seq_chr) {
        smooth_bin <- 100L
        gc_sliding <- as.numeric(
            Biostrings::letterFrequencyInSlidingView(seq_chr, smooth_bin, "GC")
        ) / smooth_bin
        left_pad <- floor((smooth_bin - 1L) / 2L)
        right_pad <- ceiling((smooth_bin - 1L) / 2L)
        gc_smoothed <- c(
            rep(gc_sliding[[1]], left_pad),
            gc_sliding,
            rep(gc_sliding[[length(gc_sliding)]], right_pad)
        )
        S4Vectors::Rle(gc_smoothed)
    })
    cov_track <- IRanges::RleList(gc_cov)
    names(cov_track) <- chr_ids

    expect_no_error(getCompartments(
        full_contacts_yeast,
        genome = cov_track,
        chromosomes = "VI"
    ))
})

test_that("insulation works", {
    hic <- HiCExperiment::contacts_yeast() |> 
        HiCExperiment::refocus('II:1-300000') |> 
        HiCExperiment::zoom(resolution = 1000)
    hic2 <- getDiamondInsulation(hic, window_size = 8000)
    
    expect_no_error(getDiamondInsulation(hic, window_size = 8000))
    expect_no_error(getBorders(hic2))
    expect_s4_class(
        hic2, 
        'HiCExperiment'
    )
    expect_s4_class(
        HiCExperiment::topologicalFeatures(hic2, 'borders'), 
        'GRanges'
    )
    expect_s4_class(
        S4Vectors::metadata(hic2)$insulation, 
        'GRanges'
    )
})

test_that("scalogram works", {
    contacts_yeast <- HiCExperiment::contacts_yeast()
    pairsFile(contacts_yeast) <- HiContactsData::HiContactsData(
        'yeast_wt', format = 'pairs.gz'
    )
    scalo1 <- scalogram(contacts_yeast['II'])
    scalo2 <- scalogram(contacts_yeast['II'], probs = c(0, 0.3, 1))
    expect_no_error(scalogram(contacts_yeast['II']))
    expect_no_error(scalogram(contacts_yeast['II'], probs = c(0, 0.3, 1)))
    expect_s3_class(plotScalogram(scalo1), 'gg')
    expect_s3_class(plotScalogram(scalo2), 'gg')

})
