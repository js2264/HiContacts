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
