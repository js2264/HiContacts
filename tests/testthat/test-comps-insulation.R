test_that("compartments works", {
    library(BSgenome.Scerevisiae.UCSC.sacCer3)
    genome <- BSgenome.Scerevisiae.UCSC.sacCer3
    GenomeInfoDb::seqlevelsStyle(genome) <- "NCBI"
    full_contacts_yeast <- HiCExperiment::full_contacts_yeast()
    comps <- getCompartments(full_contacts_yeast, genome = genome)
    comps_VI <- getCompartments(full_contacts_yeast, chromosomes = VI)
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
})
