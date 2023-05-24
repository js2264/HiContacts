test_that("Contacts works", {
    contacts_yeast <- HiCExperiment::contacts_yeast()

    expect_warning(
        Contacts(fileName(contacts_yeast), focus = 'II', resolution = 16000)
    )
    
})

test_that("v4C works", {
    contacts_yeast <- HiCExperiment::contacts_yeast() |> 
        HiCExperiment::refocus("II")
    expect_s4_class(
        virtual4C(contacts_yeast, GRanges('II:495001-505000')),
        'GRanges'
    )
    expect_s3_class({
        v4C <- virtual4C(contacts_yeast, GRanges('II:495001-505000'))
        plot4C(v4C, ggplot2::aes(x = center, y = score))
    }, 'gg')
})

test_that("cistrans works", {
    contacts_yeast <- HiCExperiment::contacts_yeast()
    expect_s3_class({
        cisTransRatio(contacts_yeast)
    }, 'tbl')
})

test_that("Ps works", {
    contacts_yeast <- HiCExperiment::contacts_yeast()
    x1 <- distanceLaw(contacts_yeast)
    x2 <- localDistanceLaw(contacts_yeast, GRanges("II:1-100000"))
    pairsFile(contacts_yeast) <- HiContactsData::HiContactsData(
        'yeast_wt', format = 'pairs.gz'
    )
    y1 <- distanceLaw(contacts_yeast, by_chr = TRUE)
    y2 <- localDistanceLaw(contacts_yeast, GRanges("II:1-100000"))
    expect_s3_class(x1, 'tbl')
    expect_s3_class(y1, 'tbl')
    expect_s3_class(x2, 'tbl')
    expect_s3_class(y2, 'tbl')
    expect_s3_class(plotPs(x1, ggplot2::aes(x = binned_distance, y = norm_p)), 'gg')
    expect_s3_class(plotPsSlope(x1, ggplot2::aes(x = binned_distance, y = norm_p)), 'gg')
})

test_that("checks and utils work", {

    expect_true(is.matrix({
        m <- matrix(data = 0, nrow = 100, ncol = 100)
        sdiag(m, k = 0) <- 3
        m
    }))

    contacts_yeast <- HiCExperiment::contacts_yeast() |> 
        HiCExperiment::refocus("II:1-30000")
    contacts_yeast_eco1 <- HiCExperiment::contacts_yeast_eco1() |> 
        HiCExperiment::refocus("II:1-30000")

    expect_true(.is_symmetrical(contacts_yeast_eco1))
    expect_true(.is_comparable(contacts_yeast, contacts_yeast_eco1))
    expect_true(.are_HiCExperiment(contacts_yeast, contacts_yeast_eco1))
    expect_true(.is_same_seqinfo(contacts_yeast, contacts_yeast_eco1))
    expect_true(.is_same_resolution(contacts_yeast, contacts_yeast_eco1))
    expect_true(.is_same_bins(contacts_yeast, contacts_yeast_eco1))
    expect_true(.is_same_regions(contacts_yeast, contacts_yeast_eco1))

})
