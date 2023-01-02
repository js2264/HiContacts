test_that("v4C works", {
    contacts_yeast <- HiCExperiment::contacts_yeast() |> 
        HiCExperiment::refocus("II")
    expect_s4_class(
        virtual4C(contacts_yeast, GRanges('II:495000-505000')),
        'GRanges'
    )
    expect_s3_class({
        v4C <- virtual4C(contacts_yeast, GRanges('II:495000-505000'))
        plot4C(v4C, ggplot2::aes(x = center, y = score))
    }, 'gg')
})

test_that("arithmetics works", {
    contacts_yeast <- HiCExperiment::contacts_yeast() |> 
        HiCExperiment::refocus('II')
    contacts_yeast_eco1 <- HiCExperiment::contacts_yeast_eco1() |> 
        HiCExperiment::refocus('II')
    full_contacts_yeast <- HiCExperiment::full_contacts_yeast() |> 
        HiCExperiment::refocus('II') |> 
        HiCExperiment::zoom(resolution = 1000)
    expect_true({
        x <- detrend(contacts_yeast)
        validObject(x)
    })
    expect_true({
        x <- autocorrelate(contacts_yeast)
        validObject(x)
    })
    expect_true({
        div_contacts <- divide(contacts_yeast_eco1, by = contacts_yeast)
        validObject(div_contacts)
    })
    expect_true({
        merged_contacts <- merge(contacts_yeast_eco1, contacts_yeast)
        validObject(merged_contacts)
    })
    expect_true({
        aggr_contacts <- aggregate(
            full_contacts_yeast, 
            targets = topologicalFeatures(full_contacts_yeast, 'centromeres')
        )
        validObject(merged_contacts)
    })
})

test_that("cistrans works", {
    full_contacts_yeast <- HiCExperiment::full_contacts_yeast()
    expect_s3_class({
        cisTransRatio(full_contacts_yeast)
    }, 'tbl')
})

test_that("plotMatrix works", {
    contacts_yeast <- HiCExperiment::contacts_yeast()
    full_contacts_yeast <- HiCExperiment::full_contacts_yeast()
    expect_s3_class(plotMatrix(contacts_yeast), 'gg')
    expect_s3_class(plotMatrix(contacts_yeast, scale = 'linear'), 'gg')
    expect_s3_class(plotMatrix(contacts_yeast, scale = 'log2'), 'gg')
    expect_s3_class(
        plotMatrix(
            contacts_yeast, 
            scale = 'log10', 
            limits = c(-1, 1), 
            cmap = bwrColors()
        ),
        'gg'
    )
    expect_s3_class(
        plotMatrix(
            full_contacts_yeast,
            scale = 'exp0.2', 
            limits = c(-1, 1), 
            cmap = bbrColors()
        ),
        'gg'
    )
    expect_s3_class(
        plotMatrix(
            full_contacts_yeast['II:1-800000'],
            borders = topologicalFeatures(full_contacts_yeast, 'centromeres'),
            scale = 'exp0.2', 
            limits = c(-1, 1), 
            cmap = bbrColors(), 
            show_grid = TRUE
        ),
        'gg'
    )
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

test_that("utils work", {

    expect_true(is.matrix({
        m <- matrix(data = 0, nrow = 100, ncol = 100)
        sdiag(m, k = 0) <- 3
        m
    }))

})
