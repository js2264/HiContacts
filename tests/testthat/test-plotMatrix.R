test_that("plotMatrix works", {
    contacts_yeast <- HiCExperiment::contacts_yeast()
    full_contacts_yeast <- HiCExperiment::full_contacts_yeast()
    cm <- Matrix::as.matrix(gi2cm(interactions(full_contacts_yeast), 'balanced'))
    y <- base::as.matrix(cm)
    expect_s3_class(plotMatrix(y), 'gg')
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
