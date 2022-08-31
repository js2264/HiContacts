test_that("contacts works", {
    expect_true({
        #fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
        #contacts_yeast <- contacts(fpath, focus = 'II', resolution = 1000)
        data(contacts_yeast)
        validObject(contacts_yeast)
    })
    expect_error(length(contacts_yeast), NA)
    expect_error(dim(contacts_yeast), NA)
    expect_error(contacts_yeast[1:10], NA)
    expect_error(contacts, NA)
})

test_that("v4C works", {
    data(contacts_yeast)
    expect_s3_class(
        virtual4C(contacts_yeast, GRanges('II:490000-510000')),
        'tbl'
    )
    expect_s3_class({
        data(contacts_yeast)
        v4C <- virtual4C(contacts_yeast, GRanges('II:490000-510000'))
        plot4C(v4C, ggplot2::aes(x = center, y = score))
    }, 'gg')
})

test_that("arithmetics works", {
    data(contacts_yeast)
    data(contacts_yeast_eco1)
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
})

test_that("cistrans works", {
    data(full_contacts_yeast)
    expect_s3_class({
        cis_trans(full_contacts_yeast)
    }, 'tbl')
})

test_that("plotMatrix works", {
    data(contacts_yeast)
    expect_s3_class(plotMatrix(contacts_yeast), 'gg')
    expect_s3_class(plotMatrix(contacts_yeast, scale = 'linear'), 'gg')
    expect_s3_class(plotMatrix(contacts_yeast, scale = 'log10'), 'gg')
    expect_s3_class(plotMatrix(contacts_yeast, scale = 'log2'), 'gg')
    expect_s3_class(
        plotMatrix(contacts_yeast, scale = 'log2', limits = c(-1, 1)),
        'gg'
    )
})

test_that("Ps works", {
    data(contacts_yeast)
    expect_s3_class(getPs(contacts_yeast), 'tbl')
    expect_s3_class(
        plotPs(getPs(contacts_yeast), aes(x = binned_distance, y = norm_p)),
        'gg'
    )
    expect_s3_class(
        plotPsSlope(getPs(contacts_yeast), aes(x = binned_distance, y = slope)),
        'gg'
    )
})
