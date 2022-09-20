test_that("contacts works", {
    contacts_yeast <- contacts_yeast()
    expect_s4_class(contacts_yeast, 'contacts')
    expect_identical(length(contacts_yeast), 74360L)
    expect_s4_class(contacts_yeast[seq_len(10)], 'contacts')
    expect_s4_class({
        sub <- c(
            rep(TRUE, length(contacts_yeast)/2), 
            rep(FALSE, length(contacts_yeast)/2)
        )
        contacts_yeast[sub]
    }, 'contacts')
    expect_s4_class({
        contacts_yeast['II:1-10000']
    }, 'contacts')
    expect_s4_class({
        contacts_yeast['II:1-10000 x II:20000-40000']
    }, 'contacts')
    expect_type(fileName(contacts_yeast), 'character')
    expect_s4_class(seqinfo(contacts_yeast), 'Seqinfo')
    expect_type(resolutions(contacts_yeast), 'integer')
    expect_equal(resolution(contacts_yeast), 1000L)
    expect_s4_class(bins(contacts_yeast), 'GRanges')
    expect_type(focus(contacts_yeast), 'character')
    expect_s4_class(interactions(contacts_yeast), 'GInteractions')
    expect_s4_class(scores(contacts_yeast), 'SimpleList')
    expect_type(scores(contacts_yeast, 1), 'double')
    expect_type(scores(contacts_yeast, 'raw'), 'double')
    expect_type(scores(contacts_yeast, 'balanced'), 'double')
    expect_s4_class(topologicalFeatures(contacts_yeast), 'SimpleList')
    expect_s4_class(topologicalFeatures(contacts_yeast, 'loops'), 'Pairs')
    expect_s4_class(topologicalFeatures(contacts_yeast, 'borders'), 'GRanges')
    expect_type(pairsFile(contacts_yeast), 'NULL')
    expect_type(anchors(contacts_yeast), 'list')
    expect_error(summary(contacts_yeast), NA)
})

test_that("checks work", {
    contacts_yeast <- contacts_yeast()
    full_contacts_yeast <- full_contacts_yeast()
    expect_true(check_resolution(contacts_yeast, 2000))
    expect_error(check_resolution(contacts_yeast, 3000))
    expect_error(is_comparable(contacts_yeast, full_contacts_yeast))
    expect_true(is_square(S4Vectors::Pairs(
        first = 'I:10000-20000', 
        second = 'I:10000-20000'
    )))
    expect_true(is_symmetrical(contacts_yeast))
    expect_true(is_symmetrical(full_contacts_yeast))
})

test_that("v4C works", {
    contacts_yeast <- contacts_yeast()
    expect_s4_class(
        virtual4C(contacts_yeast, GRanges('II:490000-510000')),
        'GRanges'
    )
    expect_s3_class({
        v4C <- virtual4C(contacts_yeast, GRanges('II:490000-510000'))
        plot4C(v4C, ggplot2::aes(x = center, y = score))
    }, 'gg')
})

test_that("arithmetics works", {
    contacts_yeast <- contacts_yeast()
    contacts_yeast_eco1 <- contacts_yeast_eco1()
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
    full_contacts_yeast <- full_contacts_yeast()
    expect_s3_class({
        cisTransRatio(full_contacts_yeast)
    }, 'tbl')
})

test_that("plotMatrix works", {
    contacts_yeast <- contacts_yeast()
    full_contacts_yeast <- full_contacts_yeast()
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
            cmap = bbrColors()
        ),
        'gg'
    )
})

test_that("Ps works", {
    contacts_yeast <- contacts_yeast()
    x <- getPs(contacts_yeast)
    pairsFile(contacts_yeast) <- HiContactsData::HiContactsData(
        'yeast_wt', format = 'pairs.gz'
    )
    y <- getPs(contacts_yeast, by_chr = TRUE)
    expect_s3_class(x, 'tbl')
    expect_s3_class(y, 'tbl')
    expect_s3_class(plotPs(x, aes(x = binned_distance, y = norm_p)), 'gg')
    expect_s3_class(plotPs(x, aes(x = binned_distance, y = norm_p)), 'gg')
    expect_s3_class(plotPsSlope(x, aes(x = binned_distance, y = norm_p)), 'gg')
})

test_that("utils works", {

    expect_type(splitCoords(GRanges('II:30-40')), 'list')
    expect_type(splitCoords('II:30-40'), 'list')

    expect_type(formatCoords('II:30-40'), 'character')
    expect_type(formatCoords('II:30-40 x II:40-50'), 'character')
    expect_type(formatCoords(GRanges('II:30-40')), 'character')

    expect_s4_class(char2pair("II:30000-50000 x II:60000-80000"), 'Pairs')

    expect_s4_class(fullContactInteractions("II", 30000, 50000, 1000), 'GInteractions')

    expect_is({
        m <- matrix(data = 0, nrow = 100, ncol = 100)
        sdiag(m, k = 0) <- 3
        m
    }, 'matrix')

    expect_equal(
        sortPairs(char2pair("II:30000-50000 x II:60000-80000")), 
        char2pair("II:30000-50000 x II:60000-80000")
    )
    expect_equal(
        sortPairs(char2pair("II:90000-100000 x II:60000-80000")), 
        char2pair("II:60000-80000 x II:90000-100000")
    )

    expect_s4_class(
        asGInteractions(
            data.frame(
                seqnames1 = 'II', 
                start1 = 23, 
                end1 = 34, 
                seqnames2 = 'II', 
                start2 = 12, 
                end2 = 54
            )
        ), 
        'GInteractions'
    )
})

test_that("parse works", {
    cool <- HiContactsData::HiContactsData(
        'yeast_wt', format = 'cool'
    )
    mcool <- HiContactsData::HiContactsData(
        'yeast_wt', format = 'mcool'
    )
    expect_s4_class({
        contacts(cool, focus = 'II:1-10000')
    }, 'contacts')
    expect_s4_class({
        contacts(cool)
    }, 'contacts')
    expect_s4_class({
        contacts(mcool, focus = 'II:1-10000', resolution = 2000)
    }, 'contacts')
    expect_s4_class({
        contacts(mcool, resolution = 16000)
    }, 'contacts')
    expect_error({
        contacts(cool, resolution = 16000)
    })
    expect_error({
        contacts(mcool)
    })
})

test_that("coerce works", {
    contacts_yeast <- contacts_yeast()
    expect_s4_class({
        as(contacts_yeast, 'GInteractions')
    }, 'GInteractions')
    expect_s4_class({
        as(contacts_yeast, 'ContactMatrix')
    }, 'ContactMatrix')
    expect_is({
        as(contacts_yeast, 'matrix')
    }, 'matrix')
})
