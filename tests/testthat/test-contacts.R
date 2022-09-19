test_that("contacts works", {
    expect_true({
        #fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
        #contacts_yeast <- contacts(fpath, focus = 'II', resolution = 1000)
        data(contacts_yeast)
        validObject(contacts_yeast)
    })
    expect_s4_class(contacts_yeast, 'contacts')
    expect_identical(length(contacts_yeast), 74360L)
    expect_s4_class(contacts_yeast[1:10], 'contacts')
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
    expect_identical(matrixType(contacts_yeast), 'sparse')
    expect_type(coolPath(contacts_yeast), 'character')
    expect_s4_class(seqinfo(contacts_yeast), 'Seqinfo')
    expect_type(resolutions(contacts_yeast), 'integer')
    expect_equal(resolution(contacts_yeast), 1000L)
    expect_s4_class(bins(contacts_yeast), 'GRanges')
    expect_type(focus(contacts_yeast), 'character')
    expect_s4_class(interactions(contacts_yeast), 'GInteractions')
    expect_s4_class(scores(contacts_yeast), 'SimpleList')
    expect_s4_class(scores(contacts_yeast, 1), 'GInteractions')
    expect_s4_class(scores(contacts_yeast, 'raw'), 'GInteractions')
    expect_type(scores(contacts_yeast)[[1]], 'integer')
    expect_type(scores(contacts_yeast)[[2]], 'double')
    expect_s4_class(features(contacts_yeast), 'SimpleList')
    expect_s4_class(features(contacts_yeast, 1), 'GRanges')
    expect_type(pairsFile(contacts_yeast), 'NULL')
    expect_type(anchors(contacts_yeast), 'list')
    expect_error(summary(contacts_yeast), NA)
})

test_that("checks work", {
    data(contacts_yeast)
    data(full_contacts_yeast)
    expect_true(check_resolution(contacts_yeast, 2000))
    expect_error(check_resolution(contacts_yeast, 3000))
    expect_error(is_comparable(contacts_yeast, full_contacts_yeast))
    expect_true(is_square(S4Vectors::Pairs(
        first = 'I:10000-20000', 
        second = 'I:10000-20000'
    )))
    expect_true(is_symmetrical(contacts_yeast))
})

test_that("v4C works", {
    data(contacts_yeast)
    expect_s4_class(
        virtual4C(contacts_yeast, GRanges('II:490000-510000')),
        'GRanges'
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
        cisTransRatio(full_contacts_yeast)
    }, 'tbl')
})

test_that("plotMatrix works", {
    data(contacts_yeast)
    data(full_contacts_yeast)
    loops <- system.file("extdata", 'S288C-loops.bedpe', package = 'HiContacts') %>% 
        rtracklayer::import() %>% 
        InteractionSet::makeGInteractionsFromGRangesPairs()
    borders <- system.file("extdata", 'S288C-borders.bed', package = 'HiContacts') %>% 
        rtracklayer::import()
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
            loops = loops,
            borders = borders,
            scale = 'exp0.2', 
            limits = c(-1, 1), 
            cmap = bbrColors()
        ),
        'gg'
    )
})

test_that("Ps works", {
    data(contacts_yeast)
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
    file <- HiContactsData::HiContactsData(
        'yeast_wt', format = 'mcool'
    )
    expect_s4_class({
        contacts(file, resolution = 16000)
    }, 'contacts')
    expect_s4_class({
        contacts(file, focus = 'II:1-10000', resolution = 2000)
    }, 'contacts')
})
