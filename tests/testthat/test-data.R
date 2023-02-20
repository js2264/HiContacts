test_that("Data works", {
    expect_s4_class(
        contacts_yeast(), 'HiCExperiment'
    )
    expect_s4_class(
        contacts_yeast_eco1(), 'HiCExperiment'
    )
})
