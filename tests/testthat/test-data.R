test_that("Data works", {
    contacts_yeast <- contacts_yeast()
    contacts_yeast_eco1 <- contacts_yeast_eco1()
    full_contacts_yeast <- full_contacts_yeast()
    expect_s4_class(
        contacts_yeast, 'HiCExperiment'
    )
    expect_s4_class(
        contacts_yeast_eco1, 'HiCExperiment'
    )
    expect_s4_class(
        full_contacts_yeast, 'HiCExperiment'
    )
})
