test_that("arithmetics works", {
    contacts_yeast <- HiCExperiment::contacts_yeast() |> 
        HiCExperiment::zoom(4000) |>
        HiCExperiment::refocus('II:1-100000')
    contacts_yeast_eco1 <- HiCExperiment::contacts_yeast_eco1() |> 
        HiCExperiment::zoom(4000) |>
        HiCExperiment::refocus('II:1-100000')
    full_contacts_yeast <- HiCExperiment::contacts_yeast() |> 
        HiCExperiment::zoom(4000) |>
        HiCExperiment::refocus('II:1-100000') |> 
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
        validObject(aggr_contacts)
    })
})
