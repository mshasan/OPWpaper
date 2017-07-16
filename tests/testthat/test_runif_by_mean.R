test_that("both means should be equal", {
    expect_equal(mean(runif_by_mean(mean = 3.0, n = 100)), 3.0)
})
