library(stringr)
context("scorePage")

data(data(Psoriasis_Etanercept_LogFC_mat))

test_that("scorePage return a list", {
    expect_is(scorePage(Psoriasis_Etanercept_LogFC_mat, type="avg"), "list")
    expect_is(scorePage(Psoriasis_Etanercept_LogFC_mat, type="max"), "list")

})
