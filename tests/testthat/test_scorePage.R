library(PAGE)
context("PAGE")

data(Psoriasis_Etanercept_LogFC)

test_that("scorePage return a list", {
    expect_is(scorePage(Psoriasis_Etanercept_LogFC, type="avg", ncore=1), "list")
    expect_is(scorePage(Psoriasis_Etanercept_LogFC, type="max", ncore=1), "list")
    
    dat_list <- lapply(data.frame(Psoriasis_Etanercept_LogFC), 
        function(x) {names(x) <- rownames(Psoriasis_Etanercept_LogFC); x} )
    expect_is(scorePage(dat_list, type="avg"), "list")
     

})
