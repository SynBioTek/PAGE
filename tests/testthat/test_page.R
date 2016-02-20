library(PAGE)
context("PAGE")

data(Psoriasis_Etanercept_LogFC)
datM <- apply(-Psoriasis_Etanercept_LogFC, 2, rank)
gene_up <- names(head(sort(Psoriasis_Etanercept_LogFC[,1], decreasing=TRUE), 100) )
gene_dn <- names(head(sort(Psoriasis_Etanercept_LogFC[,1], decreasing=FALSE), 100) )


test_that("page return a list", {
    expect_is(page1(datM, gene_up), "list")
    expect_is(page1(Psoriasis_Etanercept_LogFC, gene_up, isRank=FALSE), "list")
    
    expect_is(page2(datM, gene_up, gene_dn), "list")
    expect_is(page2(Psoriasis_Etanercept_LogFC, gene_up, gene_dn, isRank=FALSE), "list")

})

test_that("input gene check", {
    expect_error( page2(datM, paste0(gene_up, "ZZ"), gene_dn) )
})
