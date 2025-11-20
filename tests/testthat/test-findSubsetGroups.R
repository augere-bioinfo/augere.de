# library(testthat); library(augere.de); source("test-findSubsetGroups.R")

test_that("findSubsetGroups works as expected", {
    contrast.info <- processSimpleComparisons(c("A", "B"))
    expect_identical(findSubsetGroups(contrast.info), c("A", "B"))

    contrast.info <- processSimpleComparisons(list(c("A", "B"), c("B", "C", "D")))
    expect_identical(findSubsetGroups(contrast.info), c("A", "B", "C", "D"))

    contrast.info <- processSimpleComparisons(list(c("age"), c("C", "D")))
    expect_null(findSubsetGroups(contrast.info))

    contrast.info <- processSimpleComparisons(c(foo="C", bar="A", foo="D"))
    expect_identical(findSubsetGroups(contrast.info), c("A", "C", "D"))

    contrast.info <- processSimpleComparisons("foo")
    expect_null(findSubsetGroups(contrast.info))
})
