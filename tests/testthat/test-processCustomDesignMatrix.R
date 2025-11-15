# library(testthat); library(augere.de); source("test-processCustomDesignMatrix.R")

library(SummarizedExperiment)
cd <- DataFrame(treatment = rep(LETTERS[1:4], each=5), batch=factor(rep(1:5, 4)))
se <- SummarizedExperiment(colData=cd)

test_that("processCustomDesignMatrix works with a formula", {
    env <- new.env()
    env$SE <- se
    cmd <- processCustomDesignMatrix(~ batch + treatment, "SE", design.name = "FOOBAR")
    eval(parse(text=cmd), envir=env)

    expected <- model.matrix(~ batch + treatment, data=cd)
    expect_identical(env$FOOBAR, expected)
})

test_that("processCustomDesignMatrix works with a matrix", {
    env <- new.env()
    env$SE <- se
    mat <- model.matrix(~ batch + treatment, data=cd)
    cmd <- processCustomDesignMatrix(mat, "SE", design.name = "FOOBAR")
    eval(parse(text=cmd), envir=env)

    expect_identical(env$FOOBAR, mat)
})

test_that("processCustomDesignMatrix works with a function", {
    env <- new.env()
    env$SE <- se
    cmd <- processCustomDesignMatrix(function(x) {
        model.matrix(~ batch + treatment, data=colData(x))
    }, "SE", design.name = "FOOBAR")
    eval(parse(text=cmd), envir=env)

    expected <- model.matrix(~ batch + treatment, data=cd)
    expect_identical(env$FOOBAR, expected)
})
