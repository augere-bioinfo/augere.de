# library(testthat); library(augere.de); source("test-processCustomContrasts.R")

library(SummarizedExperiment)
cd <- DataFrame(treatment = rep(LETTERS[1:4], each=5), batch=factor(rep(1:5, 4)))
se <- SummarizedExperiment(colData=cd)

test_that("processCustomContrasts works for character vectors", {
    env <- new.env()
    env$SE <- se
    cmd <- processCustomDesignMatrix(~ 0 + treatment, "SE", design.name = "FOOBAR")
    eval(parse(text=cmd), envir=env)

    info <- processCustomContrasts("treatmentA - treatmentB", design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(env$WHEE, c(treatmentA = 1, treatmentB = -1, treatmentC = 0, treatmentD = 0))

    info <- processCustomContrasts(c("treatmentA - treatmentB", "treatmentC - treatmentD"), design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(colnames(env$WHEE), c("treatmentA - treatmentB", "treatmentC - treatmentD"))
})

test_that("processCustomContrasts works for functions", {
    env <- new.env()
    env$SE <- se
    cmd <- processCustomDesignMatrix(~ 0 + treatment, "SE", design.name = "FOOBAR")
    eval(parse(text=cmd), envir=env)

    info <- processCustomContrasts(function(X) {
        limma::makeContrasts(treatmentC - treatmentD, levels=X)
    }, design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(env$WHEE, c(treatmentA = 0, treatmentB = 0, treatmentC = 1, treatmentD = -1))

    info <- processCustomContrasts(function(X) {
        limma::makeContrasts(treatmentC - treatmentD, treatmentB - treatmentA, levels=X)
    }, design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(colnames(env$WHEE), c("treatmentC - treatmentD", "treatmentB - treatmentA"))
    expect_identical(env$WHEE[,1], c(treatmentA = 0, treatmentB = 0, treatmentC = 1, treatmentD = -1))
    expect_identical(env$WHEE[,2], c(treatmentA = -1, treatmentB = 1, treatmentC = 0, treatmentD = 0))

    info <- processCustomContrasts(function(X) {
        X <- limma::makeContrasts(treatmentC - treatmentD, treatmentB - treatmentA, levels=X)
        colnames(X) <- NULL
        X
    }, design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(colnames(env$WHEE), c("1", "2"))
    expect_identical(env$WHEE[,1], c(treatmentA = 0, treatmentB = 0, treatmentC = 1, treatmentD = -1))
    expect_identical(env$WHEE[,2], c(treatmentA = -1, treatmentB = 1, treatmentC = 0, treatmentD = 0))
})

test_that("processCustomContrasts works for vectors/matrices", {
    env <- new.env()
    env$SE <- se
    cmd <- processCustomDesignMatrix(~ 0 + treatment, "SE", design.name = "FOOBAR")
    eval(parse(text=cmd), envir=env)

    info <- processCustomContrasts(c(treatmentD = -1, treatmentA = 1), design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(env$WHEE, c(treatmentA = 1, treatmentB = 0, treatmentC = 0, treatmentD = -1))

    info <- processCustomContrasts(rbind(treatmentD = -1, treatmentA = 1), design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(env$WHEE, c(treatmentA = 1, treatmentB = 0, treatmentC = 0, treatmentD = -1))

    info <- processCustomContrasts(cbind(
        c(treatmentB = 1, treatmentA = -1, treatmentD = 0),
        c(treatmentB = 0, treatmentA = 1, treatmentD = -1)
    ), design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(colnames(env$WHEE), c("1", "2"))
    expect_identical(env$WHEE[,1], c(treatmentA = -1, treatmentB = 1, treatmentC = 0, treatmentD = 0))
    expect_identical(env$WHEE[,2], c(treatmentA = 1, treatmentB = 0, treatmentC = 0, treatmentD = -1))

    info <- processCustomContrasts(cbind(
        FOO=c(treatmentB = 1, treatmentA = -1, treatmentD = 0),
        BAR=c(treatmentB = 0, treatmentA = 1, treatmentD = -1)
    ), design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(info[[1]]$title, "Contrast 1")
    expect_identical(info[[1]]$type, "custom")
    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(colnames(env$WHEE), c("FOO", "BAR"))
    expect_identical(env$WHEE[,1], c(treatmentA = -1, treatmentB = 1, treatmentC = 0, treatmentD = 0))
    expect_identical(env$WHEE[,2], c(treatmentA = 1, treatmentB = 0, treatmentC = 0, treatmentD = -1))
})

test_that("processCustomContrasts works for lists", {
    env <- new.env()
    env$SE <- se
    cmd <- processCustomDesignMatrix(~ 0 + treatment, "SE", design.name = "FOOBAR")
    eval(parse(text=cmd), envir=env)

    combined <- processCustomContrasts(list(foo="treatmentA - treatmentB", bar="treatmentC - treatmentD"), design.name = "FOOBAR", contrast.name = "WHEE")
    expect_identical(combined[[1]]$title, "foo")
    combined[[1]]$title <- NULL
    expect_identical(combined[[2]]$title, "bar")
    combined[[2]]$title <- NULL

    alone1 <- processCustomContrasts("treatmentA - treatmentB", design.name = "FOOBAR", contrast.name = "WHEE")
    alone1[[1]]$title <- NULL
    expect_identical(alone1[[1]], combined[[1]])

    alone2 <- processCustomContrasts("treatmentC - treatmentD", design.name = "FOOBAR", contrast.name = "WHEE")
    alone2[[1]]$title <- NULL
    expect_identical(alone2[[1]], combined[[2]])
})
