# library(testthat); library(augere.de); source("test-processSimpleDesignMatrix.R")

library(SummarizedExperiment)
cd <- DataFrame(group = rep(LETTERS[1:4], each=5), age=1:20, batch=factor(rep(1:5, 4)))
se <- SummarizedExperiment(colData=cd)

test_that("processSimpleDesignMatrix works as expected", {
    env <- new.env()
    env$superse <- se
    cmd <- processSimpleDesignMatrix("group", NULL, NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=cmd), envir=env)
    expect_identical(nrow(env$FOOBAR), ncol(se))
    expect_identical(colnames(env$FOOBAR), sprintf("group.%s", LETTERS[1:4]))

    env <- new.env()
    env$superse <- se
    cmd <- processSimpleDesignMatrix(NULL, "age", NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=cmd), envir=env)
    expect_identical(nrow(env$FOOBAR), ncol(se))
    expect_identical(colnames(env$FOOBAR), c("(Intercept)", "covariate.age"))

    env <- new.env()
    env$superse <- se
    cmd <- processSimpleDesignMatrix("group", "age", NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=cmd), envir=env)
    expect_identical(nrow(env$FOOBAR), ncol(se))
    expect_identical(colnames(env$FOOBAR), c(sprintf("group.%s", LETTERS[1:4]), "covariate.age"))

    env <- new.env()
    env$superse <- se
    cmd <- processSimpleDesignMatrix("group", NULL, "batch", se.name="superse", design.name="FOOBAR")
    eval(parse(text=cmd), envir=env)
    expect_identical(nrow(env$FOOBAR), ncol(se))
    expect_identical(colnames(env$FOOBAR), c(sprintf("group.%s", LETTERS[1:4]), sprintf("block.batch.%s", 2:5)))
})
