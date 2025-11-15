# library(testthat); library(augere.de); source("test-processSimpleComparisons.R")

library(SummarizedExperiment)
cd <- DataFrame(group = rep(LETTERS[1:4], each=5), age=1:20, batch=factor(rep(1:5, 4)))
se <- SummarizedExperiment(colData=cd)

test_that("processSimpleComparisons works for two-group comparisons", {
    env <- new.env()
    env$superse <- se
    fun <- processSimpleDesignMatrix("group", NULL, NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=fun), envir=env)

    info <- processSimpleComparisons(c("A", "B"), design.name="FOOBAR", contrast.name="CON")
    expect_identical(length(info), 1L)

    expect_identical(info[[1]]$title, "Increase in `A` over `B`")
    expect_identical(info[[1]]$type, "versus")
    expect_identical(info[[1]]$left, "A")
    expect_identical(info[[1]]$right, "B")

    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(env$CON, c(group.A = 1, group.B = -1, group.C = 0, group.D = 0))
})

test_that("processSimpleComparisons works for covariates", {
    env <- new.env()
    env$superse <- se
    fun <- processSimpleDesignMatrix(NULL, "age", NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=fun), envir=env)

    info <- processSimpleComparisons("age", design.name="FOOBAR", contrast.name="CON")
    expect_identical(length(info), 1L)

    expect_identical(info[[1]]$title, "Effect of increasing `age`")
    expect_identical(info[[1]]$type, "covariate")
    expect_identical(info[[1]]$covariate, "age")

    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(env$CON, c(`(Intercept)` = 0, covariate.age = 1))
})

test_that("processSimpleComparisons works for multi-group comparisons", {
    env <- new.env()
    env$superse <- se
    fun <- processSimpleDesignMatrix("group", NULL, NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=fun), envir=env)

    info <- processSimpleComparisons(c("C", NA, "B", "D"), design.name="FOOBAR", contrast.name="CON")
    expect_identical(length(info), 1L)

    expect_identical(info[[1]]$title, "Increase in `C` over `(B + D)/2`")
    expect_identical(info[[1]]$type, "versus")
    expect_identical(info[[1]]$left, "C")
    expect_identical(info[[1]]$right, c("B", "D"))

    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(env$CON, c(group.A = 0, group.B = -0.5, group.C = 1, group.D = -0.5))
})

test_that("processSimpleComparisons works for ANOVAs", {
    env <- new.env()
    env$superse <- se
    fun <- processSimpleDesignMatrix("group", NULL, NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=fun), envir=env)

    info <- processSimpleComparisons(c("D", "C", "B", "A"), design.name="FOOBAR", contrast.name="CON")
    expect_identical(length(info), 1L)

    expect_identical(info[[1]]$title, "One-way ANOVA involving `D`, `C`, `B`, `A`")
    expect_identical(info[[1]]$type, "anova")
    expect_identical(info[[1]]$groups, c("D", "C", "B", "A"))

    eval(parse(text=info[[1]]$commands), envir=env)
    expect_identical(rownames(env$CON), colnames(env$FOOBAR))
    expect_identical(colnames(env$CON), c("D - A", "C - A", "B - A"))
    expect_identical(env$CON[,1], c(group.A = -1, group.B = 0, group.C = 0, group.D = 1))
    expect_identical(env$CON[,2], c(group.A = -1, group.B = 0, group.C = 1, group.D = 0))
    expect_identical(env$CON[,3], c(group.A = -1, group.B = 1, group.C = 0, group.D = 0))
})

test_that("processSimpleComparisons works for ANOVAs", {
    env <- new.env()
    env$superse <- se
    fun <- processSimpleDesignMatrix("group", NULL, NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=fun), envir=env)

    combined <- processSimpleComparisons(list(foo=c("D", "B"), bar="age"), design.name="FOOBAR", contrast.name="CON")
    expect_identical(combined[[1]]$title, "foo")
    combined[[1]]$title <- NULL
    expect_identical(combined[[2]]$title, "bar")
    combined[[2]]$title <- NULL

    alone1 <- processSimpleComparisons(c("D", "B"), design.name="FOOBAR", contrast.name="CON")
    alone1[[1]]$title <- NULL
    expect_identical(alone1[[1]], combined[[1]])

    alone2 <- processSimpleComparisons("age", design.name="FOOBAR", contrast.name="CON")
    alone2[[1]]$title <- NULL
    expect_identical(alone2[[1]], combined[[2]])
})
