# library(testthat); library(augere.de); source("test-processContrastMetadata.R")

library(SummarizedExperiment)
cd <- DataFrame(group = rep(LETTERS[1:4], each=5), age=1:20, batch=factor(rep(1:5, 4)))
se <- SummarizedExperiment(colData=cd)

test_that("processContrastMetadata works for simple comparisons", {
    env <- new.env()
    env$superse <- se
    fun <- processSimpleDesignMatrix("group", "age", NULL, se.name="superse", design.name="FOOBAR")
    eval(parse(text=fun), envir=env)

    info <- processSimpleComparisons(
        list(
            foo=c("A", "B"),
            bar=c("A", "B", "C"),
            stuff="age"
        ),
        design.name="FOOBAR",
        contrast.name="CON"
    )

    for (i in seq_along(info)) {
        cmds <- processContrastMetadata(info[[i]])
        output <- eval(parse(text=cmds), envir=env)
        expected <- info[[i]]
        expected$commands <- NULL
        expect_identical(output, expected)
    }
})

test_that("processContrastMetadata works for custom contrasts", {
    env <- new.env()
    env$superse <- se
    fun <- processCustomDesignMatrix(~0 + group, se.name="superse", design.name="FOOBAR")
    eval(parse(text=fun), envir=env)

    info <- processCustomContrasts(
        list(
            "groupA - groupB"
        ),
        design.name="FOOBAR",
        contrast.name="CON"
    )

    cmds <- processContrastMetadata(info[[1]])
    output <- eval(parse(text=cmds), envir=env)
    expected <- info[[1]]
    expected$commands <- NULL
    expect_identical(output, expected)
})
