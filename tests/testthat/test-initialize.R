# library(testthat); library(augere.de); source("test-initialize.R")

library(SummarizedExperiment)
cd <- DataFrame(group = rep(LETTERS[1:4], each=5), age=1:20, batch=factor(rep(1:5, 4)))
means <- 2^rexp(1000)
means[100] <- 0 # check that filtering behaves correctly.
mat <- matrix(rnbinom(length(means) * nrow(cd), mu=means, size=10), nrow=length(means), ncol=nrow(cd))
se <- SummarizedExperiment(list(counts=mat), colData=cd)

test_that(".initialize() works with 'default' options", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.de:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=c("A", "B"),
        covariates=NULL,
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=TRUE,
        design=NULL,
        contrasts=NULL,
        author="Chihaya Kisaragi"
    )

    # Checking that the contrasts are okay.
    expect_identical(output$contrasts[[1]]$title, "Increase in `A` over `B`")
    expect_identical(output$contrasts[[1]]$type, "versus")
    expect_identical(output$contrasts[[1]]$left, "A")
    expect_identical(output$contrasts[[0]]$right, "B")

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)

    expect_identical(env$se, se[,se$group %in% c("A", "B")])
    expect_identical(colnames(env$design), c("group.A", "group.B"))

    expect_true(inherits(env$y, "DGEList"))

    expect_lt(nrow(env$y), nrow(se)) # make sure it got some filtering.
    expect_type(env$y$genes$origin, "integer") 
    expect_false(100 %in% env$y$genes$origin) 
    expect_true(any(grepl("filterByExpr\\(.*group=", unlist(output$text))))

    expect_type(env$y$samples$norm.factors, "double") # make sure we got some normalization.
    expect_true(is.factor(env$mds.labels))
})

test_that(".initialize() works with no filtering by group", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.de:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=c("A", "B"),
        covariates=NULL,
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=FALSE,
        design=NULL,
        contrasts=NULL,
        author="Chihaya Kisaragi"
    )

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)

    expect_identical(env$se, se)
})

test_that(".initialize() works with custom sample filtering", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.de:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=c("A", "B"),
        covariates=NULL,
        block="batch",
        subset.factor="batch",
        subset.levels=c(2,4),
        subset.groups=TRUE,
        design=NULL,
        contrasts=NULL,
        author="Chihaya Kisaragi"
    )

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)

    expect_identical(env$se, se[,se$group %in% c("A", "B") & se$batch %in% c(2, 4)])
})

test_that(".initialize() works for a design matrix without groups", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.de:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups=NULL,
        comparisons=c("age"),
        covariates="age",
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=TRUE,
        design=NULL,
        contrasts=NULL,
        author="Chihaya Kisaragi"
    )

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)

    expect_identical(env$se, se) # no filtering is done on the samples.
    expect_true(any(grepl("filterByExpr\\(.*design=", unlist(output$text))))
    expect_null(env$mds.labels)
})

test_that(".find_used_groups() works as expected", {
    contrast.info <- processSimpleComparisons(c("A", "B"))
    expect_identical(augere.de:::.find_used_groups(contrast.info), c("A", "B"))

    contrast.info <- processSimpleComparisons(list(c("A", "B"), c("B", "C", "D")))
    expect_identical(augere.de:::.find_used_groups(contrast.info), c("A", "B", "C", "D"))

    contrast.info <- processSimpleComparisons(list(c("age"), c("C", "D")))
    expect_null(augere.de:::.find_used_groups(contrast.info))

    contrast.info <- processSimpleComparisons(c("C", NA, "A", "D"))
    expect_identical(augere.de:::.find_used_groups(contrast.info), c("A", "C", "D"))

    contrast.info <- processSimpleComparisons("foo")
    expect_null(augere.de:::.find_used_groups(contrast.info))
})

test_that(".initialize() works with custom matrices and contrasts", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.de:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=NULL,
        covariates=NULL,
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=TRUE,
        design=~group + batch,
        contrasts="groupC - groupD",
        author="Chihaya Kisaragi"
    )

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)

    expect_identical(env$se, se) # no filtering is done on the samples.
    expect_true(any(grepl("filterByExpr\\(.*design=", unlist(output$text))))
    expect_null(env$mds.labels)

    # Checking that the contrasts are okay.
    expect_identical(output$contrasts[[1]]$type, "custom")
})
