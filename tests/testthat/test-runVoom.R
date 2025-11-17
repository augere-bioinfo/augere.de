# library(augere.de); library(testthat); source("test-runVoom.R")

library(SummarizedExperiment)
cd <- DataFrame(group = rep(LETTERS[1:4], each=5), age=1:20, batch=factor(rep(1:5, 4)))
means <- 2^rexp(1000)
means[100] <- 0 # check that filtering behaves correctly.
mat <- matrix(rnbinom(length(means) * nrow(cd), mu=means, size=10), nrow=length(means), ncol=nrow(cd))
se <- SummarizedExperiment(list(counts=mat), colData=cd)
rownames(se) <- sprintf("GENE_%s", seq_len(nrow(se)))

test_that("runVoom works with vanilla settings", {
    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=c("C", "A"), output.dir=tmp)

    expect_s4_class(out$normalized, "SummarizedExperiment")
    expect_true("logCPM" %in% assayNames(out$normalized))
    expect_type(out$normalized$norm.factors, "double")

    filtered <- rowData(out$normalized)$retained
    expect_true(any(filtered))
    expect_true(any(!filtered))

    expect_identical(names(out$results), "Increase in `C` over `A`")
    expect_s4_class(out$results[[1]], "DFrame") 
    expect_identical(rownames(out$results[[1]]), rownames(se))
    expect_false(anyNA(out$results[[1]]$FDR[rowData(out$normalized)$retained]))
    expect_true(all(is.na(out$results[[1]]$FDR[!rowData(out$normalized)$retained])))

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "de-1", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$differential_gene_expression$contrast$left, list("C"))
    expect_identical(meta$differential_gene_expression$contrast$right, list("A"))
    expect_null(meta$differential_gene_expression$subset) # we shouldn't have any subsetting here.

    # Check that the default options are set here.
    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_false(any(grepl("treat\\(.*lfc=0.5", fname)))
    expect_true(any(grepl("voomWithQualityWeights(", fname, fixed=TRUE)))
    expect_false(any(grepl("duplicateCorrelation(", fname, fixed=TRUE)))
    expect_false(any(grepl("voom(", fname, fixed=TRUE)))
    expect_false(any(grepl("trend=TRUE", fname, fixed=TRUE)))
    expect_true(any(grepl("robust=TRUE", fname, fixed=TRUE)))
})

test_that("runVoom works with multiple comparisons", {
    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=list(foo=c("C", "A"), bar=c("B", "D")), output.dir=tmp)
    expect_identical(names(out$results), c("foo", "bar"))

    tmp1 <- tempfile()
    alone1 <- runVoom(se, group="group", comparisons=c("C", "A"), subset.group=FALSE, output.dir=tmp1)
    expect_identical(out$results$foo, alone1$results[[1]])

    tmp2 <- tempfile()
    alone2 <- runVoom(se, group="group", comparisons=c("B", "D"), subset.group=FALSE, output.dir=tmp2)
    expect_identical(out$results$bar, alone2$results[[1]])
})

test_that("runVoom works with log-fold changes", {
    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=c("C", "A"), lfc.threshold=0.5, output.dir=tmp)

    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_true(any(grepl("treat\\(.*lfc=0.5", fname)))
})

test_that("runVoom works with custom subsetting", {
    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=c("C", "A"), subset.factor="batch", subset.levels=c("1","3","5"), output.dir=tmp)

    # Subsetting works correctly.
    expect_identical(unique(out$normalized$batch), factor(c("1", "3", "5"), 1:5))

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "de-1", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$differential_gene_expression$subset, "batch IN (1,3,5)")
})

test_that("runVoom works with more annotations", {
    tmp <- tempfile()
    rowData(se)$symbol <- sprintf("FOOBAR_%s", seq_len(nrow(se)))
    out <- runVoom(se, group="group", comparisons=c("C", "A"), row.data="symbol", output.dir=tmp)

    expect_identical(out$results[[1]]$symbol, rowData(se)$symbol)
    expect_type(out$results[[1]]$FDR, "double") # make sure we didn't wipe out the other columns.

    expect_identical(rowData(out$normalized)$symbol, rowData(se)$symbol)
})

test_that("runVoom works with custom metadata", {
    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=c("C", "A"), metadata=list(X=1, Y=TRUE), output.dir=tmp)

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "de-1", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$X, 1L)
    expect_true(meta$Y)

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "normalized", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$X, 1L)
    expect_true(meta$Y)
})

test_that("runVoom works with duplicateCorrelation", {
    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=c("C", "A"), dc.block="batch", output.dir=tmp)

    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_true(any(grepl("duplicateCorrelation(", fname, fixed=TRUE)))
})

test_that("runVoom works with some of the other voom-related options", {
    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=c("C", "A"), robust=FALSE, trend=TRUE, quality=FALSE, output.dir=tmp)

    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_false(any(grepl("voomWithQualityWeights(", fname, fixed=TRUE)))
    expect_true(any(grepl("voom(", fname, fixed=TRUE)))
    expect_true(any(grepl("trend=TRUE", fname, fixed=TRUE)))
    expect_false(any(grepl("robust=TRUE", fname, fixed=TRUE)))
})

test_that("runVoom works with dry runs and no saving", {
    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=c("C", "A"), dry.run=TRUE, output.dir=tmp)
    expect_null(out)
    expect_true(file.exists(file.path(tmp, "report.Rmd")))

    tmp <- tempfile()
    out <- runVoom(se, group="group", comparisons=c("C", "A"), save.results=FALSE, output.dir=tmp)
    expect_s4_class(out$normalized, "SummarizedExperiment")
    expect_s4_class(out$results[[1]], "DFrame")
    expect_false(file.exists(file.path(tmp, "results")))
})
