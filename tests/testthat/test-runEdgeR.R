# library(augere.de); library(testthat); source("test-runEdgeR.R")

library(SummarizedExperiment)
cd <- DataFrame(group = rep(LETTERS[1:4], each=5), age=1:20, batch=factor(rep(1:5, 4)))
means <- 2^rexp(1000)
means[100] <- 0 # check that filtering behaves correctly.
mat <- matrix(rnbinom(length(means) * nrow(cd), mu=means, size=10), nrow=length(means), ncol=nrow(cd))
se <- SummarizedExperiment(list(counts=mat), colData=cd)
rownames(se) <- sprintf("GENE_%s", seq_len(nrow(se)))

test_that("runEdgeR works with vanilla settings", {
    tmp <- tempfile()
    out <- runEdgeR(se, group="group", comparisons=c("C", "A"), output.dir=tmp)

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

    roundtrip <- augere.core::readResult(file.path(tmp, "results", "de-1"))
    expect_identical(roundtrip$x, out$results[[1]])
    meta <- roundtrip$metadata
    expect_identical(meta$differential_gene_expression$contrast$left, list("C"))
    expect_identical(meta$differential_gene_expression$contrast$right, list("A"))
    expect_null(meta$differential_gene_expression$subset) # we shouldn't have any subsetting here.

    # Check that the default options are set here.
    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_false(any(grepl("glmTreat\\(.*lfc=0.5", fname)))
    expect_true(any(grepl("glmQLFTest\\(", fname)))
    expect_false(any(grepl("abundance.trend=FALSE", fname, fixed=TRUE)))
    expect_true(any(grepl("robust=TRUE", fname, fixed=TRUE)))
})

test_that("runEdgeR works with multiple comparisons", {
    tmp <- tempfile()
    out <- runEdgeR(se, group="group", comparisons=list(foo=c("C", "A"), bar=c("B", "D")), output.dir=tmp)
    expect_identical(names(out$results), c("foo", "bar"))
    expect_match(augere.core::readResult(file.path(tmp, "results", "de-1"))$meta$title, "foo")
    expect_match(augere.core::readResult(file.path(tmp, "results", "de-2"))$meta$title, "bar")

    tmp1 <- tempfile()
    alone1 <- runEdgeR(se, group="group", comparisons=c("C", "A"), subset.group=FALSE, save.results=FALSE, output.dir=tmp1)
    expect_identical(out$results$foo, alone1$results[[1]])

    tmp2 <- tempfile()
    alone2 <- runEdgeR(se, group="group", comparisons=c("B", "D"), subset.group=FALSE, save.results=FALSE, output.dir=tmp2)
    expect_identical(out$results$bar, alone2$results[[1]])
})

test_that("runEdgeR works with custom design and contrasts", {
    tmp <- tempfile()
    custom <- runEdgeR(se, design=~0 + group, contrasts="groupB - groupA", output.dir=tmp)
    meta <- augere.de::readResult(file.path(tmp, "results", "de-1"))$metadata
    expect_identical(meta$differential_gene_expression$contrast$type, "custom")

    # Making sure we get the same results with a simple design.
    tmp <- tempfile()
    comp <- runEdgeR(se, group="group", comparisons=c("B", "A"), subset.group=FALSE, save.results=FALSE, output.dir=tmp)
    expect_identical(custom$normalized, comp$normalized)
    expect_identical(unname(custom$results), unname(comp$results))
})

test_that("runEdgeR works with log-fold changes", {
    tmp <- tempfile()
    out <- runEdgeR(se, group="group", comparisons=c("C", "A"), lfc.threshold=0.5, save.results=FALSE, output.dir=tmp)

    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_true(any(grepl("glmTreat\\(.*lfc=0.5", fname)))
    expect_false(any(grepl("glmQLFTest\\(", fname)))
})

test_that("runEdgeR works with custom subsetting", {
    tmp <- tempfile()
    out <- runEdgeR(se, group="group", comparisons=c("C", "A"), subset.factor="batch", subset.levels=c("1","3","5"), output.dir=tmp)

    # Subsetting works correctly.
    expect_lt(ncol(out$normalized), ncol(se))
    expect_identical(unique(out$normalized$batch), factor(c("1", "3", "5"), 1:5))

    meta <- augere.core::readResult(file.path(tmp, "results", "de-1"))$metadata
    expect_identical(meta$differential_gene_expression$subset, "batch IN (1,3,5)")
})

test_that("runEdgeR works with more annotations", {
    tmp <- tempfile()
    rowData(se)$symbol <- sprintf("FOOBAR_%s", seq_len(nrow(se)))
    out <- runEdgeR(se, group="group", comparisons=c("C", "A"), row.data="symbol", save.results=FALSE, output.dir=tmp)

    expect_identical(out$results[[1]]$symbol, rowData(se)$symbol)
    expect_type(out$results[[1]]$FDR, "double") # make sure we didn't wipe out the other columns.

    expect_identical(rowData(out$normalized)$symbol, rowData(se)$symbol)
})

test_that("runEdgeR works with custom metadata", {
    tmp <- tempfile()
    out <- runEdgeR(se, group="group", comparisons=c("C", "A"), metadata=list(X=1, Y=TRUE), output.dir=tmp)

    meta <- augere.core::readResult(file.path(tmp, "results", "de-1"))$metadata
    expect_identical(meta$X, 1L)
    expect_true(meta$Y)

    meta <- augere.core::readResult(file.path(tmp, "results", "normalized"))$metadata
    expect_identical(meta$X, 1L)
    expect_true(meta$Y)
})

test_that("runEdgeR works with some of the other edgeR-related options", {
    tmp <- tempfile()
    out <- runEdgeR(se, group="group", comparisons=c("C", "A"), robust=FALSE, trend=FALSE, save.results=FALSE, output.dir=tmp)

    fname <- readLines(file.path(tmp, "report.Rmd"))
    expect_true(any(grepl("abundance.trend=FALSE", fname, fixed=TRUE)))
    expect_false(any(grepl("robust=TRUE", fname, fixed=TRUE)))
})

test_that("runEdgeR works with dry runs and no saving", {
    tmp <- tempfile()
    out <- runEdgeR(se, group="group", comparisons=c("C", "A"), dry.run=TRUE, output.dir=tmp)
    expect_null(out)
    expect_true(file.exists(file.path(tmp, "report.Rmd")))

    tmp <- tempfile()
    out <- runEdgeR(se, group="group", comparisons=c("C", "A"), save.results=FALSE, output.dir=tmp)
    expect_s4_class(out$normalized, "SummarizedExperiment")
    expect_s4_class(out$results[[1]], "DFrame")
    expect_false(file.exists(file.path(tmp, "results")))
})
