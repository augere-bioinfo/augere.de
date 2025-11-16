#' Differential expression analysis with edgeR
#'
#' Test for differentially expressed (DE) genes from an RNA-seq count matrix using the quasi-likelihood (QL) framework in \pkg{edgeR}. 
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object containing a count matrix where genes and samples are in rows and columns, respectively.
#' Alternatively, the output of \code{\link[augere.core]{wrapInput}} that refers to a SummarizedExperiment.
#' @param groups String specifying the \code{\link[SummarizedExperiment]{colData}(x)} column containing the grouping factor of interest,
#' see \code{\link{processSimpleDesignMatrix}} for more details.
#' This may be \code{NULL} for experimental designs with no groups, e.g., \code{covariates} only.
#' Ignored if \code{design} and \code{contrasts} are provided.
#' @param comparisons Character vector specifying two or more groups to compare from \code{groups}, or the covariate to be tested from \code{covariates}.
#' Alternatively, or a list of such character vectors specifying multiple comparisons to perform.
#' See \code{\link{processSimpleComparisons}} for more details.
#' Ignored if \code{design} and \code{contrasts} are provided.
#' @param covariates Character vector specifying the \code{\link[SummarizedExperiment]{colData}(x)} columns containing continuous covariates of interest,
#' see \code{\link{processSimpleDesignMatrix}} for more details.
#' Ignored if \code{design} and \code{contrasts} are provided.
#' @param block Character vector specifying the \code{\link[SummarizedExperiment]{colData}(x)} columns containing additional (uninteresting) blocking factors,
#' see \code{\link{processSimpleDesignMatrix}} for more details.
#' Ignored if \code{design} and \code{contrasts} are provided.
#' @param assay String or integer specifying the assay of \code{x} containing the count matrix.
#' @param row.data Character vector specifying the \code{\link[SummarizedExperiment]{rowData}(x)} columns containing extra gene annotations to include in the DE result data frames.
#' @param subset.factor String specifying the \code{\link[SummarizedExperiment]{colData}(x)} column containing the factor to use for subsetting.
#' @param subset.levels Vector containing the levels of the \code{subset.factor} to be retained.
#' @param subset.groups Boolean indicating whether to automatically subset the dataset to only those samples assigned to groups in \code{comparisons}.
#' Setting this to \code{TRUE} sacrifices some residual degrees of freedom for greater robustness against variability in irrelevant groups.
#' Ignored if \code{design} and \code{contrasts} are provided.
#' Also ignored if \code{covariates} is provided, as all samples are informative for a continuous covariate.
#' @param robust Boolean indicating whether robust empirical Bayes shrinkage should be used in \code{\link[edgeR]{glmQLFit}}. 
#' Setting this to \code{TRUE} sacrifices some precision for improved robustness against genes with extreme dispersions.
#' @param trend Boolean indicating whether to shrink the QL dispersions towards a trend fitted to the mean in \code{\link[edgeR]{glmQLFit}}.
#' Setting this to \code{TRUE} allows the analysis to adapt to difference mean-variance relationships, at the cost of some extra computational work.
#' @param lfc.threshold Number specifying a log-fold change threshold for \code{\link[edgeR]{glmTreat}}.
#' If zero, \code{\link[edgeR]{glmQLFTest}} is used instead.
#' This should only be used for contrasts that can be formulated in terms of a single coefficient.
#' @param design Matrix, function, or formula specifying the experimental design, see \code{\link{processCustomDesignMatrix}} for details.
#' If this and \code{contrasts} are specified, \code{groups}, \code{block}, \code{covariates}, \code{comparisons} and \code{subset.groups} are ignored.
#' @param contrasts String, function or matrix specifying a custom contrast, or a list of such objects; see \code{\link{processCustomContrasts}} for details.
#' If this and \code{design} are specified, \code{groups}, \code{block}, \code{covariates}, \code{comparisons} and \code{subset.groups} are ignored.
#' @param metadata Named list of additional metadata to store alongside each result.
#' @param output.dir String containing the path to an output directory in which to write the Rmarkdown file and save results.
#' @param author Character vector containg the names of the authors.
#' @param dry.run Boolean indicating whether to perform a dry run.
#' This generates the Rmarkdown report in \code{output.dir} but does not execute the analysis.
#' @param save.results Boolean indicating whether the results should be saved to file.
#'
#' @return
#' A Rmarkdown report named \code{report.Rmd} is written inside \code{output.dir} that contains the analysis commands.
#'
#' If \code{dry.run=FALSE}, a list is returned containing:
#' \itemize{
#' \item \code{results}, a list of \link[S4Vectors]{DataFrame}s of tables from all contrasts;
#' \item \code{normalized}, a \link[SummarizedExperiment]{RangedSummarizedExperiment} with normalized expression values (possibly subsetted by sample).
#' }
#' If \code{save.results=TRUE}, the results are saved in a \code{results} directory inside \code{output}.
#'
#' If \code{dry.run=TRUE}, \code{NULL} is returned.
#' Only the Rmarkdown report is saved to file.
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- loadExampleDataset()
#' 
#' tmp <- tempfile()
#' out <- runEdgeR(
#'     x, 
#'     groups="dex", 
#'     comparisons=c("trt", "untrt"),
#'     output=tmp
#' )
#'
#' list.files(tmp, recursive=TRUE)
#' out
#'
#' @export
#' @import augere.core
runEdgeR <- function(
    x,
    groups, 
    comparisons, 
    covariates=NULL,
    block=NULL, 
    subset.factor=NULL, 
    subset.levels=NULL, 
    subset.groups=TRUE,
    design=NULL, 
    contrasts=NULL, 
    robust=TRUE, 
    trend=TRUE,
    lfc.threshold=0, 
    assay=1, 
    row.data=NULL, 
    metadata=NULL,
    output.dir="edgeR", 
    author=NULL,
    dry.run=FALSE, 
    save.results=TRUE
) {
    restore.fun <- resetInputCache()
    on.exit(restore.fun(), after=FALSE, add=TRUE)

    if (is.null(author)) {
        author <- Sys.info()[["user"]]
    }

    dir.create(output.dir, showWarnings=FALSE, recursive=FALSE)
    fname <- file.path(output.dir, "report.Rmd")

    common.start <- .initialize(
        x=x,
        method="edgeR",
        assay=assay,
        groups=groups,
        comparisons=comparisons,
        covariates=covariates,
        block=block,
        subset.factor=subset.factor,
        subset.levels=subset.levels,
        subset.groups=subset.groups,
        design=design,
        contrasts=contrasts,
        author=author
    )
    writeRmd(common.start$text, file=fname)
    contrast.info <- common.start$contrasts

    template <- system.file("templates", "edgeR.Rmd", package="augere.de", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))
    replacements <- list()

    if (robust) {
        replacements$QL_OPTS <- ", robust=TRUE"
    } else {
        replacements$QL_OPTS <- ""
        parsed[["robust-eb-text"]] <- NULL
    }
    if (!trend) {
        replacements$QL_OPTS <- paste0(replacements$QL_OPTS, ", abundance.trend=FALSE")
    }

    merge.metadata <- !is.null(metadata)
    if (merge.metadata) {
        parsed[["create-common-metadata"]] <- replacePlaceholders(parsed[["create-common-metadata"]], list(COMMON_METADATA=deparseToString(metadata)))
    } else {
        parsed[["create-common-metadata"]] <- NULL
    }

    contrasts <- vector("list", length(contrast.info))
    save.names <- character(length(contrast.info))
    author.txt <- deparseToString(as.list(author))

    for (i in seq_along(contrast.info)) {
        copy <- parsed$contrast
        current <- contrast.info[[i]]
        copy[["setup-contrast"]] <- current$commands

        if (lfc.threshold == 0) {
            copy[["treat"]] <- NULL
        } else {
            copy[["regular"]] <- NULL
            copy[["treat"]] <- replacePlaceholders(copy[["treat"]], list(THRESHOLD=deparseToString(lfc.threshold)))
        }

        if (!is.null(row.data)) {
            copy[["attach-rowdata"]] <- sprintf("de.df <- cbind(SummarizedExperiment::rowData(se)[,%s,drop=FALSE], de.df)", deparseToString(row.data))
        }

        meta.cmds <- processContrastMetadata(current)
        meta.cmds[-1] <- paste0("    ", meta.cmds[-1])
        meta.cmds[1] <- paste0("    contrast=", meta.cmds[1])
        meta.cmds[length(meta.cmds)] <- paste0(meta.cmds[length(meta.cmds)], ",")
        copy[["diff-metadata"]] <- meta.cmds

        if (!is.null(subset.factor)) {
            copy[["subset-metadata"]] <- paste0(strrep(" ", 12), "subset=subset.meta,")
        }

        if (!merge.metadata) {
            copy[["merge-metadata"]] <- NULL
        } else {
            copy[["no-merge-metadata"]] <- NULL
        }

        save.name <- paste0("save-de", i)
        save.names[i] <- save.name
        contrasts[[i]] <- replacePlaceholders(
            copy,
            list(
                AUTHOR=author.txt,
                CONTRAST_NAME_SIMPLE=current$title,
                CONTRAST_NAME_DEPARSED=deparseToString(current$title),
                SAVING_CHUNK_NAME=save.name
            )
        )
    }

    parsed$contrast <- contrasts
    parsed <- replacePlaceholders(parsed, replacements)
    writeRmd(parsed, file=fname, append=TRUE)

    finish <- .finalize(assay, author, merge.metadata)
    writeRmd(finish, file=fname, append=TRUE)

    if (dry.run) {
        return(NULL)
    }

    if (save.results) {
        skip.chunks <- NULL
    } else {
        skip.chunks <- c('save-directory', save.names, "save-norm")
    }
    env <- new.env()
    compileReport(fname, env=env, skip.chunks=skip.chunks)

    list(results=env$all.results, normalized=env$norm.se)
}
