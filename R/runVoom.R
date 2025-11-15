#' Differential expression analysis with voom
#'
#' Test for differentially expressed (DE) genes from an RNA-seq count matrix using the voom algorithm from \pkg{limma}.
#'
#' @inheritParams runEdgeR
#' @param trend Boolean indicating whether variances should be shrunk towards a trend in \code{\link[limma]{eBayes}}.
#' Usually unnecessary as the observation weights already account for the mean-variance relationship. 
#' @param quality Boolean indicating whether quality weighting should be performed.
#' This reduces the influence of low-quality samples at the cost of more computational work. 
#' @param dc.block String specifying the blocking factor to use in \code{\link[limma]{duplicateCorrelation}}.
#' Typically used for uninteresting factors that cannot be used in \code{block} as they are confounded with the factors of interest. 
#' No additional blocking is performed if \code{NULL}.
#' @param robust Boolean indicating whether robust empirical Bayes shrinkage should be used in \code{\link[limma]{eBayes}}.
#' Setting this to \code{TRUE} sacrifices some precision for improved robustness against genes with extreme variances.
#' @param lfc.threshold Number specifying a threshold on the log-fold change in \code{\link[limma]{treat}}.
#' If zero, \code{\link[limma]{eBayes}} will be used instead.
#'
#' @return
#' A Rmarkdown report named \code{report.Rmd} is written inside \code{output.dir} that contains the analysis commands.
#'
#' If \code{dry.run=FALSE}, a list is returned containing:
#' \itemize{
#' \item \code{results}, a list of \link[S4Vectors]{DataFrame}s of tables from all contrasts;
#' \item \code{normalized}, a \link[SummarizedExperiment]{RangedSummarizedExperiment} with normalized expression values (possibly subsetted by sample).
#' }
#' If \code{save.all=TRUE}, the results are saved in a \code{results} directory inside \code{output}.
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
#' out <- runVoom(
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
runVoom <- function(
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
    dc.block=NULL,
    robust=TRUE, 
    trend=FALSE,
    quality=TRUE,
    lfc.threshold=0, 
    assay=1, 
    row.data=NULL, 
    metadata=NULL,
    output.dir="voom", 
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
        method="voom",
        se=x,
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

    template <- system.file("templates", "voom.Rmd", package="augere.de", mustWork=TRUE)
    parsed <- extractBlocks(readLines(template))
    replacements <- list(AUTHOR = deparseToString(as.list(author)))

    if (!is.null(dup.cor)) {
        parsed[["voom"]] <- NULL
        replacements$DUPCOR_BLOCK <- deparseToString(dup.cor)
        replacements$LM_OPTS <- ", block=dc.block, correlation=dc$consensus.correlation"
    } else {
        parsed[["duplicate-correlation"]] <- NULL
        replacements$LM_OPTS <- ""
    }

    if (quality) {
        replacements$VOOM_CMD <- "voomWithQualityWeights"
    } else {
        replacements$VOOM_CMD <- "voom"
        parsed[["quality-text"]] <- NULL
    }

    eb.args <- ""
    if (trend) {
        eb.args <- c(eb.args, "trend=TRUE")
    }

    if (robust) {
        eb.args <- c(eb.args, "robust=TRUE")
        replacements$EXTRA_EB_CAPT <- " with outliers marked in red"
    } else {
        replacements$EXTRA_EB_CAPT <- ""
        parsed[["robust-text"]] <- NULL
    }

    replacements$EB_OPTS <- paste(eb.args, collapse=", ")
    middle <- replacePlaceholders(parsed, replacements)

    merge.metadata <- !is.null(metadata)
    if (merge.metadata) {
        parsed[["create-common-metadata"]] <- replacePlaceholders(parsed[["merge-metadata-common"]], list(COMMON_METADATA=deparseToString(metadata)))
    } else {
        parsed[["create-common-metadata"]] <- NULL
    }

    contrasts <- vector("list", length(contrast.info))
    save.names <- character(length(contrast.info))
    for (i in seq_along(contrast.info)) {
        copy <- parsed$contrast
        current <- contrast.info[[i]]
        copy[["contrast-command"]] <- current$commands

        if (lfc.threshold == 0) {
            copy[["lfc-contrast"]] <- NULL
        } else {
            copy[["no-lfc-contrast"]] <- NULL
            copy[["lfc-contrast"]] <- replacePlaceholders(copy[["lfc-contrast"]], THRESHOLD=deparseToString(lfc.threshold))
        }

        if (!is.null(row.data)) {
            copy[["attach-annotation"]] <- sprintf("deres <- cbind(deres, rowData(se)[,%s,drop=FALSE])", deparseToString(row.data))
        }

        copy[["diff-metadata"]] <- formatContrastMetadata(current, groups=groups, indent=3)
        if (!is.null(subset.factor)) {
            copy[["subset-metadata"]] <- paste0(strrep(" ", 12), "subset=subset.metadata,")
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

    if (save.all) {
        skip.chunks <- NULL
    } else {
        skip.chunks <- save.names
    }
    env <- new.env()
    compileReport(fname, envir=env, skip.chunks=skip.chunks)

    list(results=env$all.results, normalized=env$processed.se)
}
