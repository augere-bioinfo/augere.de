#' Create a simple design matrix
#'
#' Construct R commands to create a simple design matrix, to be inserted into the generated Rmarkdown report.
#'
#' @param groups String specifying the \code{\link[SummarizedExperiment]{colData}} column containing the grouping factor.
#' @param covariates Character vector specifying the \code{colData} columns containing continuous covariates.
#' @param block Character vector specifying \code{colData} columns containing additional uninteresting factrors.
#' @param se.name String containing the variable name of the \link[SummarizedExperiment]{SummarizedExperiment} object.
#' @param design.name String containing the variable name of the design matrix.
#'
#' @return 
#' Character vector containing R commands that, upon evaluation, create a design matrix in the evaluation environment.
#'
#' @author Aaron Lun
#'
#' @details
#' When creating the design matrix, \code{groups}, \code{covariates} and \code{block} are treated as additive factors. 
#' \itemize{
#' \item If \code{groups} is specified, the design matrix will not have an intercept.
#' The first few columns will correspond to the levels of \code{groups} and are named by concatenating \code{groups} with the factor level.
#' If \code{groups} is not supplied, the design matrix will have an intercept in the first column.
#' \item If \code{covariates} are specified, they are represented by the columns after the per-group columns (if \code{groups} if supplied) or the intercept (otherwise).
#' \item All remaining columns will correspond to the various levels of the \code{block} factors.
#' }
#'
#' Syntactically invalid column names and levels can be used for all arguments.
#'
#' The design matrix is expected to be of full column rank, e.g., \code{groups} is not confounded with elements of \code{block}.
#'
#' @examples
#' cat(processSimpleDesignMatrix("my_groups", "age", "batch", "se"), sep="\n")
#'
#' @seealso
#' \code{\link{processCustomDesignMatrix}}, for creating more complex design matrices.
#'
#' \code{\link{processSimpleComparisons}}, to create contrasts for this design matrix.
#'
#' @export
#' @import augere.core
processSimpleDesignMatrix <- function(groups, covariates, block, se.name, design.name = "design") {
    cmds <- c(
        "model.data <- list()",
        sprintf("cd <- SummarizedExperiment::colData(%s)", se.name)
    )

    has.groups <- !is.null(groups)
    if (has.groups) {
        cmds <- c(cmds, sprintf("model.data$group. <- factor(cd[,%s])", deparseToString(groups)))
        terms <- "group."
    } else {
        terms <- character(0)
    }

    if (!is.null(covariates)) {
        for (x in seq_along(covariates)) {
            name <- paste0("covariate.", covariates[x])
            cmds <- c(cmds, sprintf("model.data$%s <- as.numeric(cd[,%s])", name, deparseToString(covariates[x])))
            terms <- c(terms, name)
        }
    }

    if (!is.null(block)) {
        for (b in seq_along(block)) {
            name <- paste0("block.", block[b], ".")
            cmds <- c(cmds, sprintf("model.data$%s <- factor(cd[,%s])", name, deparseToString(block[b])))
            terms <- c(terms, name)
        }
    }

    formula <- paste(terms, collapse=" + ")
    if (has.groups) {
        formula <- paste("~ 0 +", formula)
    } else {
        formula <- paste("~", formula)
    }

    c(cmds, sprintf("%s <- model.matrix(%s, data=model.data)", design.name, formula))
}
