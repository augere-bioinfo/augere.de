#' Create a custom design matrix
#'
#' Construct R commands to create a custom design matrix, to be inserted into the generated Rmarkdown report.
#'
#' @param design Function, formula or matrix specifying the experimental design.
#' @param design.name String containing the variable name of the design matrix.
#' @param se.name String containing the variable name of the \link[SummarizedExperiment]{SummarizedExperiment} object.
#'
#' @return 
#' Character vector containing R commands that, upon evaluation, create a design matrix in the evaluation environment.
#'
#' @author Aaron Lun
#'
#' @details
#' If \code{design} is a function, it should accept a \link[SummarizedExperiment]{SummarizedExperiment} and return the matrix.
#'
#' If \code{design} is a formula, it should use the column names of the SummarizedExperiment's \code{\link[SummarizedExperiment]{colData}}.
#' This is passed to \code{\link{model.matrix}} using \code{data=colData(se)} for a SummarizedExperiment \code{se}.
#'
#' If \code{design} is a numeric matrix, it is deparsed and used verbatim.
#' The number of rows of this matrix should be equal to the number of samples.
#'
#' The design matrix is expected to be of full column rank. 
#'
#' @examples
#' cat(processCustomDesignMatrix(~ batch + treatment, "se"))
#'
#' cat(processCustomDesignMatrix(function(x) {
#'     model.matrix(~ batch + treatment, data=colData(x))
#' }, "se"))
#' 
#' batch <- factor(rep(1:3, each=2))
#' treatment <- rep(c("Drg", "Ctrl"), 3)
#' mat <- model.matrix(~ batch + treatment)
#' cat(processCustomDesignMatrix(mat, "se"), sep="\n")
#'
#' @seealso
#' \code{\link{processSimpleDesignMatrix}}, for creating simple design matrices.
#'
#' \code{\link{processCustomContrasts}}, to create custom contrasts. 
#'
#' @export
#' @import augere.core
#' @import methods
processCustomDesignMatrix <- function(design, se.name, design.name = "design") {
    if (is.function(design)) {
        cmd <- sprintf("(%s)(%s)", deparseToString(design, indent=0), se.name)
    } else if (is(design, "formula")) {
        cmd <- sprintf("model.matrix(%s, data=colData(%s))", deparseToString(design), se.name)
    } else {
        cmd <- deparseToString(design, indent=4)
    }
    sprintf("%s <- %s", design.name, cmd)
}
