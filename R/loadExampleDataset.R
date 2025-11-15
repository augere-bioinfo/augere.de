#' Load example dataset
#'
#' Load an example RNA-seq dataset from the \pkg{airway} package.
#' This is a convenience wrapper that strips out some unnecesssary metadata. 
#' 
#' @return A \link[SummarizedExperiment]{RangedSummarizedExperiment} containing the \pkg{airway} dataset.
#'
#' @author Aaron Lun
#' 
#' @examples
#' loadExampleDataset()
#' 
#' @export
loadExampleDataset <- function() {
    if (is.null(cached$object)) {
        env <- new.env()
        data("airway", package="airway", envir=env)
        se <- env$airway
        S4Vectors::metadata(se) <- list()
        cached$object <- se
    }
    cached$object
}

cached <- new.env()
cached$object <- NULL
