#' Process contrast metadata
#'
#' Create R commands to define the contrast metadata for the results of \code{\link{runEdgeR}} and \code{\link{runVoom}}.
#' 
#' @param info List containing information for a single contrast.
#' This is typically an entry of the list returned by \code{\link{processSimpleComparisons}} or \code{\link{processCustomContrasts}}.
#'
#' @return Character vector containing arguments for constructing the \code{gpsa.differential_gene_expression} sublist of the result metadata.
#'
#' @author Aaron Lun
#' @examples
#' contrast.info <- processSimpleComparisons(list(
#'     main=c("disease", "healthy"),
#'     secondary="dosage"
#' ))
#' cat(processContrastMetadata(contrast.info[[1]]), sep="\n")
#' cat(processContrastMetadata(contrast.info[[2]]), sep="\n")
#'
#' @export
processContrastMetadata <- function(info) {
    first <- TRUE
    payload <- character(0)
    for (field in setdiff(names(info), c("commands", "title"))) {
        if (!first) {
            payload[length(payload)] <- paste0(payload[length(payload)], ",")
        } else {
            first <- FALSE
        }
        payload <- c(payload, sprintf("    %s=%s", field, deparseToString(info[[field]])))
    }
    c("list(", payload, ")")
}
