#' Subset by groups in contrasts
#'
#' List all groups involved in contrasts, to use with \code{subset.group=TRUE}.
#' This is mostly intended for use by developers of differential analysis pipelines.
#'
#' @param contrast.info List of contrast information as returned by \code{\link{processSimpleComparisons}} 
#'
#' @return Character vector of the groups involved in \code{"versus"} or \code{"anova"} contrasts.
#' If any covariate-based contrasts are present, \code{NULL} is always returned.
#'
#' @details
#' No subsetting by group is performed if covariate-based contrasts are present,
#' as all samples are potentially informative in such an analysis.
#'
#' @author Aaron Lun
#'
#' @examples
#' findSubsetGroups(processSimpleComparisons(c("disease", "healthy")))
#' findSubsetGroups(processSimpleComparisons(c("treated1", "treated2", "healthy")))
#' findSubsetGroups(processSimpleComparisons(c("age")))
#' 
#' @export
findSubsetGroups <- function(contrast.info) {
    not.grouped <- vapply(contrast.info, function(con) con$type %in% "covariate", TRUE)
    if (any(not.grouped)) {
        NULL
    } else {
        all.groups <- lapply(contrast.info, function(con) con[c("left", "right", "groups")])
        sort(unique(unlist(all.groups)))
    }
}
