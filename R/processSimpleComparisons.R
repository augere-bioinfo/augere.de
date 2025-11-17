#' Contrasts for simple comparisons
#'
#' Construct R commands to define contrast vectors/matrices for simple comparisons between groups or covariates.
#'
#' @param comparisons Character vector specifying the groups to compare or the covariates to test.
#' Non-character vectors will be coerced into character vectors.
#' Alternatively, a (possibly named) list of such vectors where each entry represents a different comparison.
#' @param design.name String containing the variable name of the design matrix.
#' @param contrast.name String containing the variable name of the contrast vector/matrix.
#'
#' @return 
#' A list containing one entry per comparison, where each entry contains:
#' \itemize{
#' \item \code{title}, the title for the comparison.
#' If \code{comparisons} is named, the corresponding (non-empty) name is used here, otherwise an appropriate title is automatically generated.
#' \item \code{type}, the type of the comparison. 
#' This is one of \code{"covariate"} (for covariates), \code{"versus"} (for comparison between two groups or two sets of groups) or \code{"anova"} (for ANOVAs).
#' \item \code{left} and \code{right}, lists of strings specifying the groups on the left (numerator) or right (denominator) of a \dQuote{versus} comparison.
#' Only present if \code{type = "versus"}.
#' \item \code{groups}, list of lists of strings specifying the groups involved in an ANOVA-like comparison.
#' Only present if \code{type = "anova"}.
#' \item \code{covariate}, list of strings specifying the covariate(s) being tested.
#' Only present if \code{type = "covariate"}.
#' \item \code{commands}, character vector of R commands that produce the desired contrast vector/matrix.
#' This assumes that the evaluation environment has a design matrix named \code{design.name}.
#' The newly defined contrast vector/matrix will be stored as \code{contrast.name} in the environment. 
#' }
#' 
#' @details
#' If a vector in \code{comparisons} is unnamed and of length 1, the sole entry is assumed to refer to a covariate in the \code{covariates} from \code{\link{processSimpleDesignMatrix}}.
#' The null hypothesis is that the coefficient of the design matrix with the same name is zero, i.e., the covariate has no effect.
#'
#' If a comparison vector is unnamed and of length 2, the entries represent the names of two groups in the specified \code{groups} factor.
#' The null hypothesis is that there is no differential expression between the two groups.
#' The log-fold change is defined as the first group (left) over the second (right). 
#' 
#' If a comparison vector is unnamed and of length 3 or greater, the null hypothesis is that all of the specified levels of \code{groups} are equal.
#' This is an ANOVA-like contrast where contrasts are formulated with respect to the last level,
#' i.e., for \code{n} coefficients, \code{n-1} log-fold changes are reported representing the differences relative to the last coefficient.
#'
#' If the comparison vector is named, all entries with the same name are assumed to represent a group of coefficients.
#' \itemize{
#' \item If there is only one unique name, all entries of the vector are assumed to refer to entries of \code{covariates}.
#' The null hypothesis is that the average of the coefficients for the specified covariates is equal to zero.
#' \item If there are exactly two unique names, these are assumed to refer to two sets of entries of \code{groups}.
#' The null hypothesis is that the average of the per-group coefficients are equal between the two sets.
#' \item If there are three or more names, these are assumed to refer to multiple sets of entries of \code{groups}.
#' The null hypothesis is that the average of the per-group coefficients are equal across all sets, i.e., an ANOVA-like comparison.
#' }
#'
#' @seealso
#' \code{\link{processSimpleDesignMatrix}}, to generate the corresponding design matrix.
#'
#' \code{\link{processCustomContrasts}}, to define more complex custom contrasts.
#'
#' @examples
#' processSimpleComparisons(c("disease", "healthy"))
#' processSimpleComparisons("dosage")
#' processSimpleComparisons(c("untreated", "treated", "healthy"))
#' processSimpleComparisons(c(treated="treatment1", treated="treatment2", healthy="healthy"))
#' processSimpleComparisons(list(
#'     main=c("disease", "healthy"),
#'     secondary="dosage"
#' ))
#'
#' @export
#' @importFrom utils head
processSimpleComparisons <- function(comparisons, design.name = "design", contrast.name = "con") { 
    if (!is.list(comparisons)) {
        comparisons <- list(comparisons)
    }
    output <- vector("list", length(comparisons))

    for (i in seq_along(comparisons)) {
        current <- comparisons[[i]]

        info <- list(title = NULL)
        if (!is.null(names(comparisons))) {
            nm <- names(comparisons)[i]
            if (nm != "") {
                info$title <- nm
            }
        }

        if (is.null(names(current))) {
            groupings <- as.list(current)
        } else {
            groupings <- split(current, names(current))
        }
        groupings <- lapply(groupings, as.character)

        ngroups <- length(groupings)
        if (ngroups == 0L) {
            stop("vectors in 'comparisons' cannot be empty")

        } else if (ngroups == 1L) {
            info$type <- "covariate"
            info$covariate <- as.list(groupings[[1]])
            if (is.null(info$title)) {
                info$title <- sprintf("Effect of increasing `%s`", .multiple_group_names(groupings[[1]]))
            }

            info$commands <- c(
                sprintf("%s <- numeric(ncol(%s))", contrast.name, design.name),
                sprintf("names(%s) <- colnames(%s)", contrast.name, design.name),
                sprintf("%s[%s] <- %s", contrast.name, deparseToString(sprintf("covariate.%s", groupings[[1]])), .multiple_group_coefs(groupings[[1]]))
            )

        } else if (ngroups == 2L) {
            info$type <- "versus"
            info$left <- as.list(groupings[[1]])
            info$right <- as.list(groupings[[2]])
            if (is.null(info$title)) {
                info$title <- sprintf("Increase in `%s` over `%s`", .multiple_group_names(groupings[[1]]), .multiple_group_names(groupings[[2]]))
            }

            info$commands <- c(
                sprintf("%s <- numeric(ncol(%s))", contrast.name, design.name),
                sprintf("names(%s) <- colnames(%s)", contrast.name, design.name),
                sprintf("%s[%s] <- %s", contrast.name, deparseToString(sprintf("group.%s", groupings[[1]])), .multiple_group_coefs(groupings[[1]])),
                sprintf("%s[%s] <- -%s", contrast.name, deparseToString(sprintf("group.%s", groupings[[2]])), .multiple_group_coefs(groupings[[2]]))
            )

        } else {
            info$type <- "anova"
            info$groups <- unname(lapply(groupings, as.list))
            all.names <- vapply(groupings, .multiple_group_names, "")
            if (is.null(info$title)) {
                concat.names <- paste(sprintf("`%s`", all.names), collapse=", ")
                info$title <- sprintf("One-way ANOVA involving %s", concat.names)
            }

            commands <- c(
                sprintf("%s <- matrix(0, ncol(%s), %s)", contrast.name, design.name, ngroups - 1),
                sprintf("rownames(%s) <- colnames(%s)", contrast.name, design.name)
            )
            for (j in seq_len(ngroups - 1)) {
                commands <- c(commands, 
                    sprintf(
                        "%s[%s,%i] <- %s",
                        contrast.name,
                        deparseToString(sprintf("group.%s", groupings[[j]])),
                        j,
                        .multiple_group_coefs(groupings[[j]])
                    )
                )
            }

            last <- groupings[[ngroups]]
            commands <- c(
                commands,
                sprintf("%s[%s,] <- -%s", contrast.name, deparseToString(sprintf("group.%s", last)), .multiple_group_coefs(last)),
                sprintf("colnames(%s) <- paste0(%s, ' - ', %s)", contrast.name, deparseToString(head(all.names, -1)), deparseToString(all.names[ngroups]))
            )
            info$commands <- commands
        }

        output[[i]] <- info
    }

    output
}

.multiple_group_names <- function(g) {
    ngroups <- length(g)
    if (ngroups == 1) {
        g
    } else {
        sprintf("(%s)/%i", paste(g, collapse=" + "), ngroups)
    }
}

.multiple_group_coefs <- function(g){
    ngroups <- length(g)
    if (ngroups == 1) {
        "1"
    } else {
        sprintf("1/%i", ngroups)
    }
}
