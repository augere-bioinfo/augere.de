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
#' \item \code{left} and \code{right}, character vectors specifying the groups on the left (numerator) or right (denominator) of a \dQuote{versus} comparison.
#' Only present if \code{type = "versus"}.
#' \item \code{groups}, character vector specifying the groups involved in an ANOVA-like comparison.
#' Only present if \code{type = "anova"}.
#' \item \code{covariate}, the name of the covariate being tested.
#' Only present if \code{type = "covariate"}.
#' \item \code{commands}, character vector of R commands that produce the desired contrast vector/matrix.
#' This assumes that the evaluation environment has a design matrix named \code{design.name}.
#' The newly defined contrast vector/matrix will be stored as \code{contrast.name} in the environment. 
#' }
#' 
#' @details
#' If a vector in \code{comparisons} is of length 1, the sole entry is assumed to refer to a covariate in the \code{covariates} from \code{\link{processSimpleDesignMatrix}}.
#' The null hypothesis is that the coefficient of the design matrix with the same name is zero, i.e., the covariate has no effect.
#'
#' If a comparison vector is of length 2, the entries represent the names of two groups in the specified \code{groups} factor.
#' The null hypothesis is that there is no differential expression between the two groups.
#' The log-fold change is defined as the first group (left) over the second (right). 
#' 
#' If a comparison vector is of length 3 or greater, the entries represent the names of groups in the specified \code{groups} factor.
#' \itemize{
#' \item If none of the strings are \code{NA}, null hypothesis is that all of the specified levels of \code{groups} are equal.
#' This is an ANOVA-like contrast where contrasts are formulated with respect to the last level,
#' i.e., for \code{n} coefficients, \code{n-1} log-fold changes are reported representing the differences relative to the last coefficient.
#' \item If any of the strings are \code{NA}, this is used to split the vector into two groups of groups.
#' The null hypothesis is that the averages of the two groups of groups are equal.
#' The log-fold change is defined as the first average over the second.
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
#' processSimpleComparisons(c("treatment1", "treatment2", NA, "healthy"))
#' processSimpleComparisons(list(
#'     main=c("disease", "healthy"),
#'     secondary="dosage"
#' ))
#'
#' @export
processSimpleComparisons <- function(comparisons, design.name = "design", contrast.name = "con") { 
    if (!is.list(comparisons)) {
        comparisons <- list(comparisons)
    }
    output <- vector("list", length(comparisons))

    for (i in seq_along(comparisons)) {
        con <- as.character(comparisons[[i]])

        info <- list(title = NULL)
        if (!is.null(names(comparisons))) {
            nm <- names(comparisons)[i]
            if (nm != "") {
                info$title <- nm
            }
        }

        ngroups <- length(con)
        if (ngroups == 0L) {
            stop("vectors in 'comparisons' cannot be empty")

        } else if (ngroups == 1L) {
            if (is.null(info$title)) {
                info$title <- sprintf("Effect of increasing `%s`", con)
            }
            info$type <- "covariate"
            info$covariate <- con 

            info$commands <- c(
                sprintf("%s <- numeric(ncol(%s))", contrast.name, design.name),
                sprintf("names(%s) <- colnames(%s)", contrast.name, design.name),
                sprintf("%s[%s] <- 1", contrast.name, deparseToString(paste0("covariate.", con)))
            )

        } else if (ngroups == 2L) {
            if (is.null(info$title)) {
                info$title <- sprintf("Increase in `%s` over `%s`", con[1], con[2])
            }

            info$type <- "versus"
            info$left <- con[1]
            info$right <- con[2]

            info$commands <- c(
                sprintf("%s <- numeric(ncol(%s))", contrast.name, design.name),
                sprintf("names(%s) <- colnames(%s)", contrast.name, design.name),
                sprintf("%s[%s] <- 1", contrast.name, deparseToString(paste0("group.", con[1]))),
                sprintf("%s[%s] <- -1", contrast.name, deparseToString(paste0("group.", con[2])))
            )

        } else if (anyNA(con)) {
            na <- is.na(con)
            if (sum(na) != 1L) {
                stop("expected no more than one NA value in the 'comparisons' vector")
            }

            separator <- which(na)
            left <- con[1:(separator - 1)]
            right <- con[(separator + 1):length(separator)]

            info$type <- "versus"
            info$left <- left
            info$right <- right

            if (is.null(info$title)) {
                info$title <- sprintf("Increase in `%s` over %s`", .multiple_group_names(left), .multiple_group_names(right))
            }
            info$commands <- c(
                sprintf("%s <- numeric(ncol(%s))", contrast.name, design.name),
                sprintf("names(%s) <- colnames(%s)", contrast.name, design.name),
                sprintf("%s[%s] <- %s", contrast.name, deparseToString(sprintf("group.%s", left)), .multiple_group_coefs(left)),
                sprintf("%s[%s] <- -%s", contrast.name, deparseToString(sprintf("group.%s", right)), .multiple_group_coefs(right))
            )

        } else {
            info$type <- "anova"
            if (is.null(info$title)) {
                info$title <- sprintf("One-way ANOVA involving %s", paste(sprintf("`%s`", con), collapse=", "))
            }
            info$groups <- con

            commands <- c(
                sprintf("%s <- matrix(0, ncol(%s), %s)", contrast.name, design.name, ngroups - 1),
                sprintf("rownames(%s) <- colnames(%s)", contrast.name ,design.name)
            )
            for (j in seq_len(ngroups - 1)) {
                commands <- c(commands, sprintf("%s[%s,%i] <- 1", contrast.name, deparseToString(paste0("group.", con[j])), j))
            }

            last <- con[ngroups]
            commands <- c(
                commands,
                sprintf("%s[%s,] <- -1", contrast.name, deparseToString(paste0("group.", last))),
                sprintf("colnames(%s) <- paste0(%s, ' - ', %s)", contrast.name, deparseToString(head(con, -1)), deparseToString(last))
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
