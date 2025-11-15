#' Custom contrasts
#'
#' Construct R commands to define custom contrast vectors/matrices.
#'
#' @param contrasts Numeric matrix or vector, a function to generate such a matrix/vector,
#' or a character vector to be passed to \code{\link[limma]{makeContrasts}}, see Details.
#' Alternatively, a list of such objects to specify multiple contrasts.
#' @param design.name String containing the variable name of the design matrix.
#' @param contrast.name String containing the variable name of the contrast vector/matrix.
#'
#' @return 
#' A list containing one entry per comparison, where each entry contains:
#' \itemize{
#' \item \code{title}, the title for the comparison.
#' If \code{comparisons} is named, the corresponding (non-empty) name is used here, otherwise an appropriate title is automatically generated.
#' \item \code{type}, the type of the comparison.
#' This is always set to \code{"custom"}.
#' \item \code{commands}, character vector of R commands that produce the desired contrast vector/matrix.
#' This assumes that the evaluation environment has a design matrix named \code{design.name}.
#' The newly defined contrast vector/matrix will be stored as \code{contrast.name} in the environment. 
#' }
#'
#' @details
#' If \code{contrasts} is a a character vector, it will be passed to \code{\link[limma]{makeContrasts}} with \code{levels=} set to the design matrix.
#' A vector of length 2 or more represents an ANOVA-like comparison.
#'
#' If \code{contrasts} is a numeric vector, it should have names that match some or all of the column names of the design matrix.
#' The values of this vector represent the entries of the contrast vector for the named coefficients; all other entries of the contrast vector are set to zero.
#' 
#' If \code{contrasts} is a numeric matrix, its row names should match some or all of the column names of the design matrix.
#' The rows of this vector represent the rows of the contrast matrix for the named coefficients; all other entries of the contrast matrix are set to zero.
#' The column names of the contrast matrix are set to those of \code{contrasts}, if available.
#' 
#' If \code{contrasts} is a function, it should accept the design matrix and return a contrast vector/matrrix.
#' This mode is useful when the deparsed design matrix is difficult to read in the Rmarkdown report.
#'
#' @examples
#' processCustomContrasts(c(coef1 = 0.5, coef2 = 0.5, coef3 = -1))
#'
#' mat <- matrix(0, 3, 2)
#' rownames(mat) <- c("grp1", "grp2", "grp3")
#' mat["grp1",] <- 1
#' mat["grp2",1] <- -1
#' mat["grp3",2] <- -1
#' processCustomContrasts(mat)
#'
#' processCustomContrasts("TREATMENT - CONTROL")
#'
#' processCustomContrasts(function(design) { 
#'     out <- numeric(ncol(design))
#'     names(out) <- colnames(design)
#'     out[c("treated", "control")] <- c(1, -1)
#'     out
#' })
#'
#' @export
#' @import augere.core 
processCustomContrasts <- function(contrasts, design.name = "design", contrast.name = "con") {
    if (!is.list(contrasts)) {
        contrasts <- list(contrasts)
    }
    output <- vector("list", length(contrasts))

    for (i in seq_along(contrasts)) {
        con <- contrasts[[i]]

        if (is.character(con)) {
            cmd <- sprintf("%s <- limma::makeContrasts(contrasts=%s, levels=%s)", contrast.name, deparseToString(con), design.name)

            # Remove dimensions if the design only has one column, so that
            # limma doesn't use the column names to rename the LFC-related
            # columns in the output table.
            cmd <- c(cmd, 
                sprintf("if (!is.null(dim(%s)) && ncol(%s) == 1) {", contrast.name, contrast.name),
                sprintf("    %s <- drop(%s)", contrast.name, contrast.name),
                        "}"
            )

        } else if (is.function(con)) {
            cmd <- sprintf("%s <- (%s)(%s)", contrast.name, deparseToString(con, indent=0), design.name)

            # Ditto to the above, but we add column names if they  aren't present already.
            cmd <- c(cmd, 
                sprintf("if (!is.null(dim(%s))) {", contrast.name),
                sprintf("    if (ncol(%s) == 1) {", contrast.name),
                sprintf("        %s <- drop(%s)", contrast.name, contrast.name),
                sprintf("    } else if (is.null(colnames(%s))) {", contrast.name),
                sprintf("        colnames(%s) <- seq_len(ncol(%s))", contrast.name, contrast.name),
                        "    }",
                        "}"
            )

        } else {
            if (!is.null(dim(con)) && ncol(con) > 1) {
                cmd <- c(
                    sprintf("%s <- local({", contrast.name),
                    sprintf("    tmp <- %s", deparseToString(con)),
                    sprintf("    out <- matrix(0, ncol(%s), ncol(tmp))", design.name), 
                    sprintf("    rownames(out) <- colnames(%s)", design.name),
                            "    out[rownames(tmp),] <- tmp"
                )
                if (is.null(colnames(con))) {
                    cmd <- c(cmd, 
                            "    colnames(out) <- seq_len(ncol(out))"
                    )
                } else {
                    cmd <- c(cmd, 
                            "    colnames(out) <- colnames(tmp)"
                    )
                }
                cmd <- c(cmd,
                            "    out",
                            "})"
                )

            } else {
                cmd <- c(
                    sprintf("%s <- local({", contrast.name),
                    sprintf("    tmp <- %s", deparseToString(drop(con))),
                    sprintf("    out <- numeric(ncol(%s))", design.name),
                    sprintf("    names(out) <- colnames(%s)", design.name),
                            "    out[names(tmp)] <- tmp", 
                            "    out",
                            "})"
                )
            }
        }

        title <- NULL
        if (!is.null(names(contrasts))) {
            name <- names(contrasts)[i]
            if (name != "") {
                title <- name
            }
        }
        if (is.null(title)) {
            title <- paste("Contrast", i)
        }

        output[[i]] <- list(title=title, type="custom", commands=cmd)
    }

    output
}
