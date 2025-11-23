#' @import augere.core
#' @importMethodsFrom alabaster.se saveObject
.finalize <- function(assay, author, merge.metadata, subset.metadata) {
    template <- system.file("templates", "finalize.Rmd", package="augere.de", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))

    if (!merge.metadata) {
        parsed[["merge-metadata"]] <- NULL
    } else {
        parsed[["no-merge-metadata"]] <- NULL
    }

    if (!subset.metadata) {
        parsed[["subset-metadata"]] <- "    summarized_experiment=setNames(list(), character(0))"
    }

    replacePlaceholders(parsed,
        list(
            ASSAY=deparseToString(assay),
            AUTHOR=deparseToString(as.list(author))
        )
    )
}

