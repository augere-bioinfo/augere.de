#' @import augere.core
#' @importMethodsFrom alabaster.se saveObject
.finalize <- function(assay, author, merge.metadata) {
    template <- system.file("templates", "finalize.Rmd", package="augere.de", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))

    if (!merge.metadata) {
        parsed[["merge-metadata"]] <- NULL
    } else {
        parsed[["no-merge-metadata"]] <- NULL
    }

    replacePlaceholders(parsed,
        list(
            ASSAY=safeDeparse(assay),
            AUTHOR=safeDeparse(as.list(author))
        )
    )
}

