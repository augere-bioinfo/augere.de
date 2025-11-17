#' @import augere.core 
.initialize <- function(
    x,
    method,
    assay,
    groups, 
    comparisons, 
    covariates,
    block,
    subset.factor,
    subset.levels,
    subset.groups,
    design,
    contrasts,
    author
) {
    template <- system.file("templates", "initialize.Rmd", package="augere.de", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))

    parsed[["create-se"]] <- processInputCommands(x, name="se")

    replacements <- list(
        METHOD=method,
        ASSAY=deparseToString(assay),
        AUTHOR=paste(sprintf("  - %s", author), collapse="\n")
    )

    if (!is.null(subset.factor)) {
        parsed[["subset-se"]] <- replacePlaceholders(
            parsed[["subset-se"]],
            list(
                SUBSET_FACTOR=deparseToString(subset.factor),
                SUBSET_LEVELS=deparseToString(subset.levels)
            )
        )
    } else {
        parsed[["subset-se"]] <- NULL
    }

    if (!is.null(design) && !is.null(contrasts)) {
        parsed[["design-matrix"]] <- processCustomDesignMatrix(design=design, se.name="se")
        contrast.info <- processCustomContrasts(contrasts)
        replacements$FILTER_OPTS <- "design=design"
        parsed[["subset-group"]] <- NULL
        parsed[["mds-grouped"]] <- NULL

    } else {
        parsed[["design-matrix"]] <- processSimpleDesignMatrix(groups=groups, block=block, covariates=covariates, se.name="se")
        contrast.info <- processSimpleComparisons(comparisons)

        if (is.null(groups)) {
            replacements$FILTER_OPTS <- "design=design"
            parsed[["subset-group"]] <- NULL
            parsed[["mds-grouped"]] <- NULL
        } else {
            replacements$FILTER_OPTS <- "group=model.data$group."
            parsed[["mds-ungrouped"]] <- NULL

            group.levels <- NULL
            if (subset.groups) {
                group.levels <- .find_used_groups(contrast.info)
            }
            if (!is.null(group.levels)) {
                parsed[["subset-group"]] <- replacePlaceholders(
                    parsed[["subset-group"]],
                    list(
                        GROUP_FACTOR=deparseToString(groups),
                        GROUP_LEVELS=deparseToString(group.levels)
                    )
                )
            } else {
                parsed[["subset-group"]] <- NULL
            }
        }
    }

    list(
        text=replacePlaceholders(parsed, replacements),
        contrasts=contrast.info
    )
}

.find_used_groups <- function(contrast.info) {
    # We don't want to throw out potentially relevant information.
    not.grouped <- vapply(contrast.info, function(con) con$type %in% "covariate", TRUE)
    if (any(not.grouped)) {
        return(NULL)
    }

    all.groups <- lapply(contrast.info, function(con) con[c("left", "right", "groups")])
    group.levels <- sort(unique(unlist(all.groups)))
}
