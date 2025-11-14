#' @import augere.core 
.initialize <- function(
    method,
    author,
    x,
    assay,
    groups, 
    comparisons, 
    covariates,
    block,
    subset.factor,
    subset.levels,
    subset.groups,
    design,
    contrasts
) {
    template <- system.file("templates", "init.Rmd", package="augere.de", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))

    parsed[["create-se"]] <- processInputCommands(x, obj="se")

    replacements <- list(
        METHOD=method,
        ASSAY=deparseToString(assay),
        AUTHOR=deparseToString(as.list(author))
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
        parsed[["subset-data"]] <- NULL
    }

    if (!is.null(design) && !is.null(contrasts)) {
        parsed[["design-matrix"]] <- processCustomDesignMatrix(design=design)
        contrast.info <- processCustomContrasts(contrasts)
        replacements$FILTER_OPTS <- "design=design"
        parsed[["subset-group"]] <- NULL

    } else {
        parsed[["design-matrix"]] <- processSimpleDesignMatrix(groups=groups, block=block, covariates=covariates)
        contrast.info <- processSimpleComparisons(comparisons)

        if (is.null(groups)) {
            replacements$FILTER_OPTS <- "design=design"
        } else {
            replacements$FILTER_ARGS <- "group=variables$group."
        }
        parsed[["mds-relabel"]] <- "labels <- variables$group."

        if (subset.groups && length(all.groups)) {
            not.grouped <- vapply(contrast.info, function(con) con$type %in% c("single", "other"), TRUE)
            group.levels <- NULL
            if (!any(not.grouped)) {
                all.groups <- lapply(contrast.info, function(con) con[c("left_groups", "right_groups", "anova_groups")])
                group.levels <- sort(unique(unlist(all.groups)))
            }

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

    list(
        text=replacePlaceholders(parsed, replacements),
        contrasts=contrast.info
    )
}
