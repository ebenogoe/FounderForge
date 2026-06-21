# ---------------------------------------------------------------------------
# Optional trait / SimParam definition
# ---------------------------------------------------------------------------

#' Build a SimParam with optional additive traits
#'
#' @param founderPop A MapPop from [build_founder_pop()].
#' @param trait_specs A list of trait definitions, each a list with elements
#'   name (character), nQtlPerChr (scalar or per-chromosome vector), mean, var,
#'   and h2. May be empty/NULL, in which case a bare SimParam is returned.
#' @param logger Optional logger.
#' @return An AlphaSimR SimParam object.
#' @noRd
define_traits <- function(founderPop, trait_specs = list(), logger = NULL) {
  logger <- as_logger(logger)
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop("Defining traits requires the 'AlphaSimR' package.", call. = FALSE)
  }
  logger$section("Defining simulation parameters")
  SP <- AlphaSimR::SimParam$new(founderPop)

  if (length(trait_specs) == 0) {
    logger$log("No traits defined; saving a bare SimParam.")
    return(SP)
  }

  h2 <- numeric(0)
  for (spec in trait_specs) {
    nm <- spec$name %||% sprintf("Trait%d", length(h2) + 1L)
    SP$addTraitA(
      nQtlPerChr = spec$nQtlPerChr,
      mean = spec$mean,
      var = spec$var
    )
    h2 <- c(h2, spec$h2)
    logger$log(sprintf(
      "Added trait '%s': nQTL/chr = %s, mean = %s, var = %s, h2 = %s.",
      nm, paste(spec$nQtlPerChr, collapse = ","),
      spec$mean, spec$var, spec$h2
    ))
  }
  if (any(!is.na(h2))) {
    SP$setVarE(h2 = h2)
    logger$log(sprintf("Set heritability h2 = c(%s).",
                       paste(round(h2, 3), collapse = ", ")))
  }
  SP
}
