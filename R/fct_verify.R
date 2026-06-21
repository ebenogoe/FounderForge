# ---------------------------------------------------------------------------
# Founder-population verification and before/after PCA validation
# ---------------------------------------------------------------------------

#' Verify a founder population preserved the expected genotype structure
#'
#' Adds a temporary SNP chip to a fresh SimParam, pulls SNP genotypes and checks
#' the dosage values. For an outbred import these must include 0, 1 and 2
#' (heterozygosity preserved). For an inbred import only 0 and 2 are expected.
#'
#' @param founderPop A MapPop.
#' @param strategy The import strategy used ("outbred"/"inbred"/"import_inbred").
#' @param n_snp_per_chr Chip size per chromosome for the check.
#' @param logger Optional logger.
#' @return list(values, pass, expected, n_seg, seg_per_chr).
#' @noRd
verify_founder <- function(founderPop, strategy = "outbred",
                           n_snp_per_chr = 100L, logger = NULL) {
  logger <- as_logger(logger)
  logger$section("Verification")

  SP_tmp <- AlphaSimR::SimParam$new(founderPop)
  seg_per_chr <- vapply(SP_tmp$segSites, function(x) x[1L], integer(1L))
  chip <- min(n_snp_per_chr, max(1L, min(seg_per_chr) - 1L))
  SP_tmp$addSnpChip(nSnpPerChr = chip)

  pop_tmp <- AlphaSimR::newPop(founderPop, simParam = SP_tmp)
  geno <- AlphaSimR::pullSnpGeno(pop_tmp, simParam = SP_tmp)
  values <- sort(unique(as.vector(geno)))

  if (identical(strategy, "outbred")) {
    pass <- all(c(0L, 1L, 2L) %in% values)
    expected <- "0, 1, 2 (heterozygosity preserved)"
  } else {
    pass <- all(values %in% c(0L, 2L)) && length(values) >= 1
    expected <- "0, 2 (fully homozygous lines)"
  }

  logger$log(sprintf("pullSnpGeno() unique values: %s",
                     paste(values, collapse = ", ")))
  logger$log(
    sprintf("Heterozygosity check: %s (expected %s).",
            if (pass) "PASSED" else "FAILED", expected),
    if (pass) "OK" else "WARN"
  )

  list(
    values = values,
    pass = pass,
    expected = expected,
    n_seg = sum(seg_per_chr),
    seg_per_chr = seg_per_chr
  )
}

#' Before/after PCA validation of an AlphaSimR import
#'
#' Computes PCA on (1) the original/filtered genotype matrix and (2) the SNP
#' genotypes pulled back out of AlphaSimR via a chip. If the two agree, the
#' import is faithful. Uses prcomp(scale. = FALSE) per the genomic convention.
#'
#' @param geno_original Integer matrix SNPs x individuals (filtered, pre-import).
#' @param founderPop A MapPop.
#' @param color_groups Optional vector (length = n individuals) for point colour.
#' @param logger Optional logger.
#' @return list(plot, variance) where plot is a ggplot/patchwork object and
#'   variance is a data.frame of variance explained for both PCAs.
#' @noRd
pca_validation <- function(geno_original, founderPop, color_groups = NULL,
                           logger = NULL) {
  logger <- as_logger(logger)
  logger$section("Before/after PCA validation")

  # PCA 1: original genotypes (individuals x SNPs)
  x1 <- t(geno_original)
  pca1 <- stats::prcomp(x1, scale. = FALSE, center = TRUE)

  # PCA 2: genotypes pulled from AlphaSimR via a chip
  SP <- AlphaSimR::SimParam$new(founderPop)
  seg_per_chr <- vapply(SP$segSites, function(x) x[1L], integer(1L))
  chip <- min(1000L, max(1L, min(seg_per_chr) - 1L))
  SP$addSnpChip(nSnpPerChr = chip)
  pop <- AlphaSimR::newPop(founderPop, simParam = SP)
  x2 <- AlphaSimR::pullSnpGeno(pop, simParam = SP)
  pca2 <- stats::prcomp(x2, scale. = FALSE, center = TRUE)

  ve <- function(p) round(100 * p$sdev^2 / sum(p$sdev^2), 2)
  ve1 <- ve(pca1)
  ve2 <- ve(pca2)

  df1 <- data.frame(PC1 = pca1$x[, 1], PC2 = pca1$x[, 2])
  df2 <- data.frame(PC1 = pca2$x[, 1], PC2 = pca2$x[, 2])
  if (!is.null(color_groups) && length(color_groups) == nrow(df1)) {
    df1$group <- color_groups
    df2$group <- color_groups
  }

  plot <- build_pca_plot(df1, df2, ve1, ve2)

  variance <- data.frame(
    PC = paste0("PC", seq_len(min(10, length(ve1)))),
    original_pct = ve1[seq_len(min(10, length(ve1)))],
    alphasimr_pct = ve2[seq_len(min(10, length(ve2)))]
  )
  logger$log(sprintf(
    "PCA validation: original PC1/PC2 = %.1f/%.1f%%, AlphaSimR PC1/PC2 = %.1f/%.1f%%.",
    ve1[1], ve1[2], ve2[1], ve2[2]
  ), "OK")

  list(plot = plot, variance = variance)
}

#' Assemble the side-by-side PCA comparison figure
#' @noRd
build_pca_plot <- function(df1, df2, ve1, ve2) {
  has_group <- "group" %in% names(df1)
  one <- function(df, title, ve) {
    aes_xy <- if (has_group) {
      ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$group)
    } else {
      ggplot2::aes(x = .data$PC1, y = .data$PC2)
    }
    ggplot2::ggplot(df, aes_xy) +
      ggplot2::geom_point(alpha = 0.7, size = 1.8) +
      ggplot2::labs(
        title = title,
        x = sprintf("PC1 (%.1f%%)", ve[1]),
        y = sprintf("PC2 (%.1f%%)", ve[2])
      ) +
      ggplot2::theme_minimal(base_size = 13)
  }
  p1 <- one(df1, "Before: original genotypes", ve1)
  p2 <- one(df2, "After: AlphaSimR import", ve2)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, ncol = 2)
  } else {
    list(before = p1, after = p2)
  }
}
