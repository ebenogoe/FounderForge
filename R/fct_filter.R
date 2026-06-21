# ---------------------------------------------------------------------------
# Filtering and imputation
# ---------------------------------------------------------------------------

#' Filter a genotype matrix by MAF, missingness, monomorphism and chr size
#'
#' Filters are applied in a fixed, documented order so results are reproducible:
#'   1. drop named individuals
#'   2. drop individuals exceeding max_miss_ind
#'   3. drop SNPs exceeding max_miss_snp
#'   4. drop monomorphic SNPs (optional)
#'   5. drop SNPs below the MAF threshold
#'   6. drop SNPs on chromosomes with fewer than min_snp_per_chr markers
#'
#' @param geno Integer matrix, SNPs x individuals.
#' @param snp_map data.frame with columns id, chr, pos (row-aligned to geno).
#' @param maf Minimum minor-allele frequency to keep a SNP (default 0.05).
#' @param max_miss_snp Maximum per-SNP missing rate (default 1 = no filter).
#' @param max_miss_ind Maximum per-individual missing rate (default 1 = no filter).
#' @param drop_mono Drop monomorphic SNPs (default TRUE).
#' @param min_snp_per_chr Minimum SNPs a chromosome must retain (default 2L;
#'   AlphaSimR requires at least 2 segregating sites per chromosome).
#' @param drop_individuals Character vector of individual ids to remove.
#' @param logger Optional logger.
#' @return list(geno, snp_map, report) where `report` is a before/after data.frame.
#' @noRd
filter_genotypes <- function(geno, snp_map,
                             maf = 0.05,
                             max_miss_snp = 1,
                             max_miss_ind = 1,
                             drop_mono = TRUE,
                             min_snp_per_chr = 2L,
                             drop_individuals = character(0),
                             logger = NULL) {
  logger <- as_logger(logger)
  logger$section("Filtering")
  steps <- list()
  record <- function(label) {
    steps[[length(steps) + 1L]] <<- data.frame(
      step = label, n_snp = nrow(geno), n_ind = ncol(geno),
      stringsAsFactors = FALSE
    )
  }
  record("input")

  # 1. drop named individuals
  if (length(drop_individuals) > 0) {
    keep <- !(colnames(geno) %in% drop_individuals)
    n_drop <- sum(!keep)
    geno <- geno[, keep, drop = FALSE]
    logger$log(sprintf("Dropped %d named individual(s).", n_drop))
    record("drop named individuals")
  }

  # 2. individual missingness
  if (max_miss_ind < 1) {
    mi <- colMeans(is.na(geno))
    keep <- mi <= max_miss_ind
    geno <- geno[, keep, drop = FALSE]
    logger$log(sprintf("Removed %d individual(s) with missing > %.2f.",
                       sum(!keep), max_miss_ind))
    record("individual missingness")
  }

  # 3. SNP missingness
  if (max_miss_snp < 1) {
    ms <- rowMeans(is.na(geno))
    keep <- ms <= max_miss_snp
    geno <- geno[keep, , drop = FALSE]
    snp_map <- snp_map[keep, , drop = FALSE]
    logger$log(sprintf("Removed %d SNP(s) with missing > %.2f.",
                       sum(!keep), max_miss_snp))
    record("SNP missingness")
  }

  # MAF/monomorphic computed on current matrix
  alt_freq <- rowMeans(geno, na.rm = TRUE) / 2
  maf_vec <- pmin(alt_freq, 1 - alt_freq)

  # 4. monomorphic
  if (isTRUE(drop_mono)) {
    keep <- !(is.na(maf_vec) | maf_vec == 0)
    geno <- geno[keep, , drop = FALSE]
    snp_map <- snp_map[keep, , drop = FALSE]
    maf_vec <- maf_vec[keep]
    logger$log(sprintf("Removed %d monomorphic SNP(s).", sum(!keep)))
    record("monomorphic")
  }

  # 5. MAF threshold
  if (maf > 0) {
    keep <- !is.na(maf_vec) & maf_vec >= maf
    geno <- geno[keep, , drop = FALSE]
    snp_map <- snp_map[keep, , drop = FALSE]
    logger$log(sprintf("Removed %d SNP(s) with MAF < %.3f.", sum(!keep), maf))
    record(sprintf("MAF >= %.3f", maf))
  }

  # 6. min SNPs per chromosome
  if (min_snp_per_chr > 1) {
    counts <- table(snp_map$chr)
    small <- names(counts)[counts < min_snp_per_chr]
    if (length(small) > 0) {
      keep <- !(as.character(snp_map$chr) %in% small)
      geno <- geno[keep, , drop = FALSE]
      snp_map <- snp_map[keep, , drop = FALSE]
      logger$log(sprintf(
        "Removed %d SNP(s) on %d chromosome(s) with < %d markers.",
        sum(!keep), length(small), min_snp_per_chr
      ), "WARN")
    }
    record(sprintf(">= %d SNPs/chr", min_snp_per_chr))
  }

  if (nrow(geno) == 0) {
    stop("All SNPs were removed by filtering. Relax the thresholds and retry.",
         call. = FALSE)
  }
  report <- do.call(rbind, steps)
  logger$log(sprintf("Filtering complete: %s SNPs x %s individuals retained.",
                     format(nrow(geno), big.mark = ","),
                     format(ncol(geno), big.mark = ",")), "OK")
  list(geno = geno, snp_map = snp_map, report = report)
}

#' Impute missing genotypes with the per-SNP rounded mean (0/1/2)
#'
#' AlphaSimR cannot accept NA in haplotype matrices. This defensive imputation
#' replaces missing calls with the rounded per-SNP mean dosage.
#'
#' @param geno Integer matrix, SNPs x individuals.
#' @param logger Optional logger.
#' @return list(geno, n_imputed)
#' @noRd
impute_missing <- function(geno, logger = NULL) {
  logger <- as_logger(logger)
  n_missing <- sum(is.na(geno))
  if (n_missing == 0) {
    logger$log("No missing genotypes; imputation not needed.")
    return(list(geno = geno, n_imputed = 0L))
  }
  row_means <- round(rowMeans(geno, na.rm = TRUE))
  row_means[is.nan(row_means)] <- 0
  na_rows <- which(rowSums(is.na(geno)) > 0)
  for (i in na_rows) {
    geno[i, is.na(geno[i, ])] <- as.integer(row_means[i])
  }
  logger$log(sprintf("Imputed %s missing cell(s) across %d SNP(s) (per-SNP rounded mean).",
                     format(n_missing, big.mark = ","), length(na_rows)), "OK")
  list(geno = geno, n_imputed = n_missing)
}
