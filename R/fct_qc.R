# ---------------------------------------------------------------------------
# Quality-control summaries
# ---------------------------------------------------------------------------

#' Compute a quality-control summary of a genotype matrix
#'
#' @param geno Integer matrix, SNPs x individuals, values 0/1/2/NA.
#' @param snp_map data.frame with columns id, chr, pos.
#' @param logger Optional logger.
#' @return A named list of summary scalars, tables and per-SNP/per-individual
#'   vectors used by the UI value boxes and plots.
#' @noRd
compute_qc_summary <- function(geno, snp_map, logger = NULL) {
  logger <- as_logger(logger)
  logger$log("Computing QC summary ...")

  n_snp <- nrow(geno)
  n_ind <- ncol(geno)

  # Allele frequency and MAF (alternate allele = dosage / 2).
  alt_freq <- rowMeans(geno, na.rm = TRUE) / 2
  maf <- pmin(alt_freq, 1 - alt_freq)
  maf[is.nan(maf)] <- NA_real_

  miss_snp <- rowMeans(is.na(geno))
  miss_ind <- colMeans(is.na(geno))

  het_snp <- rowMeans(geno == 1L, na.rm = TRUE)
  het_ind <- colMeans(geno == 1L, na.rm = TRUE)

  n_mono <- sum(maf == 0 | is.na(maf))
  snps_per_chr <- as.data.frame(table(chr = snp_map$chr),
                                responseName = "n_snps",
                                stringsAsFactors = FALSE)
  snps_per_chr$chr <- suppressWarnings(as.integer(snps_per_chr$chr))
  snps_per_chr <- snps_per_chr[order(snps_per_chr$chr), ]

  out <- list(
    n_snp = n_snp,
    n_ind = n_ind,
    n_chr = length(unique(snp_map$chr)),
    chromosomes = sort(unique(snp_map$chr)),
    snps_per_chr = snps_per_chr,
    n_monomorphic = n_mono,
    pct_missing = mean(is.na(geno)) * 100,
    mean_maf = mean(maf, na.rm = TRUE),
    mean_het = mean(het_ind, na.rm = TRUE),
    maf = maf,
    miss_snp = miss_snp,
    miss_ind = miss_ind,
    het_snp = het_snp,
    het_ind = het_ind
  )
  logger$log(sprintf(
    "QC: %s SNPs, %s individuals, %.3f%% missing, %s monomorphic, mean MAF %.3f.",
    format(n_snp, big.mark = ","), format(n_ind, big.mark = ","),
    out$pct_missing, format(n_mono, big.mark = ","), out$mean_maf
  ), "OK")
  out
}
