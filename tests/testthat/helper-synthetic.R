# Synthetic genotype data shared across tests.
make_synthetic <- function(nSNP = 60L, nInd = 24L, nChr = 3L, seed = 1L,
                           n_missing = 0L) {
  set.seed(seed)
  per <- nSNP / nChr
  chr <- rep(seq_len(nChr), each = per)
  pos <- as.numeric(unlist(lapply(seq_len(nChr), function(k) {
    sort(sample(1e6:5e6, per))
  })))
  geno <- matrix(sample(0:2, nSNP * nInd, replace = TRUE), nrow = nSNP)
  if (n_missing > 0) geno[sample(length(geno), n_missing)] <- NA_integer_
  rownames(geno) <- paste0("S", chr, "_", pos)
  colnames(geno) <- paste0("ind", seq_len(nInd))
  snp_map <- data.frame(id = rownames(geno), chr = chr, pos = pos)
  list(geno = geno, snp_map = snp_map)
}

skip_if_no_alphasimr <- function() {
  testthat::skip_if_not_installed("AlphaSimR")
}
