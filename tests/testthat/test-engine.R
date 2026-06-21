test_that("compute_qc_summary reports correct dimensions and missingness", {
  d <- make_synthetic(n_missing = 30L)
  qc <- compute_qc_summary(d$geno, d$snp_map)
  expect_equal(qc$n_snp, nrow(d$geno))
  expect_equal(qc$n_ind, ncol(d$geno))
  expect_equal(qc$n_chr, length(unique(d$snp_map$chr)))
  expect_equal(round(qc$pct_missing, 4), round(100 * 30 / length(d$geno), 4))
  expect_length(qc$maf, nrow(d$geno))
})

test_that("filter_genotypes removes low-MAF SNPs and reports steps", {
  d <- make_synthetic()
  # Force one monomorphic SNP.
  d$geno[1, ] <- 0L
  res <- filter_genotypes(d$geno, d$snp_map, maf = 0.1, drop_mono = TRUE)
  expect_lt(nrow(res$geno), nrow(d$geno))
  expect_true(all(c("step", "n_snp", "n_ind") %in% names(res$report)))
  af <- rowMeans(res$geno) / 2
  expect_true(all(pmin(af, 1 - af) >= 0.1 - 1e-9))
})

test_that("filter_genotypes errors when everything is removed", {
  d <- make_synthetic()
  expect_error(filter_genotypes(d$geno, d$snp_map, maf = 0.9),
               "All SNPs were removed")
})

test_that("impute_missing eliminates NA", {
  d <- make_synthetic(n_missing = 50L)
  imp <- impute_missing(d$geno)
  expect_equal(sum(is.na(imp$geno)), 0L)
  expect_equal(imp$n_imputed, 50L)
})

test_that("build_genmap starts each chromosome at zero and is strictly increasing", {
  d <- make_synthetic()
  ord <- order_by_map(d$geno, d$snp_map)
  gm <- build_genmap(ord$snp_map)
  for (v in gm$map_list) {
    expect_equal(unname(v[1]), 0)
    expect_true(all(diff(v) > 0))
  }
  expect_equal(nrow(gm$map_df), nrow(d$snp_map))
})

test_that("outbred build preserves heterozygosity (0,1,2 present)", {
  skip_if_no_alphasimr()
  d <- make_synthetic()
  br <- build_founder_pop(d$geno, d$snp_map, strategy = "outbred")
  expect_equal(br$founderPop@nInd, ncol(d$geno))
  v <- verify_founder(br$founderPop, "outbred")
  expect_true(v$pass)
  expect_true(all(c(0L, 1L, 2L) %in% v$values))
})

test_that("inbred build yields only homozygous dosages", {
  skip_if_no_alphasimr()
  d <- make_synthetic()
  br <- build_founder_pop(d$geno, d$snp_map, strategy = "inbred")
  v <- verify_founder(br$founderPop, "inbred")
  expect_true(all(v$values %in% c(0L, 2L)))
})

test_that("import_inbred path builds a population", {
  skip_if_no_alphasimr()
  d <- make_synthetic()
  br <- build_founder_pop(d$geno, d$snp_map, strategy = "import_inbred")
  expect_equal(br$founderPop@nInd, ncol(d$geno))
})

test_that("define_traits returns a SimParam and write_outputs creates the tree", {
  skip_if_no_alphasimr()
  d <- make_synthetic()
  br <- build_founder_pop(d$geno, d$snp_map, strategy = "outbred")
  SP <- define_traits(br$founderPop,
                      list(list(name = "Yield", nQtlPerChr = 3,
                                mean = 100, var = 25, h2 = 0.6)))
  expect_true(inherits(SP, "R6"))
  qc <- compute_qc_summary(d$geno, d$snp_map)
  out <- write_outputs(tempdir(), br, SP = SP, qc_summary = qc,
                       source_info = list(path = "synthetic", format = "numeric"))
  expect_true(dir.exists(out))
  expect_true(file.exists(file.path(out, "founder_population", "founder_pop.rds")))
  expect_true(file.exists(file.path(out, "founder_population", "genmap.rds")))
  expect_true(file.exists(file.path(out, "founder_population", "rebuild_snippet.R")))
  expect_true(length(list.files(file.path(out, "logs"))) == 1L)
})
