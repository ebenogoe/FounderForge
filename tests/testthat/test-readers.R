test_that("detect_format recognises extensions", {
  expect_equal(detect_format("a.vcf"), "vcf")
  expect_equal(detect_format("a.vcf.gz"), "vcf")
  expect_equal(detect_format("a.hmp.txt"), "hapmap")
  expect_equal(detect_format("a.raw"), "plink")
  expect_equal(detect_format("a.ped"), "plink")
  expect_equal(detect_format("a.csv"), "numeric")
  expect_equal(detect_format("a.rds"), "numeric")
})

test_that("read_vcf converts GT to dosage and drops multi-allelic sites", {
  skip_if_not_installed("vcfR")
  vcf_txt <- c(
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "s1", "s2", "s3"), collapse = "\t"),
    paste(c("1", "1000", "snp1", "A", "G", ".", ".", ".", "GT",
            "0/0", "0/1", "1/1"), collapse = "\t"),
    paste(c("1", "2000", "snp2", "C", "T", ".", ".", ".", "GT",
            "1/1", "0/0", "0|1"), collapse = "\t"),
    paste(c("2", "1500", ".", "G", "A", ".", ".", ".", "GT",
            "0/1", "1/1", "./."), collapse = "\t"),
    paste(c("2", "2500", "snp4", "T", "C,G", ".", ".", ".", "GT",
            "0/0", "0/1", "1/1"), collapse = "\t")
  )
  vf <- tempfile(fileext = ".vcf")
  writeLines(vcf_txt, vf)
  r <- read_genotypes(vf, "vcf")
  expect_equal(dim(r$geno), c(3L, 3L))  # multi-allelic dropped
  expect_equal(unname(r$geno["snp1", ]), c(0L, 1L, 2L))
  expect_true(is.na(r$geno[3, 3]))      # ./. -> NA
  expect_equal(r$snp_map$id[3], "S2_1500")  # blank ID -> fallback
})

test_that("read_hapmap handles IUPAC and two-character calls", {
  hmp <- c(
    paste(c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center",
            "protLSID", "assayLSID", "panelLSID", "QCcode",
            "L1", "L2", "L3"), collapse = "\t"),
    paste(c("m1", "A/G", "1", "100", "+", "NA", "NA", "NA", "NA", "NA", "NA",
            "A", "R", "G"), collapse = "\t"),
    paste(c("m2", "C/T", "1", "200", "+", "NA", "NA", "NA", "NA", "NA", "NA",
            "CC", "CT", "TT"), collapse = "\t"),
    paste(c("m3", "G/A", "2", "300", "+", "NA", "NA", "NA", "NA", "NA", "NA",
            "G", "N", "A"), collapse = "\t")
  )
  hf <- tempfile(fileext = ".hmp.txt")
  writeLines(hmp, hf)
  h <- read_genotypes(hf, "hapmap")
  expect_equal(dim(h$geno), c(3L, 3L))
  expect_equal(unname(h$geno["m1", ]), c(0L, 1L, 2L))
  expect_equal(unname(h$geno["m2", ]), c(0L, 1L, 2L))
  expect_true(is.na(h$geno["m3", 2]))
})

test_that("read_numeric round-trips a list(geno, snp_map) RDS", {
  d <- make_synthetic()
  f <- tempfile(fileext = ".rds")
  saveRDS(list(geno = d$geno, snp_map = d$snp_map), f)
  r <- read_genotypes(f, "numeric")
  expect_equal(dim(r$geno), dim(d$geno))
  expect_equal(r$snp_map$chr, d$snp_map$chr)
})

test_that("read_numeric auto-detects columns and excludes HapMap metadata", {
  df <- data.frame(
    `rs#` = c("S1_1", "S1_2", "S2_1"),
    alleles = c("A/G", "C/T", "G/A"),
    chrom = c(1, 1, 2),
    pos = c(10, 20, 5),
    REFERENCE_GENOME = c("A", "C", "G"),
    `IND1:lane:1` = c(0, 1, 2),
    `IND2:lane:2` = c(2, 1, 0),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  f <- tempfile(fileext = ".csv")
  data.table::fwrite(df, f)
  r <- read_genotypes(f, "numeric", opts = list(clean_sample_names = TRUE))
  expect_equal(dim(r$geno), c(3L, 2L))                 # SNP x ind, metadata dropped
  expect_equal(colnames(r$geno), c("IND1", "IND2"))    # names cleaned, REFERENCE_GENOME gone
  expect_equal(unname(r$geno["S1_1", ]), c(0L, 2L))
  expect_equal(r$snp_map$chr, c(1L, 1L, 2L))
})

test_that("bundled sample datasets all read into the internal representation", {
  for (f in c("sample_numeric.rds", "sample_numeric.csv", "sample.vcf", "sample.hmp.txt")) {
    p <- system.file("extdata", f, package = "FounderForge")
    if (!nzchar(p)) p <- file.path("..", "..", "inst", "extdata", f)
    skip_if_not(file.exists(p))
    r <- read_genotypes(p, "auto")
    expect_equal(nrow(r$geno), 300L)
    expect_equal(ncol(r$geno), 30L)
    expect_true(all(sort(unique(stats::na.omit(as.vector(r$geno)))) %in% 0:2))
  }
})

test_that("validate_geno_map rejects out-of-range codes", {
  d <- make_synthetic()
  bad <- d$geno
  bad[1, 1] <- 5L
  expect_error(validate_geno_map(bad, d$snp_map), "0/1/2")
})
