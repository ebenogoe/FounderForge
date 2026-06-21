# ---------------------------------------------------------------------------
# Genotype readers
#
# All handlers normalise to one internal representation so every downstream
# step is format-agnostic:
#
#   list(
#     geno    = <integer matrix, rows = SNPs, cols = individuals, values 0/1/2/NA>,
#     snp_map = <data.frame: id (chr), chr (int), pos (numeric, bp)>
#   )
#
# rownames(geno) == snp_map$id and the two are kept row-aligned at all times.
# ---------------------------------------------------------------------------

#' Detect genotype file format from its extension
#' @noRd
detect_format <- function(path) {
  p <- tolower(path)
  if (grepl("\\.vcf(\\.gz|\\.bgz)?$", p)) {
    return("vcf")
  }
  if (grepl("\\.hmp(\\.txt)?(\\.gz)?$", p)) {
    return("hapmap")
  }
  if (grepl("\\.(ped|bed|raw)$", p)) {
    return("plink")
  }
  if (grepl("\\.(csv|tsv|txt|rds|rdata|rda)$", p)) {
    return("numeric")
  }
  "unknown"
}

#' Read genotypes from a file into the internal representation
#'
#' @param path Path to the genotype file.
#' @param format One of "auto", "vcf", "hapmap", "numeric", "plink".
#' @param opts A named list of format-specific options (see handlers).
#' @param logger A logger from [new_run_logger()] (optional).
#' @return list(geno, snp_map). See file header for the contract.
#' @noRd
read_genotypes <- function(path, format = "auto", opts = list(), logger = NULL) {
  logger <- as_logger(logger)
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  if (identical(format, "auto")) {
    format <- detect_format(path)
    logger$log(sprintf("Auto-detected format: %s", format))
  }
  res <- switch(
    format,
    vcf     = read_vcf(path, opts, logger),
    hapmap  = read_hapmap(path, opts, logger),
    numeric = read_numeric(path, opts, logger),
    plink   = read_plink(path, opts, logger),
    stop(sprintf(
      "Unsupported or unrecognised format '%s'. Supported: vcf, hapmap, numeric, plink.",
      format
    ), call. = FALSE)
  )
  validate_geno_map(res$geno, res$snp_map)
  logger$log(sprintf(
    "Read %s SNPs x %s individuals across %s chromosome(s).",
    format(nrow(res$geno), big.mark = ","),
    format(ncol(res$geno), big.mark = ","),
    length(unique(res$snp_map$chr))
  ), "OK")
  res
}

#' Validate the internal geno/snp_map contract
#' @noRd
validate_geno_map <- function(geno, snp_map) {
  if (!is.matrix(geno)) {
    stop("Internal error: 'geno' must be a matrix.", call. = FALSE)
  }
  if (nrow(geno) != nrow(snp_map)) {
    stop(sprintf(
      "Genotype rows (%d) and SNP map rows (%d) do not match.",
      nrow(geno), nrow(snp_map)
    ), call. = FALSE)
  }
  if (!all(c("id", "chr", "pos") %in% names(snp_map))) {
    stop("SNP map must contain columns: id, chr, pos.", call. = FALSE)
  }
  bad <- setdiff(unique(as.vector(geno)), c(0L, 1L, 2L, NA_integer_))
  if (length(bad) > 0) {
    stop(sprintf(
      "Genotypes must be coded 0/1/2 (or NA). Found unexpected value(s): %s",
      paste(utils::head(bad, 5), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# VCF
# ---------------------------------------------------------------------------

#' @noRd
read_vcf <- function(path, opts = list(), logger = NULL) {
  logger <- as_logger(logger)
  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("Reading VCF requires the 'vcfR' package. Please install it.", call. = FALSE)
  }
  logger$log("Reading VCF with vcfR::read.vcfR() ...")
  vcf <- vcfR::read.vcfR(path, verbose = FALSE)

  fix <- vcf@fix
  is_biallelic <- !grepl(",", fix[, "ALT"])
  n_multi <- sum(!is_biallelic)
  if (n_multi > 0) {
    logger$log(sprintf(
      "Dropping %s multi-allelic site(s); only biallelic SNPs are supported.",
      format(n_multi, big.mark = ",")
    ), "WARN")
  }

  gt_char <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
  gt_char <- gt_char[is_biallelic, , drop = FALSE]
  fix <- fix[is_biallelic, , drop = FALSE]

  geno <- gt_to_dosage(gt_char)

  chr_raw <- fix[, "CHROM"]
  snp_map <- data.frame(
    id  = make_snp_ids(fix[, "ID"], chr_raw, fix[, "POS"]),
    chr = parse_chr(chr_raw),
    pos = suppressWarnings(as.numeric(fix[, "POS"])),
    stringsAsFactors = FALSE
  )
  rownames(geno) <- snp_map$id
  list(geno = geno, snp_map = snp_map)
}

#' Convert a character GT matrix ("0/0","0/1","1/1","./.") to a 0/1/2 dosage matrix
#' @noRd
gt_to_dosage <- function(gt_char) {
  geno <- matrix(
    NA_integer_,
    nrow = nrow(gt_char), ncol = ncol(gt_char),
    dimnames = dimnames(gt_char)
  )
  # Normalise phased "|" to unphased "/"
  g <- gsub("\\|", "/", gt_char)
  geno[g == "0/0"] <- 0L
  geno[g == "1/1"] <- 2L
  geno[g == "0/1" | g == "1/0"] <- 1L
  geno
}

# ---------------------------------------------------------------------------
# HapMap (.hmp.txt)
# ---------------------------------------------------------------------------

#' @noRd
read_hapmap <- function(path, opts = list(), logger = NULL) {
  logger <- as_logger(logger)
  logger$log("Reading HapMap (.hmp.txt) ...")
  dt <- data.table::fread(path, header = TRUE, data.table = FALSE,
                          colClasses = "character")

  std <- c("rs#", "alleles", "chrom", "pos")
  if (!all(std %in% names(dt))) {
    stop(sprintf(
      "HapMap file is missing required column(s): %s",
      paste(setdiff(std, names(dt)), collapse = ", ")
    ), call. = FALSE)
  }
  # Standard HapMap has 11 leading metadata columns; samples follow.
  meta_cols <- intersect(
    c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center",
      "protLSID", "assayLSID", "panelLSID", "QCcode"),
    names(dt)
  )
  sample_cols <- setdiff(names(dt), meta_cols)
  if (length(sample_cols) == 0) {
    stop("No sample columns detected in HapMap file.", call. = FALSE)
  }

  alleles <- strsplit(dt[["alleles"]], "/")
  calls <- as.matrix(dt[, sample_cols, drop = FALSE])
  geno <- hapmap_calls_to_dosage(calls, alleles)

  snp_map <- data.frame(
    id  = make_snp_ids(dt[["rs#"]], dt[["chrom"]], dt[["pos"]]),
    chr = parse_chr(dt[["chrom"]]),
    pos = suppressWarnings(as.numeric(dt[["pos"]])),
    stringsAsFactors = FALSE
  )
  rownames(geno) <- snp_map$id
  list(geno = geno, snp_map = snp_map)
}

#' Convert HapMap nucleotide calls to 0/1/2 (count of the minor/second allele)
#'
#' Accepts single-character IUPAC calls (A,C,G,T plus het codes R,Y,S,W,K,M) and
#' two-character calls ("AA","AG"). The first listed allele is treated as the
#' reference (dosage 0); the second as the alternate (dosage 2).
#' @noRd
hapmap_calls_to_dosage <- function(calls, alleles) {
  iupac_het <- c("R", "Y", "S", "W", "K", "M")
  na_codes <- c("N", "NN", "--", "..", "", "NA")
  geno <- matrix(NA_integer_, nrow = nrow(calls), ncol = ncol(calls),
                 dimnames = dimnames(calls))
  for (i in seq_len(nrow(calls))) {
    al <- alleles[[i]]
    if (length(al) < 2 || any(is.na(al))) next
    ref <- al[1]
    alt <- al[2]
    row <- toupper(calls[i, ])
    out <- rep(NA_integer_, length(row))
    is_na <- row %in% na_codes
    # two-character calls
    homo_ref <- row %in% c(paste0(ref, ref))
    homo_alt <- row %in% c(paste0(alt, alt))
    het2 <- row %in% c(paste0(ref, alt), paste0(alt, ref))
    # single-character calls
    homo_ref <- homo_ref | row == ref
    homo_alt <- homo_alt | row == alt
    het1 <- row %in% iupac_het
    out[homo_ref] <- 0L
    out[homo_alt] <- 2L
    out[het2 | het1] <- 1L
    out[is_na] <- NA_integer_
    geno[i, ] <- out
  }
  geno
}

# ---------------------------------------------------------------------------
# Numeric matrix (CSV / RDS / RData) of 0/1/2 dosages
# ---------------------------------------------------------------------------

#' @noRd
read_numeric <- function(path, opts = list(), logger = NULL) {
  logger <- as_logger(logger)
  logger$log("Reading numeric genotype data ...")
  obj <- load_any(path)

  # Case 1: already in the internal representation.
  if (is.list(obj) && all(c("geno", "snp_map") %in% names(obj))) {
    logger$log("Detected pre-built list(geno, snp_map).")
    geno <- as_int_matrix(obj$geno)
    return(list(geno = geno, snp_map = normalise_map(obj$snp_map)))
  }

  df <- as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)

  # Column names may be user-supplied or auto-detected from common conventions
  # (HapMap / TASSEL / HBP numeric exports).
  chr_col <- opts$chr_col %||% detect_col(names(df), c("chrom", "chr", "chromosome"))
  pos_col <- opts$pos_col %||% detect_col(names(df), c("pos", "position", "bp", "pos_bp"))
  id_col  <- opts$id_col  %||% detect_col(names(df), c("rs#", "rs", "marker", "markername", "snp", "id", "name"))

  if (!is.null(chr_col) && !is.null(pos_col)) {
    # Wide table with metadata columns + sample columns (HBP/HapMap convention).
    # Standard HapMap/TASSEL metadata columns are never samples.
    known_meta <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                    "center", "protLSID", "assayLSID", "panelLSID", "QCcode",
                    "REFERENCE_GENOME")
    meta <- unique(stats::na.omit(
      c(id_col, chr_col, pos_col, opts$drop_cols, intersect(names(df), known_meta))))
    sample_cols <- setdiff(names(df), meta)
    if (length(sample_cols) == 0) {
      stop("No sample (genotype) columns detected after removing metadata columns.",
           call. = FALSE)
    }
    # Wide table already has SNPs in rows and samples in columns: no transpose.
    geno <- as_int_matrix(df[, sample_cols, drop = FALSE])  # SNP x ind
    if (isTRUE(opts$clean_sample_names)) {
      colnames(geno) <- sub(":.*$", "", colnames(geno))
    }
    ids <- if (!is.null(id_col)) df[[id_col]] else NULL
    snp_map <- data.frame(
      id  = make_snp_ids(ids, df[[chr_col]], df[[pos_col]]),
      chr = parse_chr(df[[chr_col]]),
      pos = suppressWarnings(as.numeric(df[[pos_col]])),
      stringsAsFactors = FALSE
    )
    rownames(geno) <- snp_map$id
    logger$log(sprintf(
      "Numeric wide table: id='%s', chr='%s', pos='%s'; %d sample column(s).",
      id_col %||% "(none)", chr_col, pos_col, length(sample_cols)))
    return(list(geno = geno, snp_map = snp_map))
  }

  stop(paste0(
    "Numeric input needs a SNP map. Either supply an RDS/RData containing ",
    "list(geno, snp_map), or a table with recognisable chromosome and position ",
    "columns (e.g. 'chrom' and 'pos'), or specify chr_col / pos_col."
  ), call. = FALSE)
}

#' Find the first column whose lower-cased name matches one of `cands`
#' @noRd
detect_col <- function(nms, cands) {
  hit <- which(tolower(nms) %in% tolower(cands))
  if (length(hit)) nms[hit[1]] else NULL
}

#' Load a CSV / RDS / RData file into a single R object
#' @noRd
load_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    return(readRDS(path))
  }
  if (ext %in% c("rdata", "rda")) {
    e <- new.env()
    nms <- load(path, envir = e)
    return(get(nms[1], envir = e))
  }
  data.table::fread(path, header = TRUE, data.table = FALSE, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# PLINK (.raw from --recodeA, or .ped/.map)
# ---------------------------------------------------------------------------

#' @noRd
read_plink <- function(path, opts = list(), logger = NULL) {
  logger <- as_logger(logger)
  ext <- tolower(tools::file_ext(path))
  if (ext == "raw") {
    return(read_plink_raw(path, logger))
  }
  if (ext == "ped") {
    return(read_plink_ped(path, opts, logger))
  }
  stop(paste0(
    "PLINK support reads '.raw' files (from `plink --recodeA`) or '.ped' + ",
    "'.map' text files. Binary '.bed' is not supported in this version; ",
    "please export with `plink --recodeA` first."
  ), call. = FALSE)
}

#' Read a PLINK .raw additive-coded file (individuals in rows)
#' @noRd
read_plink_raw <- function(path, logger = NULL) {
  logger <- as_logger(logger)
  logger$log("Reading PLINK .raw (additive coding) ...")
  df <- data.table::fread(path, header = TRUE, data.table = FALSE, check.names = FALSE)
  lead <- intersect(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"), names(df))
  marker_cols <- setdiff(names(df), lead)
  ind_ids <- if ("IID" %in% names(df)) df[["IID"]] else as.character(seq_len(nrow(df)))
  geno <- t(as_int_matrix(df[, marker_cols, drop = FALSE]))  # SNP x ind
  colnames(geno) <- ind_ids
  # .raw marker names look like "rsID_A" (effect allele suffix); strip suffix.
  ids <- sub("_[ACGT0-9]+$", "", marker_cols)
  snp_map <- data.frame(
    id  = ids,
    chr = NA_integer_,
    pos = NA_real_,
    stringsAsFactors = FALSE
  )
  logger$log(paste0(
    "No map in .raw: chromosome/position unknown. Supply a .map/.bim alongside ",
    "for a real genetic map, otherwise markers are placed by order."
  ), "WARN")
  rownames(geno) <- snp_map$id
  list(geno = geno, snp_map = snp_map)
}

#' Read PLINK .ped + .map text files
#' @noRd
read_plink_ped <- function(path, opts = list(), logger = NULL) {
  logger <- as_logger(logger)
  map_path <- opts$map_path %||% sub("\\.ped$", ".map", path)
  if (!file.exists(map_path)) {
    stop(sprintf("Could not find the .map file expected at: %s", map_path),
         call. = FALSE)
  }
  logger$log("Reading PLINK .ped/.map ...")
  map <- data.table::fread(map_path, header = FALSE, data.table = FALSE)
  # .map columns: chr, snp id, genetic distance (cM), bp position
  n_snp <- nrow(map)
  ped <- data.table::fread(path, header = FALSE, data.table = FALSE,
                           colClasses = "character")
  ind_ids <- ped[[2]]
  allele_block <- as.matrix(ped[, 7:ncol(ped), drop = FALSE])

  geno <- matrix(NA_integer_, nrow = n_snp, ncol = nrow(ped))
  for (j in seq_len(n_snp)) {
    a1 <- allele_block[, 2 * j - 1]
    a2 <- allele_block[, 2 * j]
    geno[j, ] <- ped_pair_to_dosage(a1, a2)
  }
  snp_map <- data.frame(
    id  = as.character(map[[2]]),
    chr = parse_chr(map[[1]]),
    pos = suppressWarnings(as.numeric(map[[4]])),
    stringsAsFactors = FALSE
  )
  colnames(geno) <- ind_ids
  rownames(geno) <- snp_map$id
  list(geno = geno, snp_map = snp_map)
}

#' Convert a pair of allele vectors to 0/1/2 (count of the minor allele)
#' @noRd
ped_pair_to_dosage <- function(a1, a2) {
  miss <- c("0", "N", ".", "-", "")
  is_na <- a1 %in% miss | a2 %in% miss
  alleles <- sort(unique(c(a1[!is_na], a2[!is_na])))
  out <- rep(NA_integer_, length(a1))
  if (length(alleles) == 0) {
    return(out)
  }
  # Minor allele = the rarer one; default to second when tied.
  counts <- table(c(a1[!is_na], a2[!is_na]))
  minor <- names(counts)[which.min(counts)]
  out[!is_na] <- (a1[!is_na] == minor) + (a2[!is_na] == minor)
  out[is_na] <- NA_integer_
  as.integer(out)
}

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Coerce a data.frame / matrix to an integer matrix
#' @noRd
as_int_matrix <- function(x) {
  m <- as.matrix(x)
  storage.mode(m) <- "integer"
  m
}

#' Parse chromosome labels to integers where possible (strip "Chr"/"chr" prefix)
#' @noRd
parse_chr <- function(x) {
  x <- as.character(x)
  x_clean <- sub("^(chr|Chr|CHR|SBI-)", "", x)
  num <- suppressWarnings(as.integer(x_clean))
  if (all(!is.na(num) | is.na(x))) {
    return(num)
  }
  # Fall back to factor-coded integers for non-numeric chromosome names.
  as.integer(factor(x))
}

#' Build stable SNP ids, falling back to chr_pos when ids are missing/blank
#' @noRd
make_snp_ids <- function(ids, chr, pos) {
  fallback <- paste0("S", chr, "_", pos)
  if (is.null(ids)) {
    return(fallback)
  }
  ids <- as.character(ids)
  bad <- is.na(ids) | ids %in% c("", ".", "NA")
  ids[bad] <- fallback[bad]
  make.unique(ids)
}

#' Normalise a user-supplied SNP map to columns id, chr, pos
#' @noRd
normalise_map <- function(map) {
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  nm <- tolower(names(map))
  pick <- function(cands) {
    hit <- which(nm %in% cands)
    if (length(hit)) names(map)[hit[1]] else NA_character_
  }
  id_c  <- pick(c("id", "rs#", "rs", "marker", "markername", "snp", "name"))
  chr_c <- pick(c("chr", "chrom", "chromosome"))
  pos_c <- pick(c("pos", "position", "bp", "pos_bp"))
  if (is.na(chr_c) || is.na(pos_c)) {
    stop("Supplied snp_map must contain chromosome and position columns.",
         call. = FALSE)
  }
  data.frame(
    id  = if (!is.na(id_c)) as.character(map[[id_c]]) else
      make_snp_ids(NULL, map[[chr_c]], map[[pos_c]]),
    chr = parse_chr(map[[chr_c]]),
    pos = suppressWarnings(as.numeric(map[[pos_c]])),
    stringsAsFactors = FALSE
  )
}
