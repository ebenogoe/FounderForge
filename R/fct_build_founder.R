# ---------------------------------------------------------------------------
# Build an AlphaSimR founder population from a 0/1/2 genotype matrix.
#
# Primary AlphaSimR inputs produced here:
#   - genMap     : list of per-chromosome Morgan-position vectors (start at 0)
#   - haplotypes : list of per-chromosome haplotype matrices
# These are the objects newMapPop() consumes and are saved as first-class
# outputs so a user can rebuild or audit the founder population.
# ---------------------------------------------------------------------------

#' Order genotypes and map by chromosome then physical position
#' @noRd
order_by_map <- function(geno, snp_map) {
  ord <- order(snp_map$chr, snp_map$pos)
  list(geno = geno[ord, , drop = FALSE], snp_map = snp_map[ord, , drop = FALSE])
}

#' Build a per-chromosome genetic map in Morgans
#'
#' Converts physical positions to Morgans assuming a uniform recombination rate
#' (`cm_per_mb` centiMorgans per megabase; 1 cM/Mb = pos_bp / 1e8). Each
#' chromosome is shifted to start at exactly 0 Morgans (required by newMapPop);
#' this shift is constant within a chromosome and so preserves every
#' inter-marker genetic distance. Tied positions are nudged by a cumulative
#' 1e-9 Morgan epsilon to keep positions strictly increasing.
#'
#' NOTE: the uniform-recombination map is an APPROXIMATION. Supply a real
#' genetic map if one is available for your species.
#'
#' @param snp_map data.frame (id, chr, pos), already ordered by chr then pos.
#' @param cm_per_mb Assumed recombination rate (cM per Mb). Default 1.
#' @param logger Optional logger.
#' @return list(map_list, map_df) where map_list is a named list of Morgan
#'   vectors and map_df is the data.frame (markerName, chromosome, position).
#' @noRd
build_genmap <- function(snp_map, cm_per_mb = 1, logger = NULL) {
  logger <- as_logger(logger)
  chrs <- unique(snp_map$chr)
  map_list <- vector("list", length(chrs))
  names(map_list) <- as.character(chrs)
  factor_bp_to_M <- cm_per_mb / 1e8  # cM/Mb -> Morgans per bp

  for (i in seq_along(chrs)) {
    idx <- which(snp_map$chr == chrs[i])
    pos_M <- snp_map$pos[idx] * factor_bp_to_M
    pos_M <- pos_M - min(pos_M)  # start chromosome at 0 Morgans
    # enforce strictly increasing positions
    for (k in seq_along(pos_M)[-1]) {
      if (pos_M[k] <= pos_M[k - 1L]) {
        pos_M[k] <- pos_M[k - 1L] + 1e-9
      }
    }
    names(pos_M) <- snp_map$id[idx]
    map_list[[i]] <- pos_M
  }

  map_df <- data.frame(
    markerName = snp_map$id,
    chromosome = snp_map$chr,
    position   = unlist(map_list, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  logger$log(sprintf(
    "Built genetic map: %d chromosome(s), %s markers (%.0f cM/Mb, uniform).",
    length(chrs), format(nrow(snp_map), big.mark = ","), cm_per_mb
  ))
  list(map_list = map_list, map_df = map_df)
}

#' Build per-chromosome haplotype matrices for newMapPop()
#'
#' @param geno Integer matrix SNPs x individuals, ordered by chr then pos.
#' @param snp_map data.frame (id, chr, pos), same order as geno rows.
#' @param strategy "outbred" (split each genotype into two haplotypes,
#'   preserving heterozygosity) or "inbred" (one haplotype per individual,
#'   heterozygous calls coerced to the nearest homozygote).
#' @param logger Optional logger.
#' @return A named list of haplotype matrices (one per chromosome).
#' @noRd
build_haplotypes <- function(geno, snp_map, strategy = c("outbred", "inbred"),
                             logger = NULL) {
  logger <- as_logger(logger)
  strategy <- match.arg(strategy)
  chrs <- unique(snp_map$chr)
  hap_list <- vector("list", length(chrs))
  names(hap_list) <- as.character(chrs)

  if (strategy == "outbred") {
    for (i in seq_along(chrs)) {
      g <- geno[snp_map$chr == chrs[i], , drop = FALSE]  # nLoci x nInd
      hap1 <- t(g %/% 2L)        # nInd x nLoci : 0->0, 1->0, 2->1
      hap2 <- t(g - (g %/% 2L))  # nInd x nLoci : 0->0, 1->1, 2->1
      n_ind <- nrow(hap1)
      mat <- matrix(0L, nrow = 2L * n_ind, ncol = ncol(hap1))
      mat[seq(1L, 2L * n_ind, by = 2L), ] <- hap1
      mat[seq(2L, 2L * n_ind, by = 2L), ] <- hap2
      hap_list[[i]] <- mat
    }
    logger$log("Built haplotypes (outbred: 2 haplotypes/individual, het preserved).")
  } else {
    n_het <- sum(geno == 1L, na.rm = TRUE)
    if (n_het > 0) {
      logger$log(sprintf(
        "Inbred strategy: %s heterozygous call(s) coerced to nearest homozygote.",
        format(n_het, big.mark = ",")
      ), "WARN")
    }
    for (i in seq_along(chrs)) {
      g <- geno[snp_map$chr == chrs[i], , drop = FALSE]  # nLoci x nInd
      hap_list[[i]] <- t(round(g / 2))  # nInd x nLoci : 0->0, 1->0/1, 2->1
      storage.mode(hap_list[[i]]) <- "integer"
    }
    logger$log("Built haplotypes (inbred: 1 haplotype/individual).")
  }
  hap_list
}

#' Build an AlphaSimR founder population
#'
#' @param geno Integer matrix SNPs x individuals, values 0/1/2 (no NA).
#' @param snp_map data.frame (id, chr, pos), row-aligned to geno.
#' @param ploidy Ploidy level (default 2L).
#' @param strategy One of "outbred", "inbred", "import_inbred". See the upload
#'   step tooltips for the implications of each.
#' @param cm_per_mb Assumed recombination rate (cM/Mb) for the genetic map.
#' @param logger Optional logger.
#' @return list with founderPop, genMap (map_list), genMap_df, haplotypes,
#'   strategy, ind_ids, snp_map (ordered) and call_string (reproducible snippet).
#' @noRd
build_founder_pop <- function(geno, snp_map, ploidy = 2L,
                              strategy = c("outbred", "inbred", "import_inbred"),
                              cm_per_mb = 1, logger = NULL) {
  logger <- as_logger(logger)
  strategy <- match.arg(strategy)
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop("Building a founder population requires the 'AlphaSimR' package.",
         call. = FALSE)
  }
  if (anyNA(geno)) {
    stop("Genotype matrix still contains NA; run impute_missing() first.",
         call. = FALSE)
  }
  logger$section("Building AlphaSimR founder population")

  ord <- order_by_map(geno, snp_map)
  geno <- ord$geno
  snp_map <- ord$snp_map
  ind_ids <- colnames(geno)

  gm <- build_genmap(snp_map, cm_per_mb = cm_per_mb, logger = logger)

  if (strategy == "import_inbred") {
    geno_ind <- t(geno)  # individuals x SNPs, colnames = marker ids
    colnames(geno_ind) <- snp_map$id
    founderPop <- AlphaSimR::importInbredGeno(geno = geno_ind, genMap = gm$map_df)
    haplotypes <- NULL
    call_string <- paste0(
      "founderPop <- AlphaSimR::importInbredGeno(\n",
      "  geno   = geno,        # individuals x markers, 0/1/2\n",
      "  genMap = genMap_df    # data.frame(markerName, chromosome, position)\n",
      ")"
    )
  } else {
    inbred <- identical(strategy, "inbred")
    haplotypes <- build_haplotypes(
      geno, snp_map,
      strategy = if (inbred) "inbred" else "outbred",
      logger = logger
    )
    founderPop <- AlphaSimR::newMapPop(
      genMap = gm$map_list,
      haplotypes = haplotypes,
      inbred = inbred,
      ploidy = as.integer(ploidy)
    )
    call_string <- sprintf(paste0(
      "founderPop <- AlphaSimR::newMapPop(\n",
      "  genMap     = genMap,       # list of per-chromosome Morgan vectors\n",
      "  haplotypes = haplotypes,   # list of per-chromosome haplotype matrices\n",
      "  inbred     = %s,\n",
      "  ploidy     = %dL\n",
      ")"), if (inbred) "TRUE" else "FALSE", as.integer(ploidy))
  }

  logger$log(sprintf(
    "Founder population created: %d individuals, %d chromosome(s), ploidy %d.",
    founderPop@nInd, founderPop@nChr, founderPop@ploidy
  ), "OK")

  list(
    founderPop = founderPop,
    genMap = gm$map_list,
    genMap_df = gm$map_df,
    haplotypes = haplotypes,
    strategy = strategy,
    ploidy = as.integer(ploidy),
    cm_per_mb = cm_per_mb,
    ind_ids = ind_ids,
    snp_map = snp_map,
    call_string = call_string
  )
}
