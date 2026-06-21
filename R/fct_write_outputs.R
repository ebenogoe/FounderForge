# ---------------------------------------------------------------------------
# Write the organised output tree for a run
# ---------------------------------------------------------------------------

#' Write all run outputs into an organised, timestamped folder
#'
#' @param out_dir Parent folder chosen by the user.
#' @param build_result Result of [build_founder_pop()].
#' @param SP Optional SimParam from [define_traits()].
#' @param qc_summary Optional QC summary from [compute_qc_summary()].
#' @param filter_report Optional before/after data.frame from [filter_genotypes()].
#' @param verify_result Optional result of [verify_founder()].
#' @param pca_result Optional result of [pca_validation()].
#' @param source_info Named list describing the input (path, format, etc).
#' @param logger Logger whose lines become the text log.
#' @return The path to the created run folder.
#' @noRd
write_outputs <- function(out_dir, build_result, SP = NULL,
                          qc_summary = NULL, filter_report = NULL,
                          verify_result = NULL, pca_result = NULL,
                          source_info = list(), logger = NULL) {
  logger <- as_logger(logger)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(out_dir, paste0("FounderForge_run_", stamp))

  dirs <- c("founder_population", "qc", "validation", "logs")
  for (d in dirs) dir.create(file.path(run_dir, d), recursive = TRUE, showWarnings = FALSE)
  logger$section("Writing outputs")
  logger$log(sprintf("Output folder: %s", run_dir))

  fp_dir <- file.path(run_dir, "founder_population")
  saveRDS(build_result$genMap, file.path(fp_dir, "genmap.rds"))
  if (!is.null(build_result$haplotypes)) {
    saveRDS(build_result$haplotypes, file.path(fp_dir, "haplotypes.rds"))
  }
  saveRDS(build_result$genMap_df, file.path(fp_dir, "genmap_table.rds"))
  saveRDS(build_result$founderPop, file.path(fp_dir, "founder_pop.rds"))
  if (!is.null(SP)) {
    saveRDS(SP, file.path(fp_dir, "sim_param.rds"))
  }
  writeLines(make_rebuild_snippet(build_result, has_sp = !is.null(SP)),
             file.path(fp_dir, "rebuild_snippet.R"))
  logger$log("Saved genmap, haplotypes, founder_pop (+ sim_param) and rebuild snippet.")

  if (!is.null(qc_summary)) {
    write_qc_outputs(file.path(run_dir, "qc"), qc_summary, filter_report, logger)
  }

  val_dir <- file.path(run_dir, "validation")
  if (!is.null(verify_result)) {
    writeLines(format_verify(verify_result, build_result$strategy),
               file.path(val_dir, "heterozygosity_check.txt"))
  }
  if (!is.null(pca_result)) {
    save_plot(pca_result$plot, file.path(val_dir, "pca_before_after.png"),
              width = 11, height = 5)
    utils::write.csv(pca_result$variance,
                     file.path(val_dir, "pca_variance.csv"), row.names = FALSE)
    logger$log("Saved before/after PCA plot and variance table.")
  }

  # Final log with an explicit assumptions section.
  for (line in assumptions_section(build_result, SP)) logger$log(line)
  log_path <- file.path(run_dir, "logs", paste0("run_log_", stamp, ".txt"))
  writeLines(c(log_header(source_info), "", logger$lines()), log_path)
  logger$log(sprintf("Run log written: %s", log_path), "OK")

  run_dir
}

#' @noRd
write_qc_outputs <- function(qc_dir, qc, filter_report, logger) {
  summary_df <- data.frame(
    metric = c("n_snp", "n_individuals", "n_chromosomes", "pct_missing",
               "n_monomorphic", "mean_maf", "mean_heterozygosity"),
    value = c(qc$n_snp, qc$n_ind, qc$n_chr, round(qc$pct_missing, 4),
              qc$n_monomorphic, round(qc$mean_maf, 4), round(qc$mean_het, 4))
  )
  utils::write.csv(summary_df, file.path(qc_dir, "qc_summary.csv"), row.names = FALSE)
  if (!is.null(filter_report)) {
    utils::write.csv(filter_report, file.path(qc_dir, "filter_report.csv"),
                     row.names = FALSE)
  }
  save_plot(hist_plot(qc$maf, "Minor allele frequency", "MAF"),
            file.path(qc_dir, "maf_hist.png"))
  save_plot(hist_plot(qc$miss_snp, "Per-SNP missing rate", "Missing rate"),
            file.path(qc_dir, "missingness.png"))
  save_plot(snps_per_chr_plot(qc$snps_per_chr),
            file.path(qc_dir, "snps_per_chr.png"))
  logger$log("Saved QC summary CSV and plots (MAF, missingness, SNPs/chr).")
}

#' @noRd
hist_plot <- function(x, title, xlab) {
  df <- data.frame(x = x[is.finite(x)])
  ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_histogram(bins = 40, fill = "#378ADD", colour = "white") +
    ggplot2::labs(title = title, x = xlab, y = "Count") +
    ggplot2::theme_minimal(base_size = 13)
}

#' @noRd
snps_per_chr_plot <- function(spc) {
  ggplot2::ggplot(spc, ggplot2::aes(x = factor(.data$chr), y = .data$n_snps)) +
    ggplot2::geom_col(fill = "#1D9E75") +
    ggplot2::labs(title = "SNPs per chromosome", x = "Chromosome", y = "SNPs") +
    ggplot2::theme_minimal(base_size = 13)
}

#' @noRd
save_plot <- function(plot, path, width = 7, height = 5) {
  if (is.null(plot)) return(invisible(NULL))
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = 150)
}

#' @noRd
make_rebuild_snippet <- function(build_result, has_sp) {
  c(
    "# Rebuild this AlphaSimR founder population from the saved objects.",
    "# Generated by FounderForge.",
    "library(AlphaSimR)",
    "",
    "genMap     <- readRDS('genmap.rds')",
    if (!is.null(build_result$haplotypes)) "haplotypes <- readRDS('haplotypes.rds')" else
      "genMap_df  <- readRDS('genmap_table.rds')",
    "",
    "# The exact call FounderForge used:",
    build_result$call_string,
    "",
    if (has_sp) "SP <- readRDS('sim_param.rds')" else
      "SP <- SimParam$new(founderPop)",
    "pop <- newPop(founderPop, simParam = SP)"
  )
}

#' @noRd
format_verify <- function(v, strategy) {
  c(
    "Heterozygosity / structure check",
    "================================",
    sprintf("Import strategy : %s", strategy),
    sprintf("Expected dosages: %s", v$expected),
    sprintf("Observed dosages: %s", paste(v$values, collapse = ", ")),
    sprintf("Result          : %s", if (v$pass) "PASSED" else "FAILED"),
    sprintf("Segregating loci: %s", format(v$n_seg, big.mark = ","))
  )
}

#' @noRd
assumptions_section <- function(build_result, SP) {
  lines <- c(
    "",
    "ASSUMPTIONS & APPROXIMATIONS",
    "----------------------------",
    sprintf("- Genetic map: uniform recombination at %.0f cM/Mb (pos_bp/1e8). This is an approximation; supply a real map if available.",
            build_result$cm_per_mb)
  )
  if (identical(build_result$strategy, "outbred")) {
    lines <- c(lines,
      "- Haplotype phasing for heterozygous sites is arbitrary (genotypes split, not phased).")
  }
  if (!identical(build_result$strategy, "outbred")) {
    lines <- c(lines,
      "- Inbred strategy: any heterozygous calls were coerced to the nearest homozygote.")
  }
  if (!is.null(SP)) {
    lines <- c(lines,
      "- QTL positions are assigned at random by AlphaSimR.",
      "- Trait means/variances/heritabilities are user-supplied; flag any assumed values in your methods.")
  }
  lines
}

#' @noRd
log_header <- function(source_info) {
  c(
    "FounderForge run log",
    "====================",
    sprintf("Generated   : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Input file  : %s", source_info$path %||% "(unknown)"),
    sprintf("Format      : %s", source_info$format %||% "(unknown)")
  )
}
