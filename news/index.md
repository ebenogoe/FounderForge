# Changelog

## FounderForge 0.1.0

- First release.
- Point-and-click Shiny wizard to prepare marker data (VCF, HapMap,
  numeric matrices, PLINK) as an AlphaSimR founder population.
- Quality-control review, flexible filtering (MAF, missingness,
  monomorphic markers, minimum SNPs per chromosome, individual
  dropping), optional additive trait definition, and optional
  before/after PCA validation.
- Organised, timestamped outputs (`genmap.rds`, `haplotypes.rds`,
  `founder_pop.rds`, optional `sim_param.rds`, a reproducible
  `rebuild_snippet.R`, QC/validation plots and a full run log).
