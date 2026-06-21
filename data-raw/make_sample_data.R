# Generates bundled sample datasets in inst/extdata for users to test FounderForge.
set.seed(42)
nChr <- 3L; perChr <- 100L; nSNP <- nChr * perChr; nInd <- 30L

chr <- rep(seq_len(nChr), each = perChr)
pos <- as.numeric(unlist(lapply(seq_len(nChr), function(k) sort(sample(1e4:5e6, perChr)))))
ids <- paste0("S", chr, "_", pos)
ind_ids <- sprintf("LINE_%03d", seq_len(nInd))

# Outbred-style allele freqs -> 0/1/2 dosages with genuine heterozygotes
af <- runif(nSNP, 0.1, 0.9)
geno <- matrix(0L, nrow = nSNP, ncol = nInd, dimnames = list(ids, ind_ids))
for (i in seq_len(nSNP)) {
  geno[i, ] <- rbinom(nInd, 1, af[i]) + rbinom(nInd, 1, af[i])
}
nucs <- c("A", "C", "G", "T")
ref <- sample(nucs, nSNP, TRUE)
alt <- vapply(ref, function(r) sample(setdiff(nucs, r), 1), character(1))
alleles <- paste(ref, alt, sep = "/")

out <- "inst/extdata"

## 1. RDS: pre-built list(geno, snp_map) -- the simplest path
saveRDS(list(geno = geno,
             snp_map = data.frame(id = ids, chr = chr, pos = pos,
                                  stringsAsFactors = FALSE)),
        file.path(out, "sample_numeric.rds"))

## 2. CSV: HapMap-numeric wide table (HBP-style) with metadata + samples
csv <- data.frame(`rs#` = ids, alleles = alleles, chrom = chr, pos = pos,
                  check.names = FALSE, stringsAsFactors = FALSE)
csv <- cbind(csv, as.data.frame(geno, check.names = FALSE))
data.table::fwrite(csv, file.path(out, "sample_numeric.csv"))

## 3. VCF (biallelic)
gt <- ifelse(geno == 0, "0/0", ifelse(geno == 1, "0/1", "1/1"))
con <- file(file.path(out, "sample.vcf"), "w")
writeLines(c("##fileformat=VCFv4.2",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                     ind_ids), collapse = "\t")), con)
for (i in seq_len(nSNP)) {
  writeLines(paste(c(chr[i], pos[i], ids[i], ref[i], alt[i], ".", "PASS", ".", "GT",
                     gt[i, ]), collapse = "\t"), con)
}
close(con)

## 4. HapMap (.hmp.txt) nucleotide calls
two <- function(d, r, a) ifelse(d == 0, paste0(r, r), ifelse(d == 2, paste0(a, a), paste0(r, a)))
hmp_geno <- t(vapply(seq_len(nSNP), function(i) two(geno[i, ], ref[i], alt[i]), character(nInd)))
hmp <- data.frame(`rs#` = ids, alleles = alleles, chrom = chr, pos = pos,
                  strand = "+", `assembly#` = NA, center = NA, protLSID = NA,
                  assayLSID = NA, panelLSID = NA, QCcode = NA,
                  check.names = FALSE, stringsAsFactors = FALSE)
hmp <- cbind(hmp, as.data.frame(hmp_geno, stringsAsFactors = FALSE))
names(hmp)[12:(11 + nInd)] <- ind_ids
data.table::fwrite(hmp, file.path(out, "sample.hmp.txt"), sep = "\t")

cat("Wrote sample files:\n"); print(list.files(out))
cat(sprintf("Dimensions: %d SNPs x %d individuals, %d chromosomes\n", nSNP, nInd, nChr))
