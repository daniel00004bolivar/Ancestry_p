#!/usr/bin/env Rscript
# 17_summarize_mixemt.R — consolidates mixemt's mtDNA haplogroup mixtures
# (bin/*_mid.log, one per pool) into a table comparable to NGSadmix's
# autosomal results, grouped by continental macro-haplogroup.
#
# EXPLORATORY: run on reads downsampled to ~30,000 per pool (not full
# depth) so mixemt finishes in a reasonable time. Macro-haplogroup
# classification uses the 5 Native American founder lineages recognized in
# the literature (A2, B2, C1, D1, X2a) and pan-Eurasian/African lineages by
# prefix, with known exceptions (e.g. D4l2 belongs to the Asian D4 clade,
# NOT the Native American founder subclade D4h3a -- classifying purely by
# the "D" prefix would be wrong).
#
# This is a preliminary result from a single run per pool, with no
# confidence interval -- it does not replace validation with a dedicated
# haplogroup-calling tool (e.g. HaploGrep2) on non-downsampled data.

suppressPackageStartupMessages(library(data.table))

MIXEMT_DIR <- "out_global/mixemt"
POOLS <- c("POOL1", "POOL2", "HOSPITAL")
OUT_FILE <- "results/mtdna_haplogroups.tsv"
OUT_SUMMARY <- "results/mtdna_continental_summary.tsv"

# Continental macro-haplogroup classification.
# NAT = Native American founders (A2, B2, C1, D1, X2a)
# EUR = pan-European/Near Eastern clades (H, HV, I, J, K, T, U, V, W)
# AFR = clade L (L0-L6), exclusively African
# EAS = other Asian clades that are NOT Native American founders (e.g. non-h3a D4)
classify_haplogroup <- function(hap) {
    h <- gsub("\\[.*\\]", "", hap)  # strip suffixes like [1], [1][2]
    if (grepl("^D4l", h)) return("EAS")          # Asian D4, NOT Native American D4h3a
    if (grepl("^(A2|B2|C1|D1|X2a)", h)) return("NAT")
    if (grepl("^L", h)) return("AFR")
    if (grepl("^(H|HV|I|J|K|T|U|V|W)", h)) return("EUR")
    return("OTHER")
}

parse_pool <- function(pool) {
    log_file <- file.path(MIXEMT_DIR, sprintf("%s_mid.log", pool))
    if (!file.exists(log_file)) stop(sprintf("%s does not exist -- run mixemt first", log_file))
    lines <- grep("^hap[0-9]+", readLines(log_file), value = TRUE)
    dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE,
                col.names = c("hap_id", "haplogroup", "proportion", "n_reads"))
    dt[, Pool := pool]
    dt[, Continental := sapply(haplogroup, classify_haplogroup)]
    dt
}

all_haps <- rbindlist(lapply(POOLS, parse_pool))
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
fwrite(all_haps[, .(Pool, hap_id, haplogroup, proportion, n_reads, Continental)], OUT_FILE, sep = "\t")

summary_dt <- all_haps[, .(proportion = sum(proportion)), by = .(Pool, Continental)]
summary_wide <- dcast(summary_dt, Pool ~ Continental, value.var = "proportion", fill = 0)
fwrite(summary_wide, OUT_SUMMARY, sep = "\t")

cat(sprintf("Written: %s (%d haplogroups)\n", OUT_FILE, nrow(all_haps)))
cat(sprintf("Written: %s (continental summary)\n\n", OUT_SUMMARY))
print(summary_wide)

cat("\nComparison with autosomal NGSadmix K=5 (from the manuscript):\n")
cat("  POOL1    autosomal: EUR 47.2%  AMR 31.1%  AFR 16.5%\n")
cat("  POOL2    autosomal: EUR 47.9%  AMR 31.0%  AFR 16.5%\n")
cat("  HOSPITAL autosomal: EUR 36.6%  AMR 31.9%  AFR 28.2%\n")
