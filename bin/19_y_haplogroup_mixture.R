#!/usr/bin/env Rscript
# 19_y_haplogroup_mixture.R — estimates continental composition of Y
# (paternal) lineages per pool, via mixture deconvolution over diagnostic Y
# SNPs.
#
# No tool equivalent to mixemt exists for chromosome Y (mixemt is specific
# to the PhyloTree mitochondrial phylogeny). Instead, this script uses the
# same principle as 05c_continental_ancestry_v3.R (AIMs from ACGT counts)
# but solved as a constrained least-squares problem (non-negative, sum=1)
# instead of FST: Y SNPs with strong differentiation across the 5 1000G
# Phase 3 superpopulations are selected (frequencies taken directly from the
# Y VCF's INFO field, no need to extract individual genotypes), and the
# superpopulation mixture that best explains the observed allele frequency
# in each pool at those same sites is fit.
#
# Because chromosome Y is haploid and does not recombine, the "frequency"
# observed in a pool is directly the fraction of the pool's paternal
# lineages (Y chromosomes) carrying each allele -- a clean quantity, without
# the diploid pseudo-individual problem affecting NGSadmix/PCAngsd (see
# bin/README_VIGENTE.md).
#
# Requires having run beforehand, for each pool:
#   angsd -sites work/y_panel/y_markers.txt -doCounts 1 -dumpCounts 3 ...
# (see bin/README_VIGENTE.md, Y section, for the full command)

suppressPackageStartupMessages({
    library(data.table)
    library(quadprog)
})

POOLS <- c("POOL1", "POOL2", "HOSPITAL")
COUNTS_DIR <- "out_global/angsd_Y/per_pool"
FREQS_FILE <- "work/y_panel/y_freqs_sorted.txt"
OUT_FILE <- "results/y_haplogroup_mixture.tsv"

nuc_idx <- c(A = 1, C = 2, G = 3, T = 4)

ref <- fread(FREQS_FILE, header = FALSE,
             col.names = c("chr", "pos", "REF", "ALT", "AFR", "AMR", "EAS", "EUR", "SAS", "spread"))
setkey(ref, pos)

load_pool_freq <- function(pool) {
    pos_file <- file.path(COUNTS_DIR, sprintf("%s_counts.pos.gz", pool))
    cnt_file <- file.path(COUNTS_DIR, sprintf("%s_counts.counts.gz", pool))
    pos <- fread(pos_file)
    cnt <- fread(cnt_file)
    setnames(cnt, c("A", "C", "G", "T"))
    dt <- cbind(pos, cnt)
    m <- merge(dt, ref, by = "pos")
    ref_n <- mapply(function(r, a, c, g, t) c(A = a, C = c, G = g, T = t)[r], m$REF, m$A, m$C, m$G, m$T)
    alt_n <- mapply(function(r, a, c, g, t) c(A = a, C = c, G = g, T = t)[r], m$ALT, m$A, m$C, m$G, m$T)
    m[, ref_n := ref_n]
    m[, alt_n := alt_n]
    m <- m[ref_n + alt_n >= 3]  # minimum depth per site
    m[, obs_freq := alt_n / (ref_n + alt_n)]
    m[, .(pos, obs_freq, AFR, AMR, EAS, EUR, SAS)]
}

fit_mixture <- function(dt) {
    Xmat <- as.matrix(dt[, .(AFR, AMR, EAS, EUR, SAS)])
    y <- dt$obs_freq
    Dmat <- t(Xmat) %*% Xmat
    dvec <- t(Xmat) %*% y
    # Constraints: sum(q)=1 (equality), q_k >= 0 (inequality)
    Amat <- cbind(rep(1, 5), diag(5))
    bvec <- c(1, rep(0, 5))
    sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
    q <- pmax(sol$solution, 0); q <- q / sum(q)
    names(q) <- c("AFR", "AMR", "EAS", "EUR", "SAS")
    q
}

results <- list()
for (pool in POOLS) {
    dt <- load_pool_freq(pool)
    cat(sprintf("%s: %d sites with sufficient depth\n", pool, nrow(dt)))
    q <- fit_mixture(dt)
    results[[pool]] <- as.data.table(as.list(q))
    results[[pool]][, Pool := pool]
    results[[pool]][, n_sites := nrow(dt)]
}
out <- rbindlist(results)
setcolorder(out, c("Pool", "AFR", "AMR", "EAS", "EUR", "SAS", "n_sites"))
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT_FILE, sep = "\t")
cat(sprintf("\nWritten: %s\n\n", OUT_FILE))
print(out)
