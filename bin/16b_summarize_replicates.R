#!/usr/bin/env Rscript
# 16b_summarize_replicates.R — summarizes the chromosome jackknife
# (bin/16_pool_replicates.sh) into a table of jackknife mean +/- standard
# error per pool and component, in a format comparable to
# results/iadmix_summary.tsv.
#
# IMPORTANT -- label switching: NGSadmix does not guarantee that column K_i
# means the same ancestral cluster across different runs (even with the same
# seed, excluding a chromosome changes the EM optimum). Verified empirically:
# columns DO come out permuted between replicates (see bin/README_VIGENTE.md).
# Before averaging, each replicate is aligned against the full-genome
# K5.qopt (the one used in the manuscript) by trying all 120 possible
# permutations of 5 columns and keeping the one that minimizes the sum of
# squared errors against the reference.
#
# Jackknife standard error formula (delete-1-chromosome), after alignment:
#   SE = sqrt( (n-1)/n * sum( (theta_i - mean(theta_i))^2 ) )

suppressPackageStartupMessages(library(data.table))

JK_DIR   <- "out_global/angsd/pool_replicates"
IND_FILE <- "out_global/angsd/ngsadmix/individuals_full.tsv"
FULL_QOPT <- "out_global/angsd/ngsadmix/results_full/K5.qopt"
OUT_FILE <- "results/pool_jackknife_summary.tsv"

ind <- fread(IND_FILE)
pool_rows <- which(ind$super_pop == "POOL")
pool_names <- ind$sample[pool_rows]
stopifnot(length(pool_names) == 3)

qopt_files <- sprintf("%s/jk_excl_chr%d.qopt", JK_DIR, 1:22)
missing <- qopt_files[!file.exists(qopt_files)]
if (length(missing) > 0) {
    stop(sprintf("Missing %d of 22 jackknife files. Run bin/16_pool_replicates.sh first. Example missing: %s",
                 length(missing), missing[1]))
}

full <- as.matrix(fread(FULL_QOPT, header = FALSE))
K <- ncol(full)
comp_names <- paste0("K", seq_len(K))

# All permutations of 1:K (K=5 -> 120, trivial by brute force)
permute <- function(v) {
    if (length(v) <= 1) return(list(v))
    out <- list()
    for (i in seq_along(v)) {
        rest <- permute(v[-i])
        for (r in rest) out[[length(out) + 1]] <- c(v[i], r)
    }
    out
}
all_perms <- permute(seq_len(K))

align_to_reference <- function(q, ref) {
    best_sse <- Inf; best_perm <- seq_len(K)
    for (p in all_perms) {
        sse <- sum((q[, p, drop = FALSE] - ref)^2)
        if (sse < best_sse) { best_sse <- sse; best_perm <- p }
    }
    q[, best_perm, drop = FALSE]
}

# Array: [replicate, pool, component] -- already aligned to the reference
jk_array <- array(NA_real_, dim = c(22, length(pool_names), K))
for (i in seq_along(qopt_files)) {
    q <- as.matrix(fread(qopt_files[i], header = FALSE))
    q_aligned <- align_to_reference(q, full)
    jk_array[i, , ] <- q_aligned[pool_rows, ]
}

rows <- list()
for (p in seq_along(pool_names)) {
    for (k in seq_len(K)) {
        theta_i <- jk_array[, p, k]
        theta_bar <- mean(theta_i)
        se_jack <- sqrt((22 - 1) / 22 * sum((theta_i - theta_bar)^2))
        rows[[length(rows) + 1]] <- data.table(
            Pool = pool_names[p],
            Component = comp_names[k],
            Full_genome_estimate = full[pool_rows[p], k],
            Jackknife_mean = theta_bar,
            Jackknife_SE = se_jack
        )
    }
}
result <- rbindlist(rows)
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
fwrite(result, OUT_FILE, sep = "\t")

cat(sprintf("Written: %s (columns aligned by optimal permutation against the reference K5)\n\n", OUT_FILE))
print(result[, .(Pool, Component, Full_genome_estimate = round(Full_genome_estimate, 4),
                  Jackknife_mean = round(Jackknife_mean, 4), Jackknife_SE = round(Jackknife_SE, 4))])
