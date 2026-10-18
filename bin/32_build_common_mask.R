#!/usr/bin/env Rscript
# 32_build_common_mask.R — one harmonized analysis table and common SNP mask
#
# Memo items A4 (orientation), A5 (zeros), A6 (one SNP universe), C1 (coverage).
#
# For every oriented site (bin/30) present in ALL three pools (bin/31), counts
# are recoded once as derived reads (x), ancestral + derived depth (d), and
# "other" reads (neither allele; used for the error rate and QC). A derived
# count of 0 is kept as an observation.
#
# Main mask (same sites for every pool and every analysis):
#   - d >= MIN_DEPTH in all three pools
#   - d <= 99th percentile of that pool's depth (collapsed duplications/CNV)
#   - other reads <= 5% of the site's depth in every pool
#   - source MAF >= 0.05, from the equal-weight mean of the AFR/EUR/NAT proxies
#     (a property of the reference, not of any pool)
# Sensitivity scripts re-filter the same table with other thresholds.
#
# Output: out_global/revision/analysis_table.tsv.gz
#         out_global/revision/qc_summary.tsv
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))

OUT_DIR   <- "out_global/revision"
POOLS     <- c("POOL1", "POOL2", "HOSPITAL")
MIN_DEPTH <- 20L
MAX_OTHER <- 0.05
MIN_MAF   <- 0.05

sites <- fread(file.path(OUT_DIR, "sites_oriented.tsv.gz"))

tab <- copy(sites)
for (p in POOLS) {
    ct <- fread(file.path(OUT_DIR, paste0("counts_", p, ".tsv.gz")))
    ct <- merge(sites[, .(chr, pos, anc, der)], ct, by = c("chr", "pos"))
    get_base <- function(b) fifelse(b == "A", ct$A, fifelse(b == "C", ct$C, fifelse(b == "G", ct$G, ct$T)))
    ct[, `:=`(x = get_base(der), a = get_base(anc))]
    ct[, `:=`(d = x + a, other = depth - x - a)]
    setnames(ct, c("x", "d", "other", "depth"), paste0(c("x_", "d_", "o_", "depth_"), p))
    tab <- merge(tab, ct[, c("chr", "pos", paste0(c("x_", "d_", "o_", "depth_"), p)), with = FALSE],
                 by = c("chr", "pos"))
}
setorder(tab, chr, pos)
cat(sprintf("Oriented sites present in all three pools (depth >= 5): %d of %d\n", nrow(tab), nrow(sites)))

src_mean <- tab[, (f_AFR + f_EUR + f_NAT) / 3]
tab[, maf_ok := pmin(src_mean, 1 - src_mean) >= MIN_MAF]
keep <- tab$maf_ok
for (p in POOLS) {
    d <- tab[[paste0("d_", p)]]; o <- tab[[paste0("o_", p)]]; dep <- tab[[paste0("depth_", p)]]
    maxd <- quantile(d[d >= MIN_DEPTH], 0.99)
    tab[, (paste0("maxd_", p)) := maxd]
    keep <- keep & d >= MIN_DEPTH & d <= maxd & o <= MAX_OTHER * dep
}
tab[, mask_main := keep]
cat(sprintf("Main mask: %d sites (%.1f%% transversions)\n", sum(keep), 100 * mean(tab$transversion[keep])))

qc <- rbindlist(lapply(POOLS, function(p) {
    d <- tab[[paste0("d_", p)]]; x <- tab[[paste0("x_", p)]]; o <- tab[[paste0("o_", p)]]
    m <- tab$mask_main
    data.table(pool = p,
               sites_all3 = nrow(tab),
               mean_depth_all3 = mean(d), median_depth_all3 = median(d),
               sites_mask = sum(m),
               mean_depth_mask = mean(d[m]), median_depth_mask = median(d[m]),
               max_depth_cut = tab[[paste0("maxd_", p)]][1],
               error_rate = 1.5 * sum(o[m]) / sum(d[m] + o[m]),
               frac_zero_derived = mean(x[m] == 0),
               mean_derived_freq = mean(x[m] / d[m]))
}))
fwrite(qc, file.path(OUT_DIR, "qc_summary.tsv"), sep = "\t")
print(qc, digits = 4)
fwrite(tab, file.path(OUT_DIR, "analysis_table.tsv.gz"), sep = "\t")
