#!/usr/bin/env Rscript
# 35_depth_matched.R — site-by-site depth matching of the three pools
#
# Memo item A10. Instead of one global downsampling fraction, every site is
# thinned to the smallest depth among the three pools at that site, so the
# pools end up with identical per-site depth (a stricter control, since the
# pools also differ in capture profile; see registro_discrepancias.md D5/D7).
#
# Thinning reads without replacement is exactly random read downsampling at a
# site: x' ~ Hypergeometric(x derived, d - x ancestral, m = target depth).
# The resulting table has the same columns as analysis_table.tsv.gz, so
# bin/33 and bin/34 run on it unchanged with --table and --tag depthmatched.
#
# Output: out_global/revision/analysis_table_depthmatched.tsv.gz
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "out_global/revision"
POOLS   <- c("POOL1", "POOL2", "HOSPITAL")
SEED    <- 20260925

set.seed(SEED)
tab <- fread(file.path(OUT_DIR, "analysis_table.tsv.gz"))
target <- do.call(pmin, lapply(POOLS, function(p) tab[[paste0("d_", p)]]))
for (p in POOLS) {
    x <- tab[[paste0("x_", p)]]; d <- tab[[paste0("d_", p)]]
    tab[, (paste0("x_", p)) := rhyper(.N, x, d - x, target)]
    tab[, (paste0("d_", p)) := target]
}
m <- tab$mask_main
cat(sprintf("Depth-matched: %d masked sites, mean depth %.1f, median %d (seed %d)\n",
            sum(m), mean(target[m]), as.integer(median(target[m])), SEED))
fwrite(tab, file.path(OUT_DIR, "analysis_table_depthmatched.tsv.gz"), sep = "\t")
