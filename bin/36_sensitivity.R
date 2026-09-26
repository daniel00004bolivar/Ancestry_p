#!/usr/bin/env Rscript
# 36_sensitivity.R — pre-specified sensitivity analyses
#
# Memo items S7 (transversions, depth and QC filters, reference bias),
# A7 (reference-panel choice), S8 (block size).
# Every row re-runs bin/33 (ancestry) on the same analysis table with one
# change from the main analysis; bin/34 (FST) is re-run for the mask changes.
#
# Reference-bias diagnostic: mean derived-read fraction at sites where the
# GRCh37 reference base is the ancestral allele vs. where it is the derived
# allele, compared with the same quantity predicted from the fitted model.
#
# Output: out_global/revision/sensitivity_summary.tsv, reference_bias.tsv
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "out_global/revision"
POOLS   <- c("POOL1", "POOL2", "HOSPITAL")
TABLE   <- file.path(OUT_DIR, "analysis_table.tsv.gz")

tab <- fread(TABLE)
src_mean <- tab[, (f_AFR + f_EUR + f_NAT) / 3]
maf <- pmin(src_mean, 1 - src_mean)
qc_ok <- function(min_d, max_other) {
    k <- rep(TRUE, nrow(tab))
    for (p in POOLS) {
        d <- tab[[paste0("d_", p)]]; o <- tab[[paste0("o_", p)]]; dep <- tab[[paste0("depth_", p)]]
        k <- k & d >= min_d & d <= tab[[paste0("maxd_", p)]] & o <= max_other * dep
    }
    k
}
tab[, `:=`(
    mask_transversions = mask_main & transversion,
    mask_depth10       = qc_ok(10, 0.05) & maf >= 0.05,
    mask_depth40       = qc_ok(40, 0.05) & maf >= 0.05,
    mask_other2pct     = qc_ok(20, 0.02) & maf >= 0.05,
    mask_maf10         = mask_main & maf >= 0.10
)]
sens_table <- file.path(OUT_DIR, "analysis_table_sensitivity.tsv.gz")
fwrite(tab, sens_table, sep = "\t")

runs <- list(
    list(tag = "transversions", args = c("--mask", "mask_transversions")),
    list(tag = "depth10",       args = c("--mask", "mask_depth10")),
    list(tag = "depth40",       args = c("--mask", "mask_depth40")),
    list(tag = "other2pct",     args = c("--mask", "mask_other2pct")),
    list(tag = "maf10",         args = c("--mask", "mask_maf10")),
    list(tag = "ref_NATmeso",   args = c("--ref", "AFR,EUR,NATmeso")),
    list(tag = "ref_NATsa",     args = c("--ref", "AFR,EUR,NATsa")),
    list(tag = "ref_EURsw",     args = c("--ref", "AFR,EURsw,NAT")),
    list(tag = "ref_EURs",      args = c("--ref", "AFR,EURs,NAT")),
    list(tag = "ref_AFRyri",    args = c("--ref", "AFRyri,EUR,NAT")),
    list(tag = "block500",      args = c("--block", "500")),
    list(tag = "block2000",     args = c("--block", "2000"))
)
for (r in runs) {
    cat("\n=== ", r$tag, " ===\n")
    st <- system2("Rscript", c("bin/33_ancestry_counts.R", "--table", sens_table, "--tag", r$tag, r$args))
    if (st != 0) stop("bin/33 failed for ", r$tag)
    # poolfstat needs every chromosome to hold at least one full block; the
    # transversion-only set (~8,000 SNPs) is too sparse for 1000-SNP blocks.
    if (r$args[1] == "--mask")
        system2("Rscript", c("bin/34_fst_poolfstat.R", "--table", sens_table, "--tag", r$tag, r$args,
                             if (r$tag == "transversions") c("--block", "200")))
}

# ── Summary: Delta_HUN and POOL1-POOL2 under every analysis ──────────────────
tags <- c("main", "depthmatched", vapply(runs, `[[`, "", "tag"))
summ <- rbindlist(lapply(tags, function(t) {
    f <- file.path(OUT_DIR, paste0("contrasts_", t, ".tsv"))
    if (file.exists(f)) fread(f) else NULL
}), fill = TRUE)
fwrite(summ, file.path(OUT_DIR, "sensitivity_summary.tsv"), sep = "\t")

# ── Reference-bias diagnostic ────────────────────────────────────────────────
m <- tab[mask_main == TRUE]
anc_main <- fread(file.path(OUT_DIR, "ancestry_main.tsv"))
rb <- rbindlist(lapply(POOLS, function(p) {
    q <- anc_main[pool == p & parameter %like% "^q_", estimate]
    h <- as.matrix(m[, .(f_AFR, f_EUR, f_NAT)]) %*% q
    obs <- m[[paste0("x_", p)]] / m[[paste0("d_", p)]]
    rbindlist(lapply(c(TRUE, FALSE), function(ref_is_anc) {
        k <- (m$ref == m$anc) == ref_is_anc
        data.table(pool = p, ref_allele = if (ref_is_anc) "ancestral" else "derived",
                   n_sites = sum(k), observed = mean(obs[k]), expected = mean(h[k]),
                   obs_minus_exp = mean(obs[k]) - mean(h[k]))
    }))
}))
fwrite(rb, file.path(OUT_DIR, "reference_bias.tsv"), sep = "\t")
print(rb, digits = 4)
