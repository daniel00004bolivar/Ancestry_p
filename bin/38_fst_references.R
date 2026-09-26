#!/usr/bin/env Rscript
# 38_fst_references.R — pool-vs-reference differentiation (supplementary)
#
# Replaces manuscript Tables 5-6 (memo A3, A4, A6, C1 items F1-F5, R1-R5).
# One common SNP set (mask_main), one orientation (derived allele), one
# estimator, one averaging rule:
#   Hudson FST, ratio of averages (Bhatia et al. 2013):
#     num_j = (p1 - p2)^2 - v1 - v2 ,  den_j = p1 (1 - p2) + p2 (1 - p1)
#     FST   = sum(num) / sum(den)
#   with the sampling variance v of each frequency estimate corrected for how
#   it was sampled:
#     pool (N people, d reads): v = p(1-p)/(1 - c) * c,  c = 1/(2N) + 1/d - 1/(2N d)
#     reference (n chromosomes): v = p(1-p)/(n - 1)
# Uncertainty: delete-one-block jackknife, blocks of 1000 SNPs (memo S8).
# Pool-vs-pool FST with this estimator is reported as a cross-check of the
# poolfstat ANOVA estimate (bin/34).
#
# Output: out_global/revision/fst_references.tsv
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))
source("bin/lib_poolmix.R")

OUT_DIR <- "out_global/revision"
POOLS   <- c("POOL1", "POOL2", "HOSPITAL")
POOL_N  <- c(POOL1 = 50, POOL2 = 50, HOSPITAL = 200)
BLOCK   <- 1000L
REF_POPS <- c("Yoruba", "Mandenka", "BantuKenya", "French", "Basque", "Italian", "Tuscan",
              "Sardinian", "Orcadian", "Russian", "Maya", "Pima", "Colombian", "Karitiana",
              "Surui", "Han", "Japanese")
ALL_POPS <- c("Adygei","Balochi","BantuKenya","BantuSouthAfrica","Basque","Bedouin",
              "BiakaPygmy","Brahui","Burusho","Cambodian","Colombian","Dai","Daur",
              "Druze","French","Han","Han-NChina","Hazara","Hezhen","Italian",
              "Japanese","Kalash","Karitiana","Lahu","Makrani","Mandenka","Maya",
              "MbutiPygmy","Melanesian","Miao","Mongola","Mozabite","Naxi",
              "Orcadian","Oroqen","Palestinian","Papuan","Pathan","Pima","Russian",
              "San","Sardinian","She","Sindhi","Surui","Tu","Tujia","Tuscan",
              "Uygur","Xibo","Yakut","Yi","Yoruba")

tab <- fread(file.path(OUT_DIR, "analysis_table.tsv.gz"))[mask_main == TRUE]
h <- fread(cmd = "gzcat data/hgdp/hgdpGeo.txt.gz", sep = "\t", header = FALSE,
           col.names = c("bin", "chrom", "start", "end", "rsid", "anc", "der", "freqs"))
h <- h[chrom %in% paste0("chr", 1:22)][, `:=`(chr = as.integer(sub("chr", "", chrom)), pos = end)]
h <- unique(h, by = c("chr", "pos"))[tab[, .(chr, pos)], on = c("chr", "pos"), nomatch = NULL]
fr <- tstrsplit(sub(",$", "", h$freqs), ",", fixed = TRUE); names(fr) <- ALL_POPS
ref <- data.table(chr = h$chr, pos = h$pos)
for (p in REF_POPS) {
    v <- as.numeric(fr[[p]])
    ref[, (p) := 1 - v]                                        # derived-allele frequency
    ref[, (paste0("n.", p)) := round(1 / min(pmin(v[v > 0 & v < 1], 1 - v[v > 0 & v < 1])))]
}
tab <- merge(tab, ref, by = c("chr", "pos"))
setorder(tab, chr, pos)
blocks <- make_blocks(tab$chr, BLOCK)
cat(sprintf("%d sites, %d blocks\n", nrow(tab), length(unique(blocks))))

freq_var <- function(kind, name) {
    if (kind == "pool") {
        p <- tab[[paste0("x_", name)]] / tab[[paste0("d_", name)]]
        d <- tab[[paste0("d_", name)]]; twoN <- 2 * POOL_N[[name]]
        c <- 1 / twoN + 1 / d - 1 / (twoN * d)
        list(p = p, v = p * (1 - p) / (1 - c) * c)
    } else {
        p <- tab[[name]]; n <- tab[[paste0("n.", name)]][1]
        list(p = p, v = p * (1 - p) / (n - 1))
    }
}
hudson <- function(a, b) {
    num <- (a$p - b$p)^2 - a$v - b$v
    den <- a$p * (1 - b$p) + b$p * (1 - a$p)
    ok <- is.finite(num) & is.finite(den)
    list(num = num[ok], den = den[ok], blk = blocks[ok])
}
jk_ratio <- function(s) {
    ub <- unique(s$blk); g <- length(ub)
    bn <- tapply(s$num, s$blk, sum); bd <- tapply(s$den, s$blk, sum)
    est <- sum(bn) / sum(bd)
    reps <- (sum(bn) - bn) / (sum(bd) - bd)
    se <- sqrt((g - 1) / g * sum((reps - mean(reps))^2))
    c(fst = est, se = se)
}

res <- list()
for (p in POOLS) {
    a <- freq_var("pool", p)
    for (r in REF_POPS) {
        v <- jk_ratio(hudson(a, freq_var("ref", r)))
        res[[length(res) + 1]] <- data.table(pool = p, reference = r, type = "HGDP population",
                                              fst = v[["fst"]], se = v[["se"]])
    }
}
for (pr in list(c("POOL1", "POOL2"), c("POOL1", "HOSPITAL"), c("POOL2", "HOSPITAL"))) {
    v <- jk_ratio(hudson(freq_var("pool", pr[1]), freq_var("pool", pr[2])))
    res[[length(res) + 1]] <- data.table(pool = pr[1], reference = pr[2], type = "pool (cross-check of bin/34)",
                                          fst = v[["fst"]], se = v[["se"]])
}
res <- rbindlist(res)
res[, `:=`(lo95 = fst - 1.96 * se, hi95 = fst + 1.96 * se, n_sites = nrow(tab))]
fwrite(res, file.path(OUT_DIR, "fst_references.tsv"), sep = "\t")
print(res[, .(pool, reference, fst = round(fst, 5), se = round(se, 5))])
