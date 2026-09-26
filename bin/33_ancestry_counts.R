#!/usr/bin/env Rscript
# 33_ancestry_counts.R — primary three-way ancestry estimate and cohort contrasts
#
# Memo items A2 (count-based model with pool size), A9 (Delta_HUN and
# POOL1 vs POOL2 with intervals), S8 (equal-SNP blocks).
# Model and jackknife: see bin/lib_poolmix.R.
#
# Usage: Rscript bin/33_ancestry_counts.R [--table analysis_table.tsv.gz]
#          [--mask mask_main] [--tag main] [--block 1000] [--ref AFR,EUR,NAT] [--refbias]
# Other scripts (35, 36) call the same code with a different table or mask.
#
# Output: out_global/revision/ancestry_<tag>.tsv   q per pool with SE and 95% CI
#         out_global/revision/contrasts_<tag>.tsv  Delta_HUN and POOL1-POOL2
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))
source("bin/lib_poolmix.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
OUT_DIR <- "out_global/revision"
TABLE   <- get_arg("--table", file.path(OUT_DIR, "analysis_table.tsv.gz"))
MASK    <- get_arg("--mask", "mask_main")
TAG     <- get_arg("--tag", "main")
BLOCK   <- as.integer(get_arg("--block", "1000"))
REFS    <- strsplit(get_arg("--ref", "AFR,EUR,NAT"), ",")[[1]]    # proxy groups from bin/30, in AFR,EUR,NAT order
REFBIAS <- "--refbias" %in% args                                    # also fit a reference-bias term per pool

POOLS    <- c("POOL1", "POOL2", "HOSPITAL")
POOL_N   <- c(POOL1 = 50, POOL2 = 50, HOSPITAL = 200)     # people per pool (author, 25 Sep 2026)

tab <- fread(TABLE)
tab <- tab[get(MASK) == TRUE]
Fm  <- as.matrix(tab[, paste0("f_", REFS), with = FALSE])
nk  <- sapply(REFS, function(r) tab[[paste0("n_", r)]][1])
blocks <- make_blocks(tab$chr, BLOCK)
cat(sprintf("[%s] %d sites, %d jackknife blocks of ~%d SNPs\n", TAG, nrow(tab), length(unique(blocks)), BLOCK))

refder <- if (REFBIAS) tab$ref == tab$der else NULL
PARAMS <- c("q_AFR", "q_EUR", "q_NAT", "F", "ref_bias")

err <- sapply(POOLS, function(p)
    estimate_error(tab[[paste0("o_", p)]], tab[[paste0("d_", p)]] + tab[[paste0("o_", p)]]))

# One function fits all three pools on a subset, so each jackknife replicate
# drops the same block from every pool.
warm <- NULL                                   # full-data q per pool, set after the full fit
fit_all <- function(keep) {
    unlist(lapply(POOLS, function(p) {
        r <- fit_pool(tab[[paste0("x_", p)]][keep], tab[[paste0("d_", p)]][keep],
                      Fm[keep, , drop = FALSE], nk, 2 * POOL_N[[p]], err[[p]],
                      if (REFBIAS) refder[keep] else NULL,
                      starts = if (is.null(warm)) START_Q else list(warm[[p]]))
        setNames(r[PARAMS], paste0(p, ".", PARAMS))
    }))
}

t0 <- Sys.time()
full <- fit_all(rep(TRUE, nrow(tab)))
cat(sprintf("Full fit: %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
warm <- lapply(setNames(POOLS, POOLS), function(p) full[paste0(p, ".q_", SOURCES)])
jk <- jackknife(blocks, fit_all, full)

# ── Ancestry table ────────────────────────────────────────────────────────────
anc <- rbindlist(lapply(POOLS, function(p) rbindlist(lapply(PARAMS, function(v) {
    k <- paste0(p, ".", v)
    data.table(pool = p, parameter = v, estimate = jk$est[[k]], se = jk$se[[k]],
               lo95 = jk$est[[k]] - 1.96 * jk$se[[k]], hi95 = jk$est[[k]] + 1.96 * jk$se[[k]])
}))))
anc[, `:=`(n_sites = nrow(tab), n_blocks = length(unique(blocks)), error_rate = err[pool],
           refs = paste(REFS, collapse = ","), tag = TAG)]

# ── Effective pool sizes (see effective_sizes in lib_poolmix.R) ─────────────
hs <- lapply(POOLS, function(p) drop(Fm %*% jk$est[paste0(p, ".q_", SOURCES)]))
es <- effective_sizes(lapply(POOLS, function(p) tab[[paste0("x_", p)]]),
                      lapply(POOLS, function(p) tab[[paste0("d_", p)]]), hs,
                      POOL_N[POOLS], reference_pools = c(1, 2))
N_EFF <- setNames(es$n_eff, POOLS)
anc[, n_eff := N_EFF[pool]]
cat(sprintf("Effective pool size (people): %s\n", paste(sprintf("%s %.0f", POOLS, N_EFF), collapse = ", ")))

# ── Contrasts, computed inside every replicate so they keep the covariance ───
g <- nrow(jk$reps)
con <- rbindlist(lapply(paste0("q_", SOURCES), function(q) {
    d_hun <- function(v) v[[paste0("HOSPITAL.", q)]] - (v[[paste0("POOL1.", q)]] + v[[paste0("POOL2.", q)]]) / 2
    d_12  <- function(v) v[[paste0("POOL1.", q)]] - v[[paste0("POOL2.", q)]]
    rbindlist(lapply(list(Delta_HUN = d_hun, POOL1_minus_POOL2 = d_12), function(f) {
        est  <- f(jk$est)
        reps <- apply(jk$reps, 1, function(r) f(setNames(r, names(jk$est))))
        se   <- sqrt((g - 1) / g * sum((reps - mean(reps))^2))
        data.table(component = q, estimate = est, se = se, lo95 = est - 1.96 * se, hi95 = est + 1.96 * se)
    }), idcol = "contrast")
}))
# Extra variance from which people (and how much of each person's DNA) the
# pool represents: SD of individual ancestry sigma over the effective number
# of people. "_neff" columns use N_EFF; "_nominal" columns use 50/50/200.
for (sig in c(0.10, 0.15)) for (sz in c("nominal", "neff")) {
    n <- if (sz == "neff") N_EFF else POOL_N[POOLS]
    v_hun <- sig^2 / n[["HOSPITAL"]] + (sig^2 / n[["POOL1"]] + sig^2 / n[["POOL2"]]) / 4
    v_12  <- sig^2 / n[["POOL1"]] + sig^2 / n[["POOL2"]]
    extra <- fifelse(con$contrast == "Delta_HUN", v_hun, v_12)
    se_t  <- sqrt(con$se^2 + extra)
    con[, (sprintf("lo95_sd%.2f_%s", sig, sz)) := estimate - 1.96 * se_t]
    con[, (sprintf("hi95_sd%.2f_%s", sig, sz)) := estimate + 1.96 * se_t]
}
con[, `:=`(n_sites = nrow(tab), tag = TAG)]

fwrite(anc, file.path(OUT_DIR, paste0("ancestry_", TAG, ".tsv")), sep = "\t")
fwrite(con, file.path(OUT_DIR, paste0("contrasts_", TAG, ".tsv")), sep = "\t")
fwrite(as.data.table(jk$reps), file.path(OUT_DIR, paste0("jackknife_reps_", TAG, ".tsv.gz")), sep = "\t")
print(anc[, .(pool, parameter, estimate = round(estimate, 4), se = round(se, 4))])
print(con[, .(contrast, component, estimate = round(estimate, 4), lo95 = round(lo95, 4), hi95 = round(hi95, 4))])
cat(sprintf("Total time: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
