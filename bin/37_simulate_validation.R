#!/usr/bin/env Rscript
# 37_simulate_validation.R — simulation of the full pooled design
#
# Memo item A8 (and C5: a small non-sensitive synthetic dataset is written).
#
# What is matched to the real data: the masked SNP set and proxy frequencies,
# each pool's real per-site depth (the observed d of that pool at that site),
# its estimated per-base error, and its real size (50 / 50 / 200 people).
# What is simulated, per replicate and pool:
#   1. true source frequencies differ from the proxies (Balding-Nichols, F_MISMATCH)
#   2. N individuals with ancestry ~ Dirichlet around the pool's q (SD ~ 0.10)
#   3. unequal DNA contribution (log-normal weights), which shrinks the
#      effective number of people to N / (1 + CV^2). POOL1/POOL2 use CV 0.3;
#      HOSPITAL uses the CV that reproduces its estimated effective size
#      (n_eff in ancestry_main.tsv, about 14 of 200 people)
#   4. pool allele count ~ Binomial(2N_eff, weighted mean frequency)
#   5. reads ~ Binomial(d, freq with sequencing error)
# The estimand is the cohort's unweighted mean ancestry (the realised people).
#
# Scenarios: "observed" uses the main estimates of bin/33 as the truth for
# each pool; "null" gives all three pools the same ancestry (the POOL1/POOL2
# mean), so any Delta_HUN it finds is a false positive driven by depth, size
# or error differences alone.
# For speed each replicate uses a random quarter of the masked sites and
# 40 jackknife blocks; intervals are therefore wider than in the real analysis.
# Two intervals are scored: the jackknife interval alone ("lo","hi") and the
# interval that adds the between-person variance sigma^2 / N_eff used by
# bin/33 ("lo_neff","hi_neff", with sigma = IND_SD).
#
# Output: out_global/revision/simulation_replicates.tsv, simulation_summary.tsv,
#         out_global/revision/synthetic_test/ (C5)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))
source("bin/lib_poolmix.R")

OUT_DIR    <- "out_global/revision"
POOLS      <- c("POOL1", "POOL2", "HOSPITAL")
POOL_N     <- c(POOL1 = 50, POOL2 = 50, HOSPITAL = 200)
N_REPS     <- as.integer(Sys.getenv("SIM_REPS", "20"))
F_LEVELS   <- c(0.01, 0.03)
IND_SD     <- 0.10
DNA_CV_POOLS <- 0.30
SITE_FRAC  <- 0.25
N_BLOCKS   <- 40
SEED       <- 20260925

set.seed(SEED)
tab <- fread(file.path(OUT_DIR, "analysis_table.tsv.gz"))[mask_main == TRUE]
Fm_all <- as.matrix(tab[, .(f_AFR, f_EUR, f_NAT)])
nk <- c(tab$n_AFR[1], tab$n_EUR[1], tab$n_NAT[1])
err <- sapply(POOLS, function(p) estimate_error(tab[[paste0("o_", p)]], tab[[paste0("d_", p)]] + tab[[paste0("o_", p)]]))
anc_main <- fread(file.path(OUT_DIR, "ancestry_main.tsv"))
q_obs <- lapply(setNames(POOLS, POOLS), function(p) anc_main[pool == p & parameter %like% "^q_", estimate])
q_null <- (q_obs$POOL1 + q_obs$POOL2) / 2
n_eff <- sapply(POOLS, function(p) anc_main[pool == p, n_eff][1])
dna_cv <- sapply(POOLS, function(p) if (p == "HOSPITAL") sqrt(max(POOL_N[[p]] / n_eff[[p]] - 1, 0)) else DNA_CV_POOLS)
cat(sprintf("DNA contribution CV: %s\n", paste(sprintf("%s %.2f", POOLS, dna_cv), collapse = ", ")))

bn <- function(f, F) { f <- pmin(pmax(f, 1e-4), 1 - 1e-4); rbeta(length(f), f * (1 - F) / F, (1 - f) * (1 - F) / F) }
rdirichlet_sd <- function(n, q, sd) {
    a0 <- max(mean(q * (1 - q)) / sd^2 - 1, 1)           # concentration giving roughly the target SD
    g <- matrix(rgamma(n * length(q), shape = rep(q * a0, each = n)), n)
    g / rowSums(g)
}

simulate_pool <- function(q, N, P, d, e, cv) {
    qi <- rdirichlet_sd(N, q, IND_SD)
    w  <- rlnorm(N, -log(1 + cv^2) / 2, sqrt(log(1 + cv^2)))
    q_w <- colSums(qi * w) / sum(w)
    chrom_eff <- max(2, round(2 * N / (1 + cv^2)))
    h <- drop(P %*% q_w)
    c <- rbinom(length(h), chrom_eff, h)
    pi <- c / chrom_eff
    x <- rbinom(length(d), d, pi * (1 - e) + (1 - pi) * e / 3)
    list(x = x, truth = colMeans(qi))
}

run_rep <- function(scenario, Fmis, rep_id) {
    idx <- sort(sample(nrow(tab), round(SITE_FRAC * nrow(tab))))
    Fm <- Fm_all[idx, ]
    P <- cbind(bn(Fm[, 1], Fmis), bn(Fm[, 2], Fmis), bn(Fm[, 3], Fmis))
    sims <- lapply(POOLS, function(p) {
        q <- if (scenario == "null") q_null else q_obs[[p]]
        simulate_pool(q, POOL_N[[p]], P, tab[[paste0("d_", p)]][idx], err[[p]], dna_cv[[p]])
    }); names(sims) <- POOLS
    blocks <- make_blocks(tab$chr[idx], ceiling(length(idx) / N_BLOCKS))
    warm <- NULL
    fit_all <- function(keep) unlist(lapply(POOLS, function(p) {
        r <- fit_pool(sims[[p]]$x[keep], tab[[paste0("d_", p)]][idx][keep], Fm[keep, ], nk, 2 * POOL_N[[p]], err[[p]],
                      starts = if (is.null(warm)) START_Q else list(warm[[p]]))
        setNames(r[1:3], paste0(p, ".", names(r)[1:3]))
    }))
    full <- fit_all(rep(TRUE, length(idx)))
    warm <- lapply(setNames(POOLS, POOLS), function(p) full[paste0(p, ".q_", SOURCES)])
    jk <- jackknife(blocks, fit_all, full)
    g <- nrow(jk$reps)
    out <- list()
    for (k in 1:3) {
        comp <- SOURCES[k]
        for (p in POOLS) {
            key <- paste0(p, ".q_", comp)
            out[[length(out) + 1]] <- data.table(quantity = paste0("q_", comp), pool = p,
                truth = sims[[p]]$truth[k], estimate = jk$est[[key]], se = jk$se[[key]],
                extra_var = IND_SD^2 / n_eff[[p]])
        }
        dfun <- function(v) v[[paste0("HOSPITAL.q_", comp)]] - (v[[paste0("POOL1.q_", comp)]] + v[[paste0("POOL2.q_", comp)]]) / 2
        reps <- apply(jk$reps, 1, function(r) dfun(setNames(r, names(jk$est))))
        out[[length(out) + 1]] <- data.table(quantity = paste0("Delta_HUN_", comp), pool = "contrast",
            truth = sims$HOSPITAL$truth[k] - (sims$POOL1$truth[k] + sims$POOL2$truth[k]) / 2,
            estimate = dfun(jk$est), se = sqrt((g - 1) / g * sum((reps - mean(reps))^2)),
            extra_var = IND_SD^2 * (1 / n_eff[["HOSPITAL"]] + (1 / n_eff[["POOL1"]] + 1 / n_eff[["POOL2"]]) / 4))
    }
    res <- rbindlist(out)
    res[, `:=`(scenario = scenario, F_mismatch = Fmis, rep = rep_id)]
    if (scenario == "observed" && Fmis == F_LEVELS[1] && rep_id == 1) write_synthetic(idx, sims, res)
    res
}

# C5: a small synthetic dataset with its expected output (no real reads)
write_synthetic <- function(idx, sims, res) {
    dir.create(file.path(OUT_DIR, "synthetic_test"), showWarnings = FALSE)
    k <- idx[1:min(20000, length(idx))]; j <- seq_along(k)
    syn <- data.table(chr = tab$chr[k], pos = tab$pos[k], f_AFR = tab$f_AFR[k], f_EUR = tab$f_EUR[k], f_NAT = tab$f_NAT[k])
    for (p in POOLS) { syn[, (paste0("x_", p)) := sims[[p]]$x[j]]; syn[, (paste0("d_", p)) := tab[[paste0("d_", p)]][k]] }
    fwrite(syn, file.path(OUT_DIR, "synthetic_test", "synthetic_counts.tsv.gz"), sep = "\t")
    fwrite(res[pool != "contrast", .(pool, quantity, truth)], file.path(OUT_DIR, "synthetic_test", "synthetic_truth.tsv"), sep = "\t")
}

t0 <- Sys.time()
all <- list()
for (sc in c("observed", "null")) for (Fmis in F_LEVELS) for (r in seq_len(N_REPS)) {
    set.seed(SEED + 1000 * match(sc, c("observed", "null")) + 100 * match(Fmis, F_LEVELS) + r)
    all[[length(all) + 1]] <- run_rep(sc, Fmis, r)
    cat(sprintf("\r%s F=%.2f rep %d/%d  (%.1f min)", sc, Fmis, r, N_REPS, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("\n")
sim <- rbindlist(all)
sim[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se,
           lo_neff = estimate - 1.96 * sqrt(se^2 + extra_var), hi_neff = estimate + 1.96 * sqrt(se^2 + extra_var))]
fwrite(sim, file.path(OUT_DIR, "simulation_replicates.tsv"), sep = "\t")

summ <- sim[, .(reps = .N, mean_truth = mean(truth), bias = mean(estimate - truth),
                rmse = sqrt(mean((estimate - truth)^2)), mean_se = mean(se),
                coverage95_jackknife = mean(truth >= lo & truth <= hi),
                coverage95_neff = mean(truth >= lo_neff & truth <= hi_neff),
                reject_zero_jackknife = mean(lo > 0 | hi < 0),
                reject_zero_neff = mean(lo_neff > 0 | hi_neff < 0)),
            by = .(scenario, F_mismatch, quantity, pool)]
fwrite(summ, file.path(OUT_DIR, "simulation_summary.tsv"), sep = "\t")
print(summ, digits = 3)
