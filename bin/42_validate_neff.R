#!/usr/bin/env Rscript
# 42_validate_neff.R — does the effective-pool-size estimator recover the truth?
#
# The effective size reported for a pooled library (Results) is estimated from
# the residual variance of the observed read fractions around the fitted
# ancestry model. This script asks whether that estimate is right, by
# simulating pools whose effective size is known by construction and checking
# what the estimator returns.
#
# Construction. N people contribute DNA in unequal amounts: weights w_i are
# log-normal with coefficient of variation CV. A pool of N people whose
# contributions have CV carries the allele-frequency sampling variance of
#     N_eff = N / (1 + CV^2)
# equally-contributing people (the standard effective-sample-size result for a
# weighted mean: N_eff = (sum w)^2 / sum w^2, whose expectation is the above).
# CV = 0 reproduces equimolar pooling, where N_eff = N.
#
# Sweep: pool size x depth x CV, with replicates. For each simulated pool the
# script records the true N_eff, the estimate from the residual variance, and
# the estimate after calibrating out the panel-mismatch term the way the main
# analysis does (using the equimolar pools of the same run as the reference).
#
# Usage: Rscript bin/42_validate_neff.R   (SWEEP_REPS=20 by default)
# Output: out_global/revision/neff_validation.tsv  and  neff_validation_summary.tsv
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))
source("bin/lib_poolmix.R")

OUT_DIR  <- "out_global/revision"
N_SITES  <- as.integer(Sys.getenv("SWEEP_SITES", "5000"))
N_REPS   <- as.integer(Sys.getenv("SWEEP_REPS", "20"))
SEED     <- 20260925
CORES    <- as.integer(Sys.getenv("SWEEP_CORES", "8"))

POOL_N   <- c(50, 200)                 # people per pool
DEPTH_X  <- c(1, 1/3)                  # multiplier on the observed depth distribution
CV       <- c(0, 0.3, 0.6, 1.0, 2.0, 3.7)  # unevenness of DNA contribution;
                                       # CV 3.7 puts a 200-person pool at the
                                       # effective size estimated for HOSPITAL
Q_TRUE   <- c(0.15, 0.55, 0.30)        # AFR / EUR / NAT, close to the real cohorts
F_MIS    <- 0.01                       # proxy-vs-source mismatch
IND_SD   <- 0.10                       # SD of individual ancestry within a cohort

set.seed(SEED)
tab <- fread(file.path(OUT_DIR, "analysis_table.tsv.gz"))[mask_main == TRUE]
tab <- tab[sort(sample(.N, N_SITES))]
Fm  <- as.matrix(tab[, .(f_AFR, f_EUR, f_NAT)])
nk  <- c(tab$n_AFR[1], tab$n_EUR[1], tab$n_NAT[1])
depth_obs <- tab$d_POOL1
err <- estimate_error(tab$o_POOL1, tab$d_POOL1 + tab$o_POOL1)

bn <- function(f, F) { f <- pmin(pmax(f, 1e-4), 1 - 1e-4); rbeta(length(f), f * (1 - F) / F, (1 - f) * (1 - F) / F) }
rdirichlet_sd <- function(n, q, sd) {
    a0 <- max(mean(q * (1 - q)) / sd^2 - 1, 1)
    g <- matrix(rgamma(n * length(q), shape = rep(q * a0, each = n)), n)
    g / rowSums(g)
}

# One simulated pool. Returns the reads, the realised N_eff and the truth.
sim_pool <- function(N, cv, depth, P) {
    qi <- rdirichlet_sd(N, Q_TRUE, IND_SD)
    w  <- if (cv == 0) rep(1, N) else rlnorm(N, -log(1 + cv^2) / 2, sqrt(log(1 + cv^2)))
    n_eff_real <- sum(w)^2 / sum(w^2)              # realised, not just its expectation
    q_w <- colSums(qi * w) / sum(w)
    h   <- drop(P %*% q_w)
    # allele sampling with the effective number of chromosomes the weights imply
    pi  <- rbinom(length(h), max(2, round(2 * n_eff_real)), h) / max(2, round(2 * n_eff_real))
    x   <- rbinom(length(depth), depth, pi * (1 - err) + (1 - pi) * err / 3)
    list(x = x, d = depth, n_eff_true = n_eff_real, N = N, cv = cv)
}

# Raw excess variance, as a ratio of sums (see effective_sizes in lib_poolmix.R)
excess <- function(x, d, h) { v <- h * (1 - h); sum((x / d - h)^2 - v / d) / sum(v) }

run_rep <- function(rep) {
    P <- cbind(bn(Fm[, 1], F_MIS), bn(Fm[, 2], F_MIS), bn(Fm[, 3], F_MIS))
    grid <- CJ(N = POOL_N, dx = DEPTH_X, cv = CV)
    out <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
        g <- grid[i]
        depth <- pmax(5L, as.integer(round(depth_obs * g$dx)))
        s <- sim_pool(g$N, g$cv, depth, P)
        fit <- fit_pool(s$x, s$d, Fm, nk, 2 * g$N, err)
        h <- drop(Fm %*% fit[1:3])
        ex <- excess(s$x, s$d, h)
        out[[i]] <- data.table(rep = rep, N = g$N, depth_mult = g$dx, cv = g$cv,
                               n_eff_true = s$n_eff_true,
                               excess = ex,
                               n_eff_raw = 1 / (2 * max(ex, 1e-6)),
                               mean_depth = mean(depth),
                               q_AFR = fit[["q_AFR"]], q_EUR = fit[["q_EUR"]], q_NAT = fit[["q_NAT"]],
                               F_fit = fit[["F"]])
    }
    rbindlist(out)
}

t0 <- Sys.time()
res <- rbindlist(parallel::mclapply(seq_len(N_REPS), function(r) { set.seed(SEED + r); run_rep(r) }, mc.cores = CORES))
cat(sprintf("%d scenarios x %d reps in %.1f min\n", nrow(res) / N_REPS, N_REPS,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Calibration, as in the main analysis: the panel-mismatch term is shared by
# every pool of a replicate and is read off the equimolar (CV = 0) pools.
res[, shared := {
    m <- mean(excess[cv == 0] - 1 / (2 * n_eff_true[cv == 0]))
    m
}, by = .(rep, depth_mult)]
res[, n_eff_cal := 1 / (2 * pmax(excess - shared, 1e-6))]

fwrite(res, file.path(OUT_DIR, "neff_validation.tsv"), sep = "\t")

summ <- res[, .(reps = .N,
                n_eff_true = mean(n_eff_true),
                n_eff_raw = mean(n_eff_raw),
                n_eff_cal = mean(n_eff_cal),
                rel_bias_cal = mean(n_eff_cal / n_eff_true) - 1,
                sd_cal = sd(n_eff_cal),
                q_EUR_bias = mean(q_EUR) - Q_TRUE[2]),
            by = .(N, depth_mult, cv)][order(N, -depth_mult, cv)]
fwrite(summ, file.path(OUT_DIR, "neff_validation_summary.tsv"), sep = "\t")
print(summ, digits = 3)
