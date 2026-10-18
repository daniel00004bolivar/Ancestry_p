#!/usr/bin/env Rscript
# 43_neff_uncertainty.R — how precise is the effective size reported for HOSPITAL?
#
# The block jackknife gives the effective-size estimate a standard error over
# which SNPs entered the analysis. That is not the quantity a reader needs:
# the estimator also varies over which people are in the pool, how unevenly
# they contributed, and the allele sampling that follows. This script measures
# that sampling distribution directly, by simulating the study exactly as it
# was run — the same three pools, the same per-site depths, the same number of
# SNPs, the same calibration against the two equimolar pools — with HOSPITAL's
# contribution unevenness set so that its true effective size is the value
# estimated from the data.
#
# Output: out_global/revision/neff_uncertainty.tsv (one row per replicate) and
# a summary printed to stdout.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))
source("bin/lib_poolmix.R")

OUT_DIR <- "out_global/revision"
REPS    <- as.integer(Sys.getenv("UNC_REPS", "100"))
SEED    <- 20260925
CORES   <- as.integer(Sys.getenv("UNC_CORES", "8"))
POOLS   <- c("POOL1", "POOL2", "HOSPITAL")
POOL_N  <- c(POOL1 = 50, POOL2 = 50, HOSPITAL = 200)
F_MIS   <- 0.01
IND_SD  <- 0.10

tab <- fread(file.path(OUT_DIR, "analysis_table.tsv.gz"))[mask_main == TRUE]
Fm  <- as.matrix(tab[, .(f_AFR, f_EUR, f_NAT)])
nk  <- c(tab$n_AFR[1], tab$n_EUR[1], tab$n_NAT[1])
anc <- fread(file.path(OUT_DIR, "ancestry_main.tsv"))
q_obs <- lapply(setNames(POOLS, POOLS), function(p) anc[pool == p & parameter %like% "^q_", estimate])
N_EFF_OBS <- anc[pool == "HOSPITAL", n_eff][1]
err <- sapply(POOLS, function(p) estimate_error(tab[[paste0("o_", p)]], tab[[paste0("d_", p)]] + tab[[paste0("o_", p)]]))

# Unevenness that puts HOSPITAL at the observed effective size
CV_HOSP <- sqrt(max(POOL_N[["HOSPITAL"]] / N_EFF_OBS - 1, 0))
cat(sprintf("Observed effective size %.1f of %d  ->  simulated contribution CV %.2f\n",
            N_EFF_OBS, POOL_N[["HOSPITAL"]], CV_HOSP))

bn <- function(f, F) { f <- pmin(pmax(f, 1e-4), 1 - 1e-4); rbeta(length(f), f * (1 - F) / F, (1 - f) * (1 - F) / F) }
rdirichlet_sd <- function(n, q, sd) {
    a0 <- max(mean(q * (1 - q)) / sd^2 - 1, 1)
    g <- matrix(rgamma(n * length(q), shape = rep(q * a0, each = n)), n)
    g / rowSums(g)
}

one_rep <- function(rep) {
    set.seed(SEED + rep)
    P <- cbind(bn(Fm[, 1], F_MIS), bn(Fm[, 2], F_MIS), bn(Fm[, 3], F_MIS))
    xs <- list(); ds <- list(); hs <- list(); truth <- numeric()
    for (p in POOLS) {
        N  <- POOL_N[[p]]
        cv <- if (p == "HOSPITAL") CV_HOSP else 0
        qi <- rdirichlet_sd(N, q_obs[[p]], IND_SD)
        w  <- if (cv == 0) rep(1, N) else rlnorm(N, -log(1 + cv^2) / 2, sqrt(log(1 + cv^2)))
        ne <- sum(w)^2 / sum(w^2)
        h  <- drop(P %*% (colSums(qi * w) / sum(w)))
        nchr <- max(2, round(2 * ne))
        pi <- rbinom(length(h), nchr, h) / nchr
        d  <- tab[[paste0("d_", p)]]
        x  <- rbinom(length(d), d, pi * (1 - err[[p]]) + (1 - pi) * err[[p]] / 3)
        fit <- fit_pool(x, d, Fm, nk, 2 * N, err[[p]])
        xs[[p]] <- x; ds[[p]] <- d; hs[[p]] <- drop(Fm %*% fit[paste0("q_", SOURCES)])
        truth[p] <- ne
    }
    es <- effective_sizes(xs[POOLS], ds[POOLS], hs[POOLS], POOL_N[POOLS], reference_pools = c(1, 2))
    data.table(rep = rep, pool = POOLS, n_eff_true = truth[POOLS], n_eff_est = es$n_eff)
}

t0 <- Sys.time()
res <- rbindlist(parallel::mclapply(seq_len(REPS), one_rep, mc.cores = CORES))
fwrite(res, file.path(OUT_DIR, "neff_uncertainty.tsv"), sep = "\t")
cat(sprintf("%d replicates in %.1f min\n\n", REPS, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

h <- res[pool == "HOSPITAL"]
cat(sprintf("HOSPITAL, over %d replicates of the study as run:\n", nrow(h)))
cat(sprintf("  true effective size   mean %.1f  (range %.1f-%.1f)\n", mean(h$n_eff_true), min(h$n_eff_true), max(h$n_eff_true)))
cat(sprintf("  estimated             mean %.1f  SD %.1f\n", mean(h$n_eff_est), sd(h$n_eff_est)))
cat(sprintf("  relative bias         %+.1f%%\n", 100 * (mean(h$n_eff_est / h$n_eff_true) - 1)))
cat(sprintf("  2.5-97.5%% of estimates %.1f to %.1f\n", quantile(h$n_eff_est, 0.025), quantile(h$n_eff_est, 0.975)))
cat(sprintf("  estimates above 50     %.0f%%   above 100  %.0f%%\n",
            100 * mean(h$n_eff_est > 50), 100 * mean(h$n_eff_est > 100)))
