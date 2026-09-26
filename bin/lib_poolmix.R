# lib_poolmix.R — count-based admixture model for pool-seq (memo A2, A9)
#
# Shared by 33_ancestry_counts.R, 35_depth_matched.R, 36_sensitivity.R and
# 37_simulate_validation.R. Not run directly.
#
# Model, per pool and site j (derived allele, orientation fixed by bin/30):
#   h_j   = sum_k q_k f_kj                      expected derived frequency
#   mu_j  = h_j (1 - e) + (1 - h_j) e / 3       with per-base error rate e
#   x_j | d_j ~ BetaBinomial(d_j, mean mu_j, intra-class correlation rho_j)
#   rho_j = 1/(2N) + F + sum_k q_k^2 f_kj (1 - f_kj) / n_k / (h_j (1 - h_j))
# The three terms are, in order: sampling of 2N chromosomes into the pool
# (N = real number of people), drift / mismatch between the proxy panel and
# the true source populations (F, estimated), and the sampling error of the
# proxy frequencies themselves (n_k chromosomes per proxy group). Read
# sampling is the binomial part of the beta-binomial.
#
# q is parametrised with a softmax and F with a logistic, so optim works on an
# unconstrained scale. Uncertainty: delete-one-block jackknife on contiguous
# blocks of equal SNP count (memo S8), with the SAME blocks for every pool so
# contrasts between pools (Delta) carry their covariance.

suppressPackageStartupMessages(library(data.table))

SOURCES <- c("AFR", "EUR", "NAT")

.softmax <- function(z) { e <- exp(c(0, z)); e / sum(e) }

# Negative composite log-likelihood; x, d vectors; Fm matrix sites x K; nk length K.
# Optional reference-bias term (refder = TRUE where the GRCh37 base is the
# derived allele): reads carrying the non-reference allele are sampled with
# relative efficiency r = exp(par[4]), so the expected non-reference read
# fraction is a r / (a r + 1 - a) for true non-reference frequency a.
.nll <- function(par, x, d, Fm, nk, twoN, err, refder = NULL) {
    q   <- .softmax(par[1:2])
    Fd  <- plogis(par[3]) * 0.2                     # drift term capped at 0.2
    h   <- drop(Fm %*% q)
    h   <- pmin(pmax(h, 1e-6), 1 - 1e-6)
    mu  <- h * (1 - err) + (1 - h) * err / 3
    if (!is.null(refder)) {
        r  <- exp(par[4])
        a  <- ifelse(refder, 1 - mu, mu)            # non-reference allele fraction
        a  <- a * r / (a * r + 1 - a)
        mu <- ifelse(refder, 1 - a, a)
    }
    vref <- drop((Fm * (1 - Fm)) %*% (q^2 / nk)) / (h * (1 - h))
    rho <- pmin(1 / twoN + Fd + vref, 0.999)
    a   <- mu * (1 - rho) / rho
    b   <- (1 - mu) * (1 - rho) / rho
    -sum(lbeta(x + a, d - x + b) - lbeta(a, b))     # binomial coefficient omitted (constant)
}

# Fit one pool. Returns named vector q_AFR, q_EUR, q_NAT, F, ref_bias, loglik.
# ref_bias is the non-reference/reference read efficiency (1 = no bias; only
# estimated when refder is supplied).
# Several starting points, keeping the best likelihood: in noisy pools a
# single start can stall in a boundary local optimum (one source ~100%, F at
# its cap), found in the simulations of bin/37.
START_Q <- list(c(1, 1, 1), c(0.15, 0.55, 0.30), c(0.30, 0.40, 0.30), c(0.10, 0.30, 0.60), c(0.60, 0.20, 0.20))
# `starts` is a list of starting q vectors; jackknife replicates pass the
# full-data estimate as their only start.
fit_pool <- function(x, d, Fm, nk, twoN, err, refder = NULL, starts = START_Q) {
    fits <- lapply(starts, function(q0) {
        start <- c(log(q0[2:3] / q0[1]), qlogis(0.05), if (!is.null(refder)) 0)
        optim(start, .nll, x = x, d = d, Fm = Fm, nk = nk, twoN = twoN, err = err, refder = refder,
              method = "BFGS", control = list(maxit = 500, reltol = 1e-10))
    })
    o <- fits[[which.min(vapply(fits, `[[`, 0, "value"))]]
    if (o$convergence != 0) warning("optim did not converge (code ", o$convergence, ")")
    q <- .softmax(o$par[1:2])
    c(setNames(q, paste0("q_", SOURCES)), F = plogis(o$par[3]) * 0.2,
      ref_bias = if (!is.null(refder)) exp(o$par[4]) else 1,
      loglik = -o$value, converged = as.numeric(o$convergence == 0))
}

# Effective number of people in each pool, from the residual variance of the
# observed read fractions around the fitted expectation. The part shared by
# all pools (proxy mismatch) is calibrated on the reference pools, assumed to
# behave as their nominal size: excess_p = 1/(2 N_eff,p) + m.
effective_sizes <- function(xs, ds, hs, nominal, reference_pools) {
    ex <- mapply(function(x, d, h) mean(((x / d - h)^2 - h * (1 - h) / d) / (h * (1 - h))), xs, ds, hs)
    m  <- mean(ex[reference_pools] - 1 / (2 * nominal[reference_pools]))
    neff <- 1 / (2 * pmax(ex - m, 1e-6))
    neff[reference_pools] <- nominal[reference_pools]
    list(excess = ex, shared = m, n_eff = pmin(neff, nominal))
}

# Per-base error rate from reads that are neither the ancestral nor the derived
# base: they are 2 of the 3 possible wrong bases, so e = 1.5 * other / depth.
estimate_error <- function(other, depth) 1.5 * sum(other) / sum(depth)

# Equal-SNP-count contiguous blocks, never spanning two chromosomes
make_blocks <- function(chr, snps_per_block) {
    blk <- integer(length(chr)); nxt <- 0L
    for (ch in unique(chr)) {
        idx <- which(chr == ch)
        b <- (seq_along(idx) - 1L) %/% snps_per_block
        nb <- max(b) + 1L
        # fold a short last block into the previous one
        if (nb > 1 && sum(b == nb - 1L) < snps_per_block / 2) b[b == nb - 1L] <- nb - 2L
        blk[idx] <- b + nxt
        nxt <- nxt + max(b) + 1L
    }
    blk
}

# Delete-one-block jackknife. `fit_fun(keep)` returns a named numeric vector
# for the sites in logical index `keep`. Returns list(est, se, reps).
jackknife <- function(blocks, fit_fun, full = NULL, cores = 8L) {
    if (is.null(full)) full <- fit_fun(rep(TRUE, length(blocks)))
    ub <- sort(unique(blocks))
    reps <- do.call(rbind, parallel::mclapply(ub, function(b) fit_fun(blocks != b), mc.cores = cores))
    g <- length(ub)
    se <- sqrt((g - 1) / g * colSums(sweep(reps, 2, colMeans(reps))^2))
    list(est = full, se = setNames(se, names(full)), reps = reps)
}
