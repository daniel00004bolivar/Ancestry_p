#!/usr/bin/env bash
# 44_iadmix_crosscheck.sh — independent implementation check of the ancestry estimates
#
# Re-estimates the ancestry proportions of the three pools with iAdmix
# (Bansal and Libiger 2015), an independently written, peer-reviewed program
# that also accepts a pool size, and compares the result with the count-based
# model of bin/33. Same SNPs, same HGDP proxy panel, same chromosome counts;
# different likelihood (binomial, no overdispersion), different optimizer
# (L-BFGS-B), different codebase.
#
# What this does and does not check. It checks the implementation: that the
# proportions reported by bin/33 are not an artefact of how that code is
# written or optimised. It does NOT check the framework: both programs are
# supervised and share the proxy panel, so a wrong panel would mislead both.
# It also cannot speak to the effective-pool-size result, because iAdmix's
# binomial likelihood has no overdispersion term.
#
# ALLELE ORIENTATION — the trap, and how this script avoids it. iAdmix reads,
# per site, POOL_SIZE+1 genotype log-likelihoods followed by the panel
# frequencies, and pairs column k with one of the two alleles. Supplying the
# likelihoods for the other allele silently returns a plausible-looking but
# wrong answer (African-related 47% instead of 14%): no error, no warning.
# The same failure is visible in this project's earlier iAdmix runs.
#
# Rather than assume the convention, this script writes BOTH orientations and
# keeps the one iAdmix itself scores higher. The two differ by a factor of
# about 30 in log-likelihood, so the choice is not marginal, and it is made on
# iAdmix's own criterion rather than by picking whichever matches bin/33.
#
# POOL SIZE — the second trap. iAdmix's -p is the number of CHROMOSOMES, not
# of people: 100 for a 50-person pool, not 50.
#
# Output: out_global/revision/iadmix_crosscheck.tsv
set -euo pipefail
cd "$(dirname "$0")/.."

IADMIX_SRC="data/iadmix"
WORK="${TMPDIR:-/tmp}/iadmix_check.$$"
OUT="out_global/revision/iadmix_crosscheck.tsv"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

# ── Compile. The binary shipped with iAdmix is Linux x86-64; build from source.
echo "Compiling iAdmix..."
gcc -w -o "$WORK/lbfgsb.o" -c "$IADMIX_SRC/Lbfgsb.3.0/lbfgsb.c"
gcc -w "$WORK/lbfgsb.o" -lm -o "$WORK/ANCESTRY" "$IADMIX_SRC/ancestry-gll.c"

# ── Build one input file per pool, then run.
Rscript - "$WORK" <<'RS'
suppressPackageStartupMessages(library(data.table))
source("bin/lib_poolmix.R")
work <- commandArgs(trailingOnly = TRUE)[1]
tab <- fread("out_global/revision/analysis_table.tsv.gz")[mask_main == TRUE]
freqs <- as.matrix(tab[, .(f_AFR, f_EUR, f_NAT)])
CHROM <- c(POOL1 = 100, POOL2 = 100, HOSPITAL = 400)   # chromosomes, not people
for (p in names(CHROM)) {
    n <- CHROM[[p]]
    x <- tab[[paste0("x_", p)]]; d <- tab[[paste0("d_", p)]]
    e <- estimate_error(tab[[paste0("o_", p)]], d + tab[[paste0("o_", p)]])
    k  <- 0:n
    mu <- (k / n) * (1 - e) + (1 - k / n) * e / 3          # derived-allele read fraction
    LL <- outer(seq_along(x), seq_along(k),
                function(i, j) dbinom(x[i], d[i], mu[j], log = TRUE) / log(10))
    LL <- LL - apply(LL, 1, max)                            # per-site offset; the model is invariant to it
    f <- file.path(work, paste0(p, ".input"))
    for (ori in c("fwd", "rev")) {
        M <- if (ori == "fwd") LL else LL[, rev(seq_len(ncol(LL))), drop = FALSE]
        f <- file.path(work, paste0(p, ".", ori, ".input"))
        writeLines("ID AFR EUR NAT", f)                     # header names the panels only
        fwrite(as.data.table(cbind(round(M, 4), round(freqs, 6))), f,
               sep = " ", col.names = FALSE, append = TRUE)
    }
    cat(sprintf("  %s: %d SNPs, %d chromosomes\n", p, nrow(tab), n))
}
RS

echo "Running iAdmix in both allele orientations; keeping the better-scoring one..."
: > "$WORK/raw.txt"
for p in POOL1 POOL2 HOSPITAL; do
    n=100; [ "$p" = HOSPITAL ] && n=400
    best_ll=""; best_line=""; best_ori=""
    for ori in fwd rev; do
        out=$("$WORK/ANCESTRY" --input "$WORK/$p.$ori.input" -p "$n" 2>/dev/null | grep "final maxval" | tail -1)
        ll=$(echo "$out" | awk '{print $3}')
        props=$(echo "$out" | sed 's/.*ADMIX_PROP //')
        echo "    $p $ori: logL = $ll"
        if [ -z "$best_ll" ] || awk -v a="$ll" -v b="$best_ll" 'BEGIN{exit !(a>b)}'; then
            best_ll=$ll; best_line="$props"; best_ori=$ori
        fi
    done
    echo "  -> $p: keeping '$best_ori' (logL $best_ll)"
    echo "$p $best_line" >> "$WORK/raw.txt"
done

# ── Parse and compare with bin/33.
Rscript - "$WORK/raw.txt" "$OUT" <<'RS'
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
raw <- readLines(a[1])
res <- rbindlist(lapply(raw, function(l) {
    f <- strsplit(trimws(l), "\\s+")[[1]]
    g <- function(k) as.numeric(sub(paste0("^", k, ":"), "", grep(paste0("^", k, ":"), f, value = TRUE)[1]))
    data.table(pool = f[1], source = c("AFR", "EUR", "NAT"), iadmix = c(g("AFR"), g("EUR"), g("NAT")))
}))
main <- fread("out_global/revision/ancestry_main.tsv")[parameter %like% "^q_",
          .(pool, source = sub("q_", "", parameter), count_model = estimate)]
cmp <- merge(res, main, by = c("pool", "source"))
cmp[, diff_points := 100 * (iadmix - count_model)]
fwrite(cmp, a[2], sep = "\t")
print(cmp[, .(pool, source, iadmix = round(100 * iadmix, 2),
              count_model = round(100 * count_model, 2), diff_points = round(diff_points, 2))])
cat(sprintf("\nLargest absolute difference: %.2f percentage points\n", max(abs(cmp$diff_points))))
RS
echo "Wrote $OUT"
