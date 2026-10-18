#!/usr/bin/env bash
# 45_iadmix_stability.sh — is the iAdmix answer stable, given these three BAMs?
#
# The cross-check in bin/44 reports one iAdmix run per pool. This script asks
# how much that answer depends on choices made when setting the run up, which
# is what determines whether another analyst starting from the same three
# alignment files would reach the same number.
#
# Four things are varied, one at a time, on the same data:
#   1. Repeat runs of an identical input        -> is the program deterministic?
#   2. The pool size passed with -p             -> 50 / 100 / 200 / 400
#   3. Which SNPs are used                      -> five disjoint random subsets
#   4. Allele orientation                       -> already known to matter; quantified here
#
# Output: out_global/revision/iadmix_stability.tsv
set -euo pipefail
cd "$(dirname "$0")/.."

IADMIX_SRC="data/iadmix"
WORK="${TMPDIR:-/tmp}/iadmix_stab.$$"
OUT="out_global/revision/iadmix_stability.tsv"
N_SITES=10000
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

gcc -w -o "$WORK/lbfgsb.o" -c "$IADMIX_SRC/Lbfgsb.3.0/lbfgsb.c"
gcc -w "$WORK/lbfgsb.o" -lm -o "$WORK/A" "$IADMIX_SRC/ancestry-gll.c"

Rscript - "$WORK" "$N_SITES" <<'RS'
suppressPackageStartupMessages(library(data.table)); source("bin/lib_poolmix.R")
a <- commandArgs(trailingOnly = TRUE); work <- a[1]; nsit <- as.integer(a[2])
set.seed(20260925)
tab <- fread("out_global/revision/analysis_table.tsv.gz")[mask_main == TRUE]
subsets <- split(sample(nrow(tab)), rep(1:5, length.out = nrow(tab)))
subsets <- lapply(subsets, function(i) sort(i[seq_len(min(nsit, length(i)))]))

write_in <- function(p, idx, nchr, ori, path) {
    t <- tab[idx]
    x <- t[[paste0("x_", p)]]; d <- t[[paste0("d_", p)]]
    e <- estimate_error(t[[paste0("o_", p)]], d + t[[paste0("o_", p)]])
    k <- 0:nchr; mu <- (k / nchr) * (1 - e) + (1 - k / nchr) * e / 3
    LL <- outer(seq_along(x), seq_along(k), function(i, j) dbinom(x[i], d[i], mu[j], log = TRUE) / log(10))
    LL <- LL - apply(LL, 1, max)
    if (ori == "rev") LL <- LL[, rev(seq_len(ncol(LL))), drop = FALSE]
    writeLines("ID AFR EUR NAT", path)
    fwrite(as.data.table(cbind(round(LL, 4), round(as.matrix(t[, .(f_AFR, f_EUR, f_NAT)]), 6))),
           path, sep = " ", col.names = FALSE, append = TRUE)
}
# 1+4: subset 1, both orientations, nominal chromosome count
for (p in c("POOL1", "HOSPITAL")) {
    nchr <- if (p == "HOSPITAL") 400 else 100
    for (ori in c("fwd", "rev")) write_in(p, subsets[[1]], nchr, ori, file.path(work, sprintf("%s_ori%s.in", p, ori)))
    # 2: pool size
    for (n in c(50, 100, 200, 400)) write_in(p, subsets[[1]], n, "rev", file.path(work, sprintf("%s_p%d.in", p, n)))
    # 3: SNP subsets
    for (s in 1:5) write_in(p, subsets[[s]], nchr, "rev", file.path(work, sprintf("%s_sub%d.in", p, s)))
}
cat(sprintf("inputs written (%d SNPs each)\n", nsit))
RS

run() {  # run <file> <poolsize> ; prints "AFR EUR NAT logL"
    "$WORK/A" --input "$1" -p "$2" 2>/dev/null | grep "final maxval" | tail -1 \
      | sed 's/.*maxval \([^ ]*\) ADMIX_PROP AFR:\([^ ]*\) EUR:\([^ ]*\) NAT:\([^ ]*\).*/\2 \3 \4 \1/'
}

{
echo -e "pool\ttest\tsetting\tAFR\tEUR\tNAT\tlogL"
for p in POOL1 HOSPITAL; do
    nchr=100; [ "$p" = HOSPITAL ] && nchr=400
    for r in 1 2 3; do
        echo -e "$p\trepeat\trun$r\t$(run "$WORK/${p}_orirev.in" $nchr | tr ' ' '\t')"
    done
    for ori in fwd rev; do
        echo -e "$p\torientation\t$ori\t$(run "$WORK/${p}_ori${ori}.in" $nchr | tr ' ' '\t')"
    done
    for n in 50 100 200 400; do
        echo -e "$p\tpoolsize\t-p $n\t$(run "$WORK/${p}_p${n}.in" $n | tr ' ' '\t')"
    done
    for s in 1 2 3 4 5; do
        echo -e "$p\tsnp_subset\tsubset$s\t$(run "$WORK/${p}_sub${s}.in" $nchr | tr ' ' '\t')"
    done
done
} > "$OUT"

Rscript - "$OUT" <<'RS'
suppressPackageStartupMessages(library(data.table))
d <- fread(commandArgs(trailingOnly = TRUE)[1])
d[, (c("AFR","EUR","NAT")) := lapply(.SD, function(v) round(100 * v, 2)), .SDcols = c("AFR","EUR","NAT")]
print(d[, .(pool, test, setting, AFR, EUR, NAT, logL = round(logL))], row.names = FALSE)
cat("\nSpread within each test (percentage points, max - min):\n")
print(d[, .(AFR = round(diff(range(AFR)), 2), EUR = round(diff(range(EUR)), 2),
            NAT = round(diff(range(NAT)), 2)), by = .(pool, test)], row.names = FALSE)
RS
echo "Wrote $OUT"
