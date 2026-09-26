#!/usr/bin/env Rscript
# 34_fst_poolfstat.R — Pool-seq-aware differentiation between the three pools
#
# Memo items A3 (Pool-seq estimator with real pool sizes, block jackknife),
# A6 (common mask), S8 (equal-SNP blocks).
#
# Estimator: poolfstat's ANOVA FST (Hivert et al. 2018), which corrects for
# both the finite number of chromosomes (haploid pool size 2N) and the finite
# number of reads. The genome-wide value is poolfstat's ratio of averages.
# Uncertainty: poolfstat's block jackknife with NSNP_BLOCK SNPs per block.
#
# Usage: Rscript bin/34_fst_poolfstat.R [--table ...] [--mask mask_main] [--tag main]
# Output: out_global/revision/fst_pairwise_<tag>.tsv, fst_global_<tag>.tsv
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(data.table); library(poolfstat) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
OUT_DIR    <- "out_global/revision"
TABLE      <- get_arg("--table", file.path(OUT_DIR, "analysis_table.tsv.gz"))
MASK       <- get_arg("--mask", "mask_main")
TAG        <- get_arg("--tag", "main")
NSNP_BLOCK <- as.integer(get_arg("--block", "1000"))

POOLS  <- c("POOL1", "POOL2", "HOSPITAL")
POOL_N <- c(POOL1 = 50, POOL2 = 50, HOSPITAL = 200)

tab <- fread(TABLE)[get(MASK) == TRUE]

# poolfstat input in BayPass format: for each pool, reference-allele then
# alternate-allele read counts. Here "reference" = ancestral, "alt" = derived.
tmp <- tempfile("pf_")
geno <- do.call(cbind, lapply(POOLS, function(p)
    cbind(tab[[paste0("d_", p)]] - tab[[paste0("x_", p)]], tab[[paste0("x_", p)]])))
fwrite(as.data.table(geno), paste0(tmp, ".genobaypass"), sep = " ", col.names = FALSE)
writeLines(paste(2 * POOL_N[POOLS], collapse = " "), paste0(tmp, ".poolsize"))

pd <- genobaypass2pooldata(genobaypass.file = paste0(tmp, ".genobaypass"),
                           snp.pos = as.matrix(tab[, .(chr, pos)]),
                           poolsize.file = paste0(tmp, ".poolsize"),
                           poolnames = POOLS,
                           min.cov.per.pool = -1, max.cov.per.pool = 1e9,   # filters already applied by the mask
                           min.maf = -1, verbose = FALSE)
cat(sprintf("[%s] pooldata: %d SNPs, %d pools\n", TAG, pd@nsnp, pd@npools))

glob <- computeFST(pd, method = "Anova", nsnp.per.bjack.block = NSNP_BLOCK, verbose = FALSE)
pw   <- compute.pairwiseFST(pd, method = "Anova", nsnp.per.bjack.block = NSNP_BLOCK, verbose = FALSE)

res <- as.data.table(pw@values, keep.rownames = "pair")
setnames(res, gsub("[^A-Za-z0-9]+", "_", names(res)))
res[, `:=`(n_snps = pd@nsnp, block_snps = NSNP_BLOCK, tag = TAG)]
fwrite(res, file.path(OUT_DIR, paste0("fst_pairwise_", TAG, ".tsv")), sep = "\t")

g <- data.table(fst = glob$Fst[["Estimate"]], jk_mean = glob$Fst[["bjack mean"]], jk_se = glob$Fst[["bjack s.e."]],
                lo95 = glob$Fst[["CI95inf"]], hi95 = glob$Fst[["CI95sup"]], n_snps = pd@nsnp, tag = TAG)
fwrite(g, file.path(OUT_DIR, paste0("fst_global_", TAG, ".tsv")), sep = "\t")
print(res); print(g)
unlink(paste0(tmp, c(".genobaypass", ".poolsize")))
