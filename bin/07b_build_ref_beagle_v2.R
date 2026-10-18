#!/usr/bin/env Rscript
# 07b_build_ref_beagle_v2.R — granular version: converts GT genotypes from
# 110 reference 1KGP individuals (11 specific populations: IBS, CEU, YRI,
# ESN, GWD, LWK, PEL, MXL, PUR, FIN, CLM) to pseudo-Beagle format.
#
# METHODOLOGICAL WARNING: the resulting Beagle feeds NGSadmix/PCAngsd, which
# represent each pool (POOL1/POOL2/HOSPITAL: 50/50/200 real individuals) as
# ONE diploid pseudo-individual -- a structural limitation of the Beagle
# format (3 states per site: AA/Aa/aa), not configurable. See
# bin/README_VIGENTE.md and bin/16_pool_replicates.sh.
suppressPackageStartupMessages(library(data.table))

# This script is always invoked from PROJ_DIR (the calling .sh scripts do
# cd "$PROJ_DIR" first); otherwise PROJ_DIR can be set explicitly.
PROJ_DIR     <- Sys.getenv("PROJ_DIR", unset = getwd())
WORK_DIR     <- Sys.getenv("WORK_DIR", unset = file.path(PROJ_DIR, "work", "ref_beagle_v2"))
# Chromosomes to use: representative subset by default (1,6,10,13,17,20,22);
# override with the CHR_LIST env var (e.g. "1,6,10,13,17,20,22,X")
CHR_LIST     <- strsplit(Sys.getenv("CHR_LIST", unset = "1,6,10,13,17,20,22"), ",")[[1]]

SITES_FILE   <- file.path(WORK_DIR, "angsd_sites_hapmap3.txt")
GT_FILES     <- file.path(WORK_DIR, sprintf("ref_extract_v2/chr%s_gt.txt", CHR_LIST))
SAMPLES_FILE <- file.path(WORK_DIR, "ref_sample_order_v2.txt")
OUT_FILE     <- Sys.getenv("OUT_FILE", unset = "out_global/angsd/ngsadmix/ref110.beagle")

nuc_code <- c(A=0L, C=1L, G=2L, T=3L)

samples <- readLines(SAMPLES_FILE)
n_ind <- length(samples)
cat(sprintf("Individuos de referencia: %d\n", n_ind))

sites <- fread(SITES_FILE, header = FALSE, col.names = c("chr","pos","A1","A2"))
sites[, chr := as.character(chr)]
setkey(sites, chr, pos)

p_hom1 <- c(0.99, 0.005, 0.005)
p_het  <- c(0.005, 0.99, 0.005)
p_hom2 <- c(0.005, 0.005, 0.99)
p_miss <- c(1/3, 1/3, 1/3)

gt_to_probs <- function(gt, swap) {
    a <- substr(gt, 1, 1)
    b <- substr(gt, 3, 3)
    n <- length(gt)
    out <- matrix(p_miss, nrow = n, ncol = 3, byrow = TRUE)
    is00 <- a == "0" & b == "0"
    is11 <- a == "1" & b == "1"
    ishet <- (a == "0" & b == "1") | (a == "1" & b == "0")

    sel <- is00 & !swap; out[sel, ] <- matrix(p_hom1, sum(sel), 3, byrow = TRUE)
    sel <- is00 &  swap; out[sel, ] <- matrix(p_hom2, sum(sel), 3, byrow = TRUE)
    sel <- is11 & !swap; out[sel, ] <- matrix(p_hom2, sum(sel), 3, byrow = TRUE)
    sel <- is11 &  swap; out[sel, ] <- matrix(p_hom1, sum(sel), 3, byrow = TRUE)
    out[ishet, ] <- matrix(p_het, sum(ishet), 3, byrow = TRUE)
    out
}

all_rows <- list()
total_kept <- 0L
total_seen <- 0L

for (f in GT_FILES) {
    cat(sprintf("Procesando %s ...\n", f))
    dt <- fread(f, header = FALSE)
    setnames(dt, 1:4, c("chr","pos","REF","ALT"))
    dt[, chr := as.character(chr)]
    gt_cols <- 5:ncol(dt)
    total_seen <- total_seen + nrow(dt)

    m <- merge(dt, sites, by = c("chr","pos"), all = FALSE)
    if (nrow(m) == 0) next

    match_fwd <- m$REF == m$A1 & m$ALT == m$A2
    match_rev <- m$REF == m$A2 & m$ALT == m$A1
    keep <- match_fwd | match_rev
    m <- m[keep]
    swap_vec <- match_rev[keep]
    total_kept <- total_kept + nrow(m)
    if (nrow(m) == 0) next

    marker  <- paste0(m$chr, "_", m$pos)
    allele1 <- nuc_code[m$A1]
    allele2 <- nuc_code[m$A2]

    gt_mat <- as.matrix(m[, ..gt_cols])
    probs_list <- vector("list", n_ind)
    for (j in seq_len(n_ind)) {
        probs_list[[j]] <- gt_to_probs(gt_mat[, j], swap_vec)
    }
    probs_block <- do.call(cbind, probs_list)

    block <- data.table(marker = marker, allele1 = allele1, allele2 = allele2)
    block <- cbind(block, as.data.table(probs_block))
    all_rows[[f]] <- block
}

combined <- rbindlist(all_rows)
cat(sprintf("\nSitios vistos: %d | Sitios retenidos (alelos coinciden): %d (%.1f%%)\n",
            total_seen, total_kept, 100 * total_kept / total_seen))

ind_cols <- as.vector(t(outer(samples, c("p1","p2","p3"), paste0)))
setnames(combined, 4:ncol(combined), ind_cols)

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
fwrite(combined, OUT_FILE, sep = "\t", col.names = FALSE)
cat(sprintf("Escrito: %s (%d sitios, %d individuos)\n", OUT_FILE, nrow(combined), n_ind))
