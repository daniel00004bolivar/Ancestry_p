#!/usr/bin/env Rscript
# Builds a reference Beagle file (pseudo-certain genotype likelihoods from
# VCF genotypes) for a panel of 1000G individuals.
#
# Extracted from bin/12_rebuild_and_run.sh: this code used to be duplicated
# twice (panel v1 of 100 individuals, panel v2 of 110), identical except for
# paths/variables. Now it's a parameterized function invoked once per panel.
#
# Usage: Rscript lib_build_ref_beagle_panel.R <sites_file> <gt_dir> <samples_file> <out_file> <label>
#   sites_file    chr pos A1 A2 (order used by ANGSD, allele1=A1, allele2=A2)
#   gt_dir        directory with chr{1..22,X}_gt.txt (chr, pos, REF, ALT, GT...)
#   samples_file  one sample per line, same order as the GT columns
#   out_file      output Beagle (marker allele1 allele2 <3 probs per individual>)
#   label         label for log messages (e.g. "v1", "v2", "X")

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: Rscript lib_build_ref_beagle_panel.R <sites_file> <gt_dir> <samples_file> <out_file> <label> [comma-separated chr_list, default 1-22]")
}
SITES_FILE   <- args[1]
GT_DIR       <- args[2]
SAMPLES_FILE <- args[3]
OUT_FILE     <- args[4]
LABEL        <- args[5]
CHR_LIST     <- if (length(args) >= 6) strsplit(args[6], ",")[[1]] else as.character(1:22)

GT_FILES <- file.path(GT_DIR, sprintf("chr%s_gt.txt", CHR_LIST))

nuc_code <- c(A = 0L, C = 1L, G = 2L, T = 3L)
samples  <- readLines(SAMPLES_FILE)
n_ind    <- length(samples)
cat(sprintf("Panel %s: %d individuals, chromosomes: %s\n", LABEL, n_ind, paste(CHR_LIST, collapse = ",")))

sites <- fread(SITES_FILE, header = FALSE, col.names = c("chr", "pos", "A1", "A2"))
sites[, chr := as.character(chr)]
setkey(sites, chr, pos)

p_hom1 <- c(0.99, 0.005, 0.005)
p_het  <- c(0.005, 0.99, 0.005)
p_hom2 <- c(0.005, 0.005, 0.99)
p_miss <- c(1/3, 1/3, 1/3)

gt_to_probs <- function(gt, swap) {
    a <- substr(gt, 1, 1); b <- substr(gt, 3, 3); n <- length(gt)
    out <- matrix(p_miss, nrow = n, ncol = 3, byrow = TRUE)
    is00 <- a == "0" & b == "0"; is11 <- a == "1" & b == "1"
    ishet <- (a == "0" & b == "1") | (a == "1" & b == "0")
    sel <- is00 & !swap; out[sel, ] <- matrix(p_hom1, sum(sel), 3, byrow = TRUE)
    sel <- is00 &  swap; out[sel, ] <- matrix(p_hom2, sum(sel), 3, byrow = TRUE)
    sel <- is11 & !swap; out[sel, ] <- matrix(p_hom2, sum(sel), 3, byrow = TRUE)
    sel <- is11 &  swap; out[sel, ] <- matrix(p_hom1, sum(sel), 3, byrow = TRUE)
    out[ishet, ] <- matrix(p_het, sum(ishet), 3, byrow = TRUE)
    out
}

all_rows <- list(); total_kept <- 0L
for (f in GT_FILES) {
    if (!file.exists(f)) { cat(sprintf("SKIP: %s not found\n", f)); next }
    cat(sprintf("Processing %s ...\n", f))
    dt <- fread(f, header = FALSE)
    setnames(dt, 1:4, c("chr", "pos", "REF", "ALT"))
    dt[, chr := as.character(chr)]
    gt_cols <- 5:ncol(dt)
    m <- merge(dt, sites, by = c("chr", "pos"), all = FALSE)
    if (nrow(m) == 0) next
    match_fwd <- m$REF == m$A1 & m$ALT == m$A2
    match_rev <- m$REF == m$A2 & m$ALT == m$A1
    keep <- match_fwd | match_rev
    m <- m[keep]; swap_vec <- match_rev[keep]
    total_kept <- total_kept + nrow(m)
    if (nrow(m) == 0) next
    marker  <- paste0(m$chr, "_", m$pos)
    allele1 <- nuc_code[m$A1]; allele2 <- nuc_code[m$A2]
    gt_mat  <- as.matrix(m[, ..gt_cols])
    probs_list <- vector("list", n_ind)
    for (j in seq_len(n_ind)) probs_list[[j]] <- gt_to_probs(gt_mat[, j], swap_vec)
    probs_block <- do.call(cbind, probs_list)
    block <- data.table(marker = marker, allele1 = allele1, allele2 = allele2)
    block <- cbind(block, as.data.table(probs_block))
    all_rows[[f]] <- block
}
result <- rbindlist(all_rows)
fwrite(result, OUT_FILE, sep = "\t", col.names = FALSE)
cat(sprintf("Panel %s DONE: %d sites retained -> %s\n", LABEL, nrow(result), OUT_FILE))
