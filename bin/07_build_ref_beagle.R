#!/usr/bin/env Rscript
# 07_build_ref_beagle.R — convierte genotipos GT (bcftools) de 100 individuos
# de referencia 1KGP a formato pseudo-Beagle, alineado a los mismos sitios y
# codificación de alelos que el beagle de los pools (out_global/angsd/ngsadmix/pools_hapmap3.beagle.gz)
#
# Salida: out_global/angsd/ngsadmix/ref100.beagle (sin comprimir, mismo orden de marker que el panel de sitios)
suppressPackageStartupMessages(library(data.table))

SITES_FILE   <- "/tmp/angsd_sites_hapmap3.txt"          # chr pos A1 A2 (orden usado por ANGSD, allele1=A1, allele2=A2)
GT_FILES     <- sprintf("/tmp/ref_extract/chr%s_gt.txt", c(1,6,10,13,17,20,22))
SAMPLES_FILE <- "/tmp/ref_sample_order.txt"
OUT_FILE     <- "out_global/angsd/ngsadmix/ref100.beagle"

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
    # gt: vector character "0|0","0|1","1|0","1|1","."  etc
    # swap: logical vector, same length as gt (per-site allele orientation)
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
