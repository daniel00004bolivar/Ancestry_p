#!/usr/bin/env Rscript
# 30_build_oriented_sites.R — one oriented SNP universe for the revised analysis
#
# Memo items A4 (fixed allele orientation) and A7 (unadmixed source proxies).
#
# Source: UCSC hgdpGeo table (data/hgdp/hgdpGeo.txt.gz), HGDP-CEPH allele
# frequencies for 53 populations on hg19. Validated in this project:
#   - popFreqs are frequencies of the ANCESTRAL allele (cor with HapMap3
#     YRI/CEU/CHB derived-allele frequency = -0.97/-0.98/-0.99);
#   - column order is UCSC's alphabetical order of the 53 populations
#     (Yoruba matches HapMap YRI at |r| = 0.97 but CEU only at 0.73).
#
# Orientation: every site is coded ONCE as ancestral (A) / derived (D) and
# that coding is applied to every pool and every reference. Frequencies
# below are all derived-allele frequencies. Strand-ambiguous SNPs (A/T, C/G)
# are dropped; sites whose GRCh37 reference base is neither allele are dropped.
#
# Source groups (chromosome counts inferred as 1 / min non-zero frequency):
#   NAT = Maya 42 + Pima 28 + Karitiana 26 + Surui 16 + Colombian 14 = 126
#   AFR = Yoruba 42 + Mandenka 44                                  =  86
#   EUR = French 56 + Basque 48 + Italian 24 + Tuscan 14 + Sardinian 56 = 198
# Group frequency = chromosome-weighted mean of population frequencies.
#
# Output: out_global/revision/sites_oriented.tsv.gz
#   chr pos rsid ref anc der, then f_<group> n_<group> for every group in
#   GROUPS (main AFR/EUR/NAT plus the alternatives), transversion
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))

HGDP_FILE <- "data/hgdp/hgdpGeo.txt.gz"
REF_FASTA <- "data/ref/human_g1k_v37_decoy.fasta"
OUT_DIR   <- "out_global/revision"
OUT_FILE  <- file.path(OUT_DIR, "sites_oriented.tsv.gz")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

POPS <- c("Adygei","Balochi","BantuKenya","BantuSouthAfrica","Basque","Bedouin",
          "BiakaPygmy","Brahui","Burusho","Cambodian","Colombian","Dai","Daur",
          "Druze","French","Han","Han-NChina","Hazara","Hezhen","Italian",
          "Japanese","Kalash","Karitiana","Lahu","Makrani","Mandenka","Maya",
          "MbutiPygmy","Melanesian","Miao","Mongola","Mozabite","Naxi",
          "Orcadian","Oroqen","Palestinian","Papuan","Pathan","Pima","Russian",
          "San","Sardinian","She","Sindhi","Surui","Tu","Tujia","Tuscan",
          "Uygur","Xibo","Yakut","Yi","Yoruba")
GROUPS <- list(
    AFR = c(Yoruba = 42, Mandenka = 44),
    EUR = c(French = 56, Basque = 48, Italian = 24, Tuscan = 14, Sardinian = 56),
    NAT = c(Maya = 42, Pima = 28, Karitiana = 26, Surui = 16, Colombian = 14),
    # Alternative proxies for the reference-panel sensitivity (bin/36)
    NATmeso = c(Maya = 42, Pima = 28),
    NATsa   = c(Karitiana = 26, Surui = 16, Colombian = 14),
    EURsw   = c(French = 56, Basque = 48),
    EURs    = c(Italian = 24, Tuscan = 14, Sardinian = 56),
    AFRyri  = c(Yoruba = 42)
)

cat("Reading HGDP table...\n")
h <- fread(cmd = paste("gzcat", HGDP_FILE), sep = "\t", header = FALSE,
           col.names = c("bin", "chrom", "start", "end", "rsid", "anc", "der", "freqs"))
h <- h[chrom %in% paste0("chr", 1:22)]
fr <- tstrsplit(sub(",$", "", h$freqs), ",", fixed = TRUE)
stopifnot(length(fr) == length(POPS))
names(fr) <- POPS

s <- h[, .(chr = as.integer(sub("chr", "", chrom)), pos = end, rsid, anc, der)]
for (g in names(GROUPS)) {
    w <- GROUPS[[g]]
    anc_f <- Reduce(`+`, Map(function(p, n) as.numeric(fr[[p]]) * n, names(w), w)) / sum(w)
    s[, (paste0("f_", g)) := 1 - anc_f]          # derived-allele frequency
    s[, (paste0("n_", g)) := sum(w)]
}
rm(h, fr); invisible(gc())

# ── Drop strand-ambiguous and multi-mapped rows ──────────────────────────────
pair <- paste0(pmin(s$anc, s$der), pmax(s$anc, s$der))
s <- s[!pair %in% c("AT", "CG") & anc %chin% c("A","C","G","T") & der %chin% c("A","C","G","T")]
s <- unique(s, by = c("chr", "pos"))
s <- s[complete.cases(s)]

# ── GRCh37 reference base ────────────────────────────────────────────────────
reg <- tempfile(fileext = ".txt")
fwrite(s[, .(r = paste0(chr, ":", pos, "-", pos))], reg, col.names = FALSE)
cat("Looking up", nrow(s), "reference bases...\n")
fa <- system2("samtools", c("faidx", REF_FASTA, "-r", reg), stdout = TRUE)
s[, ref := toupper(fa[!startsWith(fa, ">")])]
unlink(reg)
n0 <- nrow(s)
s <- s[ref == anc | ref == der]
cat(sprintf("Reference base matches one allele: %d of %d (%.2f%%)\n", nrow(s), n0, 100 * nrow(s) / n0))

purines <- c("A", "G")
s[, transversion := (anc %chin% purines) != (der %chin% purines)]
setorder(s, chr, pos)
fcols <- as.vector(rbind(paste0("f_", names(GROUPS)), paste0("n_", names(GROUPS))))
fwrite(s[, c("chr", "pos", "rsid", "ref", "anc", "der", fcols, "transversion"), with = FALSE],
       OUT_FILE, sep = "\t")
cat(sprintf("Wrote %s: %d autosomal sites (%.1f%% transversions)\n",
            OUT_FILE, nrow(s), 100 * mean(s$transversion)))
