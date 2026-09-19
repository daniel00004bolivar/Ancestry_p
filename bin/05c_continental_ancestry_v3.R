#!/usr/bin/env Rscript
# 05c_continental_ancestry_v3.R — FST pool-seq vs 5 superpoblaciones (AFR/AMR/EAS/EUR/SAS)
#
# Combina:
#   - HapMap3: frecuencias AFR (YRI+MKK+LWK), EAS (CHB+CHD+JPT), EUR (TSI+CEU)
#   - 1KGP Phase 3: frecuencias AMR y SAS (superpoblaciones)
# Usa conteos ACGT de ANGSD (no doMaf) para frecuencias del pool.
#
# Uso:
#   Rscript bin/05c_continental_ancestry_v3.R \
#       --hapmap  data/iadmix/DATA/hapmap3.8populations.hg19.freqs \
#       --panel1k data/iadmix/DATA/hapmap3_5pops.freqs \
#       --counts  /tmp/POOL1_hapmap3_counts.txt,...
#       --samples POOL1,POOL2,HOSPITAL \
#       --out     out_global/angsd/continental
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    i <- which(args == flag); if (length(i)) args[i+1] else default
}

HAPMAP_FILE  <- get_arg("--hapmap",  "data/iadmix/DATA/hapmap3.8populations.hg19.freqs")
PANEL1K_FILE <- get_arg("--panel1k", "data/iadmix/DATA/hapmap3_5pops.freqs")
COUNTS_FILES <- strsplit(get_arg("--counts",
    "/tmp/POOL1_hapmap3_counts.txt,/tmp/POOL2_hapmap3_counts.txt,/tmp/HOSPITAL_hapmap3_counts.txt"),
    ",")[[1]]
SAMPLES      <- strsplit(get_arg("--samples", "POOL1,POOL2,HOSPITAL"), ",")[[1]]
OUT_DIR      <- get_arg("--out", "out_global/angsd/continental")
MIN_DEPTH    <- as.integer(get_arg("--min_depth",   "10"))
MIN_MAF_REF  <- as.double( get_arg("--min_maf_ref", "0.05"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Panel HapMap3: AFR, EAS, EUR ──────────────────────────────────────────
cat("Cargando HapMap3 (AFR, EAS, EUR)...\n")
hm <- fread(HAPMAP_FILE, header = TRUE)
setnames(hm, names(hm)[1], "chr")
hm[, chr := sub("^chr", "", chr)]
setnames(hm, "position", "pos")

hm[, allele1 := pmin(A1, A2)]
hm[, allele2 := pmax(A1, A2)]
for (p in c("YRI","CHB","CHD","TSI","MKK","LWK","CEU","JPT")) {
    hm[, (paste0(p,"_s")) := ifelse(A1 <= A2, get(p), 1 - get(p))]
}
hm[, AFR := (YRI_s + MKK_s + LWK_s) / 3]
hm[, EAS := (CHB_s + CHD_s + JPT_s) / 3]
hm[, EUR := (TSI_s + CEU_s) / 2]
hm_clean <- hm[, .(chr, pos, allele1, allele2, AFR, EAS, EUR)]
cat(sprintf("  HapMap3: %d SNPs\n", nrow(hm_clean)))

# ── 2. Panel 1KGP: AMR, SAS (y AFR/EAS/EUR de verificación) ─────────────────
cat("Cargando panel 1KGP 5 superpoblaciones (AMR, SAS)...\n")
kg <- fread(PANEL1K_FILE, header = TRUE)
setnames(kg, names(kg)[1], "chr")
kg[, chr := as.character(chr)]
kg[, pos := as.integer(pos)]
# Frecuencias ya estandarizadas al alelo lex-menor (hecho en 06_build_5pop_panel.sh)
# Columnas: chr pos allele1 allele2 AFR_1KGP EAS_1KGP EUR_1KGP AMR_1KGP SAS_1KGP
for (col in c("AMR_1KGP","SAS_1KGP","AFR_1KGP","EAS_1KGP","EUR_1KGP")) {
    kg[get(col) == ".", (col) := NA]
    kg[, (col) := as.numeric(get(col))]
}
cat(sprintf("  1KGP: %d SNPs con frecuencias superpoblacionales\n", nrow(kg)))

# ── 3. Combinar HapMap3 + 1KGP ───────────────────────────────────────────────
cat("Combinando paneles (merge por chr, pos, alelos)...\n")
combined <- merge(
    hm_clean,
    kg[, .(chr, pos, allele1, allele2, AMR = AMR_1KGP, SAS = SAS_1KGP)],
    by = c("chr","pos","allele1","allele2"),
    all.x = FALSE   # solo sitios en ambos paneles
)
cat(sprintf("  Combinado: %d SNPs comunes entre HapMap3 y 1KGP\n", nrow(combined)))

# ── 4. Cargar frecuencias de los pools (conteos ACGT) ────────────────────────
load_pool <- function(sname, counts_file, min_depth = 10) {
    if (!file.exists(counts_file)) { warning("No encontrado: ", counts_file); return(NULL) }
    cat(sprintf("  Cargando %s...\n", sname))
    ct <- fread(counts_file,
                col.names = c("chr","pos","depth","cA","cC","cG","cT"))
    ct[, chr := as.character(chr)]
    ct <- ct[depth >= min_depth]

    nuc <- c("A","C","G","T")
    ct[, c("cnt1","cnt2","al1","al2") := {
        m    <- cbind(cA, cC, cG, cT)
        ords <- t(apply(m, 1, order, decreasing = TRUE))
        cnt1 <- m[cbind(seq_len(nrow(m)), ords[,1])]
        cnt2 <- m[cbind(seq_len(nrow(m)), ords[,2])]
        al1  <- nuc[ords[,1]]
        al2  <- nuc[ords[,2]]
        list(cnt1, cnt2, al1, al2)
    }]
    ct <- ct[cnt2 > 0]
    ct[, allele1   := pmin(al1, al2)]
    ct[, allele2   := pmax(al1, al2)]
    ct[, freq_pool := ifelse(al1 <= al2, cnt1/(cnt1+cnt2), cnt2/(cnt1+cnt2))]
    ct[, maf_pool  := pmin(freq_pool, 1-freq_pool)]
    ct[, .(chr, pos, allele1, allele2, freq_pool, maf_pool, depth)]
}

cat("\nCargando conteos ACGT de los pools...\n")
pools <- Map(load_pool, SAMPLES, COUNTS_FILES, MoreArgs = list(min_depth = MIN_DEPTH))
names(pools) <- SAMPLES

# ── 5. FST Wright/Nei: pool vs cada continente ───────────────────────────────
wright_fst <- function(p1, p2) {
    p_bar <- (p1 + p2) / 2
    q_bar <- 1 - p_bar
    mean((p1 - p2)^2 / (4 * p_bar * q_bar + 1e-12), na.rm = TRUE)
}

continents <- c("AFR","AMR","EAS","EUR","SAS")
cont_colors <- c(AFR="#D55E00", AMR="#CC79A7", EAS="#009E73", EUR="#0072B2", SAS="#F0E442")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("FST Wright/Nei: pool-seq vs 5 superpoblaciones (HapMap3 + 1KGP)\n")
cat(sprintf("Min profundidad pool: %d | MAF ref mínima: %.2f\n", MIN_DEPTH, MIN_MAF_REF))
cat("═══════════════════════════════════════════════════════════════\n")

fst_rows    <- list()
merge_cache <- list()

for (sname in SAMPLES) {
    pd <- pools[[sname]]
    if (is.null(pd)) next

    m <- merge(pd, combined, by = c("chr","pos","allele1","allele2"))
    merge_cache[[sname]] <- m
    cat(sprintf("\n  %s: %d SNPs en panel combinado\n", sname, nrow(m)))

    for (cont in continents) {
        ref_f  <- m[[cont]]
        pool_f <- m$freq_pool
        valid  <- !is.na(ref_f) & ref_f >= MIN_MAF_REF & ref_f <= (1-MIN_MAF_REF)
        n_ok   <- sum(valid)
        if (n_ok < 100) {
            cat(sprintf("    vs %-3s: muy pocos sitios (%d)\n", cont, n_ok)); next
        }
        fst_val <- wright_fst(pool_f[valid], ref_f[valid])
        r_val   <- cor(pool_f[valid], ref_f[valid], use = "complete.obs")
        cat(sprintf("    vs %-3s: FST = %.5f  |  r = %.4f  (n = %d SNPs)\n",
                    cont, fst_val, r_val, n_ok))
        fst_rows[[paste(sname, cont)]] <- data.frame(
            pool=sname, continent=cont, n_snps=n_ok, fst=fst_val, r=r_val,
            stringsAsFactors=FALSE
        )
    }
}

fst_table <- do.call(rbind, fst_rows)
write.table(fst_table, file.path(OUT_DIR, "continental_fst_5pop.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ── 6. Ranking ────────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("Ranking por pool (de más cercano a más lejano)\n")
cat("═══════════════════════════════════════════════════════════════\n")
for (sname in SAMPLES) {
    sub <- fst_table[fst_table$pool == sname, ]
    if (nrow(sub) == 0) next
    r_order <- sub[order(sub$fst), ]
    cat(sprintf("\n  %s:\n", sname))
    for (i in seq_len(nrow(r_order))) {
        cat(sprintf("    %d. %-3s  FST = %.5f  r = %.4f\n",
                    i, r_order$continent[i], r_order$fst[i], r_order$r[i]))
    }
}

# ── 7. Heatmap FST ───────────────────────────────────────────────────────────
fst_p <- fst_table
fst_p$pool      <- factor(fst_p$pool,      levels=c("POOL1","POOL2","HOSPITAL"))
fst_p$continent <- factor(fst_p$continent, levels=rev(continents))

p_heat <- ggplot(fst_p, aes(x=pool, y=continent, fill=fst)) +
    geom_tile(color="white", linewidth=0.8) +
    geom_text(aes(label=sprintf("%.5f\n(r=%.3f)", fst, r)),
              size=3.2, lineheight=1.3) +
    scale_fill_gradientn(
        colors=c("#2166ac","#92c5de","#f7f7f7","#f4a582","#b2182b"),
        name="FST\n(menor=más similar)"
    ) +
    labs(title="Similitud genética pools vs 5 superpoblaciones",
         subtitle="FST Wright/Nei — conteos ACGT (ANGSD) vs HapMap3 + 1KGP Phase 3",
         x=NULL, y="Superpoblación de referencia") +
    theme_bw(base_size=13) +
    theme(plot.title=element_text(face="bold", size=14),
          axis.text=element_text(size=12))
ggsave(file.path(OUT_DIR, "fst_5pop_heatmap.png"), p_heat, width=7, height=5.5, dpi=180)
cat("\nFigura: fst_5pop_heatmap.png\n")

# ── 8. Barplot FST ───────────────────────────────────────────────────────────
fst_p2 <- fst_p
fst_p2$continent <- factor(fst_p2$continent, levels=continents)

p_bar <- ggplot(fst_p2, aes(x=continent, y=fst, fill=continent)) +
    geom_col(show.legend=FALSE) +
    geom_text(aes(label=sprintf("%.5f", fst)), vjust=-0.4, size=2.8) +
    facet_wrap(~pool) +
    scale_fill_manual(values=cont_colors) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
    labs(title="FST pool-seq vs 5 superpoblaciones de referencia",
         subtitle="Menor FST = mayor similitud genética",
         x="Superpoblación", y="FST medio (Wright/Nei)") +
    theme_bw(base_size=11) +
    theme(strip.text=element_text(face="bold", size=12),
          plot.title=element_text(face="bold"))
ggsave(file.path(OUT_DIR, "fst_5pop_barplot.png"), p_bar, width=10, height=4, dpi=180)
cat("Figura: fst_5pop_barplot.png\n")

# ── 9. Scatter vs AMR (el más relevante) ─────────────────────────────────────
for (sname in SAMPLES) {
    m <- merge_cache[[sname]]
    if (is.null(m) || nrow(m)==0) next
    sub   <- fst_table[fst_table$pool==sname & fst_table$continent=="AMR", ]
    if (nrow(sub)==0 || is.na(sub$fst[1])) next
    valid <- !is.na(m$AMR) & m$AMR>=MIN_MAF_REF & m$AMR<=(1-MIN_MAF_REF)
    m_ok  <- m[valid, ]
    set.seed(42)
    sub_p <- m_ok[sample(nrow(m_ok), min(20000, nrow(m_ok))), ]
    r_amr <- cor(sub_p$freq_pool, sub_p$AMR, use="complete.obs")
    p_sc <- ggplot(sub_p, aes(x=AMR, y=freq_pool)) +
        geom_point(alpha=0.07, size=0.4, color="#CC79A7") +
        geom_abline(slope=1, intercept=0, color="firebrick", linetype="dashed") +
        geom_smooth(method="lm", color="darkblue", se=FALSE, linewidth=0.8) +
        annotate("text", x=0.1, y=0.88,
                 label=sprintf("r = %.4f\nn = %d SNPs", r_amr, nrow(m_ok)),
                 hjust=0, size=4.5) +
        labs(title=sprintf("%s vs AMR (1KGP: CLM+MXL+PEL+PUR)", sname),
             subtitle="Conteos ACGT ANGSD vs frecuencias superpoblación americana 1KGP",
             x="Freq alelo1 — AMR 1KGP", y=sprintf("Freq alelo1 — %s (ANGSD)", sname)) +
        coord_fixed(xlim=c(0,1), ylim=c(0,1)) +
        theme_bw(base_size=12)
    ggsave(file.path(OUT_DIR, sprintf("scatter_%s_vs_AMR.png", sname)),
           p_sc, width=6, height=5.5, dpi=180)
}

cat("\nTabla: ", file.path(OUT_DIR, "continental_fst_5pop.tsv"), "\n")
