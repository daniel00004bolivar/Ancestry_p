#!/usr/bin/env Rscript
# 05c_continental_ancestry_v2.R — Ancestría continental via FST pool-seq vs HapMap3
#
# Versión 2: usa archivos de CONTEOS (A,C,G,T) extraídos en posiciones HapMap3,
# en lugar de los archivos doMaf (que convergen a ~0.5 para pool-seq con n=1).
#
# Frecuencia del pool = (conteo alelo más frecuente) / (conteo total en el sitio),
# orientada al alelo lex-menor del par para comparabilidad con la referencia.
#
# Uso:
#   Rscript bin/05c_continental_ancestry_v2.R \
#       --panel     data/iadmix/DATA/hapmap3.8populations.hg19.freqs \
#       --counts    /tmp/POOL1_hapmap3_counts.txt,/tmp/POOL2_hapmap3_counts.txt,/tmp/HOSPITAL_hapmap3_counts.txt \
#       --samples   POOL1,POOL2,HOSPITAL \
#       --out       out_global/angsd/continental
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    i <- which(args == flag); if (length(i)) args[i+1] else default
}

PANEL_FILE   <- get_arg("--panel",   "data/iadmix/DATA/hapmap3.8populations.hg19.freqs")
COUNTS_FILES <- strsplit(get_arg("--counts",
    "/tmp/POOL1_hapmap3_counts.txt,/tmp/POOL2_hapmap3_counts.txt,/tmp/HOSPITAL_hapmap3_counts.txt"),
    ",")[[1]]
SAMPLES      <- strsplit(get_arg("--samples", "POOL1,POOL2,HOSPITAL"), ",")[[1]]
OUT_DIR      <- get_arg("--out", "out_global/angsd/continental")
MIN_DEPTH    <- as.integer(get_arg("--min_depth", "10"))   # mínima profundidad para incluir sitio
MIN_MAF_REF  <- as.double( get_arg("--min_maf_ref", "0.05"))  # MAF mínima en referencia

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Panel HapMap3: frecuencias por continente ──────────────────────────────
cat("Cargando panel HapMap3...\n")
ref <- fread(PANEL_FILE, header = TRUE)
setnames(ref, names(ref)[1], "chr")
ref[, chr := sub("^chr", "", chr)]
setnames(ref, "position", "pos")

pop_cols <- c("YRI","CHB","CHD","TSI","MKK","LWK","CEU","JPT")

# Estandarizar: freq de A1 para el alelo lex-menor del par
ref[, allele1 := pmin(A1, A2)]
ref[, allele2 := pmax(A1, A2)]
for (p in pop_cols) {
    ref[, (paste0(p,"_s")) := ifelse(A1 <= A2, get(p), 1 - get(p))]
}
ref[, AFR := (YRI_s + MKK_s + LWK_s) / 3]
ref[, EAS := (CHB_s + CHD_s + JPT_s) / 3]
ref[, EUR := (TSI_s + CEU_s) / 2]

ref_clean <- ref[, .(chr, pos, allele1, allele2, AFR, EAS, EUR)]
setkey(ref_clean, chr, pos)
cat(sprintf("  Panel: %d SNPs\n", nrow(ref_clean)))

# ── 2. Cargar conteos del pool y computar frecuencias ─────────────────────────
# Formato archivo: chr pos totDepth totA totC totG totT
load_counts_freq <- function(sname, counts_file, min_depth = 10) {
    if (!file.exists(counts_file)) { warning("No encontrado: ", counts_file); return(NULL) }
    cat(sprintf("  Cargando conteos %s (%s)...\n", sname, counts_file))
    ct <- fread(counts_file,
                col.names = c("chr","pos","depth","cA","cC","cG","cT"))
    ct[, chr := as.character(chr)]
    ct <- ct[depth >= min_depth]
    cat(sprintf("    Sitios con depth >= %d: %d\n", min_depth, nrow(ct)))

    # Identificar los dos alelos con más conteos (los 2 más frecuentes)
    # Usar vector de alelos: A, C, G, T
    nuc <- c("A","C","G","T")
    ct[, c("cnt1","cnt2","al1","al2") := {
        # Para cada fila, ordenar los 4 conteos de mayor a menor
        m <- cbind(cA, cC, cG, cT)
        ords <- t(apply(m, 1, order, decreasing = TRUE))
        cnt_sorted <- m[cbind(seq_len(nrow(m)), ords[,1])]
        cnt2nd     <- m[cbind(seq_len(nrow(m)), ords[,2])]
        al_top     <- nuc[ords[,1]]
        al_2nd     <- nuc[ords[,2]]
        list(cnt_sorted, cnt2nd, al_top, al_2nd)
    }]

    # Solo retener sitios donde al menos 2 alelos tienen conteos (bialélicos)
    ct <- ct[cnt2 > 0]

    # Estandarizar: alelo lex-menor = allele1
    ct[, allele1 := pmin(al1, al2)]
    ct[, allele2 := pmax(al1, al2)]

    # Frecuencia del alelo lex-menor
    ct[, freq_a1 := ifelse(
        al1 <= al2,
        cnt1 / (cnt1 + cnt2),    # alelo más frecuente ES el lex-menor
        cnt2 / (cnt1 + cnt2)     # alelo lex-menor es el 2do más frecuente
    )]

    # MAF en el pool (mínimo de freq_a1 y 1-freq_a1)
    ct[, maf_pool := pmin(freq_a1, 1 - freq_a1)]

    ct[, .(chr, pos, allele1, allele2, freq_pool = freq_a1, maf_pool, depth)]
}

pools <- Map(load_counts_freq, SAMPLES, COUNTS_FILES,
             MoreArgs = list(min_depth = MIN_DEPTH))
names(pools) <- SAMPLES

# ── 3. FST Wright/Nei: pool vs continente ────────────────────────────────────
wright_fst <- function(p1, p2) {
    p_bar <- (p1 + p2) / 2
    q_bar <- 1 - p_bar
    mean((p1 - p2)^2 / (4 * p_bar * q_bar + 1e-12), na.rm = TRUE)
}

cat("\n=== FST Wright/Nei: pool vs continentes (frecuencias de conteos ACGT) ===\n")
cat(sprintf("(Min profundidad pool: %d; MAF ref mínima: %.2f)\n", MIN_DEPTH, MIN_MAF_REF))
cat("─────────────────────────────────────────────────────────────────\n")

continents <- c("AFR","EAS","EUR")
# Poblaciones individuales del panel HapMap3
pop_std_cols <- c("YRI_s","MKK_s","LWK_s","CHB_s","CHD_s","JPT_s","TSI_s","CEU_s")
pop_labels   <- c("YRI(AFR)","MKK(AFR)","LWK(AFR)","CHB(EAS)","CHD(EAS)","JPT(EAS)","TSI(EUR)","CEU(EUR)")

fst_rows   <- list()
merge_cache <- list()

for (sname in SAMPLES) {
    pd <- pools[[sname]]
    if (is.null(pd)) next

    # ref_clean no tiene las columnas pop_std; usar ref completo
    ref_full <- ref[, c("chr","pos","allele1","allele2","AFR","EAS","EUR",
                        pop_std_cols), with=FALSE]
    m <- merge(pd, ref_full, by = c("chr","pos","allele1","allele2"))
    merge_cache[[sname]] <- m
    cat(sprintf("\n  %s: %d SNPs coinciden con panel de referencia\n", sname, nrow(m)))

    # FST por continente
    for (cont in continents) {
        ref_f  <- m[[cont]]
        pool_f <- m$freq_pool
        valid  <- !is.na(ref_f) & ref_f >= MIN_MAF_REF & ref_f <= (1 - MIN_MAF_REF)
        n_ok   <- sum(valid)
        if (n_ok < 100) { cat(sprintf("    vs %s: muy pocos sitios (%d)\n", cont, n_ok)); next }
        fst_val <- wright_fst(pool_f[valid], ref_f[valid])
        r_val   <- cor(pool_f[valid], ref_f[valid], use = "complete.obs")
        cat(sprintf("    vs %-3s: FST = %.5f  |  r = %.4f  (n = %d SNPs)\n", cont, fst_val, r_val, n_ok))
        fst_rows[[paste(sname, cont)]] <- data.frame(
            pool = sname, continent = cont,
            n_snps = n_ok, fst = fst_val, r = r_val,
            stringsAsFactors = FALSE
        )
    }

    # FST por población individual
    cat(sprintf("  --- Poblaciones individuales en %s ---\n", sname))
    pop_fst_rows <- list()
    for (i in seq_along(pop_std_cols)) {
        pc <- pop_std_cols[i]; pl <- pop_labels[i]
        ref_f  <- m[[pc]]
        pool_f <- m$freq_pool
        valid  <- !is.na(ref_f) & ref_f >= MIN_MAF_REF & ref_f <= (1 - MIN_MAF_REF)
        n_ok   <- sum(valid)
        if (n_ok < 10) next
        fst_val <- wright_fst(pool_f[valid], ref_f[valid])
        r_val   <- cor(pool_f[valid], ref_f[valid], use = "complete.obs")
        cat(sprintf("    vs %-12s: FST = %.5f  r = %.4f\n", pl, fst_val, r_val))
        pop_fst_rows[[pl]] <- data.frame(pool=sname, population=pl,
            n_snps=n_ok, fst=fst_val, r=r_val, stringsAsFactors=FALSE)
    }
}

fst_table <- do.call(rbind, fst_rows)
write.table(fst_table, file.path(OUT_DIR, "continental_fst_v2.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── 4. Ranking por pool ───────────────────────────────────────────────────────
cat("\n─────────────────────────────────────────────────────────────────\n")
cat("=== Continente más cercano por pool (menor FST) ===\n")
for (sname in SAMPLES) {
    sub <- fst_table[fst_table$pool == sname, ]
    if (nrow(sub) == 0) next
    r_order <- sub[order(sub$fst), ]
    cat(sprintf("  %-10s →  1o: %-3s (FST=%.5f, r=%.4f)  2o: %-3s (FST=%.5f)  3o: %-3s (FST=%.5f)\n",
        sname,
        r_order$continent[1], r_order$fst[1], r_order$r[1],
        r_order$continent[2], r_order$fst[2],
        r_order$continent[3], r_order$fst[3]))
}

# ── 5. Figura: heatmap FST ────────────────────────────────────────────────────
fst_plot <- fst_table
fst_plot$pool      <- factor(fst_plot$pool,      levels = c("POOL1","POOL2","HOSPITAL"))
fst_plot$continent <- factor(fst_plot$continent, levels = c("AFR","EUR","EAS"))

p_heat <- ggplot(fst_plot, aes(x = continent, y = pool, fill = fst)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("FST = %.5f\nr = %.4f", fst, r)),
              size = 3.8, lineheight = 1.4) +
    scale_fill_gradientn(
        colors = c("#2166ac","#92c5de","#f7f7f7","#f4a582","#b2182b"),
        name   = "FST\n(menor = más similar)"
    ) +
    labs(title    = "Similitud genética: pools vs referencia HapMap3",
         subtitle = "FST Wright/Nei (frecuencias de conteos ACGT) vs continentes de referencia",
         x = "Continente de referencia", y = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text  = element_text(size = 12))
ggsave(file.path(OUT_DIR, "continental_fst_heatmap_v2.png"),
       p_heat, width = 7, height = 4.5, dpi = 180)

# ── 6. Figura: barras FST ─────────────────────────────────────────────────────
p_bar <- ggplot(fst_plot, aes(x = continent, y = fst, fill = continent)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.5f", fst)), vjust = -0.4, size = 3) +
    facet_wrap(~pool) +
    scale_fill_manual(values = c(AFR="#D55E00", EAS="#009E73", EUR="#0072B2")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title    = "FST pool-seq vs continentes (conteos ACGT corregidos)",
         subtitle = "Menor FST = mayor similitud; panel HapMap3 no incluye AMR",
         x = "Continente", y = "FST medio (Wright/Nei)") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 12),
          plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "continental_fst_barplot_v2.png"),
       p_bar, width = 9, height = 4, dpi = 180)

# ── 7. Scatter: pool vs continente más cercano ───────────────────────────────
for (sname in SAMPLES) {
    m <- merge_cache[[sname]]
    if (is.null(m) || nrow(m) == 0) next

    sub <- fst_table[fst_table$pool == sname, ]
    if (nrow(sub) == 0) next
    best_cont <- sub$continent[which.min(sub$fst)]

    valid <- m[[best_cont]] >= MIN_MAF_REF & m[[best_cont]] <= (1 - MIN_MAF_REF)
    m_ok  <- m[valid, ]
    set.seed(42)
    sub_plot <- m_ok[sample(nrow(m_ok), min(25000, nrow(m_ok))), ]
    r_plot   <- cor(sub_plot$freq_pool, sub_plot[[best_cont]], use = "complete.obs")

    p_sc <- ggplot(sub_plot, aes(x = .data[[best_cont]], y = freq_pool)) +
        geom_point(alpha = 0.07, size = 0.4, color = "steelblue") +
        geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "dashed") +
        geom_smooth(method = "lm", color = "darkblue", se = FALSE, linewidth = 0.8) +
        annotate("text", x = 0.1, y = 0.88,
                 label = sprintf("r = %.4f\nn = %d SNPs", r_plot, nrow(m_ok)),
                 hjust = 0, size = 4.5, color = "black") +
        labs(title    = sprintf("%s vs %s (continente más cercano)", sname, best_cont),
             subtitle = "Frecuencias de conteos ACGT (ANGSD) vs HapMap3 referencia",
             x = sprintf("Freq alelo1 — %s referencia (HapMap3)", best_cont),
             y = sprintf("Freq alelo1 — %s (ANGSD pool-seq)", sname)) +
        coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
        theme_bw(base_size = 12)
    ggsave(file.path(OUT_DIR, sprintf("scatter_v2_%s_vs_%s.png", sname, best_cont)),
           p_sc, width = 6, height = 5.5, dpi = 180)
    cat(sprintf("Figura: scatter_v2_%s_vs_%s.png\n", sname, best_cont))
}

cat("\nTabla: ", file.path(OUT_DIR, "continental_fst_v2.tsv"), "\n")
cat("\nNota: panel HapMap3 sin AMR. La ancestría amerindia (~18%% en iAdmix)\n")
cat("no tiene referencia directa — eleva FST vs todos los continentes.\n")
