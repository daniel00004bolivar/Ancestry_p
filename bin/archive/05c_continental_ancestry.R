#!/usr/bin/env Rscript
# 05c_continental_ancestry.R — Ancestría continental via FST pool-seq vs HapMap3
#
# Compara frecuencias alélicas de los pools (ANGSD doMaf) contra el panel
# de referencia HapMap3 (8 poblaciones) para identificar similitud continental.
#
# Panel: data/iadmix/DATA/hapmap3.8populations.hg19.freqs
#   AFR = YRI + MKK + LWK  (Yoruba, Maasai, Luhya)
#   EAS = CHB + CHD + JPT  (Han-Beijing, Han-Denver, Japonés)
#   EUR = TSI + CEU         (Toscana, Utah-Europeo)
#
# Método: FST Wright/Nei (varianza de frecuencias) — válido para pool-seq sin
#   conocer N. Menor FST = mayor similitud genética con ese continente.
#
# Uso:
#   Rscript bin/05c_continental_ancestry.R \
#       --panel  data/iadmix/DATA/hapmap3.8populations.hg19.freqs \
#       --mafs_dir out_global/angsd/per_pool \
#       --out    out_global/angsd/continental
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    i <- which(args == flag); if (length(i)) args[i+1] else default
}

PANEL_FILE <- get_arg("--panel",    "data/iadmix/DATA/hapmap3.8populations.hg19.freqs")
MAFS_DIR   <- get_arg("--mafs_dir", "out_global/angsd/per_pool")
OUT_DIR    <- get_arg("--out",      "out_global/angsd/continental")
MIN_MAF    <- as.double(get_arg("--min_maf", "0.05"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Cargar y estandarizar el panel HapMap3 ────────────────────────────────
cat("Cargando panel HapMap3 (8 poblaciones, ~249k SNPs)...\n")
ref <- fread(PANEL_FILE, header = TRUE)
old_chr_col <- names(ref)[1]          # "#chrom" o "chrom"
setnames(ref, old_chr_col, "chr")

# Eliminar prefijo "chr" (el panel usa hg19 con prefijo; b37 no lo tiene)
ref[, chr := sub("^chr", "", chr)]
setnames(ref, "position", "pos")

# Las columnas de frecuencia son la frecuencia de A1 en cada población
pop_cols <- c("YRI", "CHB", "CHD", "TSI", "MKK", "LWK", "CEU", "JPT")

# Estandarizar alelo: usar el lex-menor de (A1, A2) como "allele1"
# y calcular la frecuencia de ese alelo estandarizado
ref[, allele1 := pmin(A1, A2)]
ref[, allele2 := pmax(A1, A2)]

for (p in pop_cols) {
    col_std <- paste0(p, "_std")
    ref[, (col_std) := ifelse(A1 <= A2, get(p), 1.0 - get(p))]
}

# Frecuencias continentales = media de las poblaciones por continente
ref[, AFR := (YRI_std + MKK_std + LWK_std) / 3.0]
ref[, EAS := (CHB_std + CHD_std + JPT_std) / 3.0]
ref[, EUR := (TSI_std + CEU_std)            / 2.0]

ref_clean <- ref[, .(chr, pos, allele1, allele2, AFR, EAS, EUR)]
cat(sprintf("  Referencia lista: %d SNPs\n", nrow(ref_clean)))

# ── 2. Cargar frecuencias de los pools (ANGSD doMaf) ────────────────────────
load_pool <- function(sname, dir) {
    f <- file.path(dir, paste0(sname, "_maf.mafs.gz"))
    if (!file.exists(f)) { warning("No encontrado: ", f); return(NULL) }
    cat(sprintf("  Cargando %s...\n", sname))
    dt <- fread(cmd = paste("gzip -cd", shQuote(f)),
                col.names = c("chr","pos","major","minor","maf","nInd"),
                colClasses = list(character = c("chromo","major","minor")))
    dt[, chr := as.character(chr)]
    # Estandarizar: freq del alelo lex-menor
    dt[, allele1   := pmin(major, minor)]
    dt[, allele2   := pmax(major, minor)]
    dt[, freq_pool := ifelse(major <= minor, 1.0 - maf, maf)]
    dt[, maf_pool  := pmin(freq_pool, 1.0 - freq_pool)]
    dt[, .(chr, pos, allele1, allele2, freq_pool, maf_pool)]
}

samples <- c("POOL1", "POOL2", "HOSPITAL")
cat("\nCargando pools...\n")
pools <- lapply(samples, load_pool, dir = MAFS_DIR)
names(pools) <- samples

# ── 3. FST Wright/Nei: pool vs continente ────────────────────────────────────
wright_fst_vec <- function(p1, p2) {
    p_bar <- (p1 + p2) / 2.0
    q_bar <- 1.0 - p_bar
    (p1 - p2)^2 / (4.0 * p_bar * q_bar + 1e-12)
}

cat("\n=== FST Wright/Nei: pool vs continentes HapMap3 ===\n")
cat(sprintf("(MAF mínima en pool = %.2f; mínima en referencia = %.2f)\n", MIN_MAF, MIN_MAF))
cat("─────────────────────────────────────────────────────────────\n")

continents <- c("AFR", "EAS", "EUR")
fst_rows <- list()
merged_cache <- list()   # guardamos merged para los scatter plots

for (sname in samples) {
    pd <- pools[[sname]]
    if (is.null(pd)) next

    # Unir por posición y par alélico idéntico
    m <- merge(pd, ref_clean, by = c("chr", "pos", "allele1", "allele2"))
    m <- m[maf_pool >= MIN_MAF]          # filtro MAF en el pool
    merged_cache[[sname]] <- m
    cat(sprintf("\n  %s: %d SNPs compartidos con panel\n", sname, nrow(m)))

    for (cont in continents) {
        ref_f <- m[[cont]]
        pool_f <- m$freq_pool

        # Filtro MAF en referencia
        valid <- !is.na(ref_f) & ref_f >= MIN_MAF & ref_f <= (1 - MIN_MAF)
        n_ok  <- sum(valid)

        if (n_ok < 100) {
            cat(sprintf("    vs %s: muy pocos sitios (%d)\n", cont, n_ok))
            next
        }

        fst_vec <- wright_fst_vec(pool_f[valid], ref_f[valid])
        fst_val <- mean(fst_vec, na.rm = TRUE)
        r_val   <- cor(pool_f[valid], ref_f[valid], use = "complete.obs")

        cat(sprintf("    vs %-3s: FST = %.5f  |  r = %.4f  (n = %d SNPs)\n",
                    cont, fst_val, r_val, n_ok))

        fst_rows[[paste(sname, cont)]] <- data.frame(
            pool = sname, continent = cont,
            n_snps = n_ok,
            fst    = fst_val,
            r      = r_val,
            stringsAsFactors = FALSE
        )
    }
}

fst_table <- do.call(rbind, fst_rows)
write.table(fst_table, file.path(OUT_DIR, "continental_fst.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── 4. Resumen: continente más cercano por pool ───────────────────────────────
cat("\n─────────────────────────────────────────────────────────────\n")
cat("=== Continente más cercano por pool (menor FST) ===\n")
for (sname in samples) {
    sub <- fst_table[fst_table$pool == sname, ]
    if (nrow(sub) == 0) next
    best <- sub[which.min(sub$fst), ]
    rank_order <- sub[order(sub$fst), ]
    cat(sprintf(
        "  %-10s →  1o: %-3s (FST=%.5f)  2o: %-3s (FST=%.5f)  3o: %-3s (FST=%.5f)\n",
        sname,
        rank_order$continent[1], rank_order$fst[1],
        rank_order$continent[2], rank_order$fst[2],
        rank_order$continent[3], rank_order$fst[3]
    ))
}
cat("\nNota: El panel HapMap3 no incluye población AMR (amerindia/latina).\n")
cat("      La ancestría amerindia detectada por iAdmix (~18%%) eleva el FST\n")
cat("      vs todos los continentes disponibles (no hay referencia directa).\n")

# ── 5. Figura 1: Heatmap de FST ───────────────────────────────────────────────
fst_table$pool      <- factor(fst_table$pool,      levels = c("POOL1","POOL2","HOSPITAL"))
fst_table$continent <- factor(fst_table$continent, levels = c("AFR","EUR","EAS"))

p_heat <- ggplot(fst_table, aes(x = continent, y = pool, fill = fst)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("FST = %.5f\nr = %.4f", fst, r)),
              size = 3.5, lineheight = 1.3) +
    scale_fill_gradientn(
        colors  = c("#2166ac","#92c5de","#f7f7f7","#f4a582","#d6604d","#b2182b"),
        name    = "FST\n(menor = más similar)",
        limits  = c(min(fst_table$fst) * 0.9, max(fst_table$fst) * 1.1)
    ) +
    labs(
        title    = "Similitud genética: pools vs referencia HapMap3",
        subtitle = "FST Wright/Nei entre pool-seq y frecuencias continentales de referencia",
        x = "Continente de referencia",
        y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        axis.text     = element_text(size = 12)
    )
ggsave(file.path(OUT_DIR, "continental_fst_heatmap.png"),
       p_heat, width = 7, height = 4.5, dpi = 180)
cat("\nFigura: continental_fst_heatmap.png\n")

# ── 6. Figura 2: Barras de FST por pool ──────────────────────────────────────
p_bar <- ggplot(fst_table, aes(x = continent, y = fst, fill = continent)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.5f", fst)), vjust = -0.4, size = 3) +
    facet_wrap(~pool) +
    scale_fill_manual(values = c(AFR = "#D55E00", EAS = "#009E73", EUR = "#0072B2")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title    = "FST pool-seq vs continentes de referencia (HapMap3)",
        subtitle = "Menor FST = mayor similitud con ese continente",
        x = "Continente", y = "FST medio (Wright/Nei)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        strip.text  = element_text(face = "bold", size = 12),
        plot.title  = element_text(face = "bold")
    )
ggsave(file.path(OUT_DIR, "continental_fst_barplot.png"),
       p_bar, width = 9, height = 4, dpi = 180)
cat("Figura: continental_fst_barplot.png\n")

# ── 7. Figura 3: Scatter pool vs AFR (esperado como el más cercano) ───────────
for (sname in samples) {
    m <- merged_cache[[sname]]
    if (is.null(m) || nrow(m) == 0) next

    best_cont <- fst_table[fst_table$pool == sname, ]
    if (nrow(best_cont) == 0) next
    best_cont <- as.character(best_cont$continent[which.min(best_cont$fst)])

    valid <- m[[best_cont]] >= MIN_MAF & m[[best_cont]] <= (1 - MIN_MAF)
    m_ok  <- m[valid]
    set.seed(42)
    sub   <- m_ok[sample(.N, min(25000, .N))]

    r_val <- cor(sub$freq_pool, sub[[best_cont]], use = "complete.obs")
    n_tot <- nrow(m_ok)

    p_sc <- ggplot(sub, aes(x = get(best_cont), y = freq_pool)) +
        geom_point(alpha = 0.08, size = 0.4, color = "steelblue") +
        geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "dashed") +
        geom_smooth(method = "lm", color = "darkblue", se = FALSE, linewidth = 0.8) +
        annotate("text", x = 0.1, y = 0.9,
                 label = sprintf("r = %.4f\nn = %d SNPs", r_val, n_tot),
                 hjust = 0, size = 4, color = "black") +
        labs(
            title    = sprintf("%s vs %s (continente más cercano)", sname, best_cont),
            subtitle = "Frecuencias alélicas estandarizadas — ANGSD pool-seq vs HapMap3",
            x = sprintf("Frecuencia alelo1 — %s (HapMap3 referencia)", best_cont),
            y = sprintf("Frecuencia alelo1 — %s (ANGSD)", sname)
        ) +
        coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
        theme_bw(base_size = 12)

    fname <- sprintf("scatter_%s_vs_%s.png", sname, best_cont)
    ggsave(file.path(OUT_DIR, fname), p_sc, width = 6, height = 5.5, dpi = 180)
    cat(sprintf("Figura: %s\n", fname))
}

cat("\n=== Tabla guardada ===\n")
cat("  ", file.path(OUT_DIR, "continental_fst.tsv"), "\n")
cat("\nInterpretación:\n")
cat("  FST < 0.01  : diferenciación muy baja (ancestría predominante)\n")
cat("  FST 0.01-0.05: diferenciación baja (ancestría parcial)\n")
cat("  FST 0.05-0.15: diferenciación moderada (poca ancestría compartida)\n")
cat("  FST > 0.15  : diferenciación alta (ancestría distante)\n")
