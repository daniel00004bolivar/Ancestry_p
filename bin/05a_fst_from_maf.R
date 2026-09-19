#!/usr/bin/env Rscript
# 05a_fst_from_maf.R — FST Hudson-Bhatia desde archivos .mafs.gz de ANGSD
#
# Usa los archivos *_maf.mafs.gz (doMaf) en lugar de los counts completos
# para evitar problemas de memoria en análisis genoma completo.
#
# Estandarización de alelos: dado que doMajorMinor 1 asigna mayor/menor
# de forma independiente por pool, se usa el par (allele_min, allele_max)
# para definir una frecuencia consistente: freq del alelo lexicográficamente
# menor en cada sitio. Si major < minor → freq_std = 1-maf; si no → freq_std = maf.
#
# Estimador FST Bhatia et al. 2013:
#   num   = (p1 - p2)^2 - p1*(1-p1)/n1 - p2*(1-p2)/n2
#   denom = p1*(1-p2) + p2*(1-p1)
# Para pool-seq, n = 2 (diploid encoding) → corrección de muestra mínima con pools.
#
# Uso:
#   Rscript bin/05a_fst_from_maf.R \
#       --mafs_dir out_global/angsd/per_pool \
#       --out out_global/angsd/fst/hudson_fst_summary.tsv
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    i <- which(args == flag); if (length(i)) args[i+1] else default
}
MAFS_DIR <- get_arg("--mafs_dir", "out_global/angsd/per_pool")
OUT_FILE <- get_arg("--out",      "out_global/angsd/fst/hudson_fst_summary.tsv")
MIN_MAF  <- as.double(get_arg("--min_maf", "0.01"))

# ── Carga MAF y estandariza orientación de alelos ─────────────────────────────
load_maf_std <- function(sample, dir, min_maf = 0.01) {
    f <- file.path(dir, paste0(sample, "_maf.mafs.gz"))
    if (!file.exists(f)) { warning("No encontrado: ", f); return(NULL) }
    cat(sprintf("  Cargando %s...\n", basename(f)))
    # chromo, position, major, minor, knownEM, nInd
    dt <- fread(f, col.names = c("chr", "pos", "major", "minor", "maf", "nInd"),
                colClasses = list(character = c("chr","major","minor")))

    # Frecuencia estandarizada: siempre para el alelo lex-menor del par
    # Si major < minor: major es el alelo1 → freq_alelo1 = 1 - maf
    # Si major > minor: minor es el alelo1 → freq_alelo1 = maf
    dt[, allele1 := pmin(major, minor)]
    dt[, allele2 := pmax(major, minor)]
    dt[, freq1   := ifelse(major < minor, 1 - maf, maf)]

    # Filtro MAF (sobre la frecuencia estandarizada: usar el menor de freq1 y 1-freq1)
    dt[, maf_std := pmin(freq1, 1 - freq1)]
    dt <- dt[maf_std >= min_maf]

    dt[, .(chr, pos, allele1, allele2, freq1, maf_std)]
}

samples <- c("POOL1", "POOL2", "HOSPITAL")
cat("Cargando archivos MAF (frecuencias estandarizadas)...\n")
maf_list <- lapply(samples, load_maf_std, dir = MAFS_DIR, min_maf = MIN_MAF)
names(maf_list) <- samples
ok <- !sapply(maf_list, is.null)
cat(sprintf("Pools disponibles: %s\n", paste(names(ok)[ok], collapse = ", ")))

# ── FST Hudson-Bhatia pairwise ─────────────────────────────────────────────────
hudson_fst_maf <- function(dt_a, dt_b, sa, sb) {
    # Solo sitios con mismo par alélico (allele1 + allele2 iguales)
    m <- merge(
        dt_a[, .(chr, pos, allele1, allele2, p1 = freq1)],
        dt_b[, .(chr, pos, allele1, allele2, p2 = freq1)],
        by = c("chr", "pos", "allele1", "allele2")
    )
    cat(sprintf("  %s vs %s: %d sitios compartidos con mismo par alélico\n",
                sa, sb, nrow(m)))
    if (nrow(m) < 100) { warning("Muy pocos sitios para FST"); return(NULL) }

    p1 <- m$p1; p2 <- m$p2

    # FST de Wright/Nei (estimador de varianza, válido para pool-seq sin conocer N):
    #   Fst = Var(p) / (p_bar * q_bar)
    # donde Var(p) = (p1-p2)^2 / 4, p_bar = (p1+p2)/2
    p_bar <- (p1 + p2) / 2
    q_bar <- 1 - p_bar
    fst_s <- (p1 - p2)^2 / (4 * p_bar * q_bar + 1e-12)

    data.frame(
        pop1        = sa, pop2 = sb,
        n_sites     = nrow(m),
        fst_mean    = round(mean(fst_s,   na.rm = TRUE), 6),
        fst_median  = round(median(fst_s, na.rm = TRUE), 6),
        fst_sd      = round(sd(fst_s,     na.rm = TRUE), 6),
        stringsAsFactors = FALSE
    )
}

pairs <- combn(names(ok)[ok], 2, simplify = FALSE)
results <- lapply(pairs, function(p) {
    tryCatch(
        hudson_fst_maf(maf_list[[p[1]]], maf_list[[p[2]]], p[1], p[2]),
        error = function(e) { warning(e$message); NULL }
    )
})
fst_table <- do.call(rbind, Filter(Negate(is.null), results))

# ── Output ────────────────────────────────────────────────────────────────────
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
write.table(fst_table, OUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== FST Hudson-Bhatia (genoma completo, MAF estandarizado) ===\n")
print(fst_table, row.names = FALSE)
cat("\nInterpretación FST:\n")
cat("  < 0.001 : diferenciación despreciable (misma población)\n")
cat("  0.001-0.01 : diferenciación muy baja\n")
cat("  0.01-0.05 : diferenciación baja-moderada\n")
cat("  > 0.05  : diferenciación moderada-alta\n\n")
cat("Tabla guardada en:", OUT_FILE, "\n")
