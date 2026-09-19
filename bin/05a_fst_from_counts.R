#!/usr/bin/env Rscript
# 05a_fst_from_counts.R — FST Hudson desde conteos ANGSD (pool-seq)
#
# Carga los archivos .pos.gz + .counts.gz de ANGSD (-dumpCounts 3)
# que contienen conteos (totA, totC, totG, totT) por sitio por pool.
#
# Calcula FST pairwise entre POOL1, POOL2, HOSPITAL.
#
# Estimador de Hudson et al. (1992) / Bhatia et al. (2013):
#   Fst = (p1 - p2)^2 - p1*q1/(n1-1) - p2*q2/(n2-1)
#         ─────────────────────────────────────────────
#         p1*q2 + p2*q1
#
# Para pool-seq, n_i = cobertura total (lecturas); no es el número real
# de individuos, por lo que el sesgo de corrección de muestra es mínimo
# con profundidades altas (≥10x).
#
# Uso desde bash (llamado automáticamente desde 05_angsd_analysis.sh):
#   Rscript bin/05a_fst_from_counts.R \
#       --counts_dir out_global/angsd/per_pool \
#       --out out_global/angsd/fst/hudson_fst_summary.tsv
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table)
})

# ── Args ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) > 0 && length(args) >= i + 1) args[i + 1] else default
}

COUNTS_DIR <- get_arg("--counts_dir", "out_global/angsd/per_pool")
OUT_FILE   <- get_arg("--out", "out_global/angsd/fst/hudson_fst_summary.tsv")
MIN_DEPTH  <- as.integer(get_arg("--min_depth", "10"))
MIN_MAF    <- as.double(get_arg("--min_maf", "0.01"))

if (!dir.exists(COUNTS_DIR)) {
    stop("--counts_dir no existe: ", COUNTS_DIR)
}

# ── Carga conteos ANGSD (dumpCounts 3) ────────────────────────────────────────
# pos.gz: chr, pos, totDepth (3 cols)
# counts.gz: totA, totC, totG, totT (4 cols, mismas filas)
load_acgt <- function(sample_name, dir) {
    pos_f  <- file.path(dir, paste0(sample_name, "_counts.pos.gz"))
    cnt_f  <- file.path(dir, paste0(sample_name, "_counts.counts.gz"))

    if (!file.exists(pos_f) || !file.exists(cnt_f)) {
        warning("Archivos no encontrados para ", sample_name)
        return(NULL)
    }

    pos <- fread(pos_f, col.names = c("chr", "pos", "depth"))
    cnt <- fread(cnt_f, col.names = c("totA", "totC", "totG", "totT"))

    if (nrow(pos) == 0 || nrow(cnt) == 0) {
        warning("Archivos vacíos para ", sample_name)
        return(NULL)
    }

    dt <- cbind(pos, cnt)
    dt[, sample := sample_name]

    # Identificar alelo mayor y menor desde conteos
    dt[, c("allele1", "allele2", "n1", "n2") := {
        mat <- cbind(totA, totC, totG, totT)
        bases <- c("A", "C", "G", "T")
        # Ordenar por conteo descendente, tomar top 2
        lapply(seq_len(.N), function(i) {
            cts <- mat[i, ]
            ord <- order(cts, decreasing = TRUE)
            list(bases[ord[1]], bases[ord[2]], cts[ord[1]], cts[ord[2]])
        }) |> do.call(what = rbind) |> as.data.frame()
    }]

    # Simplificado: evitar el lapply que es lento para datos grandes
    # Usando operaciones vectorizadas de data.table:
    dt[, maxN := pmax(totA, totC, totG, totT)]
    dt[, minN := apply(cbind(totA, totC, totG, totT), 1, function(x) {
        x2 <- sort(x, decreasing = TRUE)
        if (length(x2) >= 2) x2[2] else 0L
    })]
    dt[, p_minor := minN / depth]   # frecuencia del segundo alelo más común

    dt
}

# Carga vectorizada más eficiente
load_acgt_fast <- function(sample_name, dir) {
    pos_f  <- file.path(dir, paste0(sample_name, "_counts.pos.gz"))
    cnt_f  <- file.path(dir, paste0(sample_name, "_counts.counts.gz"))

    if (!file.exists(pos_f) || !file.exists(cnt_f)) {
        warning("Archivos no encontrados para ", sample_name, ": ", pos_f)
        return(NULL)
    }

    pos <- fread(pos_f, col.names = c("chr", "pos", "depth"))
    cnt <- fread(cnt_f, col.names = c("totA", "totC", "totG", "totT"))

    if (nrow(pos) == 0) {
        warning("Sin datos en pos.gz para ", sample_name)
        return(NULL)
    }
    if (nrow(cnt) == 0) {
        warning("Sin datos en counts.gz para ", sample_name)
        return(NULL)
    }
    if (nrow(pos) != nrow(cnt)) {
        stop("Filas discordantes entre pos.gz y counts.gz para ", sample_name)
    }

    dt <- cbind(pos, cnt)
    dt[, sample := sample_name]

    # Frecuencia del alelo más común = major allele freq
    # Frecuencia del segundo más común = minor allele freq (p_minor)
    # Usamos sort por fila: más eficiente con apply sobre la matrix
    mat <- as.matrix(dt[, .(totA, totC, totG, totT)])
    top2 <- t(apply(mat, 1, function(x) sort(x, decreasing = TRUE)[1:2]))
    dt[, major_n := top2[, 1]]
    dt[, minor_n := top2[, 2]]
    dt[, p_minor := minor_n / depth]

    dt
}

samples <- c("POOL1", "POOL2", "HOSPITAL")
cat("Cargando conteos ANGSD...\n")
counts_list <- lapply(samples, load_acgt_fast, dir = COUNTS_DIR)
names(counts_list) <- samples

ok <- !sapply(counts_list, is.null)
if (sum(ok) < 2) {
    stop("Necesito al menos 2 archivos de conteos con datos. ",
         "¿Corriste ANGSD primero (bash bin/05_angsd_analysis.sh)?")
}

cat(sprintf("  Pools disponibles: %s\n", paste(names(ok)[ok], collapse = ", ")))

# ── Función FST Hudson pairwise ───────────────────────────────────────────────
# Estimador de Bhatia et al. 2013 (apropiado para frecuencias de pools):
#   num   = (p1 - p2)^2 - p1*q1/n1 - p2*q2/n2
#   denom = p1*q2 + p2*q1
#   Fst   = num / denom
# donde n1, n2 = profundidad de cobertura (no individuos)
hudson_fst <- function(dt_a, dt_b, sample_a, sample_b,
                       min_depth, min_maf) {
    m <- merge(
        dt_a[depth >= min_depth & p_minor >= min_maf,
             .(chr, pos, p1 = p_minor, n1 = depth)],
        dt_b[depth >= min_depth & p_minor >= min_maf,
             .(chr, pos, p2 = p_minor, n2 = depth)],
        by = c("chr", "pos")
    )

    if (nrow(m) < 100) {
        warning("Muy pocos sitios compartidos (", nrow(m), ") entre ",
                sample_a, " y ", sample_b)
        if (nrow(m) == 0) return(NULL)
    }

    p1 <- m$p1; q1 <- 1 - p1; n1 <- m$n1
    p2 <- m$p2; q2 <- 1 - p2; n2 <- m$n2

    num   <- (p1 - p2)^2 - p1*q1/(n1 - 1) - p2*q2/(n2 - 1)
    denom <- p1*q2 + p2*q1

    valid <- denom > 0
    fst_per_site <- num[valid] / denom[valid]

    data.frame(
        pop1        = sample_a,
        pop2        = sample_b,
        n_sites     = nrow(m),
        fst_mean    = mean(fst_per_site, na.rm = TRUE),
        fst_median  = median(fst_per_site, na.rm = TRUE),
        fst_sd      = sd(fst_per_site, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}

# ── Calcular FST para los 3 pares ─────────────────────────────────────────────
avail <- names(ok)[ok]
pairs <- combn(avail, 2, simplify = FALSE)
results <- vector("list", length(pairs))

for (i in seq_along(pairs)) {
    a <- pairs[[i]][1]; b <- pairs[[i]][2]
    cat(sprintf("[FST] %s vs %s ...\n", a, b))
    results[[i]] <- tryCatch(
        hudson_fst(counts_list[[a]], counts_list[[b]], a, b,
                   min_depth = MIN_DEPTH, min_maf = MIN_MAF),
        error = function(e) {
            warning("Error en FST ", a, " vs ", b, ": ", e$message)
            NULL
        }
    )
}

fst_table <- do.call(rbind, Filter(Negate(is.null), results))

if (is.null(fst_table) || nrow(fst_table) == 0) {
    cat("No se pudo calcular FST. Revisa que los archivos de conteos tengan datos.\n")
    quit(save = "no", status = 0)
}

# ── Output ────────────────────────────────────────────────────────────────────
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
write.table(fst_table, OUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== FST Hudson-Bhatia (pool-seq, ANGSD counts) ===\n")
print(fst_table, digits = 4, row.names = FALSE)
cat("\nTabla guardada en:", OUT_FILE, "\n")
cat("\nInterpretación FST (Wright 1978):\n")
cat("  FST < 0.05 : diferenciación muy baja (misma población)\n")
cat("  FST 0.05-0.15 : diferenciación moderada\n")
cat("  FST 0.15-0.25 : diferenciación grande\n")
cat("  FST > 0.25 : diferenciación muy grande\n\n")
cat("NOTA pool-seq: este FST usa cobertura (lecturas) como n_i.\n")
cat("  Con alta cobertura (≥20x) la corrección de muestra es pequeña.\n")
cat("  La referencia (--ref + realSFS) da estimados más precisos.\n")
