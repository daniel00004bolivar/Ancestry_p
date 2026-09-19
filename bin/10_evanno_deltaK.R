#!/usr/bin/env Rscript
# 10_evanno_deltaK.R — Cálculo de deltaK (Evanno et al. 2005) desde corridas NGSadmix
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Uso: Rscript 10_evanno_deltaK.R <results_dir> <out_prefix>\n",
         "  Ej: Rscript bin/10_evanno_deltaK.R out_global/angsd/ngsadmix/results_evanno_v2 out_global/angsd/ngsadmix/evanno_v2")
}

RESULTS_DIR <- args[1]
OUT_PREFIX  <- args[2]

log_files <- list.files(RESULTS_DIR, pattern = "^K[0-9]+_seed[0-9]+\\.log$", full.names = TRUE)
if (length(log_files) == 0) stop("No se encontraron archivos .log en ", RESULTS_DIR)

parse_log <- function(f) {
    lines <- readLines(f)
    last <- lines[length(lines)]
    m <- regmatches(last, regexec("best like=([^ ]+) after ([0-9]+)", last))[[1]]
    if (length(m) < 3) return(NULL)
    bn <- basename(f)
    km <- regmatches(bn, regexec("K([0-9]+)_seed([0-9]+)", bn))[[1]]
    data.table(K = as.integer(km[2]), seed = as.integer(km[3]),
               loglik = as.numeric(m[2]), iters = as.integer(m[3]))
}

dt <- rbindlist(lapply(log_files, parse_log))
setorder(dt, K, seed)

cat("\n=== Log-likelihoods por corrida ===\n")
print(dt)

# Evanno: mean y sd de L(K) por K
summary_dt <- dt[, .(mean_LK = mean(loglik), sd_LK = sd(loglik), n_runs = .N), by = K]
setorder(summary_dt, K)

# L'(K) = L(K) - L(K-1)
summary_dt[, LprimeK := mean_LK - shift(mean_LK, 1)]
# L''(K) = L'(K+1) - L'(K) = L(K+1) - 2*L(K) + L(K-1)
summary_dt[, LdoubleprimeK := shift(mean_LK, -1) - 2 * mean_LK + shift(mean_LK, 1)]
# deltaK = |L''(K)| / sd(L(K))
summary_dt[, deltaK := abs(LdoubleprimeK) / sd_LK]

cat("\n=== Tabla Evanno ===\n")
print(summary_dt, digits = 2)

tsv_file <- paste0(OUT_PREFIX, "_table.tsv")
fwrite(summary_dt, tsv_file, sep = "\t")
cat("\nTabla guardada en:", tsv_file, "\n")

# deltaK solo se puede calcular para K_min+1 ... K_max-1
evanno_plot <- summary_dt[!is.na(deltaK) & is.finite(deltaK)]

if (nrow(evanno_plot) == 0) {
    cat("AVISO: No hay suficientes valores de K para calcular deltaK (necesita al menos 3 valores de K).\n")
    quit(save = "no")
}

best_K <- evanno_plot[which.max(deltaK), K]
cat("\n>>> K óptimo (máximo deltaK):", best_K, "<<<\n")

# --- Plot 1: deltaK vs K ---
p1 <- ggplot(evanno_plot, aes(x = K, y = deltaK)) +
    geom_line(linewidth = 0.8, color = "#0072B2") +
    geom_point(size = 3, color = "#0072B2") +
    geom_point(data = evanno_plot[K == best_K], size = 5, color = "#D55E00", shape = 18) +
    scale_x_continuous(breaks = evanno_plot$K) +
    labs(title = expression(Delta * K ~ "(Evanno et al. 2005)"),
         subtitle = paste0("K óptimo = ", best_K, "  (directorio: ", basename(RESULTS_DIR), ")"),
         x = "K", y = expression(Delta * K)) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))

ggsave(paste0(OUT_PREFIX, "_deltaK.png"), p1, width = 7, height = 5, dpi = 150)
cat("Gráfico deltaK guardado en:", paste0(OUT_PREFIX, "_deltaK.png"), "\n")

# --- Plot 2: mean L(K) ± sd ---
p2 <- ggplot(summary_dt, aes(x = K, y = mean_LK)) +
    geom_ribbon(aes(ymin = mean_LK - sd_LK, ymax = mean_LK + sd_LK), alpha = 0.2, fill = "#0072B2") +
    geom_line(linewidth = 0.8, color = "#0072B2") +
    geom_point(size = 3, color = "#0072B2") +
    scale_x_continuous(breaks = summary_dt$K) +
    labs(title = "Mean log-likelihood L(K) ± SD",
         subtitle = paste0("directorio: ", basename(RESULTS_DIR)),
         x = "K", y = "Mean ln P(D|K)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))

ggsave(paste0(OUT_PREFIX, "_meanLK.png"), p2, width = 7, height = 5, dpi = 150)
cat("Gráfico mean L(K) guardado en:", paste0(OUT_PREFIX, "_meanLK.png"), "\n")

# --- Plot 3: L'(K) ---
lp_plot <- summary_dt[!is.na(LprimeK)]
p3 <- ggplot(lp_plot, aes(x = K, y = LprimeK)) +
    geom_line(linewidth = 0.8, color = "#009E73") +
    geom_point(size = 3, color = "#009E73") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_continuous(breaks = lp_plot$K) +
    labs(title = "L'(K) = L(K) - L(K-1)",
         x = "K", y = "L'(K)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))

ggsave(paste0(OUT_PREFIX, "_LprimeK.png"), p3, width = 7, height = 5, dpi = 150)
cat("Gráfico L'(K) guardado en:", paste0(OUT_PREFIX, "_LprimeK.png"), "\n")

cat("\n=== Evanno completado ===\n")
