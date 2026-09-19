#!/usr/bin/env Rscript
# 09b_plot_admixture_pca_v2.R — NGSadmix (barplots K=2..11) and PCA (PCAngsd)
# plots for the granular dataset of 11 reference populations
# (IBS, CEU, FIN, YRI, ESN, GWD, LWK, PEL, MXL, PUR, CLM).
#
# METHODOLOGICAL WARNING: the proportions/positions plotted for
# POOL1/POOL2/HOSPITAL come from treating each pool as ONE diploid
# pseudo-individual (see bin/08b_merge_beagle_ngsadmix_v2.R). They should be
# read as the pool's population average, not individual ancestry -- worth
# noting in the figure caption. See bin/README_VIGENTE.md.
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

IND_FILE <- "out_global/angsd/ngsadmix/individuals_v2.tsv"
OUT_DIR  <- "out_global/angsd/ngsadmix/results_v2"
PCA_DIR  <- "out_global/angsd/pcangsd"
POOL_NAMES <- c("POOL1","POOL2","HOSPITAL")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ind <- fread(IND_FILE)
n_ind <- nrow(ind)

# las filas de pools tienen pop=POOL1/POOL2/HOSPITAL (no "POOL"); usar super_pop=="POOL"
# para identificarlas, y sobreescribir pop a "POOL" para que agrupen/coloreen juntas.
ind[, pop := ifelse(super_pop == "POOL", "POOL", pop)]

# orden de poblaciones agrupado por continente, pools al final
pop_levels <- c("YRI","ESN","GWD","LWK",   # AFR
                 "CEU","IBS","FIN",        # EUR
                 "PEL","MXL","PUR","CLM",  # AMR
                 "POOL")

pop_colors <- c(
    YRI="#7f2704", ESN="#a63603", GWD="#d94801", LWK="#f16913",   # AFR: rojos/naranjas
    CEU="#08306b", IBS="#2171b5", FIN="#6baed6",                  # EUR: azules
    PEL="#4a1486", MXL="#88419d", PUR="#8c6bb1", CLM="#c994c7",   # AMR: morados/rosas
    POOL="black"
)

ind[, pop := factor(pop, levels = pop_levels)]

# ── 1. Admixture barplots K=2..11 ────────────────────────────────────────────
for (K in 2:11) {
    qfile <- file.path(OUT_DIR, sprintf("K%d.qopt", K))
    if (!file.exists(qfile)) { cat(sprintf("Saltando K=%d (no existe %s)\n", K, qfile)); next }
    q <- fread(qfile, header = FALSE)
    stopifnot(nrow(q) == n_ind)
    setnames(q, paste0("Cluster", seq_len(K)))
    q[, sample := ind$sample]
    q[, pop := ind$pop]

    setorder(q, pop, sample)
    q[, sample := factor(sample, levels = sample)]

    qlong <- melt(q, id.vars = c("sample","pop"),
                  measure.vars = paste0("Cluster", seq_len(K)),
                  variable.name = "cluster", value.name = "prop")

    cluster_colors <- scales::hue_pal()(K)

    # nombrar cada cluster según la población de referencia donde domina
    # (NGSadmix no etiqueta los clusters; el orden Cluster1..K es arbitrario)
    cluster_cols <- paste0("Cluster", seq_len(K))
    purity <- q[pop != "POOL", lapply(.SD, mean), by = pop, .SDcols = cluster_cols]
    best_cluster <- cluster_cols[apply(purity[, ..cluster_cols], 1, which.max)]
    names(best_cluster) <- as.character(purity$pop)
    cluster_labels <- setNames(
        sapply(cluster_cols, function(cl) {
            pops <- names(best_cluster)[best_cluster == cl]
            if (length(pops) == 0) sprintf("Otro (%s)", cl) else paste(pops, collapse = "+")
        }),
        cluster_cols
    )
    qlong[, cluster := factor(cluster, levels = cluster_cols)]

    p <- ggplot(qlong, aes(x = sample, y = prop, fill = cluster)) +
        geom_col(width = 1, color = NA) +
        scale_fill_manual(values = cluster_colors, labels = cluster_labels[cluster_cols], name = "Cluster ancestral") +
        labs(title = sprintf("NGSadmix K=%d — 3 pools + 110 individuos de referencia (11 poblaciones 1KGP)", K),
             subtitle = "Agrupado por población: YRI,ESN,GWD,LWK (AFR) | CEU,IBS,FIN (EUR) | PEL,MXL,PUR,CLM (AMR) | pools al final",
             x = NULL, y = "Proporción de ancestría") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.title = element_text(face = "bold"))

    breaks_idx <- cumsum(table(q$pop)[pop_levels])
    breaks_idx <- breaks_idx[breaks_idx > 0 & breaks_idx < n_ind]
    p <- p + geom_vline(xintercept = breaks_idx + 0.5, color = "white", linewidth = 0.8)

    ggsave(file.path(OUT_DIR, sprintf("admixture_v2_K%d.png", K)), p, width = 14, height = 4.5, dpi = 180)
    cat(sprintf("Figura: admixture_v2_K%d.png\n", K))

    pool_q <- q[pop == "POOL"]
    cat(sprintf("\n--- NGSadmix K=%d: proporciones de los pools ---\n", K))
    print(pool_q[, c("sample", paste0("Cluster", seq_len(K))), with = FALSE])

    # pureza media por población de referencia (control de calidad)
    cat(sprintf("--- NGSadmix K=%d: pureza media por población de referencia ---\n", K))
    print(q[pop != "POOL", lapply(.SD, mean), by = pop, .SDcols = paste0("Cluster", seq_len(K))])
}

# ── 2. PCA desde matriz de covarianza de PCAngsd ─────────────────────────────
cov_mat <- as.matrix(fread(file.path(PCA_DIR, "pca_combined_v2.cov"), header = FALSE))
stopifnot(nrow(cov_mat) == n_ind)
eig <- eigen(cov_mat)
pcs <- eig$vectors
varexp <- eig$values / sum(eig$values) * 100

pca_dt <- data.table(sample = ind$sample, pop = ind$pop)
for (i in 1:5) pca_dt[[paste0("PC", i)]] <- pcs[, i]

pca_dt[, is_pool := pop == "POOL"]

plot_pc <- function(pcx, pcy) {
    ggplot(pca_dt, aes(x = .data[[pcx]], y = .data[[pcy]], color = pop, size = is_pool, shape = is_pool)) +
        geom_point(alpha = 0.85) +
        ggrepel::geom_text_repel(data = pca_dt[is_pool == TRUE],
                                  aes(label = sample), color = "black", size = 3.5,
                                  fontface = "bold", show.legend = FALSE) +
        scale_color_manual(values = pop_colors, name = "Población") +
        scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2), guide = "none") +
        scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), guide = "none") +
        labs(title = sprintf("PCAngsd (v2, granular): %s vs %s", pcx, pcy),
             subtitle = "Triángulos negros = pools; puntos = individuos de referencia 1KGP (11 poblaciones)",
             x = sprintf("%s (%.1f%% var.)", pcx, varexp[as.integer(sub("PC","",pcx))]),
             y = sprintf("%s (%.1f%% var.)", pcy, varexp[as.integer(sub("PC","",pcy))])) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold"))
}

p1 <- plot_pc("PC1","PC2")
ggsave(file.path(PCA_DIR, "pca_v2_PC1_PC2.png"), p1, width = 9.5, height = 6.5, dpi = 180)
p2 <- plot_pc("PC1","PC3")
ggsave(file.path(PCA_DIR, "pca_v2_PC1_PC3.png"), p2, width = 9.5, height = 6.5, dpi = 180)
cat("Figuras: pca_v2_PC1_PC2.png, pca_v2_PC1_PC3.png\n")

fwrite(pca_dt, file.path(PCA_DIR, "pca_v2_coordinates.tsv"), sep = "\t")
cat(sprintf("\nVarianza explicada: PC1=%.1f%% PC2=%.1f%% PC3=%.1f%% PC4=%.1f%% PC5=%.1f%%\n",
            varexp[1], varexp[2], varexp[3], varexp[4], varexp[5]))

cat("\n--- Coordenadas PCA de los pools ---\n")
print(pca_dt[is_pool == TRUE, .(sample, PC1, PC2, PC3)])
