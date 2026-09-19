#!/usr/bin/env Rscript
# 09_plot_admixture_pca.R — NGSadmix (barplots K=2..5) and PCA (PCAngsd) plots
#
# METHODOLOGICAL WARNING: the proportions/positions plotted for
# POOL1/POOL2/HOSPITAL come from treating each pool as ONE diploid
# pseudo-individual (see bin/08_merge_beagle_ngsadmix.R). They should be
# read as the pool's population average, not individual ancestry -- worth
# noting in the figure caption. See bin/README_VIGENTE.md.
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

IND_FILE <- "out_global/angsd/ngsadmix/individuals.tsv"
OUT_DIR  <- "out_global/angsd/ngsadmix/results"
PCA_DIR  <- "out_global/angsd/pcangsd"
POOL_NAMES <- c("POOL1","POOL2","HOSPITAL")

ind <- fread(IND_FILE)
ind[, label := ifelse(super_pop == "POOL", sample, paste0(sample,"_",pop))]
n_ind <- nrow(ind)

cont_colors <- c(AFR="#D55E00", AMR="#CC79A7", EAS="#009E73",
                  EUR="#0072B2", SAS="#F0E442", POOL="black")

# ── 1. Admixture barplots K=2..5 ─────────────────────────────────────────────
for (K in 2:5) {
    q <- fread(file.path(OUT_DIR, sprintf("K%d.qopt", K)), header = FALSE)
    stopifnot(nrow(q) == n_ind)
    setnames(q, paste0("Cluster", seq_len(K)))
    q[, sample := ind$sample]
    q[, super_pop := ind$super_pop]
    q[, pop := ind$pop]

    # orden: agrupar por super_pop (AFR,AMR,EAS,EUR,SAS) y pools al final, destacados
    ord_levels <- c("AFR","EUR","EAS","AMR","SAS","POOL")
    q[, super_pop := factor(super_pop, levels = ord_levels)]
    setorder(q, super_pop, pop, sample)
    q[, sample := factor(sample, levels = sample)]

    qlong <- melt(q, id.vars = c("sample","super_pop","pop"),
                  measure.vars = paste0("Cluster", seq_len(K)),
                  variable.name = "cluster", value.name = "prop")

    cluster_colors <- scales::hue_pal()(K)

    # nombrar cada cluster según la superpoblación de referencia donde domina
    # (NGSadmix no etiqueta los clusters; el orden Cluster1..K es arbitrario)
    cluster_cols <- paste0("Cluster", seq_len(K))
    purity <- q[super_pop != "POOL", lapply(.SD, mean), by = super_pop, .SDcols = cluster_cols]
    best_cluster <- cluster_cols[apply(purity[, ..cluster_cols], 1, which.max)]
    names(best_cluster) <- as.character(purity$super_pop)
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
        labs(title = sprintf("NGSadmix K=%d — 3 pools + 100 individuos de referencia 1KGP", K),
             subtitle = "Agrupado por superpoblación (AFR, EUR, EAS, AMR, SAS) y pools al final",
             x = NULL, y = "Proporción de ancestría") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.title = element_text(face = "bold"))

    # líneas separadoras entre superpoblaciones
    breaks_idx <- cumsum(table(q$super_pop)[ord_levels])
    breaks_idx <- breaks_idx[breaks_idx > 0 & breaks_idx < n_ind]
    p <- p + geom_vline(xintercept = breaks_idx + 0.5, color = "white", linewidth = 0.8)

    ggsave(file.path(OUT_DIR, sprintf("admixture_K%d.png", K)), p, width = 13, height = 4.5, dpi = 180)
    cat(sprintf("Figura: admixture_K%d.png\n", K))

    # tabla específica de los 3 pools
    pool_q <- q[super_pop == "POOL"]
    cat(sprintf("\n--- NGSadmix K=%d: proporciones de los pools ---\n", K))
    print(pool_q[, c("sample", paste0("Cluster", seq_len(K))), with = FALSE])
}

# ── 2. PCA desde matriz de covarianza de PCAngsd ─────────────────────────────
cov_mat <- as.matrix(fread(file.path(PCA_DIR, "pca_combined.cov"), header = FALSE))
stopifnot(nrow(cov_mat) == n_ind)
eig <- eigen(cov_mat)
pcs <- eig$vectors
varexp <- eig$values / sum(eig$values) * 100

pca_dt <- data.table(sample = ind$sample, super_pop = ind$super_pop, pop = ind$pop)
for (i in 1:5) pca_dt[[paste0("PC", i)]] <- pcs[, i]

pca_dt[, super_pop := factor(super_pop, levels = c("AFR","EUR","EAS","AMR","SAS","POOL"))]
pca_dt[, is_pool := super_pop == "POOL"]

plot_pc <- function(pcx, pcy) {
    ggplot(pca_dt, aes(x = .data[[pcx]], y = .data[[pcy]], color = super_pop, size = is_pool, shape = is_pool)) +
        geom_point(alpha = 0.85) +
        ggrepel::geom_text_repel(data = pca_dt[is_pool == TRUE],
                                  aes(label = sample), color = "black", size = 3.5,
                                  fontface = "bold", show.legend = FALSE) +
        scale_color_manual(values = cont_colors, name = "Grupo") +
        scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2), guide = "none") +
        scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), guide = "none") +
        labs(title = sprintf("PCAngsd: %s vs %s", pcx, pcy),
             subtitle = "Triángulos negros = pools; puntos = individuos de referencia 1KGP",
             x = sprintf("%s (%.1f%% var.)", pcx, varexp[as.integer(sub("PC","",pcx))]),
             y = sprintf("%s (%.1f%% var.)", pcy, varexp[as.integer(sub("PC","",pcy))])) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold"))
}

p1 <- plot_pc("PC1","PC2")
ggsave(file.path(PCA_DIR, "pca_PC1_PC2.png"), p1, width = 8.5, height = 6.5, dpi = 180)
p2 <- plot_pc("PC1","PC3")
ggsave(file.path(PCA_DIR, "pca_PC1_PC3.png"), p2, width = 8.5, height = 6.5, dpi = 180)
cat("Figuras: pca_PC1_PC2.png, pca_PC1_PC3.png\n")

fwrite(pca_dt, file.path(PCA_DIR, "pca_coordinates.tsv"), sep = "\t")
cat(sprintf("\nVarianza explicada: PC1=%.1f%% PC2=%.1f%% PC3=%.1f%% PC4=%.1f%% PC5=%.1f%%\n",
            varexp[1], varexp[2], varexp[3], varexp[4], varexp[5]))

cat("\n--- Coordenadas PCA de los pools ---\n")
print(pca_dt[is_pool == TRUE, .(sample, PC1, PC2, PC3)])
