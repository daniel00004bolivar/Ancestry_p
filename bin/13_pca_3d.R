#!/usr/bin/env Rscript
# 13_pca_3d.R — PCA 3D (estático PNG + interactivo HTML)
suppressPackageStartupMessages({
    library(data.table)
    library(scatterplot3d)
    library(plotly)
    library(htmlwidgets)
})

args <- commandArgs(trailingOnly = TRUE)
COV_FILE <- ifelse(length(args) >= 1, args[1], "out_global/angsd/pcangsd/pca_combined.cov")
IND_FILE <- ifelse(length(args) >= 2, args[2], "out_global/angsd/ngsadmix/individuals.tsv")
OUT_PREFIX <- ifelse(length(args) >= 3, args[3], "out_global/angsd/pcangsd/pca_3d")
TITLE_SUFFIX <- ifelse(length(args) >= 4, args[4], "")
GROUP_BY <- ifelse(length(args) >= 5, args[5], "auto")  # "superpop", "pop", or "auto"

cat(sprintf("Covariance: %s\nIndividuals: %s\nOutput: %s\n", COV_FILE, IND_FILE, OUT_PREFIX))

cov_mat <- as.matrix(fread(COV_FILE, header = FALSE))
eig <- eigen(cov_mat)
pcs <- eig$vectors[, 1:3]
var_pct <- round(100 * eig$values / sum(eig$values), 1)

ind <- fread(IND_FILE)
ind[, PC1 := pcs[, 1]]
ind[, PC2 := pcs[, 2]]
ind[, PC3 := pcs[, 3]]

is_pool <- ind$super_pop == "POOL"

n_superpops <- length(unique(ind[super_pop != "POOL"]$super_pop))
n_pops <- length(unique(ind[super_pop != "POOL"]$pop))
use_pop <- (GROUP_BY == "pop") || (GROUP_BY == "auto" && n_superpops <= 3)
cat(sprintf("Superpoblaciones: %d | Poblaciones: %d | Agrupando por: %s\n",
            n_superpops, n_pops, ifelse(use_pop, "pop", "super_pop")))

if (use_pop) {
    # Panel granular: colorear por población individual
    ind[, group := ifelse(super_pop == "POOL", sample, pop)]
    pop_colors <- c(
        YRI="#D55E00", ESN="#E69F00", GWD="#F0E442", LWK="#CC79A7",
        MKK="#CC79A7",
        CEU="#0072B2", TSI="#56B4E9", IBS="#0072B2", FIN="#56B4E9",
        CHB="#009E73", JPT="#009E73", CHD="#009E73",
        MXL="#CC79A7", PEL="#882255", PUR="#AA4499", CLM="#AA4499",
        GIH="#F0E442", PJL="#F0E442",
        POOL1="#000000", POOL2="#333333", HOSPITAL="#666666"
    )
} else {
    # Panel continental: colorear por SUPERPOBLACIÓN (5 colores grandes)
    ind[, group := ifelse(super_pop == "POOL", sample, super_pop)]
    pop_colors <- c(AFR="#D55E00", AMR="#CC79A7", EAS="#009E73",
                    EUR="#0072B2", SAS="#F0E442",
                    POOL1="#000000", POOL2="#333333", HOSPITAL="#666666")
}

ind[, color := sapply(group, function(g) ifelse(g %in% names(pop_colors), pop_colors[g], "grey50"))]
ind[, pch_val := ifelse(super_pop == "POOL", 17, 16)]
ind[, size_val := ifelse(super_pop == "POOL", 2.5, 1.2)]

# ── Static 3D (PNG) ──────────────────────────────────────────────────────────
png(paste0(OUT_PREFIX, ".png"), width = 1200, height = 1000, res = 150)
par(mar = c(3, 3, 3, 1))

s3d <- scatterplot3d(
    ind[!is_pool]$PC1, ind[!is_pool]$PC2, ind[!is_pool]$PC3,
    color = ind[!is_pool]$color,
    pch = 16, cex.symbols = 1.0,
    xlab = sprintf("PC1 (%.1f%%)", var_pct[1]),
    ylab = sprintf("PC2 (%.1f%%)", var_pct[2]),
    zlab = sprintf("PC3 (%.1f%%)", var_pct[3]),
    main = paste0("PCA 3D — PCAngsd ", TITLE_SUFFIX),
    angle = 40, box = TRUE, grid = TRUE,
    type = "p"
)

s3d$points3d(
    ind[is_pool]$PC1, ind[is_pool]$PC2, ind[is_pool]$PC3,
    pch = 17, col = "black", cex = 2.5
)

pool_2d <- s3d$xyz.convert(ind[is_pool]$PC1, ind[is_pool]$PC2, ind[is_pool]$PC3)
text(pool_2d$x, pool_2d$y, labels = ind[is_pool]$sample,
     pos = 4, cex = 0.9, font = 2, col = "black")

legend_groups <- unique(ind[!is_pool]$group)
legend_colors <- sapply(legend_groups, function(g) ifelse(g %in% names(pop_colors), pop_colors[g], "grey50"))
legend("topright",
       legend = c(legend_groups, "Pools"),
       col = c(legend_colors, "black"),
       pch = c(rep(16, length(legend_groups)), 17),
       pt.cex = c(rep(1, length(legend_groups)), 2),
       cex = 0.7, bg = "white")

dev.off()
cat(sprintf("PNG guardado: %s.png\n", OUT_PREFIX))

# ── Interactive 3D (plotly HTML) ─────────────────────────────────────────────
ref_data <- ind[!is_pool]
pool_data <- ind[is_pool]

fig <- plot_ly() %>%
    add_trace(
        data = ref_data,
        x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~group,
        colors = pop_colors,
        type = "scatter3d", mode = "markers",
        marker = list(size = 4, opacity = 0.8),
        text = ~paste0(sample, " (", pop, " / ", super_pop, ")"),
        hoverinfo = "text"
    ) %>%
    add_trace(
        data = pool_data,
        x = ~PC1, y = ~PC2, z = ~PC3,
        type = "scatter3d", mode = "markers+text",
        marker = list(size = 10, color = "black", symbol = "diamond",
                      line = list(color = "white", width = 1)),
        text = ~sample, textposition = "top center",
        textfont = list(size = 12, color = "black"),
        hovertext = ~paste0(sample, "\nPC1=", round(PC1, 4),
                           "\nPC2=", round(PC2, 4),
                           "\nPC3=", round(PC3, 4)),
        hoverinfo = "text",
        name = "Pools",
        showlegend = TRUE
    ) %>%
    layout(
        title = list(text = paste0("PCA 3D interactivo — PCAngsd ", TITLE_SUFFIX)),
        scene = list(
            xaxis = list(title = sprintf("PC1 (%.1f%%)", var_pct[1])),
            yaxis = list(title = sprintf("PC2 (%.1f%%)", var_pct[2])),
            zaxis = list(title = sprintf("PC3 (%.1f%%)", var_pct[3])),
            camera = list(eye = list(x = 1.5, y = 1.5, z = 1.0))
        ),
        legend = list(x = 0.85, y = 0.95)
    )

saveWidget(fig, paste0(OUT_PREFIX, ".html"), selfcontained = TRUE)
cat(sprintf("HTML interactivo guardado: %s.html\n", OUT_PREFIX))

# ── Save coordinates ─────────────────────────────────────────────────────────
fwrite(ind[, .(sample, pop, super_pop, PC1, PC2, PC3)],
       paste0(OUT_PREFIX, "_coordinates.tsv"), sep = "\t")
cat(sprintf("Coordenadas guardadas: %s_coordinates.tsv\n", OUT_PREFIX))
