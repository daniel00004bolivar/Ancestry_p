#!/usr/bin/env Rscript
# 20_new_findings_figures.R — figures for this audit's new findings
# (jackknife, chromosome X FST, mtDNA haplogroups), in the same visual style
# as bin/15_final_figures.R (ColorBrewer Set1 palette, theme_minimal).
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

PAL5 <- c(EUR = "#377EB8", AMR = "#E41A1C", AFR = "#4DAF4A", SAS = "#FF7F00", EAS = "#984EA3")
OUTDIR <- "manuscrito_G3/manuscript_en/figures"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
POOL_LEVELS <- c("HOSPITAL", "POOL2", "POOL1")

# ══════════════════════════════════════════════════════════════════════════
# Figure A: mtDNA (maternal lineage) vs autosomal, by pool
# ══════════════════════════════════════════════════════════════════════════
mt <- fread("results/mtdna_continental_summary.tsv")
setnames(mt, c("AFR","EAS","EUR","NAT"), c("AFR","EAS","EUR","AMR"), skip_absent = TRUE)
mt_long <- melt(mt, id.vars = "Pool", variable.name = "Category", value.name = "mtDNA")

auto <- data.table(
    Pool = rep(c("POOL1","POOL2","HOSPITAL"), each = 5),
    Category = rep(c("EUR","AMR","AFR","SAS","EAS"), 3),
    Autosomal = c(0.472,0.311,0.165,0.046,0.007,
                  0.479,0.310,0.165,0.036,0.011,
                  0.366,0.319,0.282,0.023,0.011)
)

comp <- merge(auto, mt_long, by = c("Pool","Category"), all.x = TRUE)
comp[is.na(mtDNA), mtDNA := 0]
comp <- comp[Category %in% c("EUR","AMR","AFR")]
comp_long <- melt(comp, id.vars = c("Pool","Category"), measure.vars = c("Autosomal","mtDNA"),
                   variable.name = "Source", value.name = "Proportion")
comp_long[, Pool := factor(Pool, levels = POOL_LEVELS)]
comp_long[, Category := factor(Category, levels = c("EUR","AMR","AFR"),
                                labels = c("European","Native American","African"))]

pA <- ggplot(comp_long, aes(x = Category, y = Proportion * 100, fill = Category, alpha = Source)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "white", linewidth = 0.3) +
    facet_wrap(~Pool, nrow = 1) +
    scale_fill_manual(values = c("European" = unname(PAL5["EUR"]), "Native American" = unname(PAL5["AMR"]), "African" = unname(PAL5["AFR"])),
                       guide = "none") +
    scale_alpha_manual(values = c(Autosomal = 1, mtDNA = 0.45), name = NULL,
                        labels = c("Autosomal (K=5)", "Maternal (mtDNA)")) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Ancestry proportion") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
ggsave(file.path(OUTDIR, "fig_mtdna_vs_autosomal_en.png"), pA, width = 8.5, height = 4, dpi = 150)
cat("Written: fig_mtdna_vs_autosomal_en.png\n")

# ══════════════════════════════════════════════════════════════════════════
# Figure B: K=5 jackknife uncertainty
# ══════════════════════════════════════════════════════════════════════════
jk <- fread("results/pool_jackknife_summary.tsv")
jk[, Component := factor(Component, levels = c("K1","K2","K3","K4","K5"),
                          labels = c("SAS","EAS","AMR","EUR","AFR"))]
jk[, Pool := factor(Pool, levels = POOL_LEVELS)]

pB <- ggplot(jk, aes(x = Component, y = Full_genome_estimate * 100, color = Component)) +
    geom_errorbar(aes(ymin = pmax(0, (Full_genome_estimate - Jackknife_SE) * 100),
                       ymax = (Full_genome_estimate + Jackknife_SE) * 100),
                   width = 0.25, linewidth = 0.7) +
    geom_point(size = 3) +
    facet_wrap(~Pool, nrow = 1) +
    scale_color_manual(values = PAL5, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "K=5 proportion ± jackknife SE") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(OUTDIR, "fig_jackknife_uncertainty_en.png"), pB, width = 8.5, height = 4, dpi = 150)
cat("Written: fig_jackknife_uncertainty_en.png\n")

# ══════════════════════════════════════════════════════════════════════════
# Figure C: autosomal vs X FST
# ══════════════════════════════════════════════════════════════════════════
fst <- data.table(
    Pair = rep(c("POOL1-POOL2","POOL1-HOSPITAL","POOL2-HOSPITAL"), 2),
    Chromosome = rep(c("Autosomes","Chromosome X"), each = 3),
    FST = c(0.00169, 0.00135, 0.00132,
            -0.0016, 0.0276, 0.0259)
)
fst[, Pair := factor(Pair, levels = c("POOL1-POOL2","POOL1-HOSPITAL","POOL2-HOSPITAL"))]

pC <- ggplot(fst, aes(x = Pair, y = FST, fill = Chromosome)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "white", linewidth = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    scale_fill_manual(values = c("Autosomes" = "grey55", "Chromosome X" = unname(PAL5["AFR"])), name = NULL) +
    labs(x = NULL, y = expression(F[ST])) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(OUTDIR, "fig_fstX_vs_autosomal_en.png"), pC, width = 6, height = 4, dpi = 150)
cat("Written: fig_fstX_vs_autosomal_en.png\n")
