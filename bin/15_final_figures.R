#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# Paleta C (ColorBrewer Set1)
PAL5 <- c(EUR="#377EB8", AMR="#E41A1C", AFR="#4DAF4A", SAS="#FF7F00", EAS="#984EA3")
PAL5_ES <- c("Europeo"="#377EB8", "Amerindio"="#E41A1C", "Africano"="#4DAF4A",
             "Asia del Sur"="#FF7F00", "Asia Oriental"="#984EA3")
anc_labels <- c(EUR="Europeo", AMR="Amerindio", AFR="Africano",
                SAS="Asia del Sur", EAS="Asia Oriental")

OUTDIR <- "out_global/angsd/ngsadmix"

# ── Helper: identify clusters ──
id_clusters <- function(q, ind, K, by_col="super_pop") {
    q_tmp <- copy(q)
    q_tmp[, .by_var := ind[[by_col]]]
    q_tmp[, .sp := ind$super_pop]
    ref <- q_tmp[.sp != "POOL"]
    cols <- paste0("C", 1:K)
    means <- ref[, lapply(.SD, mean), by=.by_var, .SDcols=cols]
    setnames(means, ".by_var", "grp")
    cn <- character(K); used <- character(0)
    for (k in 1:K) {
        b <- means$grp[which.max(means[[paste0("C",k)]])]
        if (b %in% used) b <- paste0(b,"_",k)
        cn[k] <- b; used <- c(used, b)
    }
    cn
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. COMPOSICIÓN POOLS K=5 (barras horizontales, sin leyenda interna)
# ══════════════════════════════════════════════════════════════════════════════
ind <- fread(paste0(OUTDIR, "/individuals_full.tsv"))
q5 <- fread(paste0(OUTDIR, "/results_full/K5.qopt"), header=FALSE)
setnames(q5, paste0("C",1:5))
cn5 <- id_clusters(q5, ind, 5)
setnames(q5, paste0("C",1:5), cn5)
q5[, sample := ind$sample]; q5[, super_pop := ind$super_pop]

pools <- q5[super_pop == "POOL"]
pl <- melt(pools, id.vars=c("sample","super_pop"), measure.vars=cn5,
           variable.name="A", value.name="P")
pl[, Pct := round(P*100, 1)]
pl[, Lbl := anc_labels[as.character(A)]]
pl[, Lbl := factor(Lbl, levels=c("Europeo","Amerindio","Africano","Asia del Sur","Asia Oriental"))]
pl[, sample := factor(sample, levels=c("HOSPITAL","POOL2","POOL1"))]

p1 <- ggplot(pl, aes(x=sample, y=Pct, fill=Lbl)) +
    geom_col(width=0.6, color="white", linewidth=0.3) +
    geom_text(aes(label=ifelse(Pct>=3, paste0(Pct,"%"),"")),
              position=position_stack(vjust=0.5), size=5, fontface="bold", color="white") +
    coord_flip() +
    scale_fill_manual(values=PAL5_ES) +
    scale_y_continuous(expand=c(0,0), labels=function(x) paste0(x,"%")) +
    labs(x=NULL, y=NULL) +
    theme_minimal(base_size=15) +
    theme(legend.position="none",
          panel.grid.major.y=element_blank(), panel.grid.minor=element_blank())
ggsave(paste0(OUTDIR, "/composicion_pools_K5.png"), p1, width=12, height=4.5, dpi=150)
cat("OK: composicion_pools_K5\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. COMPARACIÓN POOLS K=5 (barras agrupadas, sin leyenda)
# ══════════════════════════════════════════════════════════════════════════════
pl2 <- copy(pl)[Pct >= 0.5]
pl2[, sample := factor(sample, levels=c("POOL1","POOL2","HOSPITAL"))]

p2 <- ggplot(pl2, aes(x=Lbl, y=Pct, fill=sample)) +
    geom_col(position=position_dodge(0.8), width=0.7) +
    geom_text(aes(label=paste0(Pct,"%")), position=position_dodge(0.8),
              vjust=-0.5, size=3.5, fontface="bold") +
    scale_fill_manual(values=c(POOL1="#377EB8", POOL2="#984EA3", HOSPITAL="#E41A1C")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.15)), labels=function(x) paste0(x,"%")) +
    labs(x=NULL, y=NULL) +
    theme_minimal(base_size=15) +
    theme(legend.position="none",
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank())
ggsave(paste0(OUTDIR, "/comparacion_pools_K5.png"), p2, width=12, height=6, dpi=150)
cat("OK: comparacion_pools_K5\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. POOLS INDIVIDUALES K=5 (faceted, sin leyenda)
# ══════════════════════════════════════════════════════════════════════════════
pl3 <- copy(pl)[Pct >= 0.5]
pl3[, sample := factor(sample, levels=c("POOL1","POOL2","HOSPITAL"))]

p3 <- ggplot(pl3, aes(x=Lbl, y=Pct, fill=Lbl)) +
    geom_col(width=0.7, show.legend=FALSE) +
    geom_text(aes(label=paste0(Pct,"%")), vjust=-0.5, size=4.5, fontface="bold") +
    facet_wrap(~sample, ncol=3) +
    scale_fill_manual(values=PAL5_ES) +
    scale_y_continuous(expand=expansion(mult=c(0,0.2)), labels=function(x) paste0(x,"%")) +
    labs(x=NULL, y=NULL) +
    theme_minimal(base_size=14) +
    theme(strip.text=element_text(face="bold", size=15),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1))
ggsave(paste0(OUTDIR, "/pools_individuales_K5.png"), p3, width=14, height=6, dpi=150)
cat("OK: pools_individuales_K5\n")

# ══════════════════════════════════════════════════════════════════════════════
# 4. COMPOSICIÓN SUB-CONTINENTAL K=6 (sin leyenda)
# ══════════════════════════════════════════════════════════════════════════════
ind2 <- fread(paste0(OUTDIR, "/individuals_v2_full.tsv"))
q6 <- fread(paste0(OUTDIR, "/results_v2_full/K6.qopt"), header=FALSE)
setnames(q6, paste0("C",1:6))
q6[, sample := ind2$sample]; q6[, super_pop := ind2$super_pop]; q6[, pop := ind2$pop]
cn6 <- id_clusters(q6, ind2, 6, "pop")
setnames(q6, paste0("C",1:6), cn6)

pools6 <- q6[super_pop == "POOL"]
p6l <- melt(pools6, id.vars=c("sample","super_pop","pop"), measure.vars=cn6,
            variable.name="Cl", value.name="P")
p6l[, Pct := round(P*100, 1)]
p6l[, sample := factor(sample, levels=c("HOSPITAL","POOL2","POOL1"))]

region_map <- c(FIN="Eur. Nórdica", ESN="Áfr. Occidental",
                IBS="Eur. Ibérica", PUR="Am. Caribe",
                PEL="Am. Andina", GWD="Áfr. Occ. (Gambia)")
p6l[, Reg := sapply(as.character(Cl), function(x) if(x %in% names(region_map)) region_map[x] else x)]
p6l[, Reg := factor(Reg, levels=c("Eur. Ibérica","Eur. Nórdica","Am. Andina","Am. Caribe","Áfr. Occidental","Áfr. Occ. (Gambia)"))]

reg_colors <- c("Eur. Ibérica"="#377EB8", "Eur. Nórdica"="#6BAED6",
                "Am. Andina"="#E41A1C", "Am. Caribe"="#FB9A99",
                "Áfr. Occidental"="#4DAF4A", "Áfr. Occ. (Gambia)"="#B2DF8A")

p4 <- ggplot(p6l, aes(x=sample, y=Pct, fill=Reg)) +
    geom_col(width=0.6, color="white", linewidth=0.3) +
    geom_text(aes(label=ifelse(Pct>=4, paste0(Pct,"%"),"")),
              position=position_stack(vjust=0.5), size=4.2, fontface="bold", color="white") +
    coord_flip() +
    scale_fill_manual(values=reg_colors) +
    scale_y_continuous(expand=c(0,0), labels=function(x) paste0(x,"%")) +
    labs(x=NULL, y=NULL) +
    theme_minimal(base_size=15) +
    theme(legend.position="none",
          panel.grid.major.y=element_blank(), panel.grid.minor=element_blank())
ggsave(paste0(OUTDIR, "/composicion_pools_K6.png"), p4, width=12, height=4.5, dpi=150)
cat("OK: composicion_pools_K6\n")

# ══════════════════════════════════════════════════════════════════════════════
# 5. ADMIXTURE BARPLOTS (sin leyenda, etiquetas de población)
# ══════════════════════════════════════════════════════════════════════════════
plot_admix_nolabel <- function(qfile, ind_file, K, outfile, is_granular=FALSE) {
    ind_d <- fread(ind_file)
    q <- fread(qfile, header=FALSE)
    setnames(q, paste0("C",1:K))
    q[, sample := ind_d$sample]; q[, super_pop := ind_d$super_pop]; q[, pop := ind_d$pop]

    by_col <- ifelse(is_granular, "pop", "super_pop")
    cn <- id_clusters(q, ind_d, K, by_col)
    setnames(q, paste0("C",1:K), cn)

    if (is_granular) {
        afr <- c("YRI","ESN","GWD","LWK")
        eur <- c("CEU","IBS","FIN")
        amr <- c("PEL","MXL","PUR","CLM")
        pord <- c(intersect(afr, unique(ind_d$pop)), intersect(eur, unique(ind_d$pop)), intersect(amr, unique(ind_d$pop)))
        ref_q <- q[super_pop != "POOL"]; ref_q[, pop := factor(pop, levels=pord)]
        setorder(ref_q, pop, sample)
    } else {
        ref_q <- q[super_pop != "POOL"]
        ref_q[, super_pop := factor(super_pop, levels=c("AFR","EUR","EAS","AMR","SAS"))]
        setorder(ref_q, super_pop, pop, sample)
    }
    pool_q <- q[super_pop == "POOL"]
    pool_q[, sample := factor(sample, levels=c("POOL1","POOL2","HOSPITAL"))]
    setorder(pool_q, sample)

    ref_q[, x_pos := .I]
    pool_q[, x_pos := nrow(ref_q) + 3 + .I]
    combined <- rbind(ref_q, pool_q, fill=TRUE)

    qlong <- melt(combined, id.vars=c("sample","super_pop","pop","x_pos"),
                  measure.vars=cn, variable.name="cluster", value.name="prop")

    grp_col <- ifelse(is_granular, "pop", "super_pop")
    ref_only <- combined[super_pop != "POOL"]
    grp_mids <- ref_only[, .(mid=mean(x_pos)), by=get(grp_col)]
    setnames(grp_mids, "get", "g")
    pool_only <- combined[super_pop == "POOL"]
    pool_mids <- data.table(g=pool_only$sample, mid=pool_only$x_pos)
    all_mids <- rbind(grp_mids, pool_mids)
    sep_x <- nrow(ref_q) + 1.5

    p <- ggplot(qlong, aes(x=x_pos, y=prop, fill=cluster)) +
        geom_col(width=0.9) +
        geom_vline(xintercept=sep_x, linetype="dashed", color="grey30", linewidth=0.4) +
        scale_x_continuous(breaks=all_mids$mid, labels=all_mids$g) +
        scale_y_continuous(expand=c(0,0)) +
        labs(x=NULL, y=NULL) +
        theme_minimal(base_size=11) +
        theme(legend.position="none",
              axis.text.x=element_text(angle=45, hjust=1, size=9, face="bold"),
              panel.grid.major.x=element_blank(), panel.grid.minor=element_blank())
    ggsave(outfile, p, width=16, height=4.5, dpi=150)
    cat(sprintf("OK: %s\n", basename(outfile)))
}

# Continental K=3, K=5
for (K in c(3, 5)) {
    plot_admix_nolabel(
        sprintf("%s/results_full/K%d.qopt", OUTDIR, K),
        paste0(OUTDIR, "/individuals_full.tsv"), K,
        sprintf("%s/results_full/admixture_full_K%d.png", OUTDIR, K), FALSE)
}

# Granular K=3, K=6
for (K in c(3, 6)) {
    plot_admix_nolabel(
        sprintf("%s/results_v2_full/K%d.qopt", OUTDIR, K),
        paste0(OUTDIR, "/individuals_v2_full.tsv"), K,
        sprintf("%s/results_v2_full/admixture_v2_full_K%d.png", OUTDIR, K), TRUE)
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. EVANNO - CORREGIDO (escala log Y, puntos grandes, K distinguibles)
# ══════════════════════════════════════════════════════════════════════════════
ev <- fread(paste0(OUTDIR, "/evanno_v2_full_table.tsv"))
ev_plot <- ev[!is.na(deltaK) & is.finite(deltaK)]

# DeltaK con escala log
p_dk <- ggplot(ev_plot, aes(x=K, y=deltaK)) +
    geom_line(linewidth=1, color="#377EB8") +
    geom_point(size=4, color="#377EB8") +
    geom_point(data=ev_plot[which.max(deltaK)], size=7, color="#E41A1C", shape=18) +
    geom_label(data=ev_plot[which.max(deltaK)],
               aes(label=paste0("K=", K, "\nΔK=", format(round(deltaK), big.mark=","))),
               vjust=-0.5, size=4, fontface="bold", fill="white", label.size=0.3) +
    scale_x_continuous(breaks=ev_plot$K) +
    scale_y_log10(labels=scales::comma) +
    labs(x="K (número de poblaciones ancestrales)", y="ΔK (escala logarítmica)") +
    theme_minimal(base_size=14) +
    theme(panel.grid.minor=element_blank())
ggsave(paste0(OUTDIR, "/evanno_v2_full_deltaK.png"), p_dk, width=8, height=5.5, dpi=150)
cat("OK: evanno_deltaK\n")

# Mean L(K) con puntos K etiquetados
p_lk <- ggplot(ev, aes(x=K, y=mean_LK)) +
    geom_ribbon(aes(ymin=mean_LK-sd_LK, ymax=mean_LK+sd_LK), alpha=0.15, fill="#377EB8") +
    geom_line(linewidth=1, color="#377EB8") +
    geom_point(size=4, color="#377EB8") +
    geom_text(aes(label=paste0("K=",K)), vjust=-1.2, size=3.5, fontface="bold") +
    scale_x_continuous(breaks=ev$K) +
    scale_y_continuous(labels=scales::comma) +
    labs(x="K (número de poblaciones ancestrales)", y="Log-verosimilitud media L(K) ± DE") +
    theme_minimal(base_size=14) +
    theme(panel.grid.minor=element_blank())
ggsave(paste0(OUTDIR, "/evanno_v2_full_meanLK.png"), p_lk, width=8, height=5.5, dpi=150)
cat("OK: evanno_meanLK\n")

cat("\n=== TODAS LAS FIGURAS FINALES GENERADAS ===\n")
