#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

ind <- fread("out_global/angsd/ngsadmix/individuals_full.tsv")
q5 <- fread("out_global/angsd/ngsadmix/results_full/K5.qopt", header=FALSE)
setnames(q5, paste0("C", 1:5))
q5[, sample := ind$sample]
q5[, super_pop := ind$super_pop]

ref <- q5[super_pop != "POOL"]
cols <- paste0("C", 1:5)
means <- ref[, lapply(.SD, mean), by=super_pop, .SDcols=cols]
cnames <- character(5)
for (k in 1:5) cnames[k] <- means$super_pop[which.max(means[[paste0("C",k)]])]
setnames(q5, cols, cnames)

pools <- q5[super_pop == "POOL"]
anc_labels <- c(EUR="Europeo", AMR="Amerindio", AFR="Africano", SAS="Asia del Sur", EAS="Asia Oriental")

pools_long <- melt(pools, id.vars=c("sample","super_pop"),
                   measure.vars=cnames, variable.name="Anc", value.name="Prop")
pools_long[, Pct := round(Prop * 100, 1)]
pools_long[, Label := anc_labels[as.character(Anc)]]
pools_long[, Label := factor(Label, levels=c("Europeo","Amerindio","Africano","Asia del Sur","Asia Oriental"))]

# ── 1. Barras horizontales apiladas ──
pools_long[, sample := factor(sample, levels=c("HOSPITAL","POOL2","POOL1"))]
p1 <- ggplot(pools_long, aes(x=sample, y=Pct, fill=Label)) +
    geom_col(width=0.6, color="white", linewidth=0.3) +
    geom_text(aes(label=ifelse(Pct >= 3, paste0(Pct,"%"), "")),
              position=position_stack(vjust=0.5), size=5, fontface="bold", color="white") +
    coord_flip() +
    scale_fill_manual(values=c("Europeo"="#4472C4","Amerindio"="#ED7D31",
                               "Africano"="#70AD47","Asia del Sur"="#FFC000","Asia Oriental"="#5B9BD5")) +
    scale_y_continuous(expand=c(0,0), labels=function(x) paste0(x,"%")) +
    labs(title="Composición ancestral de cada cohorte",
         subtitle="NGSadmix K=5 — 239,000 SNPs, 22 cromosomas", x=NULL, y=NULL, fill="Ancestría") +
    theme_minimal(base_size=15) +
    theme(plot.title=element_text(face="bold",size=18), panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(), legend.position="bottom")
ggsave("out_global/angsd/ngsadmix/composicion_pools_K5.png", p1, width=12, height=5, dpi=150)
cat("OK: composicion_pools_K5.png\n")

# ── 2. Barras agrupadas comparativas ──
pools_long2 <- copy(pools_long)[Pct >= 0.5]
pools_long2[, sample := factor(sample, levels=c("POOL1","POOL2","HOSPITAL"))]
p2 <- ggplot(pools_long2, aes(x=Label, y=Pct, fill=sample)) +
    geom_col(position=position_dodge(0.8), width=0.7) +
    geom_text(aes(label=paste0(Pct,"%")), position=position_dodge(0.8), vjust=-0.5, size=4, fontface="bold") +
    scale_fill_manual(values=c(POOL1="#4472C4", POOL2="#5B9BD5", HOSPITAL="#ED7D31")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.15)), labels=function(x) paste0(x,"%")) +
    labs(title="Comparación de ancestría entre las tres cohortes",
         subtitle="NGSadmix K=5 — 239,000 SNPs, 22 cromosomas", x=NULL, y=NULL, fill="Cohorte") +
    theme_minimal(base_size=15) +
    theme(plot.title=element_text(face="bold",size=18), panel.grid.major.x=element_blank(),
          panel.grid.minor=element_blank(), legend.position="bottom")
ggsave("out_global/angsd/ngsadmix/comparacion_pools_K5.png", p2, width=12, height=6, dpi=150)
cat("OK: comparacion_pools_K5.png\n")

# ── 3. Paneles individuales por pool ──
pools_long3 <- copy(pools_long)[Pct >= 0.5]
pools_long3[, sample := factor(sample, levels=c("POOL1","POOL2","HOSPITAL"))]
p3 <- ggplot(pools_long3, aes(x=Label, y=Pct, fill=Label)) +
    geom_col(width=0.7, show.legend=FALSE) +
    geom_text(aes(label=paste0(Pct,"%")), vjust=-0.5, size=4.5, fontface="bold") +
    facet_wrap(~sample, ncol=3) +
    scale_fill_manual(values=c("Europeo"="#4472C4","Amerindio"="#ED7D31",
                               "Africano"="#70AD47","Asia del Sur"="#FFC000","Asia Oriental"="#5B9BD5")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.2)), labels=function(x) paste0(x,"%")) +
    labs(title="Ancestría de cada cohorte por separado",
         subtitle="NGSadmix K=5 — 239,000 SNPs, 22 cromosomas", x=NULL, y=NULL) +
    theme_minimal(base_size=14) +
    theme(plot.title=element_text(face="bold",size=18), strip.text=element_text(face="bold",size=15),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1))
ggsave("out_global/angsd/ngsadmix/pools_individuales_K5.png", p3, width=14, height=6, dpi=150)
cat("OK: pools_individuales_K5.png\n")

# ── 4. K=6 granular: sub-continental por pool ──
ind2 <- fread("out_global/angsd/ngsadmix/individuals_v2_full.tsv")
q6 <- fread("out_global/angsd/ngsadmix/results_v2_full/K6.qopt", header=FALSE)
setnames(q6, paste0("C", 1:6))
q6[, sample := ind2$sample]; q6[, super_pop := ind2$super_pop]

q6[, pop := ind2$pop]
ref6 <- q6[super_pop != "POOL"]
cols6 <- paste0("C", 1:6)
means6 <- ref6[, lapply(.SD, mean), by=pop, .SDcols=cols6]
cn6 <- character(6); u6 <- character(0)
for (k in 1:6) {
    col <- paste0("C",k)
    b <- means6$pop[which.max(means6[[col]])]
    if (b %in% u6) b <- paste0(b,"(",k,")")
    cn6[k] <- b; u6 <- c(u6, b)
}
setnames(q6, cols6, cn6)

pools6 <- q6[super_pop == "POOL"]
p6_long <- melt(pools6, id.vars=c("sample","super_pop"), measure.vars=cn6,
                variable.name="Cluster", value.name="Prop")
p6_long[, Pct := round(Prop*100, 1)]
p6_long[, sample := factor(sample, levels=c("HOSPITAL","POOL2","POOL1"))]

region_map <- c(FIN="Europa Nórdica", ESN="África Occidental",
                IBS="Europa Ibérica", PUR="América Caribe",
                PEL="América Andina", GWD="África Occ. (Gambia)")
p6_long[, Region := sapply(as.character(Cluster), function(x) if(x %in% names(region_map)) region_map[x] else x)]
p6_long[, Region := factor(Region, levels=c("Europa Ibérica","Europa Nórdica",
         "América Andina","América Caribe","África Occidental","África Occ. (Gambia)"))]

p4 <- ggplot(p6_long, aes(x=sample, y=Pct, fill=Region)) +
    geom_col(width=0.6, color="white", linewidth=0.3) +
    geom_text(aes(label=ifelse(Pct >= 4, paste0(Pct,"%"), "")),
              position=position_stack(vjust=0.5), size=4.2, fontface="bold", color="white") +
    coord_flip() +
    scale_fill_manual(values=c("Europa Ibérica"="#4472C4","Europa Nórdica"="#5B9BD5",
                               "América Andina"="#ED7D31","América Caribe"="#F4B183",
                               "África Occidental"="#70AD47","África Occ. (Gambia)"="#A9D18E")) +
    scale_y_continuous(expand=c(0,0), labels=function(x) paste0(x,"%")) +
    labs(title="Composición sub-continental de cada cohorte",
         subtitle="NGSadmix K=6 — Panel granular, 11 poblaciones de referencia",
         x=NULL, y=NULL, fill="Región ancestral") +
    theme_minimal(base_size=15) +
    theme(plot.title=element_text(face="bold",size=18), panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(), legend.position="bottom")
ggsave("out_global/angsd/ngsadmix/composicion_pools_K6.png", p4, width=12, height=5, dpi=150)
cat("OK: composicion_pools_K6.png\n")

cat("\n=== Todas las figuras generadas ===\n")
