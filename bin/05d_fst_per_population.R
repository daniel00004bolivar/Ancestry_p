#!/usr/bin/env Rscript
# 05d_fst_per_population.R — FST pool-seq vs 13 poblaciones de referencia
#
# Referencias disponibles:
#   HapMap3 (8 poblaciones individuales): YRI, MKK, LWK, TSI, CEU, CHB, CHD, JPT
#   1KGP Phase 3 (5 superpoblaciones):    AFR, AMR, EAS, EUR, SAS
#
# Uso:
#   Rscript bin/05d_fst_per_population.R
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

HAPMAP_FILE  <- "data/iadmix/DATA/hapmap3.8populations.hg19.freqs"
PANEL1K_FILE <- "data/iadmix/DATA/hapmap3_5pops.freqs"
COUNTS_FILES <- c("/tmp/POOL1_hapmap3_counts.txt",
                  "/tmp/POOL2_hapmap3_counts.txt",
                  "/tmp/HOSPITAL_hapmap3_counts.txt")
SAMPLES      <- c("POOL1","POOL2","HOSPITAL")
OUT_DIR      <- "out_global/angsd/continental"
MIN_DEPTH    <- 10L
MIN_MAF_REF  <- 0.05

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Panel de referencia completo ──────────────────────────────────────────
cat("Cargando referencias...\n")

# HapMap3: 8 poblaciones individuales
hm <- fread(HAPMAP_FILE, header = TRUE)
setnames(hm, names(hm)[1], "chr")
hm[, chr := sub("^chr","",chr)]
setnames(hm, "position", "pos")
hm[, allele1 := pmin(A1,A2)]; hm[, allele2 := pmax(A1,A2)]
for (p in c("YRI","MKK","LWK","TSI","CEU","CHB","CHD","JPT"))
    hm[, (paste0(p,"_s")) := ifelse(A1<=A2, get(p), 1-get(p))]
hm_ref <- hm[, .(chr, pos, allele1, allele2,
                  YRI=YRI_s, MKK=MKK_s, LWK=LWK_s,
                  TSI=TSI_s, CEU=CEU_s,
                  CHB=CHB_s, CHD=CHD_s, JPT=JPT_s)]

# 1KGP: 5 superpoblaciones
kg <- fread(PANEL1K_FILE, header = TRUE)
setnames(kg, names(kg)[1], "chr")
kg[, chr := as.character(chr)]
kg[, pos := as.integer(pos)]
for (col in c("AFR_1KGP","AMR_1KGP","EAS_1KGP","EUR_1KGP","SAS_1KGP")) {
    kg[get(col)==".", (col) := NA]
    kg[, (col) := as.numeric(get(col))]
}
kg_ref <- kg[, .(chr, pos, allele1, allele2,
                  AFR=AFR_1KGP, AMR=AMR_1KGP, EAS=EAS_1KGP,
                  EUR=EUR_1KGP, SAS=SAS_1KGP)]

# Combinar
combined <- merge(hm_ref, kg_ref, by=c("chr","pos","allele1","allele2"))
cat(sprintf("  Panel combinado: %d SNPs, 13 referencias\n\n", nrow(combined)))

# ── 2. Cargar conteos de los pools ───────────────────────────────────────────
load_pool <- function(sname, f, min_d=10L) {
    ct <- fread(f, col.names=c("chr","pos","depth","cA","cC","cG","cT"))
    ct[, chr := as.character(chr)]
    ct <- ct[depth >= min_d]
    nuc <- c("A","C","G","T")
    ct[, c("cnt1","cnt2","al1","al2") := {
        m <- cbind(cA,cC,cG,cT)
        ords <- t(apply(m,1,order,decreasing=TRUE))
        list(m[cbind(seq_len(nrow(m)),ords[,1])],
             m[cbind(seq_len(nrow(m)),ords[,2])],
             nuc[ords[,1]], nuc[ords[,2]])
    }]
    ct <- ct[cnt2>0]
    ct[, allele1 := pmin(al1,al2)]
    ct[, allele2 := pmax(al1,al2)]
    ct[, freq_pool := ifelse(al1<=al2, cnt1/(cnt1+cnt2), cnt2/(cnt1+cnt2))]
    ct[, maf_pool  := pmin(freq_pool,1-freq_pool)]
    ct[, .(chr,pos,allele1,allele2,freq_pool,maf_pool)]
}

pools <- Map(load_pool, SAMPLES, COUNTS_FILES)
names(pools) <- SAMPLES

# ── 3. FST Wright/Nei por referencia ─────────────────────────────────────────
wright_fst <- function(p1,p2) {
    pb <- (p1+p2)/2; qb <- 1-pb
    mean((p1-p2)^2/(4*pb*qb+1e-12), na.rm=TRUE)
}

# Etiquetas con grupo continental para colorear
pop_meta <- data.frame(
    pop = c("YRI","MKK","LWK","AFR",  "TSI","CEU","EUR",  "CHB","CHD","JPT","EAS",
            "AMR",  "SAS"),
    continent = c(rep("AFR",4), rep("EUR",3), rep("EAS",4), "AMR","SAS"),
    type = c(rep("individual",3),"superpob", rep("individual",2),"superpob",
             rep("individual",3),"superpob", "superpob","superpob"),
    stringsAsFactors = FALSE
)

all_pops <- pop_meta$pop
fst_rows <- list()

cat("FST Wright/Nei por población:\n")
cat(sprintf("%-10s  %-14s  %s\n","Pool","Población","FST"))
cat(strrep("-",45),"\n")

for (sname in SAMPLES) {
    pd <- pools[[sname]]
    m  <- merge(pd, combined, by=c("chr","pos","allele1","allele2"))
    cat(sprintf("\n[%s — %d SNPs en panel]\n", sname, nrow(m)))

    for (pop in all_pops) {
        if (!pop %in% names(m)) next
        ref_f <- m[[pop]]
        pool_f <- m$freq_pool
        valid <- !is.na(ref_f) & ref_f>=MIN_MAF_REF & ref_f<=(1-MIN_MAF_REF)
        n_ok  <- sum(valid)
        if (n_ok < 100) next
        fst_val <- wright_fst(pool_f[valid], ref_f[valid])
        r_val   <- cor(pool_f[valid], ref_f[valid], use="complete.obs")
        cat(sprintf("  %-4s  FST = %.5f  r = %.4f  (n=%d)\n", pop, fst_val, r_val, n_ok))
        fst_rows[[paste(sname,pop)]] <- data.frame(
            pool=sname, pop=pop, n=n_ok, fst=fst_val, r=r_val,
            stringsAsFactors=FALSE
        )
    }
}

fst_all <- do.call(rbind, fst_rows)
fst_all <- merge(fst_all, pop_meta, by="pop")
write.table(fst_all, file.path(OUT_DIR,"fst_13pop.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ── 4. Figura: dot plot ranking ───────────────────────────────────────────────
cont_colors <- c(AFR="#D55E00", AMR="#CC79A7", EAS="#009E73", EUR="#0072B2", SAS="#F0E442")

for (sname in SAMPLES) {
    sub <- fst_all[fst_all$pool==sname, ]
    sub <- sub[order(sub$fst), ]
    sub$pop_f  <- factor(sub$pop, levels=sub$pop)
    sub$is_super <- sub$type=="superpob"

    p <- ggplot(sub, aes(x=fst, y=pop_f, color=continent, shape=is_super)) +
        geom_segment(aes(x=0, xend=fst, y=pop_f, yend=pop_f),
                     color="grey75", linewidth=0.5) +
        geom_point(aes(size=is_super)) +
        scale_color_manual(values=cont_colors, name="Continente") +
        scale_shape_manual(values=c(`TRUE`=18, `FALSE`=16),
                           labels=c(`TRUE`="superpoblación","individual"),
                           name="Tipo") +
        scale_size_manual(values=c(`TRUE`=4, `FALSE`=3), guide="none") +
        scale_x_continuous(limits=c(0, max(sub$fst)*1.15),
                           expand=expansion(mult=c(0,0.05))) +
        geom_text(aes(label=sprintf("%.5f", fst)), hjust=-0.2, size=2.8, color="black") +
        labs(title=sprintf("FST %s vs 13 referencias (HapMap3 + 1KGP)", sname),
             subtitle="Ordenadas de más cercana (menor FST) a más lejana",
             x="FST Wright/Nei (menor = más similar)", y=NULL) +
        theme_bw(base_size=12) +
        theme(plot.title=element_text(face="bold"),
              panel.grid.major.y=element_blank())
    ggsave(file.path(OUT_DIR, sprintf("fst_13pop_%s.png", sname)),
           p, width=8, height=5.5, dpi=180)
    cat(sprintf("\nFigura: fst_13pop_%s.png\n", sname))
}

# ── 5. Figura combinada: 3 pools × 13 poblaciones ────────────────────────────
# Ordenar poblaciones por FST promedio entre los 3 pools
avg_fst <- aggregate(fst~pop, data=fst_all, FUN=mean)
avg_fst <- avg_fst[order(avg_fst$fst), ]
fst_all$pop_f <- factor(fst_all$pop, levels=avg_fst$pop)
fst_all$pool  <- factor(fst_all$pool, levels=SAMPLES)

p_all <- ggplot(fst_all, aes(x=fst, y=pop_f, color=continent, shape=type)) +
    geom_point(size=3, position=position_dodge(width=0.5)) +
    scale_color_manual(values=cont_colors, name="Continente") +
    scale_shape_manual(values=c(superpob=18, individual=16), name="Tipo") +
    facet_wrap(~pool, nrow=1) +
    labs(title="FST pool-seq vs 13 referencias (HapMap3 + 1KGP Phase 3)",
         subtitle="Cada pool vs las mismas referencias, ordenadas por FST promedio",
         x="FST Wright/Nei", y=NULL) +
    theme_bw(base_size=11) +
    theme(strip.text=element_text(face="bold", size=12),
          plot.title=element_text(face="bold"),
          panel.grid.major.y=element_line(color="grey90"))
ggsave(file.path(OUT_DIR,"fst_13pop_combined.png"),
       p_all, width=13, height=5.5, dpi=180)
cat("Figura: fst_13pop_combined.png\n")
