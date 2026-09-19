#!/usr/bin/env Rscript
# 08_merge_beagle_ngsadmix.R — merges the 3-pool Beagle with the 100-
# individual 1KGP reference pseudo-Beagle, at the same sites and allele
# order, ready for NGSadmix / PCAngsd.
#
# METHODOLOGICAL WARNING: each pool is represented below as ONE diploid
# pseudo-individual (3 genotype-likelihood columns, same as a real reference
# individual) -- a structural limitation of the Beagle format, not
# configurable. NGSadmix/PCAngsd proportions should be read as the pool's
# population average, not individual ancestry. See bin/README_VIGENTE.md and
# bin/16_pool_replicates.sh.
suppressPackageStartupMessages(library(data.table))
cat("[NOTE] pools treated as diploid pseudo-individuals -- see bin/README_VIGENTE.md\n")

POOLS_BEAGLE <- "out_global/angsd/ngsadmix/pools_hapmap3.beagle.gz"
REF_BEAGLE   <- "out_global/angsd/ngsadmix/ref100.beagle"
OUT_BEAGLE   <- "out_global/angsd/ngsadmix/combined.beagle.gz"
POOL_NAMES   <- c("POOL1","POOL2","HOSPITAL")

pools <- fread(POOLS_BEAGLE, header = TRUE)
old_names <- names(pools)
ind_idx <- rep(seq_along(POOL_NAMES), each = 3)
new_ind_names <- paste0(POOL_NAMES[ind_idx], c("_p1","_p2","_p3"))
setnames(pools, 4:ncol(pools), new_ind_names)

ref <- fread(REF_BEAGLE, header = FALSE)
ref_header_names <- c("marker","allele1","allele2")
setnames(ref, 1:3, ref_header_names)

cat(sprintf("Pools: %d sitios | Referencia: %d sitios\n", nrow(pools), nrow(ref)))

m <- merge(pools, ref, by = c("marker","allele1","allele2"), suffixes = c("",""))
cat(sprintf("Sitios combinados (intersección, alelos consistentes): %d\n", nrow(m)))

dir.create(dirname(OUT_BEAGLE), showWarnings = FALSE, recursive = TRUE)
fwrite(m, OUT_BEAGLE, sep = "\t")

n_ind <- (ncol(m) - 3) / 3
cat(sprintf("Escrito: %s (%d sitios, %d individuos = 3 pools + 100 referencias)\n",
            OUT_BEAGLE, nrow(m), n_ind))

# Guardar lista ordenada de individuos con su etiqueta poblacional/continental
samples <- readLines("/tmp/ref_sample_order.txt")
panel <- fread("data/panel/integrated_call_samples_v3.20130502.ALL.panel")
pop_map <- panel[match(samples, sample), .(sample, pop, super_pop)]
ind_table <- rbind(
    data.table(sample = POOL_NAMES, pop = POOL_NAMES, super_pop = "POOL"),
    pop_map
)
fwrite(ind_table, "out_global/angsd/ngsadmix/individuals.tsv", sep = "\t")
cat("Escrito: out_global/angsd/ngsadmix/individuals.tsv\n")
