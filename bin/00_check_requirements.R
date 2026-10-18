#!/usr/bin/env Rscript
# 00_check_requirements.R — can this machine run the pipeline?
#
# Checks the command-line tools, R packages and input files the pipeline
# needs, and reports what is missing before anything long starts running.
# Exits non-zero if a requirement is missing, so run_pipeline.sh stops early.
#
# Usage: Rscript bin/00_check_requirements.R
# ─────────────────────────────────────────────────────────────────────────────

ok <- TRUE
say <- function(status, what, detail = "") {
    mark <- if (status) "  ok   " else "  MISSING "
    cat(sprintf("%s %-34s %s\n", mark, what, detail))
    if (!status) ok <<- FALSE
}

cat("\nCommand-line tools\n")
# Test what each tool must actually do, not its --version flag: BSD and GNU
# builds differ in how they report it (macOS gzip has no --version at all).
probe <- function(cmd) suppressWarnings(tryCatch(
    system(cmd, intern = TRUE, ignore.stderr = TRUE)[1], error = function(e) NA))
tmp <- tempfile(); writeLines("probe", tmp); system2("gzip", c("-f", tmp))
checks <- list(
    samtools = list(probe("samtools --version"), "not on PATH"),
    gzip     = list(probe(sprintf("gzip -dc %s.gz", tmp)), "gzip -dc does not work"),
    awk      = list(probe("echo probe | awk '{print $1}'"), "not on PATH"),
    gcc      = list(probe("gcc --version"), "not on PATH (needed only by bin/44 and bin/45)")
)
unlink(paste0(tmp, ".gz"))
for (t in names(checks)) {
    v <- checks[[t]][[1]]
    say(!is.na(v) && nzchar(v), t, if (is.na(v)) checks[[t]][[2]] else substr(v, 1, 40))
}
cat("  note   ANGSD is needed only for bin/05, which produces the read counts.\n")
cat("         If you already have out_global/angsd/per_pool/, you do not need it.\n")

cat("\nR packages\n")
need <- c(`data.table` = "1.14.0", poolfstat = "3.0.0")
for (p in names(need)) {
    has <- requireNamespace(p, quietly = TRUE)
    v   <- if (has) as.character(utils::packageVersion(p)) else NA
    good <- has && utils::compareVersion(v, need[[p]]) >= 0
    say(good, p, if (!has) sprintf("install.packages(\"%s\")", p)
                 else sprintf("%s (need >= %s)", v, need[[p]]))
}
cat(sprintf("  ok    %-34s %s\n", "R", R.version.string))
if (.Platform$OS.type == "windows")
    cat("  note   parallel::mclapply does not fork on Windows; the jackknife and\n",
        "        the simulations will run serially and take considerably longer.\n")

cat("\nReference data\n")
files <- c("data/hgdp/hgdpGeo.txt.gz"            = "bash bin/01_fetch_reference_data.sh",
           "data/ref/human_g1k_v37_decoy.fasta"  = "bash bin/01_fetch_reference_data.sh",
           "data/ref/human_g1k_v37_decoy.fasta.fai" = "samtools faidx on the reference")
for (f in names(files)) say(file.exists(f), f, if (file.exists(f)) "" else files[[f]])

cat("\nInput alignments\n")
bams <- list.files("bam", pattern = "\\.bam$", full.names = TRUE)
say(length(bams) >= 3, "bam/*.bam",
    if (length(bams)) sprintf("%d found", length(bams))
    else "controlled access; see DATA_AVAILABILITY.md")
counts <- list.files("out_global/angsd/per_pool", pattern = "_counts\\.pos\\.gz$")
if (length(counts)) cat(sprintf("  ok     %-34s %d pools already counted; bin/05 can be skipped\n",
                                "out_global/angsd/per_pool/", length(counts)))

cat(sprintf("\n%s\n\n", if (ok) "All requirements met." else
            "Some requirements are missing; see above."))
quit(status = if (ok) 0 else 1)
