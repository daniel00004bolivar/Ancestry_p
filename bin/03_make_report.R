#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(stringr)
  library(readr)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--input", type = "character", help = "Path to <prefix>.ancestry.out"),
  make_option("--outdir", type = "character", default = "reports", help = "Output directory"),
  make_option("--panel", type = "character", default = NA, help = "Panel name (e.g., N5)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || opt$input == "") {
  stop("Falta --input.\nEjemplo:\n  Rscript bin/03_make_report.R --input <prefix>.ancestry.out --outdir reports --panel N5",
       call. = FALSE)
}
if (!file.exists(opt$input)) {
  stop(paste0("No existe el archivo --input: ", opt$input), call. = FALSE)
}

outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

lines <- readLines(opt$input, warn = FALSE)

snps_n <- NA_integer_
m1 <- str_match(lines[1], "^snps in file\\s+([0-9]+)\\s+pops\\s+([0-9]+)")
if (!is.na(m1[1, 2])) snps_n <- as.integer(m1[1, 2])

idx_final <- which(str_detect(lines, "FINAL_NZ_PROPS"))
if (length(idx_final) == 0) {
  stop("No se encontró 'FINAL_NZ_PROPS' en el archivo. Verifica formato de iAdmix.", call. = FALSE)
}
final_line <- lines[max(idx_final)]

pairs <- str_extract_all(final_line, "[A-Za-z0-9]+:[0-9]*\\.?[0-9]+")[[1]]
if (length(pairs) == 0) {
  stop("No se pudieron extraer pares POP:valor desde la línea final.", call. = FALSE)
}

df <- tibble(pair = pairs) |>
  separate(pair, into = c("pop", "prop"), sep = ":", convert = TRUE) |>
  mutate(
    input = basename(opt$input),
    panel = opt$panel,
    snps_in_file = snps_n
  ) |>
  arrange(desc(prop))

out_tsv <- file.path(outdir, "tables", "admixture_props.tsv")
write_tsv(df, out_tsv)

cat("OK - wrote:\n - ", out_tsv, "\n", sep = "")
