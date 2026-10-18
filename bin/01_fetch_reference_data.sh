#!/usr/bin/env bash
# 01_fetch_reference_data.sh — download the public reference data the pipeline needs
#
# Two files, both public, neither carried in this repository because of size:
#
#   data/hgdp/hgdpGeo.txt.gz             ~53 MB   HGDP-CEPH allele frequencies
#                                                 for 53 populations on hg19,
#                                                 from the UCSC Genome Browser.
#                                                 Used by bin/30 and bin/38.
#   data/ref/human_g1k_v37_decoy.fasta   ~3 GB    GRCh37/b37 reference genome
#                                                 with decoy, from the 1000
#                                                 Genomes FTP. Used by bin/30
#                                                 (reference base per site) and
#                                                 by the alignment step.
#
# Both are skipped if already present. The reference is indexed with samtools
# faidx afterwards, which bin/30 requires.
#
# Usage: bash bin/01_fetch_reference_data.sh
set -euo pipefail
cd "$(dirname "$0")/.."

HGDP_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/hgdpGeo.txt.gz"
REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37_decoy.fasta.gz"

mkdir -p data/hgdp data/ref

if [[ -f data/hgdp/hgdpGeo.txt.gz ]]; then
    echo "data/hgdp/hgdpGeo.txt.gz already present, skipping"
else
    echo "Downloading HGDP allele frequencies (~53 MB)..."
    curl -fL --retry 3 -o data/hgdp/hgdpGeo.txt.gz "$HGDP_URL"
fi
# Sanity check: 53 comma-separated frequencies per row
cols=$(gzip -dc data/hgdp/hgdpGeo.txt.gz | head -1 | cut -f8 | tr -cd ',' | wc -c | tr -d ' ')
if [[ "$cols" -lt 52 ]]; then
    echo "ERROR: hgdpGeo.txt.gz does not look right ($cols commas in the frequency field, expected 52-53)" >&2
    exit 1
fi
echo "  ok: $cols populations per row"

if [[ -f data/ref/human_g1k_v37_decoy.fasta ]]; then
    echo "data/ref/human_g1k_v37_decoy.fasta already present, skipping"
else
    echo "Downloading GRCh37/b37 reference (~3 GB compressed; this takes a while)..."
    curl -fL --retry 3 -o data/ref/human_g1k_v37_decoy.fasta.gz "$REF_URL"
    echo "Decompressing..."
    gzip -d data/ref/human_g1k_v37_decoy.fasta.gz
fi

if [[ ! -f data/ref/human_g1k_v37_decoy.fasta.fai ]]; then
    echo "Indexing the reference..."
    samtools faidx data/ref/human_g1k_v37_decoy.fasta
fi

echo
echo "Reference data ready. Check everything else with:"
echo "  Rscript bin/00_check_requirements.R"
