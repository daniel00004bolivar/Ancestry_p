#!/bin/bash
set -e
K=$1
SEED=$2
PANEL=$3
cd /Users/danielbolivar/Downloads/Ancestria

if [ "$PANEL" == "v1" ]; then
  LIKES=out_global/angsd/ngsadmix/combined.beagle.gz
  OUTDIR=out_global/angsd/ngsadmix/results_evanno_v1
else
  LIKES=out_global/angsd/ngsadmix/combined_v2.beagle.gz
  OUTDIR=out_global/angsd/ngsadmix/results_evanno_v2
fi
mkdir -p "$OUTDIR"

tools/angsd_src/misc/NGSadmix -likes "$LIKES" -K "$K" -P 3 -seed "$SEED" -minMaf 0.01 \
  -outfiles "$OUTDIR/K${K}_seed${SEED}" \
  > "$OUTDIR/K${K}_seed${SEED}.runlog" 2>&1

echo "DONE K=$K seed=$SEED panel=$PANEL"
