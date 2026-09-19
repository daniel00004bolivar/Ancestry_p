# ===== Proyecto Ancestria (ENV) =====
# Portable: si PROJ_DIR no está definido, se asume que env.sh vive en la raíz del repo.

if [[ -z "${PROJ_DIR:-}" ]]; then
  PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  export PROJ_DIR
fi

export CONFIG="$PROJ_DIR/config"
export WORK="$PROJ_DIR/work"
export REF="$PROJ_DIR/refpanels"

# Solo necesarios si vas a reconstruir panel desde VCF
export VCF_DIR="${VCF_DIR:-$PROJ_DIR/data/vcf}"
export SAMP_DIR="${SAMP_DIR:-$PROJ_DIR/data/samples}"
export PANEL_DIR="${PANEL_DIR:-$PROJ_DIR/data/panel}"
export PANEL_FILE="${PANEL_FILE:-$PANEL_DIR/integrated_call_samples_v3.20130502.ALL.panel}"

# iAdmix
export IADMIX_DIR="${IADMIX_DIR:-$PROJ_DIR/data/iadmix}"

# tools
export BCFTOOLS="${BCFTOOLS:-bcftools}"
export SAMTOOLS="${SAMTOOLS:-samtools}"

mkdir -p \
  "$CONFIG" \
  "$WORK" "$WORK/logs" "$WORK/freq_raw" "$WORK/freq_common" \
  "$REF" \
  "$VCF_DIR" "$SAMP_DIR" "$PANEL_DIR"

