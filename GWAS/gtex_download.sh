#!/bin/bash
# GTEx v8 EUR FUSION Weights Batch Download & Auto-Extract Script
# Usage: bash gtex_download.sh [download_dir]
# Features: resume support, auto-extract, retry, idempotent

set -euo pipefail

# -------------------- Configuration --------------------
DOWNLOAD_DIR="${1:-./GTEx_v8_EUR}"
MAX_RETRIES=3
RETRY_DELAY=10
LOG_FILE="${DOWNLOAD_DIR}/download.log"

# URL list (49 tissues)
URLS=(
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Adipose_Subcutaneous.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Adipose_Visceral_Omentum.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Adrenal_Gland.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Artery_Aorta.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Artery_Coronary.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Artery_Tibial.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Amygdala.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Anterior_cingulate_cortex_BA24.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Caudate_basal_ganglia.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Cerebellar_Hemisphere.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Cerebellum.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Cortex.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Frontal_Cortex_BA9.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Hippocampus.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Hypothalamus.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Nucleus_accumbens_basal_ganglia.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Putamen_basal_ganglia.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Spinal_cord_cervical_c-1.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Brain_Substantia_nigra.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Breast_Mammary_Tissue.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Cells_Cultured_fibroblasts.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Cells_EBV-transformed_lymphocytes.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Colon_Sigmoid.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Colon_Transverse.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Esophagus_Gastroesophageal_Junction.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Esophagus_Mucosa.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Esophagus_Muscularis.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Heart_Atrial_Appendage.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Heart_Left_Ventricle.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Kidney_Cortex.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Liver.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Lung.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Minor_Salivary_Gland.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Muscle_Skeletal.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Nerve_Tibial.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Ovary.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Pancreas.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Pituitary.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Prostate.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Skin_Not_Sun_Exposed_Suprapubic.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Skin_Sun_Exposed_Lower_leg.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Small_Intestine_Terminal_Ileum.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Spleen.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Stomach.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Testis.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Thyroid.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Uterus.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Vagina.tar.gz"
"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.Whole_Blood.tar.gz"
)

# -------------------- Functions --------------------

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "$msg"
    echo "$msg" >> "$LOG_FILE"
}

download_with_resume() {
    local url="$1"
    local output="$2"
    local attempt=1

    while [[ $attempt -le $MAX_RETRIES ]]; do
        log "Downloading (attempt $attempt/$MAX_RETRIES): $(basename "$output")"

        if wget -c --progress=bar:force -O "$output" "$url" 2>&1 | tee -a "$LOG_FILE"; then
            log "Download completed: $(basename "$output")"
            return 0
        else
            log "Download failed (attempt $attempt): $(basename "$output")"
            if [[ $attempt -lt $MAX_RETRIES ]]; then
                log "Retrying in $RETRY_DELAY seconds..."
                sleep $RETRY_DELAY
            fi
        fi
        ((attempt++))
    done

    log "ERROR: All $MAX_RETRIES attempts failed for $(basename "$output")"
    return 1
}

extract_archive() {
    local archive="$1"
    local extract_dir="${archive%.tar.gz}"

    if [[ -d "$extract_dir" ]] && [[ "$(ls -A "$extract_dir" 2>/dev/null)" ]]; then
        log "Already extracted: $(basename "$extract_dir")"
        return 0
    fi

    log "Extracting: $(basename "$archive") -> $(basename "$extract_dir")/"
    mkdir -p "$extract_dir"

    if tar -xzf "$archive" -C "$extract_dir" --strip-components=0 2>>"$LOG_FILE"; then
        # Some archives have nested structure; if empty after extract, try without strip
        if [[ -z "$(ls -A "$extract_dir" 2>/dev/null)" ]]; then
            rmdir "$extract_dir" 2>/dev/null
            mkdir -p "$extract_dir"
            tar -xzf "$archive" -C "$extract_dir" 2>>"$LOG_FILE"
        fi
        log "Extraction completed: $(basename "$extract_dir")/"
        return 0
    else
        log "ERROR: Extraction failed for $(basename "$archive")"
        return 1
    fi
}

# -------------------- Main --------------------

mkdir -p "$DOWNLOAD_DIR"
touch "$LOG_FILE"

TOTAL=${#URLS[@]}
log "Starting batch download of $TOTAL GTEx v8 EUR tissue files to: $DOWNLOAD_DIR"
log "Log file: $LOG_FILE"

SUCCESS=0
FAIL=0

for i in "${!URLS[@]}"; do
    url="${URLS[$i]}"
    filename=$(basename "$url")
    filepath="${DOWNLOAD_DIR}/${filename}"
    tissue_name="${filename%.tar.gz}"

    echo ""
    log "[$((i+1))/$TOTAL] Processing: $tissue_name"

    # Skip if already extracted
    extract_dir="${filepath%.tar.gz}"
    if [[ -d "$extract_dir" ]] && [[ "$(ls -A "$extract_dir" 2>/dev/null)" ]]; then
        log "  -> Already extracted, skipping."
        SUCCESS=$((SUCCESS + 1))
        continue
    fi

    # Download if not present or incomplete
    if [[ -f "$filepath" ]]; then
        # Validate if it's a complete tar.gz (check gzip integrity)
        if tar -tzf "$filepath" >/dev/null 2>&1; then
            log "  -> Archive already downloaded and valid."
        else
            log "  -> Archive incomplete/corrupt, resuming download..."
            if ! download_with_resume "$url" "$filepath"; then
                 FAIL=$((FAIL + 1))
                continue
            fi
        fi
    else
        if ! download_with_resume "$url" "$filepath"; then
             FAIL=$((FAIL + 1))
            continue
        fi
    fi

    # Extract
    if extract_archive "$filepath"; then
        SUCCESS=$((SUCCESS + 1))
    else
         FAIL=$((FAIL + 1))
    fi

done

echo ""
log "========================================"
log "Batch processing completed!"
log "Success: $SUCCESS / $TOTAL"
log "Failed:  $FAIL / $TOTAL"
log "Output directory: $DOWNLOAD_DIR"
log "Log file: $LOG_FILE"
log "========================================"

if [[ $FAIL -gt 0 ]]; then
    exit 1
else
    exit 0
fi