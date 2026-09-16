#!/usr/bin/env bash
set -euo pipefail

deactivate 2>/dev/null || true

conda activate cytospace_v1.1.0

CYTO_PY="/mnt/DATA/home/zqy1234560915/miniconda3/envs/cytospace_v1.1.0/bin/python"
CYTO_SRC="/mnt/DATA/home/zqy1234560915/EOPC/RCTD/cytospace"

BASE="/mnt/DATA/home/zqy1234560915/EOPC/RCTD"
INPUT="${BASE}/CytoSPACE_input"

SCRNA="${INPUT}/scRNA_reference/scRNA_expression.mtx"
LABEL="${INPUT}/scRNA_reference/cell_type_labels.txt"
FRAC="${INPUT}/scRNA_reference/cell_type_fraction_estimates.txt"

cd "${CYTO_SRC}"

for SAMPLE in \
    "NEADT_11" \
    "NEADT_3" \
    "NEADT_5" \
    "NEADT_7" \
    "NEADT_8" \
    "patient 11 treatment-naive" \
    "patient 3 treatment-naive"
do
    echo "========================================"
    echo "Running CytoSPACE: ${SAMPLE}"
    echo "========================================"

    "${CYTO_PY}" -m cytospace.cytospace \
        -sp "${SCRNA}" \
        -ctp "${LABEL}" \
        -stp "${INPUT}/${SAMPLE}/ST_expression.mtx" \
        -cp "${INPUT}/${SAMPLE}/Coordinates.txt" \
        -ctfep "${FRAC}" \
        -o "${INPUT}/${SAMPLE}/cytospace_results" \
        -mcn 5 \
        -sm lap_CSPR \
        -sss \
        -nosss 5000 \
        -nop 8 \
        -se 123

    echo "Finished: ${SAMPLE}"
done
