#!/usr/bin/env bash
# Download the GSE208240 processed Perturb-seq data (Sunshine et al. 2023).
#
# Usage:  bash scripts/00_download.sh [DEST]
# Default DEST is ./data/raw
#
# Files (~1.2 GB):
#   GSE208240_CRISPRi_perturbseq_sarscov2_filtered.tar.gz
#       barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz  (CellRanger 3.1.0 filtered matrix)
#       cell_identities.csv                                (authors' per-cell guide calls)
#   GSE208240_viral_references.tar.gz
#       SARS-CoV-2 + lentivirus FASTA and GTF used to build the custom reference
set -euo pipefail

DEST="${1:-data/raw}"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208240/suppl"

mkdir -p "$DEST"
cd "$DEST"

for f in GSE208240_CRISPRi_perturbseq_sarscov2_filtered.tar.gz \
         GSE208240_viral_references.tar.gz ; do
    if [[ -f "$f" ]]; then
        echo "[skip] $f already present"
    else
        echo "[get ] $f"
        curl -fSL -O "$BASE/$f"
    fi
done

echo "[md5 ]"
md5sum ./*.tar.gz | tee md5sums.txt

echo "[untar]"
tar -xzf GSE208240_CRISPRi_perturbseq_sarscov2_filtered.tar.gz
tar -xzf GSE208240_viral_references.tar.gz

# The tarballs expand to the authors' original absolute-ish path; expose a flat symlink.
MTX_DIR="$(find . -type d -name CRISPRi_perturbseq_sarscov2_filtered | head -1)"
REF_DIR="$(find . -type d -name viral_references | head -1)"
ln -sfn "$MTX_DIR" matrix_dir
ln -sfn "$REF_DIR" viral_references

echo
echo "Done. Matrix directory: $DEST/matrix_dir"
ls -l matrix_dir/
