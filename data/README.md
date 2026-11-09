# data/

Nothing in this directory is committed. `scripts/00_download.sh` populates it:

```
data/
  raw/
    GSE208240_CRISPRi_perturbseq_sarscov2_filtered.tar.gz   1.2 GB
    GSE208240_viral_references.tar.gz                        14 KB
    md5sums.txt
    matrix_dir -> .../CRISPRi_perturbseq_sarscov2_filtered/
        barcodes.tsv.gz        90,380 barcodes
        features.tsv.gz        33,549 features (33,538 human + 9 SARS-CoV-2 + 2 lentiviral)
        matrix.mtx.gz          337,664,836 non-zero entries
        cell_identities.csv    the authors' per-cell guide calls
    viral_references -> .../viral_references/
        SCV_NC_045512.2.fasta, lenti.fasta, SARS-CoV-2_lenti.gtf
  perturbseq.h5ad     written by 01_build_anndata.py   (~0.3 GB)
  X_kept.npz          cached subset count matrix       (~1.1 GB)
```

Source: GEO [GSE208240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208240),
deposited by Sunshine et al. (2023). The analysis starts from the authors'
filtered matrix and their own guide assignments, so it is a reanalysis of
processed data, not a reprocessing from FASTQs.

`cell_identities.csv` columns used here: `cell_barcode`, `matched_library_element`,
`match_type`, `good_coverage`, `number_of_guides`, `gem_group`. Only cells with
`match_type == exact_match` and `good_coverage == True` (27,882) enter the analysis
set; barcodes with no guide detected are kept separately as the ambient-RNA
reference for infection calling.
