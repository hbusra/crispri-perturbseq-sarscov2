# Reproduce the whole analysis: download, build the object, run the notebooks.
#
#   snakemake -n                        what would run
#   snakemake --cores 4                 run it
#   snakemake --cores 4 hit_calling     run one step and whatever it needs
#   snakemake --dag | dot -Tpng > dag.png
#
# Each rule's inputs and outputs are the tables that notebook really reads and
# writes, so the graph is the actual dependency structure rather than a fixed
# 1-to-8 ordering: editing notebook 03, for example, re-runs 05, 06 and 08
# because they read its hit table, and leaves 02 and 04 alone.
#
# Notebooks run in place, so the copies in the repository keep their own output.
# nbconvert writes the .ipynb after the kernel has written the tables, which
# would leave a notebook newer than its own outputs and re-trigger the rule
# forever; `touch {output}` at the end of each rule settles that.

import os

RAW, H5AD = "data/raw", "data/perturbseq.h5ad"
TABLES, FIGURES = "results/tables", "results/figures"
MATRIX = f"{RAW}/matrix_dir/matrix.mtx.gz"

# Notebook 07 needs the organoid differential-expression tables, which are not
# public. Point ALO_LFC at a copy to include it; without it the rest still runs.
ALO_LFC = os.environ.get("ALO_LFC", "")

def T(*names):
    return [f"{TABLES}/{n}" for n in names]

def F(*names):
    return [f"{FIGURES}/{n}" for n in names]

RUN = (
    "jupyter nbconvert --to notebook --execute --inplace "
    "--ExecutePreprocessor.timeout=3600 {input.nb} > {log} 2>&1 && touch {output}"
)

CELL_STATE = f"{TABLES}/01_cell_state.tsv.gz"
HITS = f"{TABLES}/03_hits_infection.tsv"

ORGANOID = T("07_signature_overview.tsv", "07_element_signature.tsv",
             "07_top_elements.tsv") + F("organoid_signatures.png")


rule all:
    input:
        T("01_qc_per_cell_summary.tsv", "02_knockdown_summary.tsv",
          "03_concordance.tsv", "04_isg_hits.tsv", "05_underpowered_elements.tsv",
          "06_borderline_composition.tsv", "08_glmm_summary.tsv"),
        F("infection_qc.png", "knockdown.png", "concordance.png", "viral_hits.png",
          "interferon.png", "resolution.png", "infection_state.png"),
        ORGANOID if ALO_LFC else [],


rule download:
    """~1.2 GB of processed data from GEO (GSE208240)."""
    output:
        MATRIX,
        f"{RAW}/GSE208240_CRISPRi_perturbseq_sarscov2_filtered.tar.gz",
    log: "logs/00_download.log"
    shell: "bash scripts/00_download.sh {RAW} > {log} 2>&1"


rule build_anndata:
    """Parse the 337-million-entry matrix once and cache it. ~11 minutes."""
    input: MATRIX
    output:
        H5AD,
        f"{TABLES}/00_build_summary.tsv",
    log: "logs/01_build_anndata.log"
    shell:
        "python scripts/01_build_anndata.py --raw {RAW} --out data "
        "--tables {TABLES} > {log} 2>&1"


rule qc_and_infection_calling:
    input:
        H5AD,
        nb="notebook/01_qc_and_infection_calling.ipynb",
    output:
        T("01_qc_per_cell_summary.tsv", "01_qc_by_gem_group.tsv",
          "01_cells_per_element.tsv", "01_ambient_background.tsv",
          "01_classifier_agreement.tsv", "01_cell_state.tsv.gz",
          "01_state_summary.tsv", "01_host_umi_matched_depth.tsv",
          "01_depth_sensitivity.tsv", "01_infected_by_depth_decile.tsv"),
        F("infection_qc.png"),
    log: "logs/01_qc_and_infection_calling.log"
    shell: RUN


rule knockdown:
    input:
        H5AD, CELL_STATE,
        nb="notebook/02_knockdown.ipynb",
    output:
        T("02_knockdown.tsv", "02_knockdown_summary.tsv",
          "02_construct_concordance.tsv", "02_construct_concordance_summary.tsv",
          "02_element_representation.tsv"),
        F("knockdown.png"),
    log: "logs/02_knockdown.log"
    shell: RUN


rule hit_calling:
    input:
        CELL_STATE,
        nb="notebook/03_hit_calling.ipynb",
    output:
        T("03_viral_load_by_element.tsv", "03_paper_hits_descriptive.tsv",
          "03_hits_infection.tsv", "03_null_calibration.tsv", "03_hits_burden.tsv",
          "03_concordance.tsv", "03_robustness.tsv", "03_offset_sensitivity.tsv",
          "03_offset_sensitivity_summary.tsv"),
        F("concordance.png", "viral_hits.png"),
    log: "logs/03_hit_calling.log"
    shell: RUN


rule interferon:
    input:
        H5AD, CELL_STATE,
        nb="notebook/04_interferon.ipynb",
    output:
        T("04_isg_by_element.tsv", "04_isg_hits.tsv", "04_isg_paper_comparison.tsv"),
        F("interferon.png"),
    log: "logs/04_interferon.log"
    shell: RUN


rule coverage_and_resolution:
    input:
        CELL_STATE, HITS, f"{TABLES}/02_element_representation.tsv",
        nb="notebook/05_coverage_and_resolution.ipynb",
    output:
        T("05_downsample_curve.tsv", "05_detectable_effect.tsv",
          "05_underpowered_elements.tsv"),
        F("resolution.png"),
    log: "logs/05_coverage_and_resolution.log"
    shell: RUN


rule infection_state:
    input:
        H5AD, CELL_STATE, HITS,
        nb="notebook/06_infection_state.ipynb",
    output:
        T("06_state_specific_effects.tsv", "06_divergent_examples.tsv",
          "06_program_scores_by_state.tsv", "06_borderline_by_lane.tsv",
          "06_borderline_per_element.tsv", "06_borderline_composition.tsv",
          "06_borderline_host_response.tsv", "06_borderline_isg_sensitivity.tsv"),
        F("infection_state.png"),
    log: "logs/06_infection_state.log"
    shell: RUN


rule organoid_signatures:
    """Needs ALO_LFC; the organoid tables are not public."""
    input:
        H5AD, CELL_STATE,
        [ALO_LFC] if ALO_LFC else [],
        nb="notebook/07_organoid_signatures.ipynb",
    output: ORGANOID
    log: "logs/07_organoid_signatures.log"
    shell: RUN


rule sensitivity_mixed_effects:
    input:
        H5AD, CELL_STATE, HITS,
        nb="notebook/08_sensitivity_mixed_effects.ipynb",
    output: T("08_glmm_vs_mh.tsv", "08_glmm_summary.tsv")
    log: "logs/08_sensitivity_mixed_effects.log"
    shell: RUN


rule test:
    """Estimator tests: MH against statsmodels, NB recovery, null calibration."""
    shell: "python -m pytest tests -q"
