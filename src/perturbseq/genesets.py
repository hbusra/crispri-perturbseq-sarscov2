"""Gene sets used across the analysis, with provenance."""

# ----------------------------------------------------------------- gene sets
# Interferon-stimulated genes. Core type-I ISGs that are unambiguous markers of
# an autocrine/paracrine IFN response in epithelial cells. Compiled from the
# MSigDB Hallmark IFN-alpha response leading edge, restricted to genes detected
# in this dataset at 01. Script 06 reports the score with and without the
# genes that are themselves CRISPRi targets in this library, since knocking a
# gene down mechanically lowers any score that contains it.
ISG_GENES = [
    "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "MX1", "MX2",
    "OAS1", "OAS2", "OAS3", "OASL", "IFI6", "IFI27", "IFI35", "IFI44",
    "IFI44L", "IFITM1", "IFITM2", "IFITM3", "STAT1", "STAT2", "IRF7", "IRF9",
    "BST2", "RSAD2", "CMPK2", "XAF1", "HERC5", "HERC6", "EIF2AK2", "DDX58",
    "IFIH1", "PARP9", "DTX3L", "SAMD9", "SAMD9L", "EPSTI1", "LY6E", "LAP3",
    "PLSCR1", "TRIM22", "SP100", "SP110", "B2M", "PSMB8", "PSMB9", "TAP1",
    "NMI", "CXCL10", "CXCL11", "IFI16", "ZBP1", "APOL6", "ADAR",
]

# Type I/III interferons and the chemokines that mark the rare "IFN-producing"
# subset (bystander cluster L / infected cluster M in the paper).
IFN_PRODUCING = [
    "IFNB1", "IFNL1", "IFNL2", "IFNL3", "IFNA1", "IL6", "CXCL8", "CCL5",
    "TNF", "IL1B",
]

# NF-kB target genes reported as induced in infected cells (paper Fig 2E).
NFKB_GENES = [
    "NFKBIA", "NFKBIE", "NFKBIZ", "NFKB1", "NFKB2", "REL", "RELA", "RELB",
    "EGR1", "TNFAIP3", "BIRC3", "CXCL1", "CXCL2", "CXCL3", "SOD2", "ICAM1",
]

VIRAL_GENES = [
    "SCV_ORF1ab", "SCV_S", "SCV_ORF3a", "SCV_E", "SCV_M",
    "SCV_ORF6", "SCV_ORF7ab", "SCV_ORF8", "SCV_N",
]
