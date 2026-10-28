# ---------------------------------------------------------------------------
# config.R - central path configuration for the pan-cancer sncRNA atlas
#
# Every script in this repository is started from this file, so all input and
# output locations are defined in one place. No location below refers to a
# specific user account: point SNCRNA_HOME at your local copy of the data and
# everything else is derived from it.
#
#   export SNCRNA_HOME=/path/to/data/mirror
#   Rscript scripts/TCGA_sncRNA_cfPeak.R      # run from the repository root
#
# For private overrides that should not be committed, create
# `scripts/config.local.R` (git-ignored) and re-assign any variable there; it
# is sourced at the end of this file.
#
# Layout expected under $SNCRNA_HOME (all of it is data, none of it is code):
#
#   shared_reference/             gene sets, transcript tables,
#                                 genomes/, ucsc/liftover/
#   shared_utils/                 run-covarPCC.R, run-SncPccTab2GseaGmt.R, stats.R
#   shared_scripts/               run-clusterProfiler-sncGmt.R
#   biosoft/                      bedtools binaries
#   projects/WCHSU-FTC/           cfPeak/cfPeakCNN inputs, peak matrix, peak outputs
#   projects/motif-RBP-EDA/       analysis workspace: output/, figure/, tmp/, logs/
#   projects/tcga/                TCGA longRNA expression and clinical tables
#   projects/lulab_tissue_pico/   CIBERSORTx script and LM22 signature
#   smb/                          figure and table drops
# ---------------------------------------------------------------------------

# root of the local data mirror
DATA_ROOT <- Sys.getenv("SNCRNA_HOME", unset = file.path(path.expand("~"), "sncRNA-atlas-data"))

# ---- data roots -----------------------------------------------------------
REF_DIR   <- file.path(DATA_ROOT, "shared_reference")                    # gene sets, transcript tables
DATA_DIR  <- file.path(DATA_ROOT, "projects", "WCHSU-FTC")               # cfPeak/cfPeakCNN pipeline, raw inputs, peak outputs
PROJECT_DIR <- file.path(DATA_ROOT, "projects", "motif-RBP-EDA")         # analysis workspace: output/ and figure/ live here
TCGA_DIR  <- file.path(DATA_ROOT, "projects", "tcga")                    # TCGA longRNA expression and clinical tables
CIBERSORT_DIR <- file.path(DATA_ROOT, "projects", "lulab_tissue_pico")   # CIBERSORTx script and LM22 signature
SMB_DIR   <- file.path(DATA_ROOT, "smb")                                 # figure/table drops
BEDTOOLS_DIR <- file.path(DATA_ROOT, "biosoft")                          # bedtools binaries

# reference genomes and UCSC liftOver chain files
GENOME_DIR   <- file.path(REF_DIR, "genomes")
LIFTOVER_DIR <- file.path(REF_DIR, "ucsc", "liftover")

# helper scripts maintained outside this repository
SHARED_UTILS_DIR   <- file.path(DATA_ROOT, "shared_utils")               # run-covarPCC.R, run-SncPccTab2GseaGmt.R, stats.R
SHARED_SCRIPTS_DIR <- file.path(DATA_ROOT, "shared_scripts")             # run-clusterProfiler-sncGmt.R
EXSEEK_SCRIPTS_DIR <- file.path(DATA_DIR, "exSeek-dev", "scripts")      # exotic peak-calling helpers

# ---- dataset --------------------------------------------------------------
dst   <- "TCGA_small16"   # TCGA sncRNA cohort, 16 cancer types
dedup <- "_all"

# short aliases; the analysis code uses them pervasively
pre     <- DATA_DIR       # cfPeak inputs: sample tables, test sample list
outpre  <- DATA_DIR       # cfPeak outputs: count matrix, peak summary, logs
diffpre <- PROJECT_DIR    # differential expression and figure workspace

# ---- derived locations ----------------------------------------------------
OUT_DIR <- file.path(PROJECT_DIR, "output", dst)
FIG_DIR <- file.path(PROJECT_DIR, "figure")

# Some sections switch the working directory to DATA_DIR and keep writing to
# the relative path "./figure/", so both candidates are created here.
for (d in c(FIG_DIR,
            file.path(DATA_DIR, "figure"),
            OUT_DIR,
            LIFTOVER_DIR,
            file.path(PROJECT_DIR, "tmp"),
            file.path(PROJECT_DIR, "logs"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- global options -------------------------------------------------------
options(scipen = 4, stringsAsFactors = FALSE, digits = 5)
options(bedtools.path = BEDTOOLS_DIR)

# ---- local overrides (git-ignored) ---------------------------------------
if (file.exists("scripts/config.local.R")) source("scripts/config.local.R")
