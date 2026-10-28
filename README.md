# pan-cancer-sncRNA-atlas

Code repository for the pan-cancer sncRNA atlas analysis.

The atlas profiles non-canonical small RNA regions across 16 TCGA cancer types.
Peaks are called with `cfPeak` / `cfPeakCNN`, summarised into a pan-cancer peak
matrix, and analysed for differential expression, clinical and molecular
association, pathway enrichment, cross-species conservation and tissue-type
classification.

The peak calling code itself lives in the cfPeak repository
(https://github.com/HUNNNGRY/cfPeak/tree/master/scripts).


# Repository layout

```
scripts/
  config.R                                   paths, global options, output folders
  TCGA_sncRNA_util.R                         loads the peak matrix, sample table and palettes
  TCGA_sncRNA_cfPeak.R                       QC, EDA, differential expression, PCC, figures
  TCGA_longRNA_cfPeak.R                      differential expression of host mRNA
  TCGA_small_variancePartition_analysis.ipynb variance partition of the peak matrix
  gnCoor_conservation.R                      liftOver, phastCons and multi-species alignment
  random_forest_tissue_classification_cv_v2.R tissue-type classifier (k-fold CV)
  bash/
    _config.sh                               shared shell settings
    01_bed_jaccard.sh                        pairwise BED Jaccard between conditions
    02_pcc_clinical.sh                       PCC against clinical variables
    03_pcc_mrna.sh                           PCC against mRNA expression
    04_snc_pcc_gsea_gmt.sh                   long-format PCC tables for GSEA
    05_pcc_cibersort.sh                      PCC against CIBERSORTx immune fractions
    06_pathway_ora.sh                        pathway over-representation, all features
    07_pathway_ora_by_rna.sh                 pathway over-representation, per RNA class
    08_random_forest.sh                      trains and evaluates the classifier
```


# Requirements

- R (>= 4.1) with `data.table`, `ggplot2`, `ggpubr`, `pheatmap`, `ComplexHeatmap`,
  `edgeR`, `variancePartition`, `msa`, `rtracklayer`, `randomForest`, `caret`
- `bedtools` (through `bedtoolsr`)
- `clusterProfiler` for `run-clusterProfiler-sncGmt.R` (the pathway scripts)
- optional: SLURM for the per-cancer jobs submitted by `scripts/bash/*.sh`

Some helpers are maintained outside this repository and are loaded through the
`*_DIR` variables in `scripts/config.R`: `util.R` (exSeek), `sncRNA_utils.R`,
`stats.R`, `run-covarPCC.R`, `run-SncPccTab2GseaGmt.R` and
`run-clusterProfiler-sncGmt.R`. Point those variables at your own copies.


# Setup

All paths are defined once in `scripts/config.R` (R) and `scripts/bash/_config.sh`
(shell) and are derived from a single environment variable:

```bash
export SNCRNA_HOME=/path/to/data/mirror
```

The expected layout under `$SNCRNA_HOME` is documented at the top of
`scripts/config.R`. Private overrides can be placed in `scripts/config.local.R`,
which is not tracked by git.


# Running

Scripts are started from the repository root, for example:

```bash
Rscript scripts/TCGA_sncRNA_cfPeak.R
Rscript scripts/gnCoor_conservation.R
bash scripts/bash/03_pcc_mrna.sh          # cluster jobs, 16 cancer types
```

`scripts/TCGA_sncRNA_cfPeak.R` expects the peak matrix and sample table produced
by the cfPeak pipeline, and writes to `output/` and `figure/`. Sections that
were originally submitted as `sbatch --wrap` calls are kept in `scripts/bash/`
and are referenced from the script at the position where their output is used.


# Data availability

The peak matrix, sample annotation and derived tables are available from the
corresponding author on request, and from the public repositories listed in the
manuscript.


# Contact

@hunnngry (bpf19@tsinghua.org.cn)

@lulab (lulab1@tsinghua.edu.cn)


# License

MIT, see LICENSE.


# References

P. Bao et al., 2026, The Catalogue of Non-Canonical Small RNA regions Across Human Cancers. Submitted.
