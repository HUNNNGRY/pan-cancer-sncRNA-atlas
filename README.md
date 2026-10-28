# pan-cancer-sncRNA-atlas

Code for the pan-cancer sncRNA atlas: non-canonical small RNA regions across 16
TCGA cancer types, from cfPeak peak calls to differential expression, clinical
and molecular association, pathway enrichment, conservation and tissue-type
classification.

Peak calling itself is in the cfPeak repository:
<https://github.com/HUNNNGRY/cfPeak/tree/master/scripts>


## Files

```
scripts/config.R                                      paths and global options
scripts/TCGA_sncRNA_util.R                            loads peak matrix, sample table, palettes
scripts/TCGA_sncRNA_cfPeak.R                          QC, EDA, differential expression, PCC, figures
scripts/TCGA_longRNA_cfPeak.R                         host mRNA differential expression
scripts/TCGA_small_variancePartition_analysis.ipynb   variance partition
scripts/gnCoor_conservation.R                         liftOver, phastCons, alignments
scripts/random_forest_tissue_classification_cv_v2.R   tissue classifier (k-fold CV)
scripts/bash/                                         cluster steps called from TCGA_sncRNA_cfPeak.R
```


## Requirements

R (>= 4.1) with `data.table`, `ggplot2`, `ggpubr`, `pheatmap`, `ComplexHeatmap`,
`edgeR`, `variancePartition`, `msa`, `rtracklayer`, `randomForest`, `caret` and
`clusterProfiler`; `bedtools` via `bedtoolsr`; optionally SLURM.


## Setup

```bash
export SNCRNA_HOME=/path/to/data/mirror
```

All paths derive from `SNCRNA_HOME` in `scripts/config.R` and
`scripts/bash/_config.sh`; the expected layout is documented at the top of
`scripts/config.R`. Helpers kept outside this repository are configured by the
same variables. Local overrides go in `scripts/config.local.R` (git-ignored).


## Running

```bash
Rscript scripts/TCGA_sncRNA_cfPeak.R
Rscript scripts/gnCoor_conservation.R
bash scripts/bash/03_pcc_mrna.sh
```


## Contact

@hunnngry (bpf19@tsinghua.org.cn)

@lulab (lulab1@tsinghua.edu.cn)


## License

MIT, see LICENSE.


## Reference

P. Bao et al., 2026, The Catalogue of Non-Canonical Small RNA regions Across Human Cancers. Submitted.
