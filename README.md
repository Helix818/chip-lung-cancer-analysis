# R analysis code for CHIP and lung cancer study

This repository contains R scripts for the statistical analyses and descriptive figures in *Compounding Effects of Clonal Hematopoiesis and Germline Susceptibility Variants on Lung Cancer Risk*. The scripts take derived CHIP and genotype analysis tables as input. Individual-level study data are not included.

## Scripts

| File | Analysis |
|---|---|
| `R/01_descriptive_figures.R` | CHIP prevalence, variant and clone-size descriptive figures |
| `R/02_core_models.R` | CHIP–lung cancer, gene-specific, subgroup, stage and histology models |
| `R/03_snp_screening.R` | Candidate SNP screening and a comparison with an omnibus test |
| `R/04_interaction_overall.R` | Overall CHIP–SNP joint effects and interactions |
| `R/04_interaction_DNMT3A.R` | DNMT3A-specific joint effects and interactions |
| `R/04_interaction_TET2.R` | TET2-specific joint effects and interactions |
| `R/04_interaction_PPM1D.R` | PPM1D-specific joint effects and interactions |

## Inputs

Place the following authorized analysis tables in a private directory outside this repository: `CHIP_ALL.csv`, `CHIP_matched.csv`, `CH_matched.csv`, `ALL因素(分期).csv`, `SNP645人加个人信息.csv`, `CHIP_DNMT3A_SNP.csv`, `CHIP_TET2_SNP.csv`, and `CHIP_PPM1D_SNP.csv`. Keep the filenames as shown. These files can contain participant-level clinical and genotype information; do not commit them to GitHub.

## Run

The scripts target R 4.5.0 and use packages including `dplyr`, `tidyr`, `ggplot2`, `tidyverse`, `scales`, `forcats`, `mgcv`, `tibble`, `purrr`, `forestploter`, `stringr`, `readr`, `flextable`, and `officer`. Exact package versions should be recorded with the validated release.

In PowerShell, run from the repository root after setting a private input directory and a separate writable output directory:

```powershell
$env:CHIP_LC_DATA_DIR = 'PATH_TO_PRIVATE_ANALYSIS_TABLES'
$env:CHIP_LC_OUTPUT_DIR = 'PATH_TO_LOCAL_OUTPUTS'
Rscript R/02_core_models.R
Rscript R/03_snp_screening.R
Rscript R/04_interaction_overall.R
Rscript R/04_interaction_DNMT3A.R
Rscript R/04_interaction_TET2.R
Rscript R/04_interaction_PPM1D.R
Rscript R/01_descriptive_figures.R
```

Run each script in a fresh R session. Outputs are written under `CHIP_LC_OUTPUT_DIR`.

## Validation status

This is a staging copy. Several principal estimates have been independently checked against the manuscript using restricted analysis tables. The R scripts have not yet been executed in a clean R installation on the current machine, and the complete set of manuscript tables and figures has not been compared with fresh script outputs. Complete that check before tagging or citing a publication release.
