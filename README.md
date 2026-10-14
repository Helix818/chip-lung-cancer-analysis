# CHIP and lung cancer: statistical analysis code

This repository contains R code for the statistical analyses and selected figures and tables in *Compounding Effects of Clonal Hematopoiesis and Germline Susceptibility Variants on Lung Cancer Risk*. Analysis tables contain participant-level clinical and genotype information and are supplied to the scripts from a private directory; they are not part of this repository. De-identified data may be requested by qualified researchers from the corresponding author.

## Analysis scripts

| Script | Purpose |
| --- | --- |
| `R/01_descriptive_figures.R` | CHIP prevalence, variant, and clone-size descriptive figures |
| `R/02_core_models.R` | Main and gene-specific lung cancer models, subgroup models, and Figure 3B–E/G–J data and plots |
| `R/03_snp_screening.R` | Candidate SNP screen and omnibus-test audit |
| `R/04_interaction_overall.R` | Overall CHIP–SNP joint effects and interactions |
| `R/04_interaction_DNMT3A.R`, `R/04_interaction_TET2.R`, `R/04_interaction_PPM1D.R` | Original gene-specific interaction models, called by `05_gene_interactions_corrected.R` |
| `R/05_gene_interactions_corrected.R` | Derives gene-specific model inputs from each variant's VAF and runs the three original interaction models |
| `R/06_manuscript_figures.R` | Table-forest style Figure 3A, Figure 4B, and Figures S4–S6 from aggregate results |
| `R/07_supplementary_tables.R` | Corrected Table 2 replacement rows and Supplementary Tables S3–S5 from aggregate results |

## Private input tables

Place these five analysis tables in one directory outside the repository. Keep the filenames shown, or use the environment-variable overrides listed below where supported. Do not commit these files or generated participant-level files.

| File | Required content |
| --- | --- |
| `CHIP_ALL.csv` | One row per participant in the 765-person matched cohort; includes `Filename`, `TYPE`, `CHIPstatus`, and model covariates |
| `CH_matched.csv` | Variant-level calls with `Filename`, `GENE`, and each variant's own numeric `VAF` |
| `CHIP_matched.csv` | CHIP-qualified mutation rows for the Figure 3B–E mutation and clone-size panels |
| `ALL因素(分期).csv` | Stage and histology information keyed by `Filename` |
| `SNP645人加个人信息.csv` | 645-person SNP analysis cohort with `Filename`, genotype, `CHIPstatus`, outcome, and covariates |

`02_core_models.R` and `05_gene_interactions_corrected.R` also accept `CHIP_LC_COHORT_FILE` and `CHIP_LC_VARIANT_FILE` to point to the first two files at their existing locations. `02_core_models.R` additionally accepts `CHIP_LC_STAGE_FILE` and `CHIP_LC_CHIP_MUTATION_FILE`; `05_gene_interactions_corrected.R` accepts `CHIP_LC_SNP_FILE`. Other scripts use the standard filenames under `CHIP_LC_DATA_DIR`.

## Corrected classification

CHIP is defined by at least one variant with VAF ≥ 2%. Gene-specific CHIP requires a variant **in that gene** with VAF ≥ 2%; a lower-VAF variant in another gene does not make a person positive for the target gene. The Figure 3A reference comprises all people without CHIP, including those with only variants below 2% VAF. The gene-specific SNP models compare target-gene CHIP carriers with participants without CHIP; people with CHIP in other genes are excluded from those two-group models. The 645-person `CHIPstatus` column is checked against the variant-level threshold before the models run.

Figure 3G–J includes all lung cancer cases in each panel. Its `Others` category comprises records marked `benign` or `unknown` in the stage source. The `Others` slice uses `#f1f4fa`.

## Run

The scripts were run with R 4.5.0 on Windows. Install the required packages before running: `dplyr`, `tidyr`, `ggplot2`, `tidyverse`, `scales`, `forcats`, `mgcv`, `tibble`, `purrr`, `forestploter`, `stringr`, `readr`, `flextable`, `officer`, and `writexl`. `grid` is included with R.

From the repository root in PowerShell:

```powershell
$env:CHIP_LC_DATA_DIR = 'PATH_TO_PRIVATE_INPUT_DIRECTORY'
$env:CHIP_LC_OUTPUT_DIR = 'PATH_TO_OUTPUT_DIRECTORY_OUTSIDE_REPOSITORY'
Rscript R/01_descriptive_figures.R
Rscript R/02_core_models.R
Rscript R/03_snp_screening.R
Rscript R/04_interaction_overall.R
Rscript R/05_gene_interactions_corrected.R
Rscript R/06_manuscript_figures.R
Rscript R/07_supplementary_tables.R
```

Run each command in a fresh R session. `05_gene_interactions_corrected.R` creates the three `CHIP_<GENE>_SNP.csv` inputs in a temporary directory, calls the corresponding `04_interaction_<GENE>.R` scripts without changing their statistical methods, and deletes the temporary participant-level inputs when it exits. Do not separately run the three gene scripts against previously created input tables.

All persistent output goes to `CHIP_LC_OUTPUT_DIR`. The scripts write model estimates, figure files, aggregate counts, and formatted tables. The source analysis tables remain outside the repository. The forest-style PDFs are programmatic outputs; final manuscript figure layout may include editorial adjustments.

## Validation reference

The corrected 765-person cohort contains 582 participants without CHIP and 96 DNMT3A-, 35 TET2-, and 21 PPM1D-CHIP carriers. The 645-person SNP cohort contains 473 participants without CHIP and 92, 34, and 19 carriers, respectively. Figure 3A gene-specific ORs from the corrected run are 1.2732, 2.4488, and 4.2289, respectively. Figure 3G–J panel sizes are 293, 122, 18, and 27; the respective `Others` counts are 42, 17, 5, and 1. The DNMT3A × rs7178270 overall multiplicative interaction has P = 0.0179; TET2 and PPM1D have no SNP selected under the original main-table rule.

All seven run commands above completed with the stated restricted inputs and R 4.5.0. The gene-specific interaction CSVs matched the previously corrected reference outputs byte-for-byte for the joint-effect and interaction tables. This validation confirms the packaged statistical code and inputs used here; access to the private analysis tables is required to rerun it independently.
