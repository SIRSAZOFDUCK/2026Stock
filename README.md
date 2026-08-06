# Seasonal prescribing in English primary care: reproducibility release

This repository reproduces the analysis for *Seasonal trends in primary care prescribing in England: a longitudinal ecological study*. It contains the canonical analysis, the exact R dependency lock, source provenance, two small version-specific public reference inputs, and the sealed 1 MB validation reference required by the sensitivity and diagnostic stages.

The monthly prescribing archives and registered-patient files are **not distributed in this repository**. They remain available from their public publishers and must be downloaded independently. Generated analysis outputs are also ignored, apart from the sealed validation reference. No manuscript or project-management files are included.

## Repository contents

- `analysis_main.R` — short authoritative coordinator that validates the runtime and dispatches each analysis stage.
- `R/02_import_qc.R` — frozen-source import, aggregation, recoding and Stage 2 quality control.
- `R/03_primary_analysis.R` — eligibility, model fitting, characterisation and primary result tables.
- `R/05_sensitivities.R` — the four prespecified trend, HAC, working-day and threshold sensitivity analyses.
- `R/06_diagnostics.R` — targeted diagnostics, cohort flow and missingness accounting.
- `R/07_reporting.R` — publication figures and the prespecified 2021 window diagnostics.
- `verify_inputs.R` — verifies every acquired source against the frozen size and SHA-256 manifest before analysis.
- `renv.lock`, `.Rprofile`, `renv/activate.R`, `renv/settings.json`, `DESCRIPTION` — the exact R 4.6.1 project environment.
- `reproducibility/input_manifest.csv` — 82-source provenance record with publisher, period, canonical URL, licence, schema, expected filename, size and SHA-256.
- `bnf_code_current_202505_version_88.csv` — the May 2025 BNF hierarchy lookup used in the analysis.
- `reproducibility/source_snapshots/bank-holidays_2026-08-05.json` — the dated GOV.UK England-and-Wales bank-holiday snapshot used for the working-day sensitivity.
- `Outputs/canonical/qc/stage4/stage4_snapshot*` — a small sealed numerical reference used only to verify that the primary class and all-eligible-drug analysis is unchanged before the sensitivity stages run. It is not primary data.

The allowlist-style `.gitignore` excludes everything else, including primary data, generated results, manuscripts, local package libraries and temporary files.

### Release-only implementation changes

The accepted monolithic script was split into the modules above so that the live workflow can be reviewed and maintained by stage. The historical Stage 3 and Stage 4 transition runners and Word-export staging were removed from the public runner: the former are superseded by the bundled sealed Stage 4 reference and the latter is presentation-only. The default source locations are repository-relative `data/epd/` and `data/list_size/`, and SHA-256 recording accepts either `shasum` or `sha256sum`. These implementation changes do not alter source selection, processing, models, inference, scientific outputs or release gates. Each run records hashes for the coordinator and every module in `analysis_source_manifest.csv`.

## Study scope

The primary analysis covers January 2022 through December 2025. A descriptive window assessment additionally uses January–December 2021. The source manifest therefore records:

- 60 monthly English Prescribing Dataset (EPD) archives, January 2021–December 2025;
- 20 quarterly *Patients Registered at a GP Practice* files, January 2021–October 2025;
- one May 2025 BNF hierarchy lookup; and
- one dated England-and-Wales bank-holiday snapshot.

The class analysis is the primary inferential family. The all-eligible-drug analysis is secondary and exploratory. Benjamini–Hochberg false-discovery-rate control is applied separately within the 220 eligible classes and 974 eligible drugs; no hierarchy-wide error-rate claim is made.

## Source acquisition

### English Prescribing Dataset

- Publisher: NHS Business Services Authority
- Dataset: English Prescribing Dataset with SNOMED code
- Source: https://opendata.nhsbsa.net/dataset/english-prescribing-dataset-epd-with-snomed-code
- Licence: Open Government Licence v3.0

Download the 60 monthly ZIP archives from January 2021 through December 2025. Preserve the publisher filenames:

```text
EPD_SNOMED_202101.ZIP
...
EPD_SNOMED_202512.ZIP
```

Place them in `data/epd/`, or point `STOCK2026_EPD_DIR` to another directory.

### Registered-patient denominators

- Publisher: NHS England
- Dataset: Patients Registered at a GP Practice
- Source: https://digital.nhs.uk/data-and-information/publications/statistical/patients-registered-at-a-gp-practice
- Licence: Open Government Licence v3.0

Download the all-practice CSV for January, April, July and October in each year from 2021 through 2025. Preserve the filenames:

```text
gp-reg-pat-prac-all_202101.csv
gp-reg-pat-prac-all_202104.csv
...
gp-reg-pat-prac-all_202510.csv
```

Place them in `data/list_size/`, or point `STOCK2026_LIST_SIZE_DIR` to another directory.

The complete URLs, expected sizes, SHA-256 hashes, schemas and file roles are in `reproducibility/input_manifest.csv`. The manifest is the authority if a publisher page lists several similarly named resources.

## System requirements

- R 4.6.1 exactly. The analysis stops under another R version.
- A POSIX-like command line tested on macOS. Linux should also work when the named utilities are present.
- `bsdtar` from libarchive, required because the frozen EPD archives use ZIP method 6.3.
- `shasum` or `sha256sum` for provenance verification.
- Internet access for the one-time `renv` package restore.
- Sufficient compute for a large national dataset. The 60 EPD archives contain more than 300 GB of uncompressed CSV data. Monthly import checkpoints allow an interrupted import to resume without rereading completed months.

## Restore the R environment

From the repository root, install/use R 4.6.1 and run:

```sh
Rscript -e 'source("renv/activate.R"); renv::restore(prompt = FALSE)'
Rscript -e 'source("renv/activate.R"); renv::status()'
```

The accepted lockfile SHA-256 is:

```text
85a15983d9fc8ce84b93a0212734ca068e5f7cddc2fe459f725e5fdf7eb750b5
```

Do not update packages before reproduction. A changed R version or lockfile defines a different computational environment.

## Verify acquired inputs

If the data use the default `data/epd/` and `data/list_size/` locations:

```sh
Rscript --vanilla verify_inputs.R
```

For data stored elsewhere:

```sh
STOCK2026_EPD_DIR="/path/to/EPD" \
STOCK2026_LIST_SIZE_DIR="/path/to/GP_ListSize" \
Rscript --vanilla verify_inputs.R
```

The verifier checks all 82 required files, including the bundled BNF and calendar snapshots. It hashes only files whose byte size first matches the manifest. Do not start the expensive analysis if verification fails.

To check only the two bundled version-specific reference inputs while reviewing the repository:

```sh
Rscript --vanilla verify_inputs.R --bundled-only
```

## Check setup without importing data

After restoring packages and verifying the sources:

```sh
STOCK2026_SETUP_ONLY=true Rscript analysis_main.R
```

This checks the R version, active `renv` project, package interfaces, all declared input paths, `bsdtar`, and the output location without importing the monthly EPD records.

## Reproduce the complete analysis

Run the accepted stages in this order from the repository root:

```sh
STOCK2026_RUN_STAGE=stage5_trend Rscript analysis_main.R
STOCK2026_RUN_STAGE=stage5_hac Rscript analysis_main.R
STOCK2026_RUN_STAGE=stage5_working_days Rscript analysis_main.R
STOCK2026_RUN_STAGE=stage5_threshold Rscript analysis_main.R
STOCK2026_RUN_STAGE=stage6 Rscript analysis_main.R
STOCK2026_RUN_STAGE=full Rscript analysis_main.R
```

When source data are outside the default directories, prefix every command with the two source-directory environment variables shown above. `STOCK2026_OUTPUT_DIR` may be used for a single run, but the complete ordered workflow should use the default `Outputs/canonical` location because later stages verify earlier snapshots there.

The first command performs the raw import and creates per-month checkpoints. Subsequent commands reuse those checkpoints, refit the primary models, and add the uniform-HAC, working-day, threshold and targeted-diagnostic results. The final `full` route creates the complete R-generated reporting outputs. Word-table rendering is intentionally absent because it is presentation-only and is not required to reproduce any analysis result.

Every stage stops on a failed analytical or provenance gate. Warnings, package versions, session information, source flow, missingness, exclusions, model routes and stage manifests are written beneath `Outputs/canonical/qc/`.

## Expected reconciliation

A successful reproduction should recover:

| Quantity | Classes | Drugs |
|---|---:|---:|
| Observed after restriction to BNF chapters 1–14 | 344 | 2,155 |
| Eligible | 220 | 974 |
| Seasonality detected at 5% within-family BH FDR | 125 | 391 |
| Meaningfully seasonal | 30 | 88 |

The analysis should also recover 4,641,000,754 prescription items within BNF chapters 1–14. The machine-readable stage completion and QC files are authoritative; do not weaken a failed gate merely to obtain these headline totals.

## Data handling and Git safety

Primary data belong only in ignored local directories or outside the repository. Before any commit, verify the tracked inventory with:

```sh
git status --short
git ls-files
```

Do not use `git add -f` for ignored data or outputs. Avoid `git clean -fdx` when primary data are stored under `data/`, because it deletes ignored files. The repository intentionally contains no manuscript text beyond this README.

## Provenance and licensing

The source manifest records the exact frozen inputs used for the accepted run, including checksums and the absence of contemporaneous download dates where that information was not available. It does not claim that a current publisher download is identical until its hash has been verified.

The public source data and bundled public reference files are released under the Open Government Licence v3.0. No patient-level data were used or linked. An explicit software licence has not been assigned to the analysis code; add one before public release if downstream reuse should be licensed.
