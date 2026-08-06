#### SEASONALITY IN PRIMARY CARE PRESCRIBING IN ENGLAND PRIMARY CARE 
#### Nadine Stock, Gillian Carr, Islam Omar, Saran Shantikumar, July 2026

## Coordinator: runtime, inputs, provenance, and stage dispatch ----------------

.run_analysis <- function() {

# The canonical runtime is deliberately fixed. renv restores packages, not R.
canonical_r_version <- "4.6.1"
if (!identical(as.character(getRversion()), canonical_r_version) || nzchar(R.version$status)) {
  stop(
    "The canonical analysis requires the stable R ", canonical_r_version,
    " release; found ", R.version.string,
    if (nzchar(R.version$status)) paste0(" (status: ", R.version$status, ")") else "",
    ". Install/use the recorded R release, activate Analysis/renv, and rerun."
  )
}

# Resolve the repository root from the executed file, with getwd() as the
# interactive fallback. Modules use this path without changing directories.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
analysis_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

# Modules define cohesive pipeline stages; all run in this function's explicit
# context so the established object names and stage interface remain stable.
module_files <- c(
  "R/01_import_qc.R",
  "R/02_primary_analysis.R",
  "R/03_sensitivities.R",
  "R/04_diagnostics.R",
  "R/05_reporting.R"
)
for (module_file in module_files) {
  source(file.path(analysis_dir, module_file), local = environment())
}

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv is not available. From the repository root, run renv::restore() before the analysis.")
}
active_project <- tryCatch(renv::project(), error = function(e) "")
if (!nzchar(active_project) ||
    !identical(normalizePath(active_project, mustWork = TRUE), analysis_dir)) {
  stop("The repository renv project is not active. Start R/Rscript from the repository root and rerun.")
}

# Authoritative analysis constants. Later sensitivity specifications are added
# in their own implementation stage rather than silently changing these values.
study_start                  <- as.Date("2022-01-01")
study_end                    <- as.Date("2025-12-01")
window_start                 <- as.Date("2021-01-01")
window_end                   <- as.Date("2021-12-01")
eligibility_min_items_year   <- 1000L
trend_spline_df              <- 3L
harmonic_period_months       <- c(12L, 6L)
fdr_alpha                    <- 0.05
diagnostic_alpha             <- 0.05
amplitude_lci_threshold      <- 1.10
stl_strength_threshold       <- 0.50
coefficient_draw_seed        <- 1L
coefficient_draw_count       <- 2000L
plot_jitter_seed             <- 1L

required_packages <- c(
  "data.table", "dplyr", "stringr", "tidyr", "conflicted", "purrr",
  "lubridate", "readr", "tibble", "MASS", "AER", "sandwich", "lmtest",
  "car", "ggplot2", "ggrepel", "jsonlite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Packages are missing from the active renv library: ",
    paste(missing_packages, collapse = ", "),
    ". Run renv::restore() from Analysis; analysis_main.R never installs packages."
  )
}

# Explicit declarations make direct dependencies visible to renv::dependencies().
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(conflicted)
  library(purrr)
  library(lubridate)
  library(readr)
  library(tibble)
  library(MASS)
  library(AER)
  library(sandwich)
  library(lmtest)
  library(car)
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
})

conflicted::conflict_prefer("select", "dplyr", quiet = TRUE)
conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)
conflicted::conflict_prefer("first", "dplyr", quiet = TRUE)

# Machine-specific source locations may be overridden without editing code.
epd_dir <- normalizePath(
  Sys.getenv(
    "STOCK2026_EPD_DIR",
    unset = file.path(analysis_dir, "data", "epd")
  ),
  mustWork = FALSE
)
ls_dir <- normalizePath(
  Sys.getenv(
    "STOCK2026_LIST_SIZE_DIR",
    unset = file.path(analysis_dir, "data", "list_size")
  ),
  mustWork = FALSE
)
bnf_path <- normalizePath(
  Sys.getenv(
    "STOCK2026_BNF_PATH",
    unset = file.path(analysis_dir, "bnf_code_current_202505_version_88.csv")
  ),
  mustWork = FALSE
)
calendar_path <- file.path(
  analysis_dir, "reproducibility", "source_snapshots",
  "bank-holidays_2026-08-05.json"
)
out_dir <- normalizePath(
  Sys.getenv(
    "STOCK2026_OUTPUT_DIR",
    unset = file.path(analysis_dir, "Outputs", "canonical")
  ),
  mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- out_dir

# Each implementation stage has an explicit execution route. Runtime metadata
# are written inside that stage's QC directory so a later run cannot silently
# invalidate the provenance record of an already completed stage.
setup_only <- tolower(Sys.getenv("STOCK2026_SETUP_ONLY", unset = "false")) %in%
  c("1", "true", "yes")
run_stage <- tolower(Sys.getenv("STOCK2026_RUN_STAGE", unset = "full"))
valid_run_stages <- c(
  "full", "stage2", "stage5_trend", "stage5_hac",
  "stage5_working_days", "stage5_threshold", "stage6"
)
if (!run_stage %in% valid_run_stages) {
  stop("STOCK2026_RUN_STAGE must be one of: ", paste(valid_run_stages, collapse = ", "), ".")
}
provenance_stage <- if (setup_only) "stage1" else run_stage

expected_ym <- as.integer(format(seq(study_start, study_end, by = "month"), "%Y%m"))
window_ym <- as.integer(format(seq(window_start, window_end, by = "month"), "%Y%m"))
epd_files <- file.path(epd_dir, sprintf("EPD_SNOMED_%d.ZIP", expected_ym))
epd_2021_files <- file.path(epd_dir, sprintf("EPD_SNOMED_%d.ZIP", window_ym))
list_size_ym <- as.integer(format(seq(study_start, study_end, by = "3 months"), "%Y%m"))
list_size_files <- file.path(ls_dir, sprintf("gp-reg-pat-prac-all_%d.csv", list_size_ym))
list_size_2021_ym <- as.integer(format(seq(window_start, window_end, by = "3 months"), "%Y%m"))
list_size_2021_files <- file.path(
  ls_dir, sprintf("gp-reg-pat-prac-all_%d.csv", list_size_2021_ym)
)

# Validate all frozen inputs before any source archive is opened.
declared_inputs <- c(bnf_path, calendar_path, epd_files, epd_2021_files,
                     list_size_files, list_size_2021_files)
missing_inputs <- declared_inputs[!file.exists(declared_inputs)]
if (length(missing_inputs)) {
  stop(
    "Required frozen input file(s) are missing:\n",
    paste(missing_inputs, collapse = "\n"),
    "\nSet STOCK2026_EPD_DIR, STOCK2026_LIST_SIZE_DIR or STOCK2026_BNF_PATH if the mounted paths differ."
  )
}

archive_executable <- Sys.which("bsdtar")
if (!nzchar(archive_executable)) {
  stop("The system 'bsdtar' executable is required to read the ZIP 6.3-compressed frozen EPD archives.")
}
# Stream an EPD archive into memory without extracting it to disk.
read_epd_archive <- function(path, select = NULL, nrows = Inf, colClasses = NULL) {
  data.table::fread(
    cmd = paste(shQuote(archive_executable), "-xOf", shQuote(path)),
    select = select,
    nrows = nrows,
    colClasses = colClasses
  )
}
read_epd_header <- function(path) {
  # fread(cmd=) materialises command output before reading it. Limit this command
  # to one line so schema inspection does not expand each 6--8 GB CSV twice.
  data.table::fread(
    cmd = paste(shQuote(archive_executable), "-xOf", shQuote(path), "2>/dev/null | head -n 1"),
    nrows = 0L
  )
}

# Promote completed checkpoints atomically so partial files cannot be reused.
atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  data.table::fwrite(x, tmp)
  if (!file.rename(tmp, path)) stop("Could not promote completed checkpoint: ", path)
  invisible(path)
}

atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Could not promote completed checkpoint: ", path)
  invisible(path)
}

is_blank <- function(x) is.na(x) | trimws(as.character(x)) == ""

# Use either common system implementation to calculate SHA-256 provenance.
sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  shasum <- Sys.which("shasum")
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(shasum)) {
    ans <- system2(shasum, c("-a", "256", shQuote(path)), stdout = TRUE)
  } else if (nzchar(sha256sum)) {
    ans <- system2(sha256sum, shQuote(path), stdout = TRUE)
  } else {
    return(NA_character_)
  }
  sub("[[:space:]].*$", "", ans[1])
}

# Record the exact runtime and source-code hashes for this execution route.
repro_dir <- file.path(out_dir, "qc", provenance_stage, "reproducibility")
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)
script_path <- file.path(analysis_dir, "analysis_main.R")
lock_path <- file.path(analysis_dir, "renv.lock")
analysis_source_paths <- file.path(analysis_dir, c("analysis_main.R", module_files))
analysis_source_manifest <- data.frame(
  relative_path = c("analysis_main.R", module_files),
  sha256 = vapply(analysis_source_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
package_versions <- data.frame(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(package_versions, file.path(repro_dir, "package_versions.csv"), row.names = FALSE)
utils::write.csv(
  analysis_source_manifest,
  file.path(repro_dir, "analysis_source_manifest.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(repro_dir, "sessionInfo.txt"))
writeLines(c(
  paste("run_started", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "="),
  paste("run_stage", provenance_stage, sep = "="),
  paste("R.version.string", R.version.string, sep = "="),
  paste("R.status", R.version$status, sep = "="),
  paste("R.platform", R.version$platform, sep = "="),
  paste("R.arch", R.version$arch, sep = "="),
  paste("R.executable", file.path(R.home("bin"), "R"), sep = "="),
  paste("operating_system", paste(Sys.info(), collapse = "; "), sep = "="),
  paste("analysis_script_sha256", sha256_file(script_path), sep = "="),
  paste("analysis_source_files", nrow(analysis_source_manifest), sep = "="),
  paste("renv_lock_sha256", sha256_file(lock_path), sep = "="),
  paste("coefficient_draw_seed", coefficient_draw_seed, sep = "="),
  paste("coefficient_draw_count", coefficient_draw_count, sep = "="),
  paste("plot_jitter_seed", plot_jitter_seed, sep = "="),
  paste("epd_dir", epd_dir, sep = "="),
  paste("list_size_dir", ls_dir, sep = "="),
  paste("bnf_path", bnf_path, sep = "="),
  paste("calendar_path", calendar_path, sep = "="),
  paste("archive_reader", archive_executable, sep = "="),
  paste("archive_reader_version", system2(archive_executable, "--version", stdout = TRUE)[1], sep = "="),
  paste("output_dir", out_dir, sep = "=")
), file.path(repro_dir, "run_manifest.txt"))

warnings_seen <- character()
on.exit({
  warning_table <- data.frame(
    warning = if (length(warnings_seen)) unique(warnings_seen) else character(),
    stringsAsFactors = FALSE
  )
  utils::write.csv(warning_table, file.path(repro_dir, "warnings.csv"), row.names = FALSE)
}, add = TRUE)

if (setup_only) {
  message(
    "Stage 1 setup check passed: ", R.version.string,
    "; active renv project ", analysis_dir,
    "; ", length(declared_inputs), " declared frozen inputs present."
  )
  return(invisible(TRUE))
}

withCallingHandlers({

# Import/QC may stop cleanly after a requested month or the Stage 2 route.
if (!run_import_qc(environment())) return(invisible(TRUE))
run_primary_analysis(environment())

# Sensitivity and diagnostic routes stop after writing their own gated outputs.
if (run_stage == "stage5_trend") {
  run_trend_sensitivity(environment())
  return(invisible(TRUE))
}
if (run_stage == "stage5_hac") {
  run_hac_sensitivity(environment())
  return(invisible(TRUE))
}
if (run_stage == "stage5_working_days") {
  run_working_days_sensitivity(environment())
  return(invisible(TRUE))
}
if (run_stage == "stage5_threshold") {
  run_threshold_sensitivity(environment())
  return(invisible(TRUE))
}
if (run_stage == "stage6") {
  run_diagnostics(environment())
  return(invisible(TRUE))
}

# Only the full route continues into publication and window-assessment figures.
run_reporting(environment())

}, warning = function(w) {
  warning_text <- conditionMessage(w)
  warnings_seen <<- c(warnings_seen, warning_text)
  message("WARNING: ", warning_text)
  invokeRestart("muffleWarning")
})

invisible(TRUE)
}

.run_analysis()
