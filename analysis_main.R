#### SEASONALITY IN PRIMARY CARE PRESCRIBING IN ENGLAND PRIMARY CARE 
#### Nadine Stock, Gillian Carr, Islam Omar, Saran Shantikumar, July 2026

### 1. Set up ----------

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

# Resolve the Analysis project root from the executed file, with getwd() as an
# interactive fallback. No later section changes the working directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
analysis_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv is not available. From Analysis, run renv::restore() before the analysis.")
}
active_project <- tryCatch(renv::project(), error = function(e) "")
if (!nzchar(active_project) ||
    !identical(normalizePath(active_project, mustWork = TRUE), analysis_dir)) {
  stop("The Analysis renv project is not active. Start R/Rscript from the Analysis directory and rerun.")
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
  "full", "stage2", "stage4", "stage5_trend", "stage5_hac",
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

repro_dir <- file.path(out_dir, "qc", provenance_stage, "reproducibility")
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)
script_path <- file.path(analysis_dir, "analysis_main.R")
lock_path <- file.path(analysis_dir, "renv.lock")
package_versions <- data.frame(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(package_versions, file.path(repro_dir, "package_versions.csv"), row.names = FALSE)
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


### 2. Get data -------------

## BNF drug class lookup table

# Retrived from https://opendata.nhsbsa.net/dataset/bnf-code-information-current-year (May 2025 version)

bnf_ref <- fread(bnf_path, colClasses = "character") %>% # Read in and force character to avoid dropping leading zeros
  select(BNF_SECTION,
         BNF_SECTION_CODE,
         BNF_PARAGRAPH,
         BNF_PARAGRAPH_CODE) %>% # Keep required columns
  rename(BNF_SECTION_NAME = BNF_SECTION,
         BNF_CLASS_NAME = BNF_PARAGRAPH,
         BNF_CLASS_CODE = BNF_PARAGRAPH_CODE) %>% # Rename
  distinct() %>% # Remove duplicate rows
  rename_with(tolower) # lower case column headers
bnf_ref <- as.data.table(bnf_ref)
if (anyDuplicated(bnf_ref$bnf_class_code)) stop("BNF reference has duplicate class-code rows.")


## EPD data. The 60 frozen archives are read once. Per-month aggregates and QC
## files are checkpoints, so an interrupted run resumes without rereading
## completed 6--8 GB uncompressed months.

stage2_dir <- file.path(out_dir, "qc", "stage2")
epd_qc_month_dir <- file.path(stage2_dir, "epd_monthly")
cache_dir_2021 <- file.path(data_dir, "cache_2021")
for (d in c(stage2_dir, epd_qc_month_dir, cache_dir_2021)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

valid_chapters <- sprintf("%02d", 1:14)
all_epd_ym <- c(expected_ym, window_ym)
all_epd_files <- c(epd_files, epd_2021_files)
requested_months_raw <- trimws(Sys.getenv("STOCK2026_IMPORT_MONTHS", unset = ""))
partial_epd_run <- nzchar(requested_months_raw)
if (partial_epd_run) {
  requested_months <- suppressWarnings(as.integer(
    trimws(strsplit(requested_months_raw, ",", fixed = TRUE)[[1]])
  ))
  if (anyNA(requested_months) || !all(requested_months %in% all_epd_ym)) {
    stop("STOCK2026_IMPORT_MONTHS must be a comma-separated subset of 202101--202512 available to this study.")
  }
  import_index <- match(requested_months, all_epd_ym)
} else {
  import_index <- seq_along(all_epd_ym)
}

epd_paths_for_month <- function(ym) {
  is_main <- ym %in% expected_ym
  list(
    qc = file.path(epd_qc_month_dir, sprintf("epd_file_qc_%d.csv", ym)),
    mapping = file.path(epd_qc_month_dir, sprintf("epd_drug_mapping_%d.csv", ym)),
    duplicates = file.path(epd_qc_month_dir, sprintf("epd_candidate_duplicates_%d.csv", ym)),
    unmatched = file.path(epd_qc_month_dir, sprintf("epd_unmatched_lookup_%d.csv", ym)),
    drug = if (is_main) file.path(out_dir, sprintf("data_bydrug_%d.csv", ym)) else NA_character_,
    class = if (is_main) {
      file.path(out_dir, sprintf("data_byclass_%d.csv", ym))
    } else {
      file.path(cache_dir_2021, sprintf("data_byclass_2021_%d.csv", ym))
    }
  )
}

for (ii in seq_along(import_index)) {
  i <- import_index[ii]
  ym <- all_epd_ym[i]
  epd_file <- all_epd_files[i]
  role <- if (ym %in% expected_ym) "main_2022_2025" else "descriptive_2021"
  paths <- epd_paths_for_month(ym)
  required_checkpoints <- c(paths$qc, paths$mapping, paths$duplicates,
                            paths$unmatched, paths$class)
  if (role == "main_2022_2025") required_checkpoints <- c(required_checkpoints, paths$drug)
  qc_passed <- FALSE
  if (file.exists(paths$qc)) {
    prior_qc <- tryCatch(fread(paths$qc), error = function(e) NULL)
    qc_passed <- !is.null(prior_qc) && nrow(prior_qc) == 1L &&
      identical(prior_qc$qc_status[1], "PASS")
  }
  if (all(file.exists(required_checkpoints)) && qc_passed) {
    message("Skipping EPD ", ym, " (complete Stage 2 checkpoint)")
    next
  }

  message("Processing EPD ", ym, " (", ii, " of ", length(import_index), ")")
  started <- proc.time()[["elapsed"]]
  hdr <- names(read_epd_header(epd_file))

  legacy_fields <- c(
    "YEAR_MONTH", "PRACTICE_CODE", "CHEMICAL_SUBSTANCE_BNF_DESCR",
    "BNF_CHEMICAL_SUBSTANCE", "BNF_CODE", "BNF_CHAPTER_PLUS_CODE",
    "ITEMS", "SNOMED_CODE"
  )
  current_fields <- c(
    "YEAR_MONTH", "PRACTICE_CODE", "BNF_CHEMICAL_SUBSTANCE",
    "BNF_CHEMICAL_SUBSTANCE_CODE", "BNF_PRESENTATION_CODE",
    "BNF_CHAPTER_PLUS_CODE", "ITEMS", "SNOMED_CODE"
  )
  is_legacy <- all(legacy_fields %in% hdr) && !"BNF_PRESENTATION_CODE" %in% hdr
  is_current <- all(current_fields %in% hdr) && "BNF_PRESENTATION_CODE" %in% hdr
  if (is_legacy == is_current) stop("Unrecognised or ambiguous EPD schema in ", basename(epd_file))
  schema_variant <- if (is_legacy) "legacy_through_202502" else "snomed_from_202503"
  expected_schema <- if (ym <= 202502) "legacy_through_202502" else "snomed_from_202503"
  if (!identical(schema_variant, expected_schema)) {
    stop("EPD schema does not match the declared boundary for ", ym, ".")
  }

  source_fields <- if (is_legacy) legacy_fields else current_fields
  data <- read_epd_archive(
    epd_file,
    select = source_fields,
    colClasses = list(character = setdiff(source_fields, "ITEMS"))
  )
  canonical_fields <- c(
    "year_month_raw", "practice_code", "bnf_drug_name_raw",
    "bnf_chemical_code_raw", "bnf_presentation_code_raw",
    "bnf_chapter_plus_raw", "items_raw", "snomed_code_raw"
  )
  setnames(data, source_fields, canonical_fields)
  for (nm in setdiff(canonical_fields, "items_raw")) {
    set(data, j = nm, value = trimws(as.character(data[[nm]])))
  }

  items_source_type <- typeof(data$items_raw)
  item_was_blank <- is_blank(data$items_raw)
  items_numeric <- suppressWarnings(as.numeric(data$items_raw))
  item_nonnumeric <- sum(!item_was_blank & is.na(items_numeric))
  set(data, j = "items_raw", value = items_numeric)
  data[, year_month_normalised := gsub("-", "", year_month_raw, fixed = TRUE)]
  data[, `:=`(
    bnf_chapter_code = substr(bnf_presentation_code_raw, 1, 2),
    bnf_class_code = substr(bnf_presentation_code_raw, 1, 6),
    bnf_drug_code = substr(bnf_presentation_code_raw, 1, 9),
    chapter_plus_code = substr(bnf_chapter_plus_raw, 1, 2),
    bnf_chapter_name = substring(bnf_chapter_plus_raw, 5)
  )]

  candidate_key <- c("year_month_normalised", "practice_code", "bnf_presentation_code_raw")
  extended_key <- c(candidate_key, "bnf_chemical_code_raw", "snomed_code_raw")
  candidate_duplicate_excess <- sum(duplicated(data, by = candidate_key))
  extended_duplicate_excess <- sum(duplicated(data, by = extended_key))
  duplicate_groups <- 0L
  duplicate_max_records <- 1L
  duplicate_detail <- data.table(
    year_month = character(), practice_code = character(),
    bnf_presentation_code_raw = character(), records = integer(),
    n_snomed_codes = integer(), n_chemical_codes = integer(), items = numeric()
  )
  if (candidate_duplicate_excess > 0L) {
    dup_mask <- duplicated(data, by = candidate_key) |
      duplicated(data, by = candidate_key, fromLast = TRUE)
    duplicate_detail <- data[dup_mask, .(
      records = .N,
      n_snomed_codes = uniqueN(snomed_code_raw),
      n_chemical_codes = uniqueN(bnf_chemical_code_raw),
      items = sum(items_raw, na.rm = TRUE)
    ), by = .(
      year_month = year_month_normalised,
      practice_code,
      bnf_presentation_code_raw
    )][order(-records, practice_code, bnf_presentation_code_raw)]
    duplicate_groups <- nrow(duplicate_detail)
    duplicate_max_records <- max(duplicate_detail$records)
  }

  month_mismatch_rows <- data[year_month_normalised != sprintf("%d", ym) |
                                is_blank(year_month_normalised), .N]
  missing_practice_rows <- data[is_blank(practice_code), .N]
  missing_presentation_rows <- data[is_blank(bnf_presentation_code_raw), .N]
  short_presentation_rows <- data[!is_blank(bnf_presentation_code_raw) &
                                    nchar(bnf_presentation_code_raw) < 9L, .N]
  missing_drug_name_rows <- data[is_blank(bnf_drug_name_raw), .N]
  missing_chemical_code_rows <- data[is_blank(bnf_chemical_code_raw), .N]
  missing_snomed_rows <- data[is_blank(snomed_code_raw), .N]
  missing_chapter_plus_rows <- data[is_blank(bnf_chapter_plus_raw), .N]
  chapter_prefix_mismatch_rows <- data[
    !is_blank(bnf_chapter_plus_raw) & bnf_chapter_code != chapter_plus_code, .N
  ]
  chemical_code_mismatch_all_rows <- data[
    !is_blank(bnf_chemical_code_raw) & bnf_drug_code != bnf_chemical_code_raw, .N
  ]
  chemical_code_mismatch_all_items <- data[
    !is_blank(bnf_chemical_code_raw) & bnf_drug_code != bnf_chemical_code_raw,
    sum(items_raw, na.rm = TRUE)
  ]
  chemical_code_mismatch_in_scope_rows <- data[
    bnf_chapter_code %chin% valid_chapters & !is_blank(bnf_chemical_code_raw) &
      bnf_drug_code != bnf_chemical_code_raw, .N
  ]
  chemical_code_mismatch_in_scope_items <- data[
    bnf_chapter_code %chin% valid_chapters & !is_blank(bnf_chemical_code_raw) &
      bnf_drug_code != bnf_chemical_code_raw, sum(items_raw, na.rm = TRUE)
  ]
  missing_items_rows <- sum(is.na(data$items_raw))
  negative_items_rows <- data[!is.na(items_raw) & items_raw < 0, .N]
  noninteger_items_rows <- data[
    !is.na(items_raw) & abs(items_raw - round(items_raw)) > sqrt(.Machine$double.eps), .N
  ]
  zero_items_rows <- data[!is.na(items_raw) & items_raw == 0, .N]

  qc_row <- data.table(
    source_id = sprintf("EPD_%d", ym), year_month = ym, analytical_role = role,
    source_file = epd_file, schema_variant = schema_variant,
    expected_schema = expected_schema, raw_row_count = nrow(data),
    raw_column_count = length(hdr), source_columns = paste(hdr, collapse = ";"),
    items_source_type = items_source_type,
    candidate_key = paste(candidate_key, collapse = " + "),
    candidate_duplicate_excess_rows = candidate_duplicate_excess,
    candidate_duplicate_groups = duplicate_groups,
    candidate_duplicate_max_records = duplicate_max_records,
    extended_key = paste(extended_key, collapse = " + "),
    extended_duplicate_excess_rows = extended_duplicate_excess,
    month_mismatch_rows = month_mismatch_rows,
    missing_practice_rows = missing_practice_rows,
    missing_presentation_rows = missing_presentation_rows,
    short_presentation_rows = short_presentation_rows,
    missing_drug_name_rows = missing_drug_name_rows,
    missing_chemical_code_rows = missing_chemical_code_rows,
    missing_snomed_rows = missing_snomed_rows,
    missing_chapter_plus_rows = missing_chapter_plus_rows,
    chapter_prefix_mismatch_rows = chapter_prefix_mismatch_rows,
    chemical_code_mismatch_all_rows = chemical_code_mismatch_all_rows,
    chemical_code_mismatch_all_items = chemical_code_mismatch_all_items,
    chemical_code_mismatch_in_scope_rows = chemical_code_mismatch_in_scope_rows,
    chemical_code_mismatch_in_scope_items = chemical_code_mismatch_in_scope_items,
    missing_items_rows = missing_items_rows,
    nonnumeric_items_rows = item_nonnumeric,
    negative_items_rows = negative_items_rows,
    noninteger_items_rows = noninteger_items_rows,
    zero_items_rows = zero_items_rows
  )

  fatal_counts <- c(
    month_mismatch_rows = month_mismatch_rows,
    missing_practice_rows = missing_practice_rows,
    missing_presentation_rows = missing_presentation_rows,
    short_presentation_rows = short_presentation_rows,
    missing_drug_name_rows = missing_drug_name_rows,
    missing_chemical_code_rows = missing_chemical_code_rows,
    chapter_prefix_mismatch_rows = chapter_prefix_mismatch_rows,
    chemical_code_mismatch_in_scope_rows = chemical_code_mismatch_in_scope_rows,
    missing_items_rows = missing_items_rows,
    nonnumeric_items_rows = item_nonnumeric,
    negative_items_rows = negative_items_rows,
    noninteger_items_rows = noninteger_items_rows
  )
  failed_names <- names(fatal_counts)[fatal_counts > 0]
  if (length(failed_names)) {
    qc_row[, `:=`(
      qc_status = "FAIL",
      failure_reason = paste(paste0(failed_names, "=", fatal_counts[failed_names]), collapse = "; "),
      elapsed_seconds = proc.time()[["elapsed"]] - started
    )]
    atomic_fwrite(qc_row, paths$qc)
    stop("EPD ", ym, " failed source validation: ", qc_row$failure_reason)
  }

  raw_total_items <- data[, sum(items_raw)]
  retained <- data[bnf_chapter_code %chin% valid_chapters]
  retained_items <- retained[, sum(items_raw)]
  removed_rows <- nrow(data) - nrow(retained)
  removed_items <- raw_total_items - retained_items

  mapping <- retained[, .(
    raw_records = .N,
    items = sum(items_raw)
  ), by = .(
    bnf_drug_code, bnf_chemical_code_raw, bnf_drug_name_raw
  )]
  mapping[, year_month := as.integer(ym)]
  setcolorder(mapping, c(
    "year_month", "bnf_drug_code", "bnf_chemical_code_raw",
    "bnf_drug_name_raw", "raw_records", "items"
  ))
  mapping[, chemical_code_match := bnf_drug_code == bnf_chemical_code_raw]

  name_choice <- mapping[, .(name_items = sum(items)),
                         by = .(bnf_drug_code, bnf_drug_name_raw)]
  setorder(name_choice, bnf_drug_code, -name_items, bnf_drug_name_raw)
  name_choice <- name_choice[, .SD[1L], by = bnf_drug_code][,
    .(bnf_drug_code, bnf_drug_name = bnf_drug_name_raw)
  ]
  chapter_choice <- retained[, .(chapter_items = sum(items_raw)),
                             by = .(bnf_chapter_code, bnf_chapter_name)]
  setorder(chapter_choice, bnf_chapter_code, -chapter_items, bnf_chapter_name)
  chapter_choice <- chapter_choice[, .SD[1L], by = bnf_chapter_code][,
    .(bnf_chapter_code, bnf_chapter_name)
  ]

  drug_agg <- retained[, .(items = sum(items_raw)), by = .(
    bnf_chapter_code, bnf_class_code, bnf_drug_code
  )]
  drug_agg[, year_month := as.integer(ym)]
  drug_agg <- merge(drug_agg, chapter_choice, by = "bnf_chapter_code", all.x = TRUE, sort = FALSE)
  drug_agg <- merge(drug_agg, bnf_ref, by = "bnf_class_code", all.x = TRUE, sort = FALSE)
  drug_agg <- merge(drug_agg, name_choice, by = "bnf_drug_code", all.x = TRUE, sort = FALSE)
  setcolorder(drug_agg, c(
    "year_month", "bnf_chapter_code", "bnf_chapter_name",
    "bnf_section_code", "bnf_section_name", "bnf_class_code",
    "bnf_class_name", "bnf_drug_code", "bnf_drug_name", "items"
  ))
  setorder(drug_agg, year_month, bnf_drug_code)
  drug_duplicate_rows <- sum(duplicated(drug_agg, by = c("year_month", "bnf_drug_code")))

  class_agg <- drug_agg[, .(items = sum(items)), by = .(
    year_month, bnf_chapter_code, bnf_chapter_name,
    bnf_section_code, bnf_section_name, bnf_class_code, bnf_class_name
  )]
  setorder(class_agg, year_month, bnf_class_code)
  class_duplicate_rows <- sum(duplicated(class_agg, by = c("year_month", "bnf_class_code")))
  class_drug_item_difference <- class_agg[, sum(items)] - drug_agg[, sum(items)]

  unmatched <- class_agg[is.na(bnf_class_name), .(
    year_month, bnf_chapter_code, bnf_class_code, items
  )]
  unmatched_items <- unmatched[, sum(items)]
  unmatched_rows <- if (nrow(unmatched)) {
    retained[bnf_class_code %chin% unmatched$bnf_class_code, .N]
  } else 0L

  qc_row[, `:=`(
    raw_total_items = raw_total_items,
    retained_chapter_01_14_rows = nrow(retained),
    retained_chapter_01_14_items = retained_items,
    removed_chapter_rows = removed_rows,
    removed_chapter_items = removed_items,
    removed_chapter_item_share = if (raw_total_items > 0) removed_items / raw_total_items else NA_real_,
    unmatched_lookup_raw_rows = unmatched_rows,
    unmatched_lookup_items = unmatched_items,
    unmatched_lookup_item_share = if (retained_items > 0) unmatched_items / retained_items else NA_real_,
    national_drug_rows = nrow(drug_agg),
    national_class_rows = nrow(class_agg),
    drug_month_duplicate_rows = drug_duplicate_rows,
    class_month_duplicate_rows = class_duplicate_rows,
    class_drug_item_difference = class_drug_item_difference,
    qc_status = if (drug_duplicate_rows == 0L && class_duplicate_rows == 0L &&
                     class_drug_item_difference == 0) "PASS" else "FAIL",
    failure_reason = if (drug_duplicate_rows == 0L && class_duplicate_rows == 0L &&
                          class_drug_item_difference == 0) "" else
      "post-aggregation duplicate or class/drug reconciliation failure",
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )]

  atomic_fwrite(mapping, paths$mapping)
  atomic_fwrite(head(duplicate_detail, 100L), paths$duplicates)
  atomic_fwrite(unmatched, paths$unmatched)
  if (qc_row$qc_status != "PASS") {
    atomic_fwrite(qc_row, paths$qc)
    stop("EPD ", ym, " failed post-aggregation validation.")
  }
  if (role == "main_2022_2025") atomic_fwrite(drug_agg, paths$drug)
  atomic_fwrite(class_agg, paths$class)
  atomic_fwrite(qc_row, paths$qc)
  message(
    "Completed EPD ", ym, ": ", format(nrow(data), big.mark = ","),
    " raw rows in ", round(qc_row$elapsed_seconds, 1), " seconds"
  )
  rm(data, retained, mapping, name_choice, chapter_choice,
     drug_agg, class_agg, unmatched, duplicate_detail)
  invisible(gc())
}

if (partial_epd_run) {
  message("Requested EPD checkpoint month(s) completed; stopping before cross-month Stage 2 checks.")
  return(invisible(TRUE))
}

epd_qc_paths <- file.path(epd_qc_month_dir, sprintf("epd_file_qc_%d.csv", all_epd_ym))
epd_mapping_paths <- file.path(epd_qc_month_dir, sprintf("epd_drug_mapping_%d.csv", all_epd_ym))
if (any(!file.exists(c(epd_qc_paths, epd_mapping_paths)))) {
  stop("Not all 60 EPD Stage 2 checkpoints are present.")
}
epd_file_qc <- rbindlist(lapply(epd_qc_paths, fread), use.names = TRUE, fill = TRUE)
setorder(epd_file_qc, year_month)
atomic_fwrite(epd_file_qc, file.path(stage2_dir, "epd_file_qc.csv"))

epd_unmatched_paths <- file.path(
  epd_qc_month_dir, sprintf("epd_unmatched_lookup_%d.csv", all_epd_ym)
)
epd_unmatched_lookup_qc <- rbindlist(lapply(epd_unmatched_paths, function(f) {
  fread(f, colClasses = list(character = c("bnf_chapter_code", "bnf_class_code")))
}), use.names = TRUE, fill = TRUE)
setorder(epd_unmatched_lookup_qc, year_month, bnf_class_code)
unmatched_detail_items <- epd_unmatched_lookup_qc[, sum(items)]
unmatched_summary_items <- epd_file_qc[, sum(unmatched_lookup_items)]
if (unmatched_detail_items != unmatched_summary_items) {
  stop("Per-code unmatched BNF lookup detail does not reconcile with the per-file QC totals.")
}
atomic_fwrite(
  epd_unmatched_lookup_qc,
  file.path(stage2_dir, "epd_unmatched_lookup_qc.csv")
)

mapping_all <- rbindlist(lapply(epd_mapping_paths, function(f) {
  fread(f, colClasses = list(character = c(
    "bnf_drug_code", "bnf_chemical_code_raw", "bnf_drug_name_raw"
  )))
}), use.names = TRUE)
mapping_all[, normalised_drug_name := bnf_drug_name_raw |>
  str_to_lower() |>
  str_replace_all("[[:punct:]]", " ") |>
  str_replace_all("\\s+", " ") |>
  str_trim()]

drug_code_name_qc <- mapping_all[, .(
  first_month = min(year_month), last_month = max(year_month),
  n_months = uniqueN(year_month),
  n_raw_chemical_codes = uniqueN(bnf_chemical_code_raw),
  n_raw_names = uniqueN(bnf_drug_name_raw),
  raw_chemical_codes = paste(sort(unique(bnf_chemical_code_raw)), collapse = " | "),
  raw_names = paste(sort(unique(bnf_drug_name_raw)), collapse = " | "),
  raw_records = sum(raw_records), items = sum(items),
  chemical_mismatch_items = sum(items[!chemical_code_match])
), by = bnf_drug_code][order(-items)]
atomic_fwrite(drug_code_name_qc, file.path(stage2_dir, "epd_drug_code_name_qc.csv"))

chemical_code_qc <- mapping_all[, .(
  first_month = min(year_month), last_month = max(year_month),
  n_derived_drug_codes = uniqueN(bnf_drug_code),
  n_raw_names = uniqueN(bnf_drug_name_raw),
  derived_drug_codes = paste(sort(unique(bnf_drug_code)), collapse = " | "),
  raw_names = paste(sort(unique(bnf_drug_name_raw)), collapse = " | "),
  raw_records = sum(raw_records), items = sum(items)
), by = bnf_chemical_code_raw][order(-items)]
atomic_fwrite(chemical_code_qc, file.path(stage2_dir, "epd_chemical_code_qc.csv"))

name_code_qc <- mapping_all[, .(
  first_month = min(year_month), last_month = max(year_month),
  n_derived_drug_codes = uniqueN(bnf_drug_code),
  derived_drug_codes = paste(sort(unique(bnf_drug_code)), collapse = " | "),
  raw_names = paste(sort(unique(bnf_drug_name_raw)), collapse = " | "),
  raw_records = sum(raw_records), items = sum(items)
), by = normalised_drug_name][n_derived_drug_codes > 1L][order(-items)]
atomic_fwrite(name_code_qc, file.path(stage2_dir, "epd_normalised_name_code_qc.csv"))


## List size data

# Explicitly selected quarterly files from the England all-practice extracts.
all_list_size_files <- c(list_size_files, list_size_2021_files)
all_list_size_ym <- c(list_size_ym, list_size_2021_ym)
list_size_source_qc <- vector("list", length(all_list_size_files))
list_size_observed <- vector("list", length(all_list_size_files))

for (i in seq_along(all_list_size_files)) {
  f <- all_list_size_files[i]
  ym <- all_list_size_ym[i]
  hdr <- names(fread(f, nrows = 0L))
  required <- c(
    "PUBLICATION", "EXTRACT_DATE", "TYPE", "CODE", "SEX", "AGE",
    "NUMBER_OF_PATIENTS"
  )
  missing_columns <- setdiff(required, hdr)
  if (length(missing_columns)) {
    stop(basename(f), ": denominator column(s) missing: ", paste(missing_columns, collapse = ", "))
  }
  schema_variant <- if ("CCG_CODE" %in% hdr) {
    "legacy_CCG_fields"
  } else if ("SUB_ICB_LOC_CODE" %in% hdr) {
    "transition_SubICB_LOC_fields"
  } else if ("SUB_ICB_LOCATION_CODE" %in% hdr) {
    "current_SubICB_fields"
  } else {
    stop(basename(f), ": unrecognised denominator schema.")
  }

  d <- fread(f, colClasses = list(character = setdiff(required, "NUMBER_OF_PATIENTS")))
  patient_source_type <- typeof(d$NUMBER_OF_PATIENTS)
  patients <- suppressWarnings(as.numeric(d$NUMBER_OF_PATIENTS))
  nonnumeric_patients <- sum(!is_blank(d$NUMBER_OF_PATIENTS) & is.na(patients))
  set(d, j = "NUMBER_OF_PATIENTS", value = patients)
  extract_values <- unique(d$EXTRACT_DATE)
  parsed_extract <- as.Date(
    extract_values,
    tryFormats = c("%Y-%m-%d", "%d%b%Y", "%d-%b-%y", "%d-%b-%Y")
  )
  expected_extract <- as.Date(paste0(ym, "01"), "%Y%m%d")

  duplicate_practice_codes <- sum(duplicated(d$CODE))
  missing_practice_codes <- sum(is_blank(d$CODE))
  missing_patient_counts <- sum(is.na(d$NUMBER_OF_PATIENTS))
  negative_patient_counts <- sum(d$NUMBER_OF_PATIENTS < 0, na.rm = TRUE)
  noninteger_patient_counts <- sum(
    abs(d$NUMBER_OF_PATIENTS - round(d$NUMBER_OF_PATIENTS)) > sqrt(.Machine$double.eps),
    na.rm = TRUE
  )
  publication_ok <- identical(sort(unique(d$PUBLICATION)), "GP_PRAC_PAT_LIST")
  type_ok <- identical(sort(unique(d$TYPE)), "GP")
  sex_ok <- identical(sort(unique(d$SEX)), "ALL")
  age_ok <- identical(sort(unique(d$AGE)), "ALL")
  extract_date_ok <- length(parsed_extract) == 1L && !is.na(parsed_extract) &&
    identical(parsed_extract, expected_extract)
  national_total <- sum(d$NUMBER_OF_PATIENTS, na.rm = TRUE)
  plausible_total <- national_total >= 55e6 && national_total <= 70e6

  fail_counts <- c(
    duplicate_practice_codes = duplicate_practice_codes,
    missing_practice_codes = missing_practice_codes,
    missing_patient_counts = missing_patient_counts,
    nonnumeric_patient_counts = nonnumeric_patients,
    negative_patient_counts = negative_patient_counts,
    noninteger_patient_counts = noninteger_patient_counts
  )
  failed <- names(fail_counts)[fail_counts > 0]
  metadata_failures <- c(
    if (!publication_ok) "publication_scope",
    if (!type_ok) "practice_type",
    if (!sex_ok) "sex_scope",
    if (!age_ok) "age_scope",
    if (!extract_date_ok) "extract_date",
    if (!plausible_total) "national_total_range"
  )
  failed_count_labels <- if (length(failed)) {
    paste0(failed, "=", fail_counts[failed])
  } else {
    character()
  }
  failure_reason <- paste(c(failed_count_labels, metadata_failures), collapse = "; ")
  status <- if (nzchar(failure_reason)) "FAIL" else "PASS"

  list_size_source_qc[[i]] <- data.table(
    source_id = sprintf("LIST_SIZE_%d", ym), year_month = ym,
    analytical_role = if (ym %in% list_size_ym) "main_2022_2025" else "descriptive_2021",
    source_file = f, schema_variant = schema_variant, source_rows = nrow(d),
    source_columns = paste(hdr, collapse = ";"), patient_source_type = patient_source_type,
    extract_date_raw = paste(extract_values, collapse = " | "),
    extract_date = if (length(parsed_extract)) as.character(parsed_extract[1]) else NA_character_,
    publication_values = paste(sort(unique(d$PUBLICATION)), collapse = " | "),
    type_values = paste(sort(unique(d$TYPE)), collapse = " | "),
    sex_values = paste(sort(unique(d$SEX)), collapse = " | "),
    age_values = paste(sort(unique(d$AGE)), collapse = " | "),
    coverage_scope = "England all-practice publication; no row-level country field",
    duplicate_practice_codes = duplicate_practice_codes,
    missing_practice_codes = missing_practice_codes,
    missing_patient_counts = missing_patient_counts,
    nonnumeric_patient_counts = nonnumeric_patients,
    negative_patient_counts = negative_patient_counts,
    noninteger_patient_counts = noninteger_patient_counts,
    national_list_size = national_total, plausible_total = plausible_total,
    qc_status = status, failure_reason = failure_reason
  )
  list_size_observed[[i]] <- data.table(year_month = ym, list_size = national_total)
  if (status != "PASS") stop(basename(f), " failed denominator validation: ", failure_reason)
}

list_size_source_qc <- rbindlist(list_size_source_qc)
setorder(list_size_source_qc, year_month)
list_size_source_qc[, quarterly_change_pct := 100 * (national_list_size / shift(national_list_size) - 1),
                    by = analytical_role]
list_size_source_qc[, discontinuity_flag := !is.na(quarterly_change_pct) &
                      abs(quarterly_change_pct) > 2]
atomic_fwrite(list_size_source_qc, file.path(stage2_dir, "list_size_source_qc.csv"))
if (list_size_source_qc[discontinuity_flag == TRUE, .N] > 0L) {
  stop("A denominator quarterly change exceeded the pre-declared 2% review threshold.")
}

list.size <- rbindlist(list_size_observed)[year_month %in% list_size_ym]
setorder(list.size, year_month)
list.size[, `:=`(
  date = as.Date(paste0(year_month, "01"), "%Y%m%d"),
  list_size_source_month = year_month
)]
full <- data.table(date = seq(study_start, study_end, by = "month"))
full[, year_month := as.integer(format(date, "%Y%m"))]
list_size <- merge(
  full,
  list.size[, .(date, list_size, list_size_source_month)],
  by = "date", all.x = TRUE, sort = TRUE
)
list_size[, `:=`(
  list_size = nafill(list_size, type = "locf"),
  list_size_source_month = nafill(list_size_source_month, type = "locf")
)]
list_size[, list_size_carried_forward := year_month != list_size_source_month]
list_size <- list_size[, .(
  year_month, list_size, list_size_source_month, list_size_carried_forward
)]
if (nrow(list_size) != 48L || anyNA(list_size) || anyDuplicated(list_size$year_month)) {
  stop("The carried-forward denominator does not form one complete 48-month series.")
}
atomic_fwrite(list_size, file.path(out_dir, "listsize.csv"))


### 3. Process data -----------

# `data_dir` and `expected_ym` were declared once in the setup section.

fails <- character(0)   # validation failures collected, reported together

# Ensure year month is in the same format across dataframes
normalise_ym <- function(x, fname) {
  if (is.numeric(x)) return(as.integer(x))
  x <- trimws(as.character(x))
  if (all(grepl("^\\d{6}$", x))) return(as.integer(x))
  if (all(grepl("^\\d{4}-\\d{2}$", x)))
    return(as.integer(paste0(substr(x, 1, 4), substr(x, 6, 7))))
  stop(sprintf("%s: unrecognised year_month format (e.g. %s).",
               fname, paste(head(unique(x), 3), collapse = ", ")))
}

## Collate the monthly files for drug and drug class

read_monthly_set <- function(prefix) {
  paths <- file.path(data_dir, sprintf("%s_%d.csv", prefix, expected_ym))
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(sprintf("%d expected file(s) not found in %s:\n%s",
                 length(missing), data_dir,
                 paste(basename(missing), collapse = "\n")))
  }
  req <- c("year_month", "bnf_chapter_code", "bnf_class_code", "items")
  if (prefix == "data_bydrug") req <- c(req, "bnf_drug_code")
  cc <- list(character = c("bnf_chapter_code", "bnf_section_code", "bnf_class_code"))
  if (prefix == "data_bydrug") {
    cc$character <- c(cc$character, "bnf_drug_code")
  }
  
  read_one <- function(i) {
    f  <- paths[i]
    dt <- fread(f, colClasses = cc)
    miss_col <- setdiff(req, names(dt))
    if (length(miss_col))
      stop(sprintf("%s: required column(s) missing: %s",
                   basename(f), paste(miss_col, collapse = ", ")))
    dt[, year_month := normalise_ym(year_month, basename(f))]
    u <- unique(dt$year_month)
    if (length(u) != 1L || u != expected_ym[i])
      stop(sprintf("%s: contains year_month %s; expected %d",
                   basename(f), paste(head(u, 3), collapse = ", "), expected_ym[i]))
    ## harmonise key column types so batch differences cannot propagate
    dt[, `:=`(bnf_chapter_code = sprintf("%02d", as.integer(bnf_chapter_code)),
              bnf_section_code = as.character(bnf_section_code),
              bnf_class_code   = sprintf("%06d", as.integer(bnf_class_code)),
              items            = as.numeric(items))]
    dt
  }
  ## match by name: batches may differ in column order as well as formats
  rbindlist(lapply(seq_along(paths), read_one), use.names = TRUE)
}

class_monthly <- read_monthly_set("data_byclass")
drug_monthly  <- read_monthly_set("data_bydrug")
setorder(class_monthly, year_month, bnf_class_code)
setorder(drug_monthly,  year_month, bnf_drug_code)

## Validate prescribing frames
if (anyDuplicated(class_monthly, by = c("year_month", "bnf_class_code"))) {
  fails <- c(fails, "class_monthly: duplicate year_month x class rows.")
}
if (anyDuplicated(drug_monthly, by = c("year_month", "bnf_drug_code"))) {
  fails <- c(fails, "drug_monthly: duplicate year_month x drug rows.")
}

if (class_monthly[is.na(items) | items < 1 | is.na(bnf_class_code), .N] > 0) {
  fails <- c(fails, "class_monthly: missing/non-positive items or missing class codes.")
}
if (drug_monthly[is.na(items) | items < 1 | is.na(bnf_class_code), .N] > 0) {
  fails <- c(fails, "drug_monthly: missing/non-positive items or missing class codes.")
}

valid_chapters <- sprintf("%02d", 1:14)
if (!all(class_monthly$bnf_chapter_code %in% valid_chapters) ||
    !all(drug_monthly$bnf_chapter_code %in% valid_chapters)) {
  fails <- c(fails, "chapters outside 1-14 present - import filter has slipped.")
}

## drug code prefix must reproduce the class code 
n_prefix_bad <- drug_monthly[substr(bnf_drug_code, 1, 6) != bnf_class_code, .N]
if (n_prefix_bad > 0)
  fails <- c(fails, sprintf("drug_monthly: %d rows where substr(drug code, 1, 6) != class code.",
                            n_prefix_bad))

## Validate the denominator

listsize <- fread(file.path(data_dir, "listsize.csv"))
listsize[, year_month := normalise_ym(year_month, "listsize.csv")]

if (!identical(listsize[, sort(unique(year_month))], expected_ym)) {
  fails <- c(fails, "listsize: months do not match the 2022-2025 window exactly.")
}
if (nrow(listsize) != length(expected_ym)) {
  fails <- c(fails, "listsize: duplicate months present.")
}
if (listsize[is.na(list_size) | list_size < 55e6 | list_size > 70e6, .N] > 0) {
  fails <- c(fails, "listsize: missing or implausible values (expected ~55-70 million).")
}

## Cross-check: class-level vs drug-level totals per month - totals should agree 

qc_month <- merge(
  class_monthly[, .(class_items = sum(items)), by = year_month],
  drug_monthly[,  .(drug_items  = sum(items)), by = year_month],
  by = "year_month"
)
qc_month <- merge(qc_month, listsize, by = "year_month")
qc_month[, items_diff := class_items - drug_items]

if (qc_month[items_diff != 0, .N] > 0) {
  warning(sprintf("class vs drug item totals differ in %d month(s) - see outputs/qc/input_qc_by_month.csv",
                  qc_month[items_diff != 0, .N]))
}

if (length(fails)) stop(paste(c("Input validation failed:", fails), collapse = "\n  - "))


### 4. Create shared covariate data frame ---------

# Frozen England-and-Wales bank-holiday calendar and monthly working-day count.
calendar_json <- jsonlite::fromJSON(calendar_path)
calendar_events <- as.data.table(calendar_json$`england-and-wales`$events)
calendar_events[, bank_holiday_date := as.Date(date)]
calendar_events <- calendar_events[
  bank_holiday_date >= study_start &
    bank_holiday_date <= as.Date("2025-12-31")
]
setnames(calendar_events, "title", "bank_holiday_title")
calendar_events <- calendar_events[, .(
  bank_holiday_date, bank_holiday_title, notes, bunting
)]
if (nrow(calendar_events) != 35L || anyNA(calendar_events$bank_holiday_date) ||
    anyDuplicated(calendar_events[, .(bank_holiday_date, bank_holiday_title)])) {
  stop("The frozen 2022-2025 England-and-Wales bank-holiday snapshot failed validation.")
}
atomic_fwrite(calendar_events, file.path(stage2_dir, "bank_holidays_2022_2025.csv"))

calendar_days <- data.table(date = seq(study_start, as.Date("2025-12-31"), by = "day"))
calendar_days[, `:=`(
  year_month = as.integer(format(date, "%Y%m")),
  weekday_number = as.POSIXlt(date)$wday,
  is_bank_holiday = date %in% calendar_events$bank_holiday_date
)]
calendar_days[, is_working_day := weekday_number %in% 1:5 & !is_bank_holiday]
working_days_monthly <- calendar_days[, .(
  working_days = sum(is_working_day),
  weekday_days = sum(weekday_number %in% 1:5),
  bank_holidays_on_weekdays = sum(is_bank_holiday & weekday_number %in% 1:5)
), by = year_month]
setorder(working_days_monthly, year_month)
if (nrow(working_days_monthly) != 48L || any(working_days_monthly$working_days <= 0L) ||
    anyDuplicated(working_days_monthly$year_month)) {
  stop("The working-day calendar does not form one valid 48-month series.")
}
working_day_spot_checks <- data.table(
  year_month = c(202206L, 202209L, 202305L),
  expected_working_days = c(20L, 21L, 20L),
  rationale = c(
    "Spring and Platinum Jubilee holidays",
    "State Funeral of Queen Elizabeth II",
    "Early May, Coronation and Spring holidays"
  )
)
working_day_spot_checks <- merge(
  working_day_spot_checks,
  working_days_monthly[, .(year_month, observed_working_days = working_days)],
  by = "year_month", all.x = TRUE
)
working_day_spot_checks[, passed := observed_working_days == expected_working_days]
atomic_fwrite(working_day_spot_checks, file.path(stage2_dir, "working_day_spot_checks.csv"))
if (!all(working_day_spot_checks$passed)) stop("A special-holiday working-day spot check failed.")
atomic_fwrite(working_days_monthly, file.path(stage2_dir, "working_days_2022_2025.csv"))

covar <- data.table(year_month = expected_ym)
covar[, month_date := as.Date(paste0(year_month, "01"), "%Y%m%d")]
covar[, t := seq_len(.N)]
covar[, days_in_month := as.integer(lubridate::days_in_month(month_date))]

# Fourier terms: annual (12-month) and first harmonic (6-month) periods
covar[, `:=`(
  sin12 = sin(2 * pi * t / harmonic_period_months[1]),
  cos12 = cos(2 * pi * t / harmonic_period_months[1]),
  sin6  = sin(2 * pi * t / harmonic_period_months[2]),
  cos6  = cos(2 * pi * t / harmonic_period_months[2])
)]

# Spline basis stored as fixed columns so every series shares identical knots
S <- splines::ns(covar$t, df = trend_spline_df)
covar[, `:=`(trend1 = S[, 1], trend2 = S[, 2], trend3 = S[, 3])]

# Denominator and combined offset (sum of logs; no large product formed)
covar <- merge(covar, listsize, by = "year_month", sort = TRUE)
covar[, offset_log_patient_days := log(list_size) + log(days_in_month)]
covar <- merge(covar, working_days_monthly[, .(year_month, working_days)],
               by = "year_month", sort = TRUE)
covar[, offset_log_patient_working_days := log(list_size) + log(working_days)]

## Validate the frame 
stopifnot(
  nrow(covar) == 48,
  !anyNA(covar),
  sum(covar$days_in_month) == 1461,                      # incl. 29 Feb 2024
  covar[year_month == 202402, days_in_month] == 29,
  all(diff(covar$t) == 1)
)

## Stage 2 flow accounting and updated source manifest -------------------------------

expected_all_epd_ym <- sort(c(expected_ym, window_ym))
if (!identical(epd_file_qc$year_month, expected_all_epd_ym) ||
    any(epd_file_qc$qc_status != "PASS")) {
  stop("The 60-file EPD QC table is incomplete, duplicated, or contains a failed file.")
}
if (epd_file_qc[schema_variant == "legacy_through_202502", .N] != 50L ||
    epd_file_qc[schema_variant == "snomed_from_202503", .N] != 10L) {
  stop("The observed EPD schema counts do not match the declared February/March 2025 boundary.")
}

epd_flow <- data.table(
  stage = c(
    "Frozen EPD source records", "Records retained in BNF chapters 01-14",
    "Records excluded outside BNF chapters 01-14", "Frozen EPD source items",
    "Items retained in BNF chapters 01-14", "Items excluded outside BNF chapters 01-14",
    "National drug-month rows", "National class-month rows"
  ),
  records_or_rows = c(
    sum(epd_file_qc$raw_row_count), sum(epd_file_qc$retained_chapter_01_14_rows),
    sum(epd_file_qc$removed_chapter_rows), NA_real_, NA_real_, NA_real_,
    sum(epd_file_qc$national_drug_rows), sum(epd_file_qc$national_class_rows)
  ),
  items = c(
    NA_real_, NA_real_, NA_real_, sum(epd_file_qc$raw_total_items),
    sum(epd_file_qc$retained_chapter_01_14_items), sum(epd_file_qc$removed_chapter_items),
    sum(epd_file_qc$retained_chapter_01_14_items),
    sum(epd_file_qc$retained_chapter_01_14_items)
  ),
  scope = c(rep("2021 descriptive plus 2022-2025 main", 6L),
            rep("monthly aggregates across all 60 months", 2L))
)
atomic_fwrite(epd_flow, file.path(stage2_dir, "epd_flow_summary.csv"))

stage0_manifest_path <- file.path(analysis_dir, "reproducibility", "input_manifest.csv")
if (!file.exists(stage0_manifest_path)) stop("Stage 0 input_manifest.csv is missing.")
stage2_manifest <- fread(stage0_manifest_path)
stage2_manifest[, `:=`(stage2_qc_status = NA_character_, stage2_qc_file = NA_character_)]

epd_update_rows <- match(epd_file_qc$source_id, stage2_manifest$source_id)
if (anyNA(epd_update_rows)) stop("An EPD source is absent from the Stage 0 input manifest.")
stage2_manifest[epd_update_rows, `:=`(
  schema_variant = epd_file_qc$schema_variant,
  row_count = epd_file_qc$raw_row_count,
  row_count_status = "verified during canonical Stage 2 import",
  column_count = epd_file_qc$raw_column_count,
  columns = epd_file_qc$source_columns,
  declared_grain = paste0(
    "published EPD record; candidate practice/presentation/month key is non-unique; ",
    "all source records are summed to unique national drug-month and class-month rows"
  ),
  candidate_key = epd_file_qc$candidate_key,
  duplicate_key_count = epd_file_qc$candidate_duplicate_excess_rows,
  missing_key_count = epd_file_qc$month_mismatch_rows +
    epd_file_qc$missing_practice_rows + epd_file_qc$missing_presentation_rows,
  measure_total = epd_file_qc$raw_total_items,
  stage2_qc_status = epd_file_qc$qc_status,
  stage2_qc_file = file.path(
    "qc", "stage2", "epd_monthly", sprintf("epd_file_qc_%d.csv", epd_file_qc$year_month)
  ),
  notes = paste0(
    "Stage 2 retained ", epd_file_qc$retained_chapter_01_14_rows,
    " rows and ", epd_file_qc$retained_chapter_01_14_items,
    " items in BNF chapters 01-14; excluded ", epd_file_qc$removed_chapter_rows,
    " rows and ", epd_file_qc$removed_chapter_items,
    " items outside scope. In-scope chemical-code mismatches: ",
    epd_file_qc$chemical_code_mismatch_in_scope_rows, "."
  )
)]

list_update_rows <- match(list_size_source_qc$source_id, stage2_manifest$source_id)
if (anyNA(list_update_rows)) stop("A list-size source is absent from the Stage 0 input manifest.")
stage2_manifest[list_update_rows, `:=`(
  schema_variant = list_size_source_qc$schema_variant,
  row_count = list_size_source_qc$source_rows,
  row_count_status = "verified during canonical Stage 2 import",
  declared_grain = "one England all-practice denominator row per GP practice and quarterly source month",
  candidate_key = "EXTRACT_DATE + CODE",
  duplicate_key_count = list_size_source_qc$duplicate_practice_codes,
  missing_key_count = list_size_source_qc$missing_practice_codes,
  measure_total = list_size_source_qc$national_list_size,
  stage2_qc_status = list_size_source_qc$qc_status,
  stage2_qc_file = file.path("qc", "stage2", "list_size_source_qc.csv"),
  notes = paste0(
    list_size_source_qc$coverage_scope, "; national total ",
    list_size_source_qc$national_list_size, "."
  )
)]

bnf_manifest_row <- grep("^BNF_REFERENCE", stage2_manifest$source_id)
calendar_manifest_row <- grep("^BANK_HOLIDAYS", stage2_manifest$source_id)
if (length(bnf_manifest_row) != 1L || length(calendar_manifest_row) != 1L) {
  stop("BNF or bank-holiday source is not uniquely identified in the Stage 0 manifest.")
}
stage2_manifest[bnf_manifest_row, `:=`(
  stage2_qc_status = "PASS",
  stage2_qc_file = file.path("qc", "stage2", "epd_file_qc.csv"),
  notes = paste0(notes, " Stage 2 verified one descriptor row per BNF class code.")
)]
stage2_manifest[calendar_manifest_row, `:=`(
  stage2_qc_status = "PASS",
  stage2_qc_file = file.path("qc", "stage2", "working_days_2022_2025.csv"),
  notes = paste0(notes, " Stage 2 reproduced 48 positive monthly working-day counts and passed special-holiday checks.")
)]
setorder(stage2_manifest, source_id)
atomic_fwrite(stage2_manifest, file.path(stage2_dir, "input_manifest_stage2.csv"))

## Save pre-recode objects and quality-control outcomes -------------------------------

qc_month <- merge(qc_month, working_days_monthly[, .(year_month, working_days)],
                  by = "year_month", all.x = TRUE, sort = TRUE)
if (nrow(qc_month) != 48L || anyNA(qc_month$working_days)) {
  stop("Monthly prescribing, denominator, and working-day data did not join one-to-one.")
}
atomic_save_rds(class_monthly, file.path(data_dir, "class_monthly_pre_recode.rds"))
atomic_save_rds(drug_monthly,  file.path(data_dir, "drug_monthly_pre_recode.rds"))
atomic_save_rds(covar,         file.path(data_dir, "covariate_frame.rds"))
atomic_fwrite(qc_month,        file.path(data_dir, "input_qc_by_month.csv"))

cat(sprintf(paste0(
  "Complete.\n",
  "  Months:            %d (%d - %d)\n",
  "  Classes observed:  %d\n",
  "  Drugs observed:    %d\n",
  "  Total items:       %s\n",
  "  Listsize:          %s - %s (%d distinct values)\n",
  "  Saved: pre-recode panels, covariate_frame.rds, input_qc_by_month.csv\n"),
  nrow(covar), min(covar$year_month), max(covar$year_month),
  class_monthly[, uniqueN(bnf_class_code)],
  drug_monthly[, uniqueN(bnf_drug_code)],
  format(class_monthly[, sum(items)], big.mark = ","),
  format(min(covar$list_size), big.mark = ","),
  format(max(covar$list_size), big.mark = ","),
  uniqueN(covar$list_size)))

### 5. Reconcile recodes - in case BNF codes change between 2022-2025 for the same drug

## normalise substance name so drift in spellings collapse 
normalise_name <- function(x) {
  x |>
    str_to_lower() |>
    str_remove("\\s*\\([^()]*\\)\\s*$") |>   # drop a single trailing "(qualifier)"
    str_replace_all("[[:punct:]]", " ") |>
    str_replace_all("\\s+", " ") |>
    str_trim()
}

# Pick the code/name columns for a level
.level_cols <- function(level) {
  switch(level,
         drug  = list(code = "bnf_drug_code",  name = "bnf_drug_name"),
         class = list(code = "bnf_class_code", name = "bnf_class_name"),
         stop("level must be 'drug' or 'class'"))
}

## Detect clusters that share a normalised name
detect_recodes <- function(x, level = c("drug", "class"), n_months_req = 48L) {
  level <- match.arg(level)
  cols  <- .level_cols(level)
  x <- x |> rename(.code = all_of(cols$code), .name = all_of(cols$name))
  
  code_activity <- x |>
    group_by(.code) |>
    summarise(name           = .name[which.max(year_month)],
              months_present = n_distinct(year_month),
              total_items    = sum(items),
              active         = list(sort(unique(year_month))),
              .groups = "drop") |>
    mutate(nname = normalise_name(name))
  
  code_activity |>
    group_by(nname) |>
    filter(n() > 1) |>
    summarise(
      n_codes            = n(),
      codes              = paste(.code, collapse = " + "),
      example_name       = dplyr::first(name),
      combined_items     = sum(total_items),
      best_single_months = max(months_present),
      union_months       = n_distinct(unlist(active)),
      max_simultaneous   = max(as.integer(table(unlist(active)))),
      .groups = "drop") |>
    mutate(level         = level,
           clean_recode  = max_simultaneous == 1L,
           recovers      = union_months > best_single_months,
           heals_to_full = union_months == n_months_req) |>
    filter(recovers) |>
    arrange(desc(combined_items))
}

apply_recode_crosswalk <- function(x, crosswalk, level = c("drug", "class")) {
  level <- match.arg(level)
  stopifnot(all(c("from_code", "to_code") %in% names(crosswalk)))
  
  if (level == "class") {
    class_desc <- x |>
      arrange(bnf_class_code, desc(year_month)) |>
      group_by(bnf_class_code) |>
      slice(1L) |>
      ungroup() |>
      select(bnf_class_code, bnf_chapter_code, bnf_chapter_name,
             bnf_section_code, bnf_section_name, bnf_class_name)
    out <- x |>
      left_join(crosswalk, by = c("bnf_class_code" = "from_code")) |>
      mutate(bnf_class_code = coalesce(to_code, bnf_class_code)) |>
      select(year_month, bnf_class_code, items) |>
      group_by(year_month, bnf_class_code) |>
      summarise(items = sum(items), .groups = "drop") |>
      left_join(class_desc, by = "bnf_class_code")
    if (anyNA(out$bnf_class_name))
      warning("apply_recode_crosswalk: a canonical class code has no descriptor - check crosswalk targets.")
    return(out |>
             select(year_month, bnf_chapter_code, bnf_chapter_name, bnf_section_code,
                    bnf_section_name, bnf_class_code, bnf_class_name, items))
  }
  
  ## drug level: remap the 9-char code, re-derive class from its prefix
  class_desc <- x |>
    arrange(bnf_class_code, desc(year_month)) |>
    group_by(bnf_class_code) |>
    slice(1L) |>
    ungroup() |>
    select(bnf_class_code, bnf_chapter_code, bnf_chapter_name,
           bnf_section_code, bnf_section_name, bnf_class_name)
  out <- x |>
    left_join(crosswalk, by = c("bnf_drug_code" = "from_code")) |>
    mutate(bnf_drug_code  = coalesce(to_code, bnf_drug_code),
           bnf_class_code = substr(bnf_drug_code, 1, 6)) |>
    group_by(bnf_drug_code) |>
    mutate(bnf_drug_name = bnf_drug_name[which.max(year_month)]) |>
    ungroup() |>
    group_by(year_month, bnf_class_code, bnf_drug_code, bnf_drug_name) |>
    summarise(items = sum(items), .groups = "drop") |>
    left_join(class_desc, by = "bnf_class_code")
  if (anyNA(out$bnf_class_name))
    warning("apply_recode_crosswalk: a canonical class code has no descriptor - check crosswalk targets.")
  out |>
    select(year_month, bnf_chapter_code, bnf_chapter_name, bnf_section_code,
           bnf_section_name, bnf_class_code, bnf_class_name,
           bnf_drug_code, bnf_drug_name, items)
}

# Make list of candidates to re-include
if (exists("drug_monthly")) {
  cand_drug <- detect_recodes(drug_monthly, "drug")
  message("cand_drug: ", nrow(cand_drug), " candidate cluster(s).")
  if (exists("out_dir")) readr::write_csv(cand_drug, file.path(out_dir, "recode_candidates_drug.csv"))
}
if (exists("class_monthly")) {
  cand_class <- detect_recodes(class_monthly, "class")
  message("cand_class: ", nrow(cand_class), " candidate cluster(s).")
  if (exists("out_dir")) readr::write_csv(cand_class, file.path(out_dir, "recode_candidates_class.csv"))
}

# Review candidates - no candidate classes, only beclomet is candidate for drugs to reinclude
cand_drug
cand_class

cand_drug |> filter(heals_to_full, combined_items >= 4000)   # keeps the only candidate we need

# Select which of the codes we want to keep for this
cand_drug |> slice(1) |> pull(codes)     # the two full codes for the line below
drug_monthly |> filter(bnf_drug_code %in% c("0301011AB","0302000AA")) |>
  group_by(bnf_drug_code) |>
  summarise(first = min(year_month), last = max(year_month), items = sum(items), .groups = "drop")

# Unify codes so same code is used for the above candidate drug
xwalk_drug <- tibble::tribble(
  ~from_code,   ~to_code,
  "0301011AB",  "0302000AA"   # Trimbow: reconciled to the current (2025) corticosteroid code
)

accepted_candidate <- cand_drug |>
  filter(str_detect(codes, "0301011AB") & str_detect(codes, "0302000AA"))
if (nrow(accepted_candidate) != 1L ||
    !accepted_candidate$clean_recode || !accepted_candidate$heals_to_full) {
  stop("The pre-declared Trimbow recode is not reproduced as one clean, full-window candidate.")
}

drug_monthly_before_recode <- as.data.table(copy(drug_monthly))
recode_source_evidence <- drug_monthly_before_recode[
  bnf_drug_code %chin% c("0301011AB", "0302000AA"),
  .(
    first_month = min(year_month), last_month = max(year_month),
    months_present = uniqueN(year_month), items = sum(items),
    drug_names = paste(sort(unique(bnf_drug_name)), collapse = " | ")
  ), by = bnf_drug_code
]
if (nrow(recode_source_evidence) != 2L) {
  stop("Both source and target codes for the accepted recode must be observed.")
}
recode_month_overlap <- intersect(
  drug_monthly_before_recode[bnf_drug_code == "0301011AB", year_month],
  drug_monthly_before_recode[bnf_drug_code == "0302000AA", year_month]
)
if (length(recode_month_overlap)) {
  stop("The accepted recode codes overlap in month and cannot be combined without review.")
}

pre_recode_drug_items <- drug_monthly_before_recode[, sum(items)]
pre_recode_class_items <- as.data.table(class_monthly)[, sum(items)]
drug_monthly <- apply_recode_crosswalk(drug_monthly, xwalk_drug, "drug")

# Rederive class_monthly
class_monthly <- drug_monthly |>
  group_by(year_month, bnf_chapter_code, bnf_chapter_name, bnf_section_code,
           bnf_section_name, bnf_class_code, bnf_class_name) |>
  summarise(items = sum(items), .groups = "drop")

drug_monthly_after_recode <- as.data.table(drug_monthly)
class_monthly_after_recode <- as.data.table(class_monthly)
post_recode_drug_items <- drug_monthly_after_recode[, sum(items)]
post_recode_class_items <- class_monthly_after_recode[, sum(items)]
post_recode_drug_duplicates <- sum(duplicated(
  drug_monthly_after_recode, by = c("year_month", "bnf_drug_code")
))
post_recode_class_duplicates <- sum(duplicated(
  class_monthly_after_recode, by = c("year_month", "bnf_class_code")
))
post_recode_prefix_mismatches <- drug_monthly_after_recode[
  substr(bnf_drug_code, 1L, 6L) != bnf_class_code, .N
]
post_recode_from_code_rows <- drug_monthly_after_recode[bnf_drug_code == "0301011AB", .N]
post_recode_target_months <- drug_monthly_after_recode[
  bnf_drug_code == "0302000AA", uniqueN(year_month)
]
post_recode_monthly_reconciliation <- merge(
  drug_monthly_after_recode[, .(drug_items = sum(items)), by = year_month],
  class_monthly_after_recode[, .(class_items = sum(items)), by = year_month],
  by = "year_month", all = TRUE
)
post_recode_monthly_reconciliation[, item_difference := class_items - drug_items]

recode_pass <- pre_recode_drug_items == pre_recode_class_items &&
  post_recode_drug_items == pre_recode_drug_items &&
  post_recode_class_items == pre_recode_class_items &&
  post_recode_drug_duplicates == 0L && post_recode_class_duplicates == 0L &&
  post_recode_prefix_mismatches == 0L && post_recode_from_code_rows == 0L &&
  post_recode_target_months == 48L &&
  all(post_recode_monthly_reconciliation$item_difference == 0)

recode_qc <- data.table(
  from_code = "0301011AB", to_code = "0302000AA",
  evidence_status = "clean non-overlapping normalised-name candidate that restores all 48 months",
  overlap_months_before = length(recode_month_overlap),
  target_months_after = post_recode_target_months,
  from_code_rows_after = post_recode_from_code_rows,
  pre_recode_drug_items = pre_recode_drug_items,
  post_recode_drug_items = post_recode_drug_items,
  pre_recode_class_items = pre_recode_class_items,
  post_recode_class_items = post_recode_class_items,
  post_recode_drug_duplicate_rows = post_recode_drug_duplicates,
  post_recode_class_duplicate_rows = post_recode_class_duplicates,
  post_recode_prefix_mismatch_rows = post_recode_prefix_mismatches,
  maximum_monthly_class_drug_item_difference = max(abs(post_recode_monthly_reconciliation$item_difference)),
  qc_status = if (recode_pass) "PASS" else "FAIL"
)
recode_crosswalk_record <- merge(
  as.data.table(xwalk_drug), recode_source_evidence,
  by.x = "from_code", by.y = "bnf_drug_code", all.x = TRUE
)
recode_crosswalk_record[, `:=`(
  decision = "accepted",
  rationale = paste0(
    "Same normalised medicine name; no simultaneous months; union restores the 48-month series; ",
    "canonical target is the current 2025 code."
  )
)]
atomic_fwrite(recode_source_evidence, file.path(stage2_dir, "accepted_recode_source_evidence.csv"))
atomic_fwrite(recode_crosswalk_record, file.path(stage2_dir, "accepted_recode_crosswalk.csv"))
atomic_fwrite(post_recode_monthly_reconciliation,
              file.path(stage2_dir, "post_recode_class_drug_reconciliation.csv"))
atomic_fwrite(recode_qc, file.path(stage2_dir, "recode_application_qc.csv"))
if (!recode_pass) stop("The accepted drug-code recode failed a preservation or uniqueness check.")

atomic_save_rds(class_monthly, file.path(data_dir, "class_monthly.rds"))
atomic_save_rds(drug_monthly, file.path(data_dir, "drug_monthly.rds"))

## Final Stage 2 gate ---------------------------------------------------------------

stage2_check <- function(check_id, passed, observed, requirement) {
  data.table(
    check_id = check_id,
    qc_status = if (isTRUE(passed)) "PASS" else "FAIL",
    observed = as.character(observed),
    requirement = requirement
  )
}
invalid_item_rows <- epd_file_qc[, sum(
  missing_items_rows + nonnumeric_items_rows + negative_items_rows + noninteger_items_rows
)]
missing_epd_key_rows <- epd_file_qc[, sum(
  month_mismatch_rows + missing_practice_rows + missing_presentation_rows
)]
stage2_qc_summary <- rbindlist(list(
  stage2_check("EPD_FILE_COVERAGE", nrow(epd_file_qc) == 60L,
               paste(nrow(epd_file_qc), "files"), "60 frozen monthly archives"),
  stage2_check("EPD_MONTH_COVERAGE", identical(epd_file_qc$year_month, expected_all_epd_ym),
               paste(min(epd_file_qc$year_month), max(epd_file_qc$year_month), sep = "-"),
               "exactly 202101-202512 with 2021 descriptive only"),
  stage2_check("EPD_SCHEMA_BOUNDARY",
               epd_file_qc[schema_variant == "legacy_through_202502", .N] == 50L &&
                 epd_file_qc[schema_variant == "snomed_from_202503", .N] == 10L,
               paste0("legacy=", epd_file_qc[schema_variant == "legacy_through_202502", .N],
                      "; current=", epd_file_qc[schema_variant == "snomed_from_202503", .N]),
               "legacy through 202502; current from 202503"),
  stage2_check("EPD_REQUIRED_KEYS", missing_epd_key_rows == 0,
               missing_epd_key_rows, "no missing/mismatched month, practice, or presentation key"),
  stage2_check("EPD_ITEMS_VALID", invalid_item_rows == 0,
               invalid_item_rows, "no missing, nonnumeric, negative, or noninteger ITEMS"),
  stage2_check("EPD_IN_SCOPE_CHEMICAL_CODE",
               sum(epd_file_qc$chemical_code_mismatch_in_scope_rows) == 0,
               sum(epd_file_qc$chemical_code_mismatch_in_scope_rows),
               "presentation-derived drug code agrees with supplied chemical code in chapters 01-14"),
  stage2_check(
    "EPD_BNF_LOOKUP_ACCOUNTING",
    unmatched_detail_items == unmatched_summary_items,
    paste0(
      format(unmatched_summary_items, big.mark = ","), " retained unmatched items; ",
      uniqueN(epd_unmatched_lookup_qc$bnf_class_code), " class codes; ",
      uniqueN(epd_unmatched_lookup_qc$year_month), " months"
    ),
    "all gaps against the May 2025 lookup are retained, quantified by code/month, and reconcile to source QC"
  ),
  stage2_check("EPD_NATIONAL_GRAIN",
               sum(epd_file_qc$drug_month_duplicate_rows + epd_file_qc$class_month_duplicate_rows) == 0,
               sum(epd_file_qc$drug_month_duplicate_rows + epd_file_qc$class_month_duplicate_rows),
               "unique drug-month and class-month aggregates"),
  stage2_check("EPD_CLASS_DRUG_RECONCILIATION",
               all(epd_file_qc$class_drug_item_difference == 0),
               max(abs(epd_file_qc$class_drug_item_difference)),
               "exact item equality within every source month"),
  stage2_check("DENOMINATOR_SOURCE_FILES",
               nrow(list_size_source_qc) == 20L && all(list_size_source_qc$qc_status == "PASS"),
               paste(nrow(list_size_source_qc), "files"), "20 quarterly England all-practice files"),
  stage2_check("DENOMINATOR_DISCONTINUITY",
               !any(list_size_source_qc$discontinuity_flag),
               max(abs(list_size_source_qc$quarterly_change_pct), na.rm = TRUE),
               "no within-role quarterly change exceeds 2%"),
  stage2_check("DENOMINATOR_MONTHLY_JOIN",
               nrow(list_size) == 48L && !anyNA(list_size) && !anyDuplicated(list_size$year_month),
               paste(nrow(list_size), "monthly rows"), "one complete 48-month LOCF series with source indicator"),
  stage2_check("WORKING_DAY_CALENDAR",
               nrow(calendar_events) == 35L && nrow(working_days_monthly) == 48L &&
                 all(working_days_monthly$working_days > 0L),
               paste0(nrow(calendar_events), " events; ", nrow(working_days_monthly), " months"),
               "frozen England/Wales holidays and 48 positive month counts"),
  stage2_check("SPECIAL_HOLIDAY_SPOT_CHECKS", all(working_day_spot_checks$passed),
               paste(sum(working_day_spot_checks$passed), "of", nrow(working_day_spot_checks)),
               "all pre-declared special-holiday months agree"),
  stage2_check("ACCEPTED_RECODE", recode_pass, recode_qc$qc_status,
               "clean candidate; no overlap/double count; totals and hierarchy preserved"),
  stage2_check("FINAL_MAIN_PANEL",
               identical(sort(unique(drug_monthly_after_recode$year_month)), expected_ym) &&
                 identical(sort(unique(class_monthly_after_recode$year_month)), expected_ym),
               paste(uniqueN(drug_monthly_after_recode$year_month), "drug months;",
                     uniqueN(class_monthly_after_recode$year_month), "class months"),
               "exactly 48 main-analysis months at both levels")
), use.names = TRUE)
atomic_fwrite(stage2_qc_summary, file.path(stage2_dir, "stage2_qc_summary.csv"))
if (any(stage2_qc_summary$qc_status != "PASS")) {
  stop("Stage 2 final QC gate failed; see qc/stage2/stage2_qc_summary.csv.")
}

if (run_stage == "stage2") {
  stage2_completion <- data.table(
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stage = "Stage 2 import and source validation",
    status = "PASS",
    epd_files = nrow(epd_file_qc), denominator_files = nrow(list_size_source_qc),
    raw_epd_rows = sum(epd_file_qc$raw_row_count),
    retained_items = sum(epd_file_qc$retained_chapter_01_14_items),
    main_months = uniqueN(drug_monthly_after_recode$year_month),
    checks_passed = sum(stage2_qc_summary$qc_status == "PASS"),
    checks_total = nrow(stage2_qc_summary)
  )
  atomic_fwrite(stage2_completion, file.path(stage2_dir, "stage2_completion.csv"))

  # A Stage 2 manifest covers Stage 2 artefacts only. Downstream outputs may
  # already exist when Stage 2 is revalidated and must not enter this snapshot.
  manifest_exclusions <- c(
    file.path(stage2_dir, "stage2_output_manifest.csv"),
    file.path(repro_dir, "warnings.csv")
  )
  stage2_output_files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  stage2_relative_paths <- substring(stage2_output_files, nchar(out_dir) + 2L)
  downstream_output <-
    grepl("^(eligibility_|appendix1_|appendix_|screen_|model_failures|characterisation_|results_)",
          stage2_relative_paths) |
    grepl("^(figures|results)/", stage2_relative_paths) |
    grepl("^reproducibility/", stage2_relative_paths) |
    grepl("^qc/(stage1|stage3|full)/", stage2_relative_paths)
  keep_stage2 <- file.info(stage2_output_files)$isdir %in% FALSE &
    !stage2_output_files %in% manifest_exclusions &
    !grepl("\\.tmp$", stage2_output_files) &
    !downstream_output
  stage2_output_files <- stage2_output_files[keep_stage2]
  stage2_output_manifest <- data.table(
    relative_path = substring(stage2_output_files, nchar(out_dir) + 2L),
    size_bytes = unname(file.info(stage2_output_files)$size),
    sha256 = vapply(stage2_output_files, sha256_file, character(1))
  )
  setorder(stage2_output_manifest, relative_path)
  atomic_fwrite(stage2_output_manifest, file.path(stage2_dir, "stage2_output_manifest.csv"))

  message(
    "Stage 2 complete: all ", nrow(stage2_qc_summary),
    " final checks passed; stopping before eligibility and modelling."
  )
  return(invisible(TRUE))
}

### 6. Eligibility and coverage

# Eligible if prescribed in every month AND at least 1000 items dispensed nationally each calendar year
if (!exists("data_dir")) stop("data_dir not set - run the setup/import section first.")
if (!exists("expected_ym")) stop("expected_ym was not created by canonical setup.")
if (!exists("class_monthly")) class_monthly <- readRDS(file.path(data_dir, "class_monthly.rds"))
if (!exists("drug_monthly"))  drug_monthly  <- readRDS(file.path(data_dir, "drug_monthly.rds"))

n_months_req <- length(expected_ym)                    # 48
n_years_req  <- n_distinct(expected_ym %/% 100L)       # 4
min_items_yr <- eligibility_min_items_year

## Eligibility flags per drug class

elig_class <- class_monthly |>
  group_by(bnf_class_code, bnf_class_name) |>
  summarise(n_months = n_distinct(year_month), total_items = sum(items), .groups = "drop") |>
  left_join(
    class_monthly |>
      group_by(bnf_class_code, year = year_month %/% 100L) |>
      summarise(items_yr = sum(items), .groups = "drop") |>
      group_by(bnf_class_code) |>
      summarise(n_years = n(), min_year_items = min(items_yr), .groups = "drop"),
    by = "bnf_class_code") |>
  mutate(rule_every_month = n_months == n_months_req,
         rule_min_volume  = n_years == n_years_req & min_year_items >= min_items_yr,
         eligible         = rule_every_month & rule_min_volume)


# Eligibility flags per drug )note that drug names for the same code may have changed, so aggregate by code, and use more recent name)

drug_meta <- drug_monthly |>
  group_by(bnf_drug_code) |>
  summarise(bnf_class_code = bnf_class_code[which.max(year_month)],
            bnf_class_name = bnf_class_name[which.max(year_month)],
            bnf_drug_name  = bnf_drug_name[which.max(year_month)],
            .groups = "drop")

elig_drug <- drug_meta |>
  left_join(
    drug_monthly |>
      group_by(bnf_drug_code) |>
      summarise(n_months = n_distinct(year_month), total_items = sum(items), .groups = "drop"),
    by = "bnf_drug_code") |>
  left_join(
    drug_monthly |>
      group_by(bnf_drug_code, year = year_month %/% 100L) |>
      summarise(items_yr = sum(items), .groups = "drop") |>
      group_by(bnf_drug_code) |>
      summarise(n_years = n(), min_year_items = min(items_yr), .groups = "drop"),
    by = "bnf_drug_code") |>
  mutate(rule_every_month = n_months == n_months_req,
         rule_min_volume  = n_years == n_years_req & min_year_items >= min_items_yr,
         eligible         = rule_every_month & rule_min_volume)

## one row per series (guards against grouping anomalies)
stopifnot(
  nrow(elig_class) == n_distinct(class_monthly$bnf_class_code),
  nrow(elig_drug)  == n_distinct(drug_monthly$bnf_drug_code)
)

## List exclusions and summarise coverage

coverage_one <- function(tab, level_name) {
  eligible_flag <- tab$eligible %in% TRUE
  tibble(
    level = level_name,
    total_series = nrow(tab),
    eligible = sum(eligible_flag),
    excluded = sum(!eligible_flag),
    excluded_item_share = sum(tab$total_items[!eligible_flag]) / sum(tab$total_items)
  )
}
coverage <- bind_rows(
  coverage_one(elig_class, "class"),
  coverage_one(elig_drug, "drug")
)

add_reason <- function(tab) {
  tab |>
    filter(!eligible) |>
    mutate(reason = case_when(
      !rule_every_month & !rule_min_volume ~ "not every month; <1000 items/year",
      !rule_every_month                    ~ "not every month",
      TRUE                                 ~ "<1000 items/year")) |>
    arrange(desc(total_items))
}
excl_class <- add_reason(elig_class)
excl_drug  <- add_reason(elig_drug)

# Summarise reasons for exclusion
reason_breakdown <- function(excl, elig, level_name) {
  lvls <- c("<1000 items/year", "not every month", "not every month; <1000 items/year")
  n_total <- nrow(elig)
  
  tab <- excl |>
    count(reason, name = "n") |>
    right_join(tibble(reason = lvls), by = "reason") |>   # keep empty categories, fix order
    mutate(n = coalesce(n, 0L),
           pct_of_excluded = 100 * n / sum(n),
           pct_of_all      = 100 * n / n_total) |>
    transmute(reason, n,
              pct_of_excluded = round(pct_of_excluded, 1),
              pct_of_all      = round(pct_of_all, 1))
  
  cat(sprintf("\n%s: %d excluded of %d (%.1f%% of all %s)\n",
              level_name, sum(tab$n), n_total, 100 * sum(tab$n) / n_total, level_name))
  print(tab)
  invisible(tab)
}

excl_drug_summary  <- reason_breakdown(excl_drug,  elig_drug,  "drugs")
excl_class_summary <- reason_breakdown(excl_class, elig_class, "classes")

## Filter data for analysis sets 

class_monthly_elig <- class_monthly |>
  semi_join(filter(elig_class, eligible), by = "bnf_class_code")

drug_monthly_elig <- drug_monthly |>
  semi_join(filter(elig_drug, eligible), by = "bnf_drug_code")

## Save and report

saveRDS(elig_class,         file.path(data_dir, "eligibility_class.rds"))
saveRDS(elig_drug,          file.path(data_dir, "eligibility_drug.rds"))
saveRDS(class_monthly_elig, file.path(data_dir, "class_monthly_eligible.rds"))
saveRDS(drug_monthly_elig,  file.path(data_dir, "drug_monthly_eligible.rds"))
fwrite(coverage,           file.path(data_dir, "appendix1_coverage_summary.csv"))
fwrite(excl_class,         file.path(data_dir, "appendix1_exclusions_class.csv"))
fwrite(excl_drug,          file.path(data_dir, "appendix1_exclusions_drug.csv"))
fwrite(excl_class_summary, file.path(data_dir, "appendix1_exclusions_class_summary.csv"))
fwrite(excl_drug_summary,  file.path(data_dir, "appendix1_exclusions_drug_summary.csv"))

cc <- filter(coverage, level == "class")
cd <- filter(coverage, level == "drug")
cat(sprintf(paste0(
  "Complete.\n",
  "  Classes: %d eligible / %d total (excluded item share %.3f%%)\n",
  "  Drugs:   %d eligible / %d total (excluded item share %.3f%%)\n",
  "  Saved: eligibility tables, filtered analysis sets, Appendix 1 CSVs\n"),
  cc$eligible, cc$total_series, 100 * cc$excluded_item_share,
  cd$eligible, cd$total_series, 100 * cd$excluded_item_share))

## Check high-volume exclusions in case an exclusion is due to a BNF code change
if (nrow(excl_drug)) {
  cat("  Largest excluded drugs:\n")
  print(excl_drug |> slice_head(n = 20) |>
          select(bnf_drug_code, bnf_drug_name, total_items, reason))
}
if (nrow(excl_class)) {
  cat("  Largest excluded classes:\n")
  print(excl_class |> slice_head(n = 20) |>
          select(bnf_class_code, bnf_class_name, total_items, reason))
}

### 7. Model fitting Poisson GLM (fixed spline + Fourier + patient-days offset) -------
# Returns a one-row data frame - seasonality p value, distribution, route, diagnostics, harmonic coefficients, convergence flags

# Packages
.needed <- c("MASS", "AER", "sandwich", "lmtest", "car")
.missing <- .needed[!vapply(.needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("install required package(s): ", paste(.missing, collapse = ", "))

# formulas reference columns in the per-series modelling frame (offset included in-formula so Poisson and negative binomial treated identically)
.f_full <- items ~ trend1 + trend2 + trend3 + sin12 + cos12 + sin6 + cos6 + offset(off)
.f_red  <- items ~ trend1 + trend2 + trend3 + offset(off)
.harmonics <- c("sin12", "cos12", "sin6", "cos6")

# HAC route: joint Wald on the four harmonics under Newey-West covariance
.do_hac <- function(fp, bw_cap = 12) {
  bw     <- tryCatch(sandwich::bwNeweyWest(fp, prewhite = 1L), error = function(e) NA_real_)
  bw_use <- if (is.na(bw)) bw_cap else min(bw, bw_cap)
  capped <- !is.na(bw) && bw > bw_cap
  V  <- suppressWarnings(
    sandwich::NeweyWest(fp, lag = bw_use, prewhite = 1L, adjust = TRUE))
  wt <- car::linearHypothesis(fp, paste0(.harmonics, " = 0"), vcov. = V, test = "Chisq")
  list(p = wt[["Pr(>Chisq)"]][2], bw = bw_use, capped = capped)
}

# NB fit to treat non-convergence as failure, and will fall back to HAC
.nb_fit <- function(form, d) {
  bad <- FALSE
  f <- withCallingHandlers(
    tryCatch(MASS::glm.nb(form, data = d), error = function(e) NULL),
    warning = function(w) {
      if (grepl("alternation|iteration limit|theta|converge",
                conditionMessage(w), ignore.case = TRUE)) bad <<- TRUE
      invokeRestart("muffleWarning")
    })
  if (is.null(f) || bad || !isTRUE(f$converged)) return(NULL)
  f
}

# fit series function

fit_test_series <- function(series, covar,
                            alpha_disp = diagnostic_alpha,
                            alpha_lb = diagnostic_alpha,
                            lb_lag = harmonic_period_months[1],
                            f_full = .f_full,
                            f_red = .f_red,
                            offset_col = "offset_log_patient_days") {
  res <- data.frame(
    p_value = NA_real_, distribution = NA_character_, route = NA_character_,
    disp_ratio = NA_real_, disp_p = NA_real_, lb_p = NA_real_,
    nw_lag = NA_real_, hac_capped = NA, theta = NA_real_,
    b_sin12 = NA_real_, b_cos12 = NA_real_, b_sin6 = NA_real_, b_cos6 = NA_real_,
    converged = FALSE, note = NA_character_, stringsAsFactors = FALSE)
  
  tryCatch({
    d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
    d <- d[order(d$t), ]
    if (nrow(d) != nrow(covar) || anyNA(d$items))
      stop("series does not cover all 48 months")
    if (!offset_col %in% names(d)) stop("offset column not found: ", offset_col)
    d$off <- d[[offset_col]]
    
    fp <- glm(f_full, family = poisson, data = d)
    if (!isTRUE(fp$converged)) stop("Poisson full model did not converge")
    
    pear       <- residuals(fp, type = "pearson")
    disp_ratio <- sum(pear^2) / df.residual(fp)
    disp_p <- tryCatch(AER::dispersiontest(fp)$p.value, error = function(e) NA_real_)
    lb_p   <- tryCatch(Box.test(pear, lag = lb_lag, type = "Ljung-Box")$p.value,
                       error = function(e) NA_real_)
    res$disp_ratio <- disp_ratio; res$disp_p <- disp_p; res$lb_p <- lb_p
    
    overdisp <- !is.na(disp_p) && disp_p < alpha_disp
    autocorr <- !is.na(lb_p)   && lb_p   < alpha_lb
    
    if (autocorr) {
      h <- .do_hac(fp)
      b <- coef(fp)
      res$p_value <- h$p; res$distribution <- "poisson"; res$route <- "HAC-Wald"
      res$nw_lag <- h$bw; res$hac_capped <- h$capped
      
    } else if (overdisp) {
      fnb_full <- .nb_fit(f_full, d)
      fnb_red  <- .nb_fit(f_red,  d)
      if (is.null(fnb_full) || is.null(fnb_red)) {   # NB failed/not converged -> HAC
        h <- .do_hac(fp); b <- coef(fp)
        res$p_value <- h$p; res$distribution <- "poisson"; res$route <- "HAC-Wald(NBfail)"
        res$nw_lag <- h$bw; res$hac_capped <- h$capped
      } else {
        res$p_value <- lmtest::lrtest(fnb_red, fnb_full)[["Pr(>Chisq)"]][2]
        res$distribution <- "negbin"; res$route <- "NB-LRT"; res$theta <- fnb_full$theta
        b <- coef(fnb_full)
      }
      
    } else {
      fr <- glm(f_red, family = poisson, data = d)
      res$p_value <- lmtest::lrtest(fr, fp)[["Pr(>Chisq)"]][2]
      res$distribution <- "poisson"; res$route <- "Poisson-LRT"
      b <- coef(fp)
    }
    
    res$b_sin12 <- unname(b["sin12"]); res$b_cos12 <- unname(b["cos12"])
    res$b_sin6  <- unname(b["sin6"]);  res$b_cos6  <- unname(b["cos6"])
    res$converged <- TRUE
  }, error = function(e) res$note <<- conditionMessage(e))
  
  res
}


### 8. Run fit over every eligible drug class, applying BH at 5% FDR, reports route/distrbution properties ---------

stopifnot(exists("fit_test_series"), exists("covar"), exists("class_monthly_elig"))

screen_class <- class_monthly_elig |>
  select(bnf_class_code, bnf_class_name, year_month, items) |>
  group_by(bnf_class_code, bnf_class_name) |>
  group_modify(~ fit_test_series(.x, covar)) |>
  ungroup()

# BH across the class family; a class is significant at FDR 5%
stage3_legacy_inference <- FALSE
if (stage3_legacy_inference) {
  screen_class <- screen_class |>
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           significant = !is.na(p_adj) & p_adj < fdr_alpha) |>
    arrange(p_adj)
} else {
  screen_class <- screen_class |>
    mutate(class_q_bh = p.adjust(p_value, method = "BH"),
           class_significant = !is.na(class_q_bh) & class_q_bh < fdr_alpha,
           inference_scope = "primary_inferential",
           multiplicity_family = "all_eligible_classes") |>
    arrange(class_q_bh)
}

# proportions the methods commit to reporting
route_summary_class <- screen_class |>
  count(distribution, route, name = "n") |>
  mutate(pct = round(100 * n / sum(n), 1))

n_fail_class <- sum(!screen_class$converged)

saveRDS(screen_class,        file.path(data_dir, "screen_class.rds"))
fwrite(screen_class,         file.path(data_dir, "screen_class.csv"))
fwrite(route_summary_class,  file.path(data_dir, "screen_class_route_summary.csv"))

cat(sprintf("Class screen: %d classes | %d significant (BH 5%%) | %d non-converged\n",
            nrow(screen_class),
            if (stage3_legacy_inference) sum(screen_class$significant)
            else sum(screen_class$class_significant),
            n_fail_class))
print(route_summary_class)

# Interpretation
# Current outputs use class_q_bh/class_significant. The unqualified historical
# p_adj/significant names exist only in the Stage 3 reproduction branch.
# distribution = the distribution used for the modelling fit
# route = the inferential path the diagnostics selected: NB-LRT (overdispersion, no autocorrelation), HAC-Wald (residual autocorrelation), Poisson-LRT (neither overdispersion or autocorrelation)
# disp_ratio - dispersion diagnostic; disp_p - p value <0.05 if overdispersion
# lb_p - Ljung-Box p value for autocorrelation - p < 0.05 significant autocirrelation
# nw-lag - Newey-West bandwidth chosen
# theta - NB dispersion parameter
# sin/cos - harmonic coefficients - will convert to peak months later

### 9. Run fit over every eligible drug, applying BH at 5% FDR, reports route/distribution properties
# Current inference applies BH once across all eligible drugs. Parent-class
# significance is attached afterward as descriptive context only.

stopifnot(exists("fit_test_series"), exists("covar"),
          exists("drug_monthly_elig"), exists("screen_class"))

# fit every eligible drug
screen_drug <- drug_monthly_elig |>
  select(bnf_class_code, bnf_drug_code, bnf_drug_name, year_month, items) |>
  group_by(bnf_class_code, bnf_drug_code, bnf_drug_name) |>
  group_modify(~ fit_test_series(.x, covar)) |>
  ungroup()

if (stage3_legacy_inference) {
  sig_classes <- screen_class |> filter(significant) |> pull(bnf_class_code)
  screen_drug <- screen_drug |>
    mutate(parent_class_sig = bnf_class_code %in% sig_classes)

  # Historical selected-parent family, retained only for Stage 3 reproduction.
  sig_family <- screen_drug |>
    filter(parent_class_sig) |>
    mutate(p_adj = p.adjust(p_value, "BH")) |>
    select(bnf_drug_code, p_adj)

  screen_drug <- screen_drug |>
    left_join(sig_family, by = "bnf_drug_code") |>
    mutate(significant = parent_class_sig & !is.na(p_adj) & p_adj < fdr_alpha,
           p_adj_all   = p.adjust(p_value, "BH"),
           sig_all     = !is.na(p_adj_all) & p_adj_all < fdr_alpha) |>
    left_join(screen_class |> select(bnf_class_code, bnf_class_name, class_p_adj = p_adj),
              by = "bnf_class_code") |>
    arrange(!parent_class_sig, p_adj, p_adj_all)
} else {
  # Stage 4 authority: one BH family across every eligible drug. Parent-class
  # status is attached only after both complete testing families are defined.
  screen_drug <- screen_drug |>
    mutate(drug_all_q_bh = p.adjust(p_value, method = "BH"),
           drug_significant = !is.na(drug_all_q_bh) & drug_all_q_bh < fdr_alpha,
           inference_scope = "secondary_exploratory",
           multiplicity_family = "all_eligible_drugs") |>
    left_join(
      screen_class |>
        select(bnf_class_code, bnf_class_name, class_q_bh, class_significant),
      by = "bnf_class_code"
    ) |>
    mutate(parent_class_significant = class_significant)

  # Retained only as an explicitly named legacy audit field. It is not used for
  # characterisation, tables, figures, narrative counts or main conclusions.
  conditional_legacy <- screen_drug |>
    filter(parent_class_significant) |>
    transmute(
      bnf_drug_code,
      conditional_drug_q_legacy = p.adjust(p_value, method = "BH")
    )
  screen_drug <- screen_drug |>
    left_join(conditional_legacy, by = "bnf_drug_code") |>
    mutate(
      conditional_drug_significant_legacy = parent_class_significant &
        !is.na(conditional_drug_q_legacy) &
        conditional_drug_q_legacy < fdr_alpha
    ) |>
    arrange(drug_all_q_bh, bnf_drug_code)
}

route_summary_drug <- screen_drug |>
  count(distribution, route, name = "n") |>
  mutate(pct = round(100 * n / sum(n), 1))

n_fail_drug <- sum(!screen_drug$converged)

saveRDS(screen_drug,       file.path(data_dir, "screen_drug.rds"))
fwrite(screen_drug,        file.path(data_dir, "screen_drug.csv"))
fwrite(route_summary_drug, file.path(data_dir, "screen_drug_route_summary.csv"))

cat(sprintf(paste0(
  "Drug analysis: %d eligible drugs fitted | %d non-converged\n",
  if (stage3_legacy_inference)
    "  PRIMARY (in significant classes): %d of %d significant (pooled BH 5%%)\n"
  else
    "  LEGACY AUDIT (selected parent classes): %d of %d significant\n",
  "  SECONDARY EXPLORATORY (all drugs): %d of %d significant (BH 5%%)\n",
  "    of which in NON-significant classes: %d\n"),
  nrow(screen_drug), n_fail_drug,
  if (stage3_legacy_inference) sum(screen_drug$significant)
  else sum(screen_drug$conditional_drug_significant_legacy),
  if (stage3_legacy_inference) sum(screen_drug$parent_class_sig)
  else sum(screen_drug$parent_class_significant),
  if (stage3_legacy_inference) sum(screen_drug$sig_all)
  else sum(screen_drug$drug_significant),
  nrow(screen_drug),
  if (stage3_legacy_inference)
    sum(screen_drug$sig_all & !screen_drug$parent_class_sig)
  else
    sum(screen_drug$drug_significant & !screen_drug$parent_class_significant)))
print(route_summary_drug)

# Model-level failures remain visible in the screening tables and are also
# collected in one release-facing file.
model_failures <- bind_rows(
  screen_class |>
    filter(!converged) |>
    transmute(level = "class", code = as.character(bnf_class_code), note),
  screen_drug |>
    filter(!converged) |>
    transmute(level = "drug", code = as.character(bnf_drug_code), note)
)
fwrite(model_failures, file.path(data_dir, "model_failures.csv"))

### 10. Characterise significant seasonality ------
# peak:trough ratio (with 95% CI), calendar months of max and min,
# modality (one or two cycles of seasonality)
stopifnot(exists(".nb_fit"), exists(".f_full"), exists("covar"),
          exists("screen_class"), exists("screen_drug"),
          exists("class_monthly_elig"), exists("drug_monthly_elig"))

.f_1h <- items ~ trend1 + trend2 + trend3 + sin12 + cos12 + offset(off)
.harm <- c("sin12", "cos12", "sin6", "cos6")

# 12-month harmonic basis (columns ordered sin12, cos12, sin6, cos6) reused by
# both the point curve and the bootstrap
.harm_basis <- cbind(
  sin(2 * pi * (1:12) / harmonic_period_months[1]),
  cos(2 * pi * (1:12) / harmonic_period_months[1]),
  sin(2 * pi * (1:12) / harmonic_period_months[2]),
  cos(2 * pi * (1:12) / harmonic_period_months[2])
)

# seasonal curve over calendar months 1 to 12 from the four coefficients
# (t = 1 is January, so month m maps directly onto the harmonic arguments)
.seasonal_curve <- function(b_sin12, b_cos12, b_sin6, b_cos6, m = 1:12) {
  b_sin12 * sin(2 * pi * m / harmonic_period_months[1]) +
    b_cos12 * cos(2 * pi * m / harmonic_period_months[1]) +
    b_sin6 * sin(2 * pi * m / harmonic_period_months[2]) +
    b_cos6 * cos(2 * pi * m / harmonic_period_months[2])
}
# count local maxima of the 12-month curve, treating it as circular
.n_peaks <- function(s) {
  sum(s > c(s[12], s[-12]) & s > c(s[-1], s[1]))
}

# capped Newey-West HAC covariance (same bandwidth rule as the fit-test engine)
.hac_vcov <- function(fp, bw_cap = 12) {
  bw     <- tryCatch(sandwich::bwNeweyWest(fp, prewhite = 1L), error = function(e) NA_real_)
  bw_use <- if (is.na(bw)) bw_cap else min(bw, bw_cap)
  suppressWarnings(sandwich::NeweyWest(fp, lag = bw_use, prewhite = 1L, adjust = TRUE))
}

# 95% CI for the peak-to-trough ratio by parametric bootstrap of the four
# harmonic coefficients, drawn from the route-appropriate covariance (HAC for
# autocorrelated series, model-based otherwise). Non-smoothness of max/min and
# peak-month instability are handled naturally by recomputing the ratio per draw.
.ptr_ci <- function(f_full, route,
                    seed = coefficient_draw_seed,
                    B = coefficient_draw_count) {
  cf <- coef(f_full)
  if (!all(.harm %in% names(cf))) return(c(NA_real_, NA_real_))
  b <- cf[.harm]
  V <- tryCatch(if (startsWith(route, "HAC")) .hac_vcov(f_full) else vcov(f_full),
                error = function(e) NULL)
  if (is.null(V)) return(c(NA_real_, NA_real_))
  V4 <- V[.harm, .harm, drop = FALSE]; V4 <- (V4 + t(V4)) / 2   # symmetrise
  tryCatch({
    set.seed(seed)                                   # order-independent reproducibility
    draws <- MASS::mvrnorm(B, mu = b, Sigma = V4)
    S  <- draws %*% t(.harm_basis)                   # B x 12 seasonal curves
    mx <- S[cbind(seq_len(B), max.col( S, "first"))]
    mn <- S[cbind(seq_len(B), max.col(-S, "first"))]
    unname(quantile(exp(mx - mn), c(0.025, 0.975), na.rm = TRUE))
  }, error = function(e) c(NA_real_, NA_real_))
}

# seasonal reproducibility: does the within-year pattern repeat across years?
# Detrend (trend-only fit), then average the pairwise correlation between each
# year's 12-month profile. Genuine seasonality -> high (~0.8-1.0); a structural
# break or one-off level shift -> low (~0), because the pattern does not recur.
# Scale-invariant, so smooth growth/decline in a truly seasonal series is not
# penalised.
.seasonal_reproducibility <- function(series, covar) {
  tryCatch({
    d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
    d <- d[order(d$t), ]; d$off <- d$offset_log_patient_days
    ft <- glm(items ~ trend1 + trend2 + trend3 + offset(off), poisson, data = d)
    d$lr   <- log(d$items) - predict(ft, type = "link")   # detrended log-residual
    d$year <- d$year_month %/% 100L; d$mon <- d$year_month %% 100L
    years <- sort(unique(d$year))
    mat <- sapply(years, function(y) { dy <- d[d$year == y, ]; dy$lr[order(dy$mon)] })
    if (!is.matrix(mat) || ncol(mat) < 2) return(NA_real_)
    cc <- suppressWarnings(cor(mat))
    mean(cc[lower.tri(cc)], na.rm = TRUE)
  }, error = function(e) NA_real_)
}

# STL seasonal & trend strength (Wang, Smith & Hyndman): decompose the log-rate
# with a FLEXIBLE loess trend, then measure the variance share of each component.
# Because the trend is locally adaptive, a steep or accelerating trend is removed
# properly, so a trend-dominated series with no real cycle scores low on seasonal
# strength even though the rigid-spline reproducibility metric can be fooled by
# leftover trend curvature. Seasonal strength stays high only when a genuine
# recurring within-year cycle remains after flexible detrending.
.stl_strength <- function(series, covar,
                          offset_col = "offset_log_patient_days") {
  tryCatch({
    d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
    d <- d[order(d$t), ]
    if (!offset_col %in% names(d)) stop("offset column not found: ", offset_col)
    x <- ts(log(d$items) - d[[offset_col]], frequency = 12)
    fit <- stl(x, s.window = "periodic", robust = TRUE)
    cmp <- fit$time.series
    rem <- cmp[, "remainder"]; sea <- cmp[, "seasonal"]; tr <- cmp[, "trend"]
    c(seasonal = max(0, 1 - var(rem) / var(sea + rem)),
      trend    = max(0, 1 - var(rem) / var(tr  + rem)))
  }, error = function(e) c(seasonal = NA_real_, trend = NA_real_))
}

characterise_one <- function(scr_row, series, covar) {
  s   <- .seasonal_curve(scr_row$b_sin12, scr_row$b_cos12,
                         scr_row$b_sin6,  scr_row$b_cos6)
  out <- tibble(
    peak_trough_ratio = exp(max(s) - min(s)),
    ptr_lci = NA_real_, ptr_uci = NA_real_,
    peak_month        = month.abb[which.max(s)],
    trough_month      = month.abb[which.min(s)],
    amp_annual        = sqrt(scr_row$b_sin12^2 + scr_row$b_cos12^2),
    amp_semiannual    = sqrt(scr_row$b_sin6^2  + scr_row$b_cos6^2),
    n_peaks           = .n_peaks(s),
    seasonal_reproducibility = NA_real_,
    stl_seasonal_strength = NA_real_, stl_trend_strength = NA_real_,
    aic_1h = NA_real_, aic_2h = NA_real_, modality = NA_character_)
  
  # refit 1- and 2-harmonic models under the chosen distribution (for AIC
  # modality) and reuse the 2-harmonic fit for the peak:trough CI
  d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
  d <- d[order(d$t), ]; d$off <- d$offset_log_patient_days
  fits <- if (identical(scr_row$distribution, "negbin")) {
    f2 <- .nb_fit(.f_full, d); f1 <- .nb_fit(.f_1h, d)
    if (is.null(f1) || is.null(f2)) NULL else list(f1 = f1, f2 = f2)
  } else NULL
  if (is.null(fits))   # poisson series, or NB refit failed -> Poisson AICs
    fits <- list(f1 = glm(.f_1h,  poisson, data = d),
                 f2 = glm(.f_full, poisson, data = d))
  
  out$aic_1h <- AIC(fits$f1); out$aic_2h <- AIC(fits$f2)
  out$modality <- if (out$aic_1h <= out$aic_2h) "unimodal (1 harmonic preferred)"
  else if (out$n_peaks >= 2)    "bimodal (2 harmonics, 2 peaks)"
  else                          "unimodal, non-sinusoidal (2 harmonics, 1 peak)"
  
  ci <- .ptr_ci(fits$f2, scr_row$route)
  out$ptr_lci <- ci[1]; out$ptr_uci <- ci[2]
  out$seasonal_reproducibility <- .seasonal_reproducibility(series, covar)
  st <- .stl_strength(series, covar)
  out$stl_seasonal_strength <- unname(st["seasonal"])
  out$stl_trend_strength    <- unname(st["trend"])
  out
}

characterise_level <- function(scr, monthly, id_cols, keep_flag) {
  sig <- scr |> filter({{keep_flag}})
  if (nrow(sig) == 0) return(tibble())
  sig |>
    group_by(across(all_of(id_cols))) |>
    group_modify(function(row, key) {
      ser <- monthly |> semi_join(key, by = id_cols)
      characterise_one(row, ser, covar)
    }) |>
    ungroup() |>
    left_join(scr, by = id_cols) |>
    arrange(desc(peak_trough_ratio))
}

if (stage3_legacy_inference) {
  char_class <- characterise_level(screen_class, class_monthly_elig,
                                   c("bnf_class_code"), significant)
  char_drug  <- characterise_level(screen_drug, drug_monthly_elig,
                                   c("bnf_drug_code"), significant | sig_all)
} else {
  char_class <- characterise_level(screen_class, class_monthly_elig,
                                   c("bnf_class_code"), class_significant)
  char_drug  <- characterise_level(screen_drug, drug_monthly_elig,
                                   c("bnf_drug_code"), drug_significant)
}

saveRDS(char_class, file.path(data_dir, "characterisation_class.rds"))
saveRDS(char_drug,  file.path(data_dir, "characterisation_drug.rds"))
fwrite(char_class,  file.path(data_dir, "characterisation_class.csv"))
fwrite(char_drug,   file.path(data_dir, "characterisation_drug.csv"))

cat(sprintf(paste0(
  "Characterisation complete.\n",
  "  Classes: %d characterised | peak months: %s\n",
  "  Drugs:   %d characterised (significant or exploratory-significant)\n"),
  nrow(char_class),
  paste(head(names(sort(table(char_class$peak_month), decreasing = TRUE)), 3), collapse = ", "),
  nrow(char_drug)))
cat("\nTop 10 classes by amplitude (peak:trough ratio with 95% CI):\n")
print(char_class |> slice_head(n = 10) |>
        mutate(ptr = sprintf("%.2f (%.2f-%.2f)", peak_trough_ratio, ptr_lci, ptr_uci)) |>
        select(bnf_class_name, ptr, peak_month, modality,
               all_of(if (stage3_legacy_inference) "p_adj" else "class_q_bh")))


### 11. Reporting --------

suppressMessages(library(ggplot2))
stopifnot(exists("char_class"), exists("char_drug"), exists(".seasonal_curve"),
          exists(".f_full"), exists(".nb_fit"), exists("covar"),
          exists("class_monthly_elig"))

# "Meaningful" seasonality requires BOTH: an appreciable, precisely-estimated
# amplitude (lower 95% CI of the peak:trough ratio at/above the threshold) AND a
# genuine recurring within-year cycle that survives flexible detrending (STL
# seasonal strength). STL uses a locally-adaptive trend, so it excludes BOTH
# structural breaks and trend-dominated series (e.g. a steeply rising drug whose
# apparent seasonality is trend leakage) - the two failure modes the harmonic
# model can misread. seasonal_reproducibility is retained as a reported
# diagnostic (it agrees on breaks but can be fooled by trend curvature).
# NOTE: calibrated to real data - appreciable-amplitude classes split cleanly
# into an excluded floor (breaks/trend-dominated, seasonal strength <= 0.36) and
# a genuine group (>= 0.66), with an empty band between. 0.50 sits mid-gap, so
# the cut is robust anywhere in 0.40-0.60. stl_trend_strength is NOT a cutoff
# (it is high for genuine strongly-trending seasonal classes too) but is kept as
# a descriptive column for characterising trend-dominated series in the text.
meaningful_threshold <- amplitude_lci_threshold

fig_dir <- file.path(data_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

## Results table - classes -----------------------------------------------

if (stage3_legacy_inference) {
  results_class <- char_class |>
    transmute(bnf_class_code, bnf_class_name,
              peak_trough_ratio = round(peak_trough_ratio, 3),
              ptr_lci = round(ptr_lci, 3), ptr_uci = round(ptr_uci, 3),
              peak_month, trough_month, modality, n_peaks,
              seasonal_reproducibility = round(seasonal_reproducibility, 3),
              stl_seasonal_strength = round(stl_seasonal_strength, 3),
              stl_trend_strength    = round(stl_trend_strength, 3),
              distribution, route, hac_capped, p_adj,
              meaningful = ptr_lci >= meaningful_threshold &
                stl_seasonal_strength >= stl_strength_threshold) |>
    arrange(desc(meaningful), desc(peak_trough_ratio))
} else {
  results_class <- char_class |>
    transmute(bnf_class_code, bnf_class_name,
              inference_scope = "primary_inferential",
              multiplicity_family = "all_eligible_classes",
              peak_trough_ratio = round(peak_trough_ratio, 3),
              ptr_lci = round(ptr_lci, 3), ptr_uci = round(ptr_uci, 3),
              peak_month, trough_month, modality, n_peaks,
              seasonal_reproducibility = round(seasonal_reproducibility, 3),
              stl_seasonal_strength = round(stl_seasonal_strength, 3),
              stl_trend_strength    = round(stl_trend_strength, 3),
              distribution, route, hac_capped, class_q_bh, class_significant,
              meaningful = ptr_lci >= meaningful_threshold &
                stl_seasonal_strength >= stl_strength_threshold) |>
    arrange(desc(meaningful), desc(peak_trough_ratio))
}

## Results table - drugs -------------------------------------------------

if (stage3_legacy_inference) {
  results_drug <- char_drug |>
    transmute(bnf_drug_code, bnf_drug_name, bnf_class_name,
              peak_trough_ratio = round(peak_trough_ratio, 3),
              ptr_lci = round(ptr_lci, 3), ptr_uci = round(ptr_uci, 3),
              peak_month, trough_month, modality,
              seasonal_reproducibility = round(seasonal_reproducibility, 3),
              stl_seasonal_strength = round(stl_seasonal_strength, 3),
              stl_trend_strength    = round(stl_trend_strength, 3),
              distribution, route,
              p_adj_primary = p_adj, significant_primary = significant,
              p_adj_all, sig_all, parent_class_sig,
              meaningful = ptr_lci >= meaningful_threshold &
                stl_seasonal_strength >= stl_strength_threshold) |>
    arrange(desc(meaningful), desc(peak_trough_ratio))
} else {
  results_drug <- char_drug |>
    transmute(bnf_drug_code, bnf_drug_name, bnf_class_name,
              inference_scope = "secondary_exploratory",
              multiplicity_family = "all_eligible_drugs",
              peak_trough_ratio = round(peak_trough_ratio, 3),
              ptr_lci = round(ptr_lci, 3), ptr_uci = round(ptr_uci, 3),
              peak_month, trough_month, modality,
              seasonal_reproducibility = round(seasonal_reproducibility, 3),
              stl_seasonal_strength = round(stl_seasonal_strength, 3),
              stl_trend_strength    = round(stl_trend_strength, 3),
              distribution, route, drug_all_q_bh, drug_significant,
              parent_class_significant, class_q_bh,
              conditional_drug_q_legacy,
              conditional_drug_significant_legacy,
              meaningful = ptr_lci >= meaningful_threshold &
                stl_seasonal_strength >= stl_strength_threshold) |>
    arrange(desc(meaningful), desc(peak_trough_ratio))
}

## Figures ---------------------------------------------------------------

# observed + fitted monthly rate (per 1000 registered patients) for one class
.fit_frame <- function(code) {
  ser  <- class_monthly_elig |> filter(bnf_class_code == code) |> select(year_month, items)
  d    <- covar |> left_join(ser, by = "year_month") |> arrange(t)
  d$off <- d$offset_log_patient_days
  dist <- char_class$distribution[char_class$bnf_class_code == code][1]
  fit  <- if (identical(dist, "negbin")) .nb_fit(.f_full, d) else NULL
  if (is.null(fit)) fit <- glm(.f_full, poisson, data = d)
  nm <- char_class$bnf_class_name[char_class$bnf_class_code == code][1]
  d |> transmute(class = nm, month_date,
                 observed = items / list_size * 1000,
                 fitted   = as.numeric(fitted(fit)) / list_size * 1000)
}

# Exemplars: strongest MEANINGFUL seasonality (so structural-break and trivial
# classes are excluded by construction), excluding the extreme vaccine/antiviral
# series whose amplitude would flatten the shared axes
exemplar_class_codes <- results_class |>
  filter(meaningful, peak_trough_ratio < 10) |>
  slice_head(n = 6) |> pull(bnf_class_code)

if (length(exemplar_class_codes)) {
  fit_df <- bind_rows(lapply(exemplar_class_codes, .fit_frame)) %>%
    arrange(class)
  fit_df$class <- factor(fit_df$class,
                         levels = results_class |> filter(bnf_class_code %in% exemplar_class_codes) |>
                           arrange(desc(peak_trough_ratio)) |> pull(bnf_class_name))
  fit_df$class <- factor(fit_df$class,
                         levels = sort(unique(as.character(fit_df$class))))
  
  p_fit <- ggplot(fit_df, aes(month_date)) +
    geom_point(aes(y = observed), size = 0.9, colour = "grey45") +
    geom_line(aes(y = fitted), colour = "#B2182B", linewidth = 0.7) +
    facet_wrap(~ class, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Items per 1000 registered patients",
         title = "Observed monthly prescribing with fitted trend and seasonal model") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, "fig_observed_fitted.png"), p_fit,
         width = 9, height = 8, dpi = 300)
  
  # seasonal multiplicative factor by calendar month (from harmonic coefficients)
  shape_df <- char_class |>
    filter(bnf_class_code %in% exemplar_class_codes) |>
    group_by(bnf_class_name) |>
    reframe(month = 1:12,
            factor = { s <- .seasonal_curve(b_sin12, b_cos12, b_sin6, b_cos6, 1:12)
            f <- exp(s); f / mean(f) })
  p_shape <- ggplot(shape_df, aes(month, factor)) +
    geom_hline(yintercept = 1, colour = "grey70", linetype = 2) +
    geom_line(colour = "#2166AC", linewidth = 0.7) +
    facet_wrap(~ bnf_class_name, ncol = 2) +
    scale_x_continuous(breaks = c(1,4,7,10), labels = month.abb[c(1,4,7,10)]) +
    labs(x = NULL, y = "Seasonal factor (relative to annual mean)",
         title = "Fitted seasonal shape by calendar month") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, "fig_seasonal_shape.png"), p_shape,
         width = 9, height = 8, dpi = 300)
}

## Appendix assembly -----------------------------------------------------

# distribution / route proportions across both levels
appendix_routes <- bind_rows(
  screen_class |> count(distribution, route, name = "n") |> mutate(level = "class"),
  screen_drug  |> count(distribution, route, name = "n") |> mutate(level = "drug")) |>
  group_by(level) |> mutate(pct = round(100 * n / sum(n), 1)) |> ungroup() |>
  select(level, distribution, route, n, pct)

# recode crosswalk (if reconciliation was applied this session)
if (exists("xwalk_drug"))  fwrite(xwalk_drug,  file.path(data_dir, "appendix_recode_crosswalk_drug.csv"))
if (exists("xwalk_class")) fwrite(xwalk_class, file.path(data_dir, "appendix_recode_crosswalk_class.csv"))

## Save ------------------------------------------------------------------

fwrite(results_class,   file.path(data_dir, "results_class.csv"))
fwrite(results_drug,    file.path(data_dir, "results_drug.csv"))
fwrite(appendix_routes, file.path(data_dir, "appendix_route_proportions.csv"))
saveRDS(results_class,  file.path(data_dir, "results_class.rds"))
saveRDS(results_drug,   file.path(data_dir, "results_drug.rds"))
if (!stage3_legacy_inference) {
  inference_scope_metadata <- data.table(
    level = c("class", "drug", "hierarchy"),
    analysis_role = c("primary_inferential", "secondary_exploratory", "not_applicable"),
    multiplicity_family = c(
      "all_eligible_classes", "all_eligible_drugs", "no_joint_hierarchical_family"
    ),
    family_size = c(nrow(screen_class), nrow(screen_drug), NA_integer_),
    adjusted_value = c("class_q_bh", "drug_all_q_bh", NA_character_),
    discovery_flag = c("class_significant", "drug_significant", NA_character_),
    interpretation = c(
      "BH FDR 5% within the complete eligible-class family.",
      "BH FDR 5% within the complete eligible-drug family; hypothesis-generating.",
      "Neither family is claimed to provide 5% FDR control across the whole hierarchy."
    )
  )
  atomic_fwrite(inference_scope_metadata, file.path(data_dir, "inference_scope_metadata.csv"))
}

cat(sprintf(paste0(
  "Reporting outputs complete.\n",
  "  Meaningful = peak:trough lower CI >= %.2f AND STL seasonal strength >= %.2f\n",
  "  Classes: %d significant | %d meaningful\n",
  "  Drugs:   %d characterised | %d meaningful\n",
  "  Figures: %d exemplar classes -> %s\n",
  "  Tables:  results_class.csv, results_drug.csv, appendix_route_proportions.csv\n"),
  meaningful_threshold, stl_strength_threshold,
  nrow(results_class), sum(results_class$meaningful),
  nrow(results_drug), sum(results_drug$meaningful),
  length(exemplar_class_codes), fig_dir))
cat("\nMain-text classes (meaningful seasonality, by amplitude):\n")
print(results_class |> filter(meaningful) |>
        mutate(ptr = sprintf("%.2f (%.2f-%.2f)", peak_trough_ratio, ptr_lci, ptr_uci)) |>
        select(bnf_class_name, ptr, peak_month, modality, stl_seasonal_strength) |>
        slice_head(n = 20))


### 12. Publication outputs (tables and figures) ------------------------------
# Writes publication-ready tables (meaningful column names, only the columns
# needed) and figures into <data_dir>/results. Main-text and appendix material
# are kept separate. Appendix includes observed+fitted panels for EVERY eligible
# class and drug, paginated by BNF chapter, for visual inspection.
#
# Requires Sections 4-11 objects: covar, class_monthly_elig, drug_monthly_elig,
# elig_class, elig_drug, screen_class, screen_drug, char_class, char_drug,
# results_class, results_drug, coverage, excl_class, excl_drug, appendix_routes,
# and the helpers .f_full, .nb_fit, .seasonal_curve. xwalk_drug is optional.

suppressMessages(library(ggplot2))
stopifnot(exists("results_class"), exists("results_drug"), exists("char_class"),
          exists("covar"), exists("class_monthly_elig"), exists("drug_monthly_elig"),
          exists(".f_full"), exists(".seasonal_curve"))
has_repel <- requireNamespace("ggrepel", quietly = TRUE)

res_dir <- file.path(data_dir, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
for (d in c(res_dir, tab_dir, fig_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- formatting helpers -----------------------------------------------------

fmt_p  <- function(p) ifelse(is.na(p), NA_character_,
                             ifelse(p < 0.001, "<0.001", formatC(signif(p, 3), format = "g")))
fmt_ci <- function(est, lo, hi) sprintf("%.2f (%.2f\u2013%.2f)", est, lo, hi)

relabel_modality <- function(m) dplyr::case_when(
  grepl("^bimodal", m)                 ~ "Bimodal",
  grepl("non-sinusoidal", m)           ~ "Unimodal (non-sinusoidal)",
  grepl("^unimodal", m)                ~ "Unimodal",
  TRUE                                 ~ m)
relabel_route <- function(r) dplyr::case_when(
  grepl("^HAC", r)          ~ "Poisson + Newey\u2013West HAC (Wald)",
  r == "NB-LRT"             ~ "Negative binomial (LRT)",
  r == "Poisson-LRT"        ~ "Poisson (LRT)",
  TRUE                      ~ r)

# BNF class (paragraph) codes are stored internally as integers, so a leading
# zero (chapters 1-9) is lost and long codes could take scientific form. Render
# them zero-padded and DOTTED (e.g. "03.01.01" = chapter.section.paragraph):
# the dots make the value intrinsically non-numeric, so it survives re-reading
# in Excel or fread without being coerced back to a number. Drug (chemical
# substance) codes contain letters so are already safe as character; the few
# all-digit ones (e.g. some nutrition codes) are kept verbatim - read any
# re-imported table with colClasses = "character" to preserve them exactly.
bnf_class_dotted <- function(code) {
  s <- sprintf("%06d", as.integer(code))
  sub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1.\\2.\\3", s)
}

# BNF chapter lookup from the eligible monthly frames
chap_class <- class_monthly_elig |>
  distinct(bnf_class_code, bnf_chapter_code, bnf_chapter_name)
chap_drug <- drug_monthly_elig |>
  distinct(bnf_drug_code, bnf_chapter_code, bnf_chapter_name)

## ---- MAIN TEXT: Table 1, meaningful-seasonal classes ------------------------

if (stage3_legacy_inference) {
  main_table_classes <- results_class |>
    filter(meaningful) |>
    transmute(
      `BNF class code`                = bnf_class_dotted(bnf_class_code),
      `BNF class`                     = bnf_class_name,
      `Peak-to-trough ratio (95% CI)` = fmt_ci(peak_trough_ratio, ptr_lci, ptr_uci),
      `Peak month`                    = peak_month,
      `Trough month`                  = trough_month,
      `Seasonal shape`                = relabel_modality(modality),
      `Seasonal strength`             = sprintf("%.2f", stl_seasonal_strength),
      `Adjusted p`                    = fmt_p(p_adj))
} else {
  main_table_classes <- results_class |>
    filter(meaningful) |>
    transmute(
      `BNF class code`                = bnf_class_dotted(bnf_class_code),
      `BNF class`                     = bnf_class_name,
      `Peak-to-trough ratio (95% CI)` = fmt_ci(peak_trough_ratio, ptr_lci, ptr_uci),
      `Peak month`                    = peak_month,
      `Trough month`                  = trough_month,
      `Seasonal shape`                = relabel_modality(modality),
      `Seasonal strength`             = sprintf("%.2f", stl_seasonal_strength),
      `Class q (BH; all eligible classes)` = fmt_p(class_q_bh))
}
fwrite(main_table_classes, file.path(tab_dir, "table1_meaningful_classes.csv"))

## ---- MAIN TEXT: Table 2, meaningful-seasonal drugs -------------------------
# Parallel to Table 1 but drug-level. Stage 4 uses only the complete all-drug
# family; parent-class significance is descriptive context, not a selection rule.
if (stage3_legacy_inference) {
  main_table_drugs <- results_drug |>
    filter(meaningful) |>
    transmute(
      `BNF class code`                = sub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1.\\2.\\3",
                                            substr(as.character(bnf_drug_code), 1, 6)),
      `BNF class`                     = bnf_class_name,
      `BNF drug code`                 = as.character(bnf_drug_code),
      `Drug (chemical substance)`     = bnf_drug_name,
      `Peak-to-trough ratio (95% CI)` = fmt_ci(peak_trough_ratio, ptr_lci, ptr_uci),
      `Peak month`                    = peak_month,
      `Trough month`                  = trough_month,
      `Seasonal shape`                = relabel_modality(modality),
      `Seasonal strength`             = sprintf("%.2f", stl_seasonal_strength),
      `In significant class`          = ifelse(parent_class_sig, "Yes", "No"),
      `Adjusted p`                    = fmt_p(dplyr::coalesce(p_adj_primary, p_adj_all)))
} else {
  main_table_drugs <- results_drug |>
    filter(meaningful) |>
    transmute(
      `BNF class code`                = sub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1.\\2.\\3",
                                            substr(as.character(bnf_drug_code), 1, 6)),
      `BNF class`                     = bnf_class_name,
      `BNF drug code`                 = as.character(bnf_drug_code),
      `Drug (chemical substance)`     = bnf_drug_name,
      `Peak-to-trough ratio (95% CI)` = fmt_ci(peak_trough_ratio, ptr_lci, ptr_uci),
      `Peak month`                    = peak_month,
      `Trough month`                  = trough_month,
      `Seasonal shape`                = relabel_modality(modality),
      `Seasonal strength`             = sprintf("%.2f", stl_seasonal_strength),
      `Parent class significant`      = ifelse(parent_class_significant, "Yes", "No"),
      `Drug q (BH; all eligible drugs)` = fmt_p(drug_all_q_bh))
}
fwrite(main_table_drugs, file.path(tab_dir, "table2_meaningful_drugs.csv"))

## ---- MAIN TEXT: analytic accounting (funnel) -------------------------------

if (stage3_legacy_inference) {
  sig_class_n <- sum(screen_class$significant, na.rm = TRUE)
  sig_drug_n <- sum(screen_drug$significant, na.rm = TRUE)
  accounting <- tibble::tibble(
    Stage = c("Eligible series",
              "Statistically significant (BH FDR 5%)",
              "Significant but trivial amplitude (ratio < 1.05)",
              "Meaningful seasonality (amplitude CI \u2265 threshold and seasonal strength \u2265 threshold)"),
    Classes = c(sum(elig_class$eligible), sig_class_n,
                sum(results_class$peak_trough_ratio < 1.05, na.rm = TRUE),
                sum(results_class$meaningful, na.rm = TRUE)),
    Drugs = c(sum(elig_drug$eligible), sig_drug_n,
              sum(results_drug$peak_trough_ratio < 1.05, na.rm = TRUE),
              sum(results_drug$meaningful, na.rm = TRUE)))
} else {
  sig_class_n <- sum(screen_class$class_significant, na.rm = TRUE)
  sig_drug_n <- sum(screen_drug$drug_significant, na.rm = TRUE)
  accounting <- tibble::tibble(
    Stage = c("Eligible series",
              "Statistically significant within stated complete family (BH FDR 5%)",
              "Significant but trivial amplitude (ratio < 1.05)",
              "Meaningful seasonality (amplitude CI \u2265 threshold and seasonal strength \u2265 threshold)"),
    Classes = c(sum(elig_class$eligible), sig_class_n,
                sum(results_class$peak_trough_ratio < 1.05, na.rm = TRUE),
                sum(results_class$meaningful, na.rm = TRUE)),
    Drugs = c(sum(elig_drug$eligible), sig_drug_n,
              sum(results_drug$peak_trough_ratio < 1.05, na.rm = TRUE),
              sum(results_drug$meaningful, na.rm = TRUE)),
    `Class inference scope` = rep("Primary; all eligible classes", 4),
    `Drug inference scope` = rep("Secondary exploratory; all eligible drugs", 4))
}
fwrite(accounting, file.path(tab_dir, "table_accounting.csv"))

## ---- APPENDIX: full class and drug results ---------------------------------

if (stage3_legacy_inference) {
  appendix_table_classes <- results_class |>
    left_join(chap_class, by = "bnf_class_code") |>
    transmute(
      `BNF chapter` = bnf_chapter_name, `BNF class code` = bnf_class_dotted(bnf_class_code),
      `BNF class` = bnf_class_name, `Peak-to-trough ratio` = sprintf("%.3f", peak_trough_ratio),
      `95% CI lower` = sprintf("%.3f", ptr_lci), `95% CI upper` = sprintf("%.3f", ptr_uci),
      `Peak month` = peak_month, `Trough month` = trough_month,
      `Seasonal shape` = relabel_modality(modality),
      `Seasonal strength (STL)` = sprintf("%.3f", stl_seasonal_strength),
      `Trend strength (STL)` = sprintf("%.3f", stl_trend_strength),
      `Cross-year reproducibility` = sprintf("%.3f", seasonal_reproducibility),
      `Distribution` = tools::toTitleCase(distribution),
      `Inference method` = relabel_route(route),
      `HAC bandwidth capped` = ifelse(is.na(hac_capped), "", ifelse(hac_capped, "Yes", "No")),
      `Adjusted p` = fmt_p(p_adj), `Meaningful seasonality` = ifelse(meaningful, "Yes", "No"))
} else {
  appendix_table_classes <- results_class |>
    left_join(chap_class, by = "bnf_class_code") |>
    transmute(
      `BNF chapter` = bnf_chapter_name, `BNF class code` = bnf_class_dotted(bnf_class_code),
      `BNF class` = bnf_class_name, `Inference scope` = "Primary inferential",
      `Peak-to-trough ratio` = sprintf("%.3f", peak_trough_ratio),
      `95% CI lower` = sprintf("%.3f", ptr_lci), `95% CI upper` = sprintf("%.3f", ptr_uci),
      `Peak month` = peak_month, `Trough month` = trough_month,
      `Seasonal shape` = relabel_modality(modality),
      `Seasonal strength (STL)` = sprintf("%.3f", stl_seasonal_strength),
      `Trend strength (STL)` = sprintf("%.3f", stl_trend_strength),
      `Cross-year reproducibility` = sprintf("%.3f", seasonal_reproducibility),
      `Distribution` = tools::toTitleCase(distribution),
      `Inference method` = relabel_route(route),
      `HAC bandwidth capped` = ifelse(is.na(hac_capped), "", ifelse(hac_capped, "Yes", "No")),
      `Class q (BH; all eligible classes)` = fmt_p(class_q_bh),
      `Meaningful seasonality` = ifelse(meaningful, "Yes", "No"))
}
fwrite(appendix_table_classes, file.path(tab_dir, "appendixA1_all_classes.csv"))

if (stage3_legacy_inference) {
  appendix_table_drugs <- results_drug |>
    left_join(chap_drug, by = "bnf_drug_code") |>
    transmute(
      `BNF chapter` = bnf_chapter_name,
      `BNF class code` = sub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1.\\2.\\3",
                             substr(as.character(bnf_drug_code), 1, 6)),
      `BNF class` = bnf_class_name, `BNF drug code` = as.character(bnf_drug_code),
      `Drug (chemical substance)` = bnf_drug_name,
      `Peak-to-trough ratio` = sprintf("%.3f", peak_trough_ratio),
      `95% CI lower` = sprintf("%.3f", ptr_lci), `95% CI upper` = sprintf("%.3f", ptr_uci),
      `Peak month` = peak_month, `Trough month` = trough_month,
      `Seasonal shape` = relabel_modality(modality),
      `Seasonal strength (STL)` = sprintf("%.3f", stl_seasonal_strength),
      `Trend strength (STL)` = sprintf("%.3f", stl_trend_strength),
      `Distribution` = tools::toTitleCase(distribution), `Inference method` = relabel_route(route),
      `In significant class` = ifelse(parent_class_sig, "Yes", "No"),
      `Significant (within class)` = ifelse(significant_primary, "Yes", "No"),
      `Adjusted p (within class)` = fmt_p(p_adj_primary),
      `Significant (all-drug scan)` = ifelse(sig_all, "Yes", "No"),
      `Adjusted p (all-drug scan)` = fmt_p(p_adj_all),
      `Meaningful seasonality` = ifelse(meaningful, "Yes", "No"))
} else {
  appendix_table_drugs <- results_drug |>
    left_join(chap_drug, by = "bnf_drug_code") |>
    transmute(
      `BNF chapter` = bnf_chapter_name,
      `BNF class code` = sub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1.\\2.\\3",
                             substr(as.character(bnf_drug_code), 1, 6)),
      `BNF class` = bnf_class_name, `BNF drug code` = as.character(bnf_drug_code),
      `Drug (chemical substance)` = bnf_drug_name,
      `Inference scope` = "Secondary exploratory",
      `Peak-to-trough ratio` = sprintf("%.3f", peak_trough_ratio),
      `95% CI lower` = sprintf("%.3f", ptr_lci), `95% CI upper` = sprintf("%.3f", ptr_uci),
      `Peak month` = peak_month, `Trough month` = trough_month,
      `Seasonal shape` = relabel_modality(modality),
      `Seasonal strength (STL)` = sprintf("%.3f", stl_seasonal_strength),
      `Trend strength (STL)` = sprintf("%.3f", stl_trend_strength),
      `Distribution` = tools::toTitleCase(distribution), `Inference method` = relabel_route(route),
      `Parent class significant` = ifelse(parent_class_significant, "Yes", "No"),
      `Drug q (BH; all eligible drugs)` = fmt_p(drug_all_q_bh),
      `Meaningful seasonality` = ifelse(meaningful, "Yes", "No"))
}
fwrite(appendix_table_drugs, file.path(tab_dir, "appendixA2_all_drugs.csv"))

# coverage, exclusions, crosswalk, route proportions
cov_one <- function(tab, lvl) {
  elig <- tab$eligible %in% TRUE                    # NA / non-TRUE counted as excluded
  tibble(level = lvl,
         total_series        = nrow(tab),
         eligible            = sum(elig),
         excluded            = sum(!elig),
         excluded_item_share = sum(tab$total_items[!elig]) / sum(tab$total_items))
}
coverage <- bind_rows(cov_one(elig_class, "class"), cov_one(elig_drug, "drug"))

if (exists("coverage"))   fwrite(coverage,   file.path(tab_dir, "appendix_coverage_summary.csv"))
if (exists("excl_class") && nrow(excl_class))
  fwrite(excl_class |> mutate(bnf_class_code = bnf_class_dotted(bnf_class_code)),
         file.path(tab_dir, "appendix_exclusions_class.csv"))
if (exists("excl_drug") && nrow(excl_drug))
  fwrite(excl_drug |> mutate(bnf_drug_code = as.character(bnf_drug_code)),
         file.path(tab_dir, "appendix_exclusions_drug.csv"))
if (exists("xwalk_drug")) fwrite(xwalk_drug, file.path(tab_dir, "appendix_recode_crosswalk.csv"))
if (exists("appendix_routes")) fwrite(appendix_routes, file.path(tab_dir, "appendix_route_proportions.csv"))

## ---- STAGE 3: reproduce and seal the pre-change baseline -------------------
# Stage 3 deliberately stops after the publication tables. The extensive
# all-series plotting and separate 2021 window diagnostics below do not alter
# the cohort or inference and are not needed to establish this baseline.
if (run_stage == "stage3") {
  stage3_dir <- file.path(out_dir, "qc", "stage3")
  snapshot_dir <- file.path(stage3_dir, "baseline_snapshot")
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  baseline_dir <- normalizePath(
    Sys.getenv(
      "STOCK2026_BASELINE_DIR",
      unset = file.path(analysis_dir, "Outputs")
    ),
    mustWork = TRUE
  )
  if (identical(normalizePath(out_dir, mustWork = TRUE), baseline_dir)) {
    stop("Stage 3 baseline directory must differ from the canonical output directory.")
  }

  stage3_checks <- list()
  field_comparisons <- list()
  add_stage3_check <- function(check_id, scope, pass, expected, observed, details = "") {
    stage3_checks[[length(stage3_checks) + 1L]] <<- data.table(
      check_id = check_id,
      scope = scope,
      pass = isTRUE(pass),
      expected = as.character(expected),
      observed = as.character(observed),
      details = as.character(details)
    )
  }

  key_string <- function(x, key) {
    do.call(
      paste,
      c(lapply(x[, ..key], function(z) ifelse(is.na(z), "<NA>", as.character(z))),
        sep = "\r")
    )
  }

  normalise_comparison_codes <- function(x) {
    if ("bnf_class_code" %in% names(x)) {
      z <- trimws(as.character(x$bnf_class_code))
      numeric_code <- grepl("^[0-9]{1,6}$", z)
      z[numeric_code] <- sprintf("%06d", as.integer(z[numeric_code]))
      x[, bnf_class_code := z]
    }
    x
  }

  compare_frames <- function(scope, current, reference, key, fields = NULL,
                             rel_tol = 1e-10, abs_tol = 1e-14,
                             allowed_added_keys = character(),
                             allowed_removed_keys = character()) {
    cur <- as.data.table(copy(current))
    ref <- as.data.table(copy(reference))
    cur <- normalise_comparison_codes(cur)
    ref <- normalise_comparison_codes(ref)
    missing_key <- setdiff(key, intersect(names(cur), names(ref)))
    if (length(missing_key)) {
      add_stage3_check(
        paste0(scope, "_keys"), scope, FALSE,
        paste(key, collapse = ";"),
        paste(intersect(names(cur), names(ref)), collapse = ";"),
        paste("Missing comparison key(s):", paste(missing_key, collapse = ", "))
      )
      return(invisible(FALSE))
    }

    cur_key <- key_string(cur, key)
    ref_key <- key_string(ref, key)
    missing_keys <- setdiff(ref_key, cur_key)
    added_keys <- setdiff(cur_key, ref_key)
    keys_ok <- !anyDuplicated(cur_key) && !anyDuplicated(ref_key) &&
      setequal(missing_keys, allowed_removed_keys) &&
      setequal(added_keys, allowed_added_keys)
    add_stage3_check(
      paste0(scope, "_keys"), scope, keys_ok,
      sprintf("%d reference keys; allowed removed=%s; allowed added=%s",
              length(ref_key), paste(allowed_removed_keys, collapse = ";"),
              paste(allowed_added_keys, collapse = ";")),
      sprintf("%d current keys; removed=%s; added=%s", length(cur_key),
              paste(missing_keys, collapse = ";"), paste(added_keys, collapse = ";")),
      if (keys_ok && !length(c(missing_keys, added_keys))) "Key membership agrees exactly."
      else if (keys_ok) "The pre-declared boundary membership difference is present and no other key differs."
      else
        sprintf("Missing=%d; added=%d; current duplicates=%d; reference duplicates=%d",
                length(missing_keys), length(added_keys),
                anyDuplicated(cur_key), anyDuplicated(ref_key))
    )
    if (!keys_ok) return(invisible(FALSE))

    common_keys <- intersect(cur_key, ref_key)
    cur <- cur[match(common_keys, key_string(cur, key))]
    ref <- ref[match(common_keys, key_string(ref, key))]
    setorderv(cur, key)
    setorderv(ref, key)
    if (is.null(fields)) fields <- setdiff(union(names(cur), names(ref)), key)
    for (field in fields) {
      present <- field %in% names(cur) && field %in% names(ref)
      if (!present) {
        field_comparisons[[length(field_comparisons) + 1L]] <<- data.table(
          scope = scope, field = field, n_compared = 0L, n_mismatch = NA_integer_,
          max_abs_diff = NA_real_, exact_equal = FALSE, within_tolerance = FALSE,
          details = "Field missing from current or reference output."
        )
        next
      }

      a <- cur[[field]]
      b <- ref[[field]]
      exact_equal <- identical(a, b)
      if (is.numeric(a) && is.numeric(b)) {
        both_na <- is.na(a) & is.na(b)
        same_inf <- is.infinite(a) & is.infinite(b) & sign(a) == sign(b)
        finite_pair <- is.finite(a) & is.finite(b)
        delta <- rep(Inf, length(a))
        delta[finite_pair] <- abs(a[finite_pair] - b[finite_pair])
        equal <- both_na | same_inf |
          (finite_pair & delta <= abs_tol + rel_tol * pmax(abs(b), 1))
        max_abs_diff <- if (any(finite_pair)) max(delta[finite_pair]) else NA_real_
      } else {
        both_na <- is.na(a) & is.na(b)
        equal <- both_na | (!is.na(a) & !is.na(b) & as.character(a) == as.character(b))
        max_abs_diff <- NA_real_
      }
      field_comparisons[[length(field_comparisons) + 1L]] <<- data.table(
        scope = scope, field = field, n_compared = length(a),
        n_mismatch = sum(!equal), max_abs_diff = max_abs_diff,
        exact_equal = exact_equal, within_tolerance = all(equal),
        details = if (exact_equal) "Exact vector agreement." else if (all(equal))
          "Values agree; storage type/attributes or sub-tolerance numeric representation differs."
        else "One or more values differ beyond the declared tolerance."
      )
    }
    invisible(TRUE)
  }

  total_items <- sum(class_monthly$items)
  add_stage3_check("total_items", "accounting", total_items == 4641000754,
                   "4641000754", sprintf("%.0f", total_items))
  add_stage3_check("observed_classes", "accounting", nrow(elig_class) == 344L,
                   "344", nrow(elig_class))
  add_stage3_check("eligible_classes", "accounting", sum(elig_class$eligible) == 220L,
                   "220", sum(elig_class$eligible))
  add_stage3_check("excluded_classes", "accounting", sum(!elig_class$eligible) == 124L,
                   "124", sum(!elig_class$eligible))
  add_stage3_check("observed_drugs", "accounting", nrow(elig_drug) == 2155L,
                   "2155", nrow(elig_drug))
  add_stage3_check("eligible_drugs", "accounting", sum(elig_drug$eligible) == 974L,
                   "974", sum(elig_drug$eligible))
  add_stage3_check("excluded_drugs", "accounting", sum(!elig_drug$eligible) == 1181L,
                   "1181", sum(!elig_drug$eligible))

  class_share <- coverage$excluded_item_share[coverage$level == "class"]
  drug_share <- coverage$excluded_item_share[coverage$level == "drug"]
  add_stage3_check("excluded_class_item_share", "accounting",
                   abs(class_share - 0.00585129338248757) < 1e-14,
                   "0.00585129338248757", format(class_share, digits = 17))
  add_stage3_check("excluded_drug_item_share", "accounting",
                   abs(drug_share - 0.0076596621901777) < 1e-14,
                   "0.0076596621901777", format(drug_share, digits = 17))
  sole_month_class <- sum(!elig_class$rule_every_month & elig_class$rule_min_volume)
  sole_month_drug <- sum(!elig_drug$rule_every_month & elig_drug$rule_min_volume)
  add_stage3_check("sole_every_month_class_exclusion", "accounting",
                   sole_month_class == 0L, "0", sole_month_class)
  add_stage3_check("sole_every_month_drug_exclusion", "accounting",
                   sole_month_drug == 1L, "1", sole_month_drug)

  add_stage3_check("class_models_complete", "models",
                   nrow(screen_class) == 220L && all(screen_class$converged) &&
                     all(is.na(screen_class$note) | !nzchar(screen_class$note)),
                   "220 fitted; 0 failures",
                   sprintf("%d fitted; %d failures", nrow(screen_class), sum(!screen_class$converged)))
  add_stage3_check("drug_models_complete", "models",
                   nrow(screen_drug) == 974L && all(screen_drug$converged) &&
                     all(is.na(screen_drug$note) | !nzchar(screen_drug$note)),
                   "974 fitted; 0 failures",
                   sprintf("%d fitted; %d failures", nrow(screen_drug), sum(!screen_drug$converged)))
  add_stage3_check("class_discoveries", "models", sum(screen_class$significant) == 125L,
                   "125", sum(screen_class$significant))
  add_stage3_check("meaningful_classes", "models", sum(results_class$meaningful) == 30L,
                   "30", sum(results_class$meaningful))

  route_string <- function(x) {
    x |> count(route, name = "n") |> arrange(route) |>
      transmute(value = paste0(route, "=", n)) |> pull(value) |> paste(collapse = ";")
  }
  class_routes_observed <- route_string(screen_class)
  drug_routes_observed <- route_string(screen_drug)
  add_stage3_check("class_model_routes", "models",
                   identical(class_routes_observed,
                             "HAC-Wald=165;NB-LRT=46;Poisson-LRT=9"),
                   "HAC-Wald=165;NB-LRT=46;Poisson-LRT=9", class_routes_observed)
  add_stage3_check("drug_model_routes", "models",
                   identical(drug_routes_observed,
                             "HAC-Wald=672;HAC-Wald(NBfail)=2;NB-LRT=208;Poisson-LRT=92"),
                   "canonical frozen-environment routes: HAC-Wald=672; HAC-Wald(NBfail)=2; NB-LRT=208; Poisson-LRT=92",
                   drug_routes_observed,
                   "The sealed baseline had NB-LRT=206 and NB-failure HAC=4; the four series and scientific effect are checked below.")

  required_reference_files <- c(
    "appendix1_coverage_summary.csv",
    "eligibility_class.rds", "eligibility_drug.rds",
    "screen_class.csv", "screen_drug.csv",
    "characterisation_class.csv", "characterisation_drug.csv",
    "results_class.csv", "results_drug.csv",
    file.path("results", "tables", "table1_meaningful_classes.csv"),
    file.path("results", "tables", "table2_meaningful_drugs.csv"),
    file.path("results", "tables", "table_accounting.csv"),
    file.path("results", "tables", "appendixA1_all_classes.csv"),
    file.path("results", "tables", "appendixA2_all_drugs.csv"),
    file.path("results", "tables", "appendix_coverage_summary.csv"),
    file.path("results", "tables", "appendix_route_proportions.csv")
  )
  missing_reference <- required_reference_files[
    !file.exists(file.path(baseline_dir, required_reference_files))
  ]
  add_stage3_check("reference_outputs_present", "reference", !length(missing_reference),
                   sprintf("%d sealed files", length(required_reference_files)),
                   sprintf("%d present", length(required_reference_files) - length(missing_reference)),
                   paste(missing_reference, collapse = "; "))
  if (length(missing_reference)) {
    atomic_fwrite(rbindlist(stage3_checks), file.path(stage3_dir, "stage3_qc_summary.csv"))
    stop("Stage 3 reference output(s) missing: ", paste(missing_reference, collapse = ", "))
  }

  compare_frames("eligibility_class", elig_class,
                 readRDS(file.path(baseline_dir, "eligibility_class.rds")),
                 "bnf_class_code")
  compare_frames("eligibility_drug", elig_drug,
                 readRDS(file.path(baseline_dir, "eligibility_drug.rds")),
                 "bnf_drug_code")
  compare_frames("screen_class", screen_class,
                 fread(file.path(baseline_dir, "screen_class.csv")),
                 "bnf_class_code")
  compare_frames("screen_drug", screen_drug,
                 fread(file.path(baseline_dir, "screen_drug.csv"),
                       colClasses = list(character = "bnf_drug_code")),
                 "bnf_drug_code")
  compare_frames("characterisation_class", char_class,
                 fread(file.path(baseline_dir, "characterisation_class.csv")),
                 "bnf_class_code")
  compare_frames("characterisation_drug", char_drug,
                 fread(file.path(baseline_dir, "characterisation_drug.csv"),
                       colClasses = list(character = "bnf_drug_code")),
                 "bnf_drug_code")
  compare_frames("results_class", results_class,
                 fread(file.path(baseline_dir, "results_class.csv")),
                 "bnf_class_code")
  compare_frames("results_drug", results_drug,
                 fread(file.path(baseline_dir, "results_drug.csv"),
                       colClasses = list(character = "bnf_drug_code")),
                 "bnf_drug_code")

  compare_character_table <- function(scope, relative_path, key,
                                      allowed_added_keys = character(),
                                      allowed_removed_keys = character()) {
    compare_frames(
      scope,
      fread(file.path(out_dir, relative_path), colClasses = "character"),
      fread(file.path(baseline_dir, relative_path), colClasses = "character"),
      key,
      rel_tol = 0,
      abs_tol = 0,
      allowed_added_keys = allowed_added_keys,
      allowed_removed_keys = allowed_removed_keys
    )
  }
  compare_character_table("table1_classes", file.path("results", "tables", "table1_meaningful_classes.csv"),
                          "BNF class code")
  compare_character_table("table2_drugs", file.path("results", "tables", "table2_meaningful_drugs.csv"),
                          "BNF drug code", allowed_added_keys = "1307000Q0")
  compare_character_table("table_accounting", file.path("results", "tables", "table_accounting.csv"),
                          "Stage")
  compare_character_table("appendix_classes", file.path("results", "tables", "appendixA1_all_classes.csv"),
                          "BNF class code")
  compare_character_table("appendix_drugs", file.path("results", "tables", "appendixA2_all_drugs.csv"),
                          "BNF drug code")
  compare_character_table("appendix_routes", file.path("results", "tables", "appendix_route_proportions.csv"),
                          c("level", "distribution", "route"))

  # The sealed early coverage CSV contains the known calculation bug. Require
  # the repaired early output to agree with the later, correct sealed table.
  compare_frames(
    "coverage_repair",
    fread(file.path(out_dir, "appendix1_coverage_summary.csv")),
    fread(file.path(baseline_dir, "results", "tables", "appendix_coverage_summary.csv")),
    "level"
  )
  old_coverage <- fread(file.path(baseline_dir, "appendix1_coverage_summary.csv"))
  bug_repaired <- all(old_coverage$excluded == 0L) &&
    identical(as.integer(coverage$excluded), c(124L, 1181L))
  add_stage3_check(
    "known_coverage_bug_repaired", "mechanical_repair", bug_repaired,
    "sealed early file incorrectly has 0/0 exclusions; canonical file has 124/1181",
    paste0("sealed=", paste(old_coverage$excluded, collapse = "/"),
           "; canonical=", paste(coverage$excluded, collapse = "/")),
    "This is the pre-declared mechanical correction and does not alter cohort membership."
  )

  # Scientific stability checks distinguish harmless representation/runtime
  # differences from a change in the study conclusions. The historical
  # environment was not recorded; the canonical R/renv environment is fixed.
  ref_screen_class <- normalise_comparison_codes(
    fread(file.path(baseline_dir, "screen_class.csv"))
  )
  new_screen_class <- normalise_comparison_codes(as.data.table(copy(screen_class)))
  class_compare <- merge(
    new_screen_class, ref_screen_class,
    by = "bnf_class_code", suffixes = c("_new", "_old")
  )
  class_discovery_same <- setequal(
    class_compare[significant_new == TRUE, bnf_class_code],
    class_compare[significant_old == TRUE, bnf_class_code]
  )
  class_route_same <- all(class_compare$route_new == class_compare$route_old)
  class_p_max_delta <- max(abs(class_compare$p_value_new - class_compare$p_value_old),
                           na.rm = TRUE)
  add_stage3_check(
    "class_inference_stable", "scientific_stability",
    class_discovery_same && class_route_same && class_p_max_delta < 2e-8,
    "identical discovery identities/routes; maximum raw-p difference <2e-8",
    sprintf("discoveries identical=%s; routes identical=%s; max raw-p difference=%.3g",
            class_discovery_same, class_route_same, class_p_max_delta),
    "Small numerical differences are attributable to the unrecorded historical runtime."
  )

  ref_screen_drug <- normalise_comparison_codes(fread(
    file.path(baseline_dir, "screen_drug.csv"),
    colClasses = list(character = "bnf_drug_code")
  ))
  new_screen_drug <- normalise_comparison_codes(as.data.table(copy(screen_drug)))
  drug_compare <- merge(
    new_screen_drug, ref_screen_drug,
    by = "bnf_drug_code", suffixes = c("_new", "_old")
  )
  primary_drug_same <- setequal(
    drug_compare[significant_new == TRUE, bnf_drug_code],
    drug_compare[significant_old == TRUE, bnf_drug_code]
  )
  all_drug_same <- setequal(
    drug_compare[sig_all_new == TRUE, bnf_drug_code],
    drug_compare[sig_all_old == TRUE, bnf_drug_code]
  )
  route_change_codes <- sort(drug_compare[route_new != route_old, bnf_drug_code])
  expected_route_change_codes <- sort(c("0409010AA", "0704010M0", "0704020P0", "0802020G0"))
  route_changes_non_significant <- drug_compare[
    bnf_drug_code %in% route_change_codes,
    all(!significant_new & !significant_old & !sig_all_new & !sig_all_old)
  ]
  unchanged_route_p_delta <- drug_compare[route_new == route_old,
    max(abs(p_value_new - p_value_old), na.rm = TRUE)]
  add_stage3_check(
    "drug_discovery_sets_stable", "scientific_stability",
    primary_drug_same && all_drug_same,
    "identical conditional and all-drug discovery identities",
    sprintf("conditional identical=%s; all-drug identical=%s",
            primary_drug_same, all_drug_same)
  )
  add_stage3_check(
    "drug_route_differences_explained", "scientific_stability",
    setequal(route_change_codes, expected_route_change_codes) &&
      route_changes_non_significant && unchanged_route_p_delta < 2e-8,
    paste0("four known route changes, all non-significant; unchanged-route max raw-p ",
           "difference <2e-8"),
    sprintf("codes=%s; all non-significant=%s; unchanged-route max raw-p difference=%.3g",
            paste(route_change_codes, collapse = ";"), route_changes_non_significant,
            unchanged_route_p_delta),
    "Three historical NB failures now converge and one current NB fit falls back to HAC; no inferential flag changes."
  )

  ref_results_class <- normalise_comparison_codes(
    fread(file.path(baseline_dir, "results_class.csv"))
  )
  new_results_class <- normalise_comparison_codes(as.data.table(copy(results_class)))
  class_result_compare <- merge(
    new_results_class, ref_results_class,
    by = "bnf_class_code", suffixes = c("_new", "_old")
  )
  class_meaningful_same <- setequal(
    class_result_compare[meaningful_new == TRUE, bnf_class_code],
    class_result_compare[meaningful_old == TRUE, bnf_class_code]
  )
  class_point_same <-
    max(abs(class_result_compare$peak_trough_ratio_new -
              class_result_compare$peak_trough_ratio_old), na.rm = TRUE) < 1e-9 &&
    all(class_result_compare$peak_month_new == class_result_compare$peak_month_old) &&
    all(class_result_compare$trough_month_new == class_result_compare$trough_month_old) &&
    all(class_result_compare$modality_new == class_result_compare$modality_old)
  add_stage3_check(
    "class_characterisation_stable", "scientific_stability",
    class_meaningful_same && class_point_same,
    "identical meaningful set, point amplitudes, timing and modality",
    sprintf("meaningful identical=%s; characterisation identical=%s",
            class_meaningful_same, class_point_same)
  )

  ref_results_drug <- fread(
    file.path(baseline_dir, "results_drug.csv"),
    colClasses = list(character = "bnf_drug_code")
  )
  new_results_drug <- as.data.table(copy(results_drug))
  drug_result_compare <- merge(
    new_results_drug, ref_results_drug,
    by = "bnf_drug_code", suffixes = c("_new", "_old")
  )
  meaningful_added <- sort(drug_result_compare[
    meaningful_new == TRUE & meaningful_old == FALSE, bnf_drug_code])
  meaningful_removed <- sort(drug_result_compare[
    meaningful_new == FALSE & meaningful_old == TRUE, bnf_drug_code])
  silver <- drug_result_compare[bnf_drug_code == "1307000Q0"]
  silver_boundary_change <- nrow(silver) == 1L &&
    silver$ptr_lci_old < amplitude_lci_threshold &&
    silver$ptr_lci_new >= amplitude_lci_threshold &&
    silver$stl_seasonal_strength_new >= stl_strength_threshold
  drug_point_same <-
    max(abs(drug_result_compare$peak_trough_ratio_new -
              drug_result_compare$peak_trough_ratio_old), na.rm = TRUE) < 1e-9 &&
    all(drug_result_compare$peak_month_new == drug_result_compare$peak_month_old) &&
    all(drug_result_compare$trough_month_new == drug_result_compare$trough_month_old) &&
    all(drug_result_compare$modality_new == drug_result_compare$modality_old)
  add_stage3_check(
    "drug_characterisation_stable", "scientific_stability",
    drug_point_same && identical(meaningful_added, "1307000Q0") &&
      !length(meaningful_removed) && silver_boundary_change,
    "identical point characterisation; only Silver nitrate crosses the 1.10 lower-CI boundary",
    sprintf(paste0("point characterisation identical=%s; added=%s; removed=%s; ",
                   "Silver nitrate lower CI %.3f -> %.3f, strength %.3f"),
            drug_point_same, paste(meaningful_added, collapse = ";"),
            paste(meaningful_removed, collapse = ";"), silver$ptr_lci_old,
            silver$ptr_lci_new, silver$stl_seasonal_strength_new),
    "Drug findings are exploratory; retain the canonical result and emphasise continuous estimates and threshold sensitivity."
  )

  ref_char_class <- normalise_comparison_codes(
    fread(file.path(baseline_dir, "characterisation_class.csv"))
  )
  new_char_class <- normalise_comparison_codes(as.data.table(copy(char_class)))
  class_ci_compare <- merge(
    new_char_class[, .(bnf_class_code, ptr_lci_new = ptr_lci, ptr_uci_new = ptr_uci)],
    ref_char_class[, .(bnf_class_code, ptr_lci_old = ptr_lci, ptr_uci_old = ptr_uci)],
    by = "bnf_class_code"
  )
  class_ci_max_relative <- max(
    abs(class_ci_compare$ptr_lci_new / class_ci_compare$ptr_lci_old - 1),
    abs(class_ci_compare$ptr_uci_new / class_ci_compare$ptr_uci_old - 1),
    na.rm = TRUE
  )
  ref_char_drug <- fread(
    file.path(baseline_dir, "characterisation_drug.csv"),
    colClasses = list(character = "bnf_drug_code")
  )
  drug_ci_compare <- merge(
    as.data.table(char_drug)[, .(bnf_drug_code, ptr_lci_new = ptr_lci, ptr_uci_new = ptr_uci)],
    ref_char_drug[, .(bnf_drug_code, ptr_lci_old = ptr_lci, ptr_uci_old = ptr_uci)],
    by = "bnf_drug_code"
  )
  drug_ci_max_relative <- max(
    abs(drug_ci_compare$ptr_lci_new / drug_ci_compare$ptr_lci_old - 1),
    abs(drug_ci_compare$ptr_uci_new / drug_ci_compare$ptr_uci_old - 1),
    na.rm = TRUE
  )
  add_stage3_check(
    "bootstrap_ci_runtime_difference_bounded", "scientific_stability",
    class_ci_max_relative < 0.01 && drug_ci_max_relative < 0.06,
    "maximum relative CI change <1% for classes and <6% for drugs",
    sprintf("class maximum=%.3f%%; drug maximum=%.3f%%",
            100 * class_ci_max_relative, 100 * drug_ci_max_relative),
    "The historical package/RNG environment was not recorded; canonical intervals are deterministic under the frozen environment."
  )

  ref_elig_drug <- normalise_comparison_codes(as.data.table(
    readRDS(file.path(baseline_dir, "eligibility_drug.rds"))
  ))
  new_elig_drug <- normalise_comparison_codes(as.data.table(copy(elig_drug)))
  name_compare <- merge(
    new_elig_drug[, .(bnf_drug_code, name_new = bnf_drug_name, eligible_new = eligible)],
    ref_elig_drug[, .(bnf_drug_code, name_old = bnf_drug_name, eligible_old = eligible)],
    by = "bnf_drug_code"
  )
  name_difference_codes <- sort(name_compare[name_new != name_old, bnf_drug_code])
  expected_name_difference_codes <- sort(c("0914011B0", "0914011E0", "0914011P0", "0914081A0"))
  whitespace_only <- name_compare[bnf_drug_code %in% name_difference_codes,
    all(trimws(name_new) == trimws(name_old) & !eligible_new & !eligible_old)]
  add_stage3_check(
    "excluded_name_whitespace_cleaned", "mechanical_repair",
    setequal(name_difference_codes, expected_name_difference_codes) && whitespace_only,
    "four excluded drug labels differ only by removed trailing whitespace",
    sprintf("codes=%s; whitespace-only and excluded=%s",
            paste(name_difference_codes, collapse = ";"), whitespace_only)
  )

  ref_elig_class <- normalise_comparison_codes(as.data.table(
    readRDS(file.path(baseline_dir, "eligibility_class.rds"))
  ))
  new_elig_class <- normalise_comparison_codes(as.data.table(copy(elig_class)))
  lookup_class_names <- merge(
    new_elig_class[, .(bnf_class_code, name_new = bnf_class_name, eligible_new = eligible)],
    ref_elig_class[, .(bnf_class_code, name_old = bnf_class_name, eligible_old = eligible)],
    by = "bnf_class_code"
  )
  lookup_class_codes <- sort(lookup_class_names[
    is.na(name_old) & !is.na(name_new) & name_new == "", bnf_class_code])
  lookup_drug_codes <- sort(name_compare[
    bnf_drug_code == "0914081A0" & is.na(ref_elig_drug$bnf_class_name[
      match(bnf_drug_code, ref_elig_drug$bnf_drug_code)]) &
      new_elig_drug$bnf_class_name[match(bnf_drug_code, new_elig_drug$bnf_drug_code)] == "",
    bnf_drug_code
  ])
  add_stage3_check(
    "unmatched_lookup_label_representation", "mechanical_repair",
    setequal(lookup_class_codes, c("091304", "091408")) &&
      identical(lookup_drug_codes, "0914081A0") &&
      lookup_class_names[bnf_class_code %in% lookup_class_codes,
                         all(!eligible_new & !eligible_old)],
    "historical NA versus canonical blank label only for retained May-2025 lookup gaps 091304/091408",
    sprintf("class codes=%s; drug code=%s",
            paste(lookup_class_codes, collapse = ";"), paste(lookup_drug_codes, collapse = ";")),
    "Both class series and the affected drug are below eligibility thresholds; item accounting is unchanged."
  )

  explained_differences <- data.table(
    difference_id = c(
      "coverage_summary", "identifier_storage", "excluded_name_whitespace",
      "bnf_lookup_label_representation", "class_runtime_numerics",
      "drug_route_selection", "bootstrap_intervals"
    ),
    cause = c(
      "Pre-declared sequential-summary bug in the historical early coverage output.",
      "Historical CSV/RDS reads lost leading zeroes on numeric-looking six-digit class codes.",
      "Canonical import trims four trailing spaces in excluded nutrition-product labels.",
      "The two class codes absent from the May 2025 BNF snapshot are blank canonically and NA historically.",
      "Historical R/package environment was not recorded; canonical environment is R 4.6.1 with renv.",
      "Negative-binomial convergence differs for four non-significant drug series under the frozen environment.",
      "Parametric-bootstrap draws/covariances differ from the unrecorded historical environment."
    ),
    extent = c(
      "Historical early output 0/0 exclusions; canonical 124/1181 and correct volume shares.",
      "Representation only; normalized class-code membership agrees.",
      paste(name_difference_codes, collapse = ";"),
      paste(lookup_class_codes, collapse = ";"),
      sprintf("Maximum class raw-p difference %.3g; routes and discoveries identical.", class_p_max_delta),
      paste(route_change_codes, collapse = ";"),
      sprintf("Maximum relative CI change %.3f%% class and %.3f%% drug; Silver nitrate lower CI %.3f to %.3f.",
              100 * class_ci_max_relative, 100 * drug_ci_max_relative,
              silver$ptr_lci_old, silver$ptr_lci_new)
    ),
    scientific_effect = c(
      "None: cohort membership and item totals are unchanged.",
      "None.", "None: all four series are excluded.",
      "None: retained items reconcile and the affected series are not eligible.",
      "None: 125 class discoveries and 30 meaningful classes are unchanged.",
      "None: no conditional or all-drug discovery flag changes.",
      "One exploratory boundary classification: Silver nitrate; class findings are unchanged."
    ),
    disposition = c(
      "Corrected in canonical code and checked against the correct sealed publication table.",
      "Normalize only for baseline comparison; retain character codes in canonical outputs.",
      "Retain cleaned canonical labels.",
      "Retain the documented blank canonical labels; do not impute names from a different reference date.",
      "Retain canonical values from the frozen current environment.",
      "Retain diagnostic routing produced by the frozen current environment.",
      "Retain canonical intervals; report continuous estimates and threshold sensitivity."
    )
  )
  atomic_fwrite(explained_differences, file.path(stage3_dir, "explained_baseline_differences.csv"))

  field_comparison <- rbindlist(field_comparisons, use.names = TRUE, fill = TRUE)
  allowed_mismatch_fields <- list(
    eligibility_class = c("bnf_class_name"),
    eligibility_drug = c("bnf_class_name", "bnf_drug_name"),
    screen_class = c("p_value", "disp_ratio", "disp_p", "lb_p", "nw_lag", "theta",
                     "b_sin12", "b_cos12", "b_sin6", "b_cos6", "p_adj"),
    screen_drug = c("p_value", "distribution", "route", "disp_ratio", "disp_p", "lb_p",
                    "nw_lag", "hac_capped", "theta", "b_sin12", "b_cos12", "b_sin6",
                    "b_cos6", "p_adj", "p_adj_all", "class_p_adj"),
    characterisation_class = c("ptr_lci", "ptr_uci", "aic_2h", "p_value", "p_adj"),
    characterisation_drug = c("ptr_lci", "ptr_uci", "theta", "p_value", "p_adj",
                              "p_adj_all", "class_p_adj"),
    results_class = c("ptr_lci", "ptr_uci", "p_adj"),
    results_drug = c("ptr_lci", "ptr_uci", "p_adj_primary", "p_adj_all", "meaningful"),
    table1_classes = c("Peak-to-trough ratio (95% CI)"),
    table2_drugs = c("Peak-to-trough ratio (95% CI)"),
    table_accounting = c("Drugs"),
    appendix_classes = c("95% CI lower", "95% CI upper"),
    appendix_drugs = c("95% CI lower", "95% CI upper", "Meaningful seasonality"),
    appendix_routes = c("n", "pct"),
    coverage_repair = character()
  )
  for (scope_name in unique(field_comparison$scope)) {
    z <- field_comparison[scope == scope_name]
    mismatch_fields <- z[within_tolerance == FALSE, field]
    allowed_fields <- allowed_mismatch_fields[[scope_name]]
    if (is.null(allowed_fields)) allowed_fields <- character()
    unexpected_fields <- setdiff(mismatch_fields, allowed_fields)
    add_stage3_check(
      paste0(scope_name, "_values"), scope_name,
      !length(unexpected_fields),
      sprintf("no unexplained field difference across %d fields", nrow(z)),
      sprintf("different fields=%s; unexpected fields=%s",
              paste(mismatch_fields, collapse = ";"),
              paste(unexpected_fields, collapse = ";")),
      sprintf("%d fields are byte/type exact; allowed differences are evaluated by explicit scientific-stability checks",
              sum(z$exact_equal))
    )
  }

  stage3_qc_summary <- rbindlist(stage3_checks, use.names = TRUE, fill = TRUE)
  setorder(stage3_qc_summary, scope, check_id)
  atomic_fwrite(field_comparison, file.path(stage3_dir, "baseline_field_comparison.csv"))
  atomic_fwrite(stage3_qc_summary, file.path(stage3_dir, "stage3_qc_summary.csv"))

  reference_manifest <- data.table(
    relative_path = required_reference_files,
    bytes = as.numeric(file.info(file.path(baseline_dir, required_reference_files))$size),
    sha256 = vapply(file.path(baseline_dir, required_reference_files), sha256_file, character(1))
  )
  atomic_fwrite(reference_manifest, file.path(stage3_dir, "sealed_reference_manifest.csv"))

  if (any(!stage3_qc_summary$pass)) {
    failed <- stage3_qc_summary[pass == FALSE, check_id]
    stop("Stage 3 completion gate failed: ", paste(failed, collapse = ", "),
         ". See ", file.path(stage3_dir, "stage3_qc_summary.csv"), ".")
  }

  snapshot_files <- c(
    "appendix1_coverage_summary.csv", "appendix1_exclusions_class.csv",
    "appendix1_exclusions_drug.csv", "screen_class.csv", "screen_drug.csv",
    "screen_class_route_summary.csv", "screen_drug_route_summary.csv",
    "model_failures.csv", "characterisation_class.csv", "characterisation_drug.csv",
    "results_class.csv", "results_drug.csv", "appendix_route_proportions.csv",
    file.path("results", "tables", "table1_meaningful_classes.csv"),
    file.path("results", "tables", "table2_meaningful_drugs.csv"),
    file.path("results", "tables", "table_accounting.csv"),
    file.path("results", "tables", "appendixA1_all_classes.csv"),
    file.path("results", "tables", "appendixA2_all_drugs.csv"),
    file.path("results", "tables", "appendix_coverage_summary.csv"),
    file.path("results", "tables", "appendix_route_proportions.csv")
  )
  for (relative_path in snapshot_files) {
    source_path <- file.path(out_dir, relative_path)
    target_path <- file.path(snapshot_dir, relative_path)
    dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
    temporary_path <- paste0(target_path, ".tmp")
    if (!file.copy(source_path, temporary_path, overwrite = TRUE) ||
        !file.rename(temporary_path, target_path)) {
      stop("Could not seal Stage 3 snapshot file: ", relative_path)
    }
  }
  snapshot_manifest <- data.table(
    relative_path = snapshot_files,
    bytes = as.numeric(file.info(file.path(snapshot_dir, snapshot_files))$size),
    sha256 = vapply(file.path(snapshot_dir, snapshot_files), sha256_file, character(1))
  )
  atomic_fwrite(snapshot_manifest, file.path(stage3_dir, "stage3_snapshot_manifest.csv"))

  stage3_completion <- data.table(
    stage = "stage3",
    status = "PASS",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    checks_passed = sum(stage3_qc_summary$pass),
    checks_total = nrow(stage3_qc_summary),
    compared_fields = nrow(field_comparison),
    fields_beyond_tolerance = sum(!field_comparison$within_tolerance),
    exact_field_vectors = sum(field_comparison$exact_equal),
    snapshot_files = nrow(snapshot_manifest),
    analysis_script_sha256 = sha256_file(script_path),
    renv_lock_sha256 = sha256_file(lock_path)
  )
  atomic_fwrite(stage3_completion, file.path(stage3_dir, "stage3_completion.csv"))
  message(
    "Stage 3 complete: ", nrow(stage3_qc_summary),
    " checks passed; ", nrow(field_comparison),
    " baseline fields compared; stopping before all-series plotting and 2021 diagnostics."
  )
  return(invisible(TRUE))
}

## ---- STAGE 4: authoritative all-eligible-drug BH family --------------------
if (run_stage == "stage4") {
  stage4_dir <- file.path(out_dir, "qc", "stage4")
  snapshot_dir <- file.path(stage4_dir, "stage4_snapshot")
  stage3_dir <- file.path(out_dir, "qc", "stage3")
  stage3_snapshot_dir <- file.path(stage3_dir, "baseline_snapshot")
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  stage4_checks <- list()
  add_stage4_check <- function(check_id, pass, expected, observed, details = "") {
    stage4_checks[[length(stage4_checks) + 1L]] <<- data.table(
      check_id = check_id, pass = isTRUE(pass), expected = as.character(expected),
      observed = as.character(observed), details = as.character(details)
    )
  }
  pad_class_code <- function(x) {
    z <- trimws(as.character(x))
    numeric_code <- grepl("^[0-9]{1,6}$", z)
    z[numeric_code] <- sprintf("%06d", as.integer(z[numeric_code]))
    z
  }
  flag_codes <- function(x, flag, code) sort(as.character(x[[code]][x[[flag]] %in% TRUE]))

  required_stage3_files <- c(
    "screen_class.csv", "screen_drug.csv", "characterisation_class.csv",
    "characterisation_drug.csv", "results_class.csv", "results_drug.csv"
  )
  missing_stage3 <- required_stage3_files[
    !file.exists(file.path(stage3_snapshot_dir, required_stage3_files))
  ]
  add_stage4_check(
    "stage3_reference_present", !length(missing_stage3),
    paste(length(required_stage3_files), "Stage 3 baseline files"),
    paste(length(required_stage3_files) - length(missing_stage3), "present"),
    paste(missing_stage3, collapse = ";")
  )
  if (length(missing_stage3)) stop("Stage 4 requires the completed Stage 3 snapshot.")

  stage3_manifest_path <- file.path(stage3_dir, "stage3_snapshot_manifest.csv")
  stage3_manifest <- fread(stage3_manifest_path)
  stage3_manifest_paths <- file.path(stage3_snapshot_dir, stage3_manifest$relative_path)
  stage3_hash_now <- vapply(stage3_manifest_paths, sha256_file, character(1))
  stage3_snapshot_intact <- all(file.exists(stage3_manifest_paths)) &&
    all(as.numeric(file.info(stage3_manifest_paths)$size) == stage3_manifest$bytes) &&
    identical(unname(stage3_hash_now), stage3_manifest$sha256)
  add_stage4_check(
    "stage3_snapshot_unchanged", stage3_snapshot_intact,
    paste(nrow(stage3_manifest), "unchanged snapshot files"),
    paste(sum(stage3_hash_now == stage3_manifest$sha256), "matching hashes")
  )

  add_stage4_check("class_family_size", nrow(screen_class) == 220L, "220", nrow(screen_class))
  add_stage4_check("class_discoveries", sum(screen_class$class_significant) == 125L,
                   "125", sum(screen_class$class_significant))
  add_stage4_check("drug_family_size", nrow(screen_drug) == 974L, "974", nrow(screen_drug))
  add_stage4_check("drug_discoveries", sum(screen_drug$drug_significant) == 391L,
                   "391", sum(screen_drug$drug_significant))
  add_stage4_check(
    "drug_discoveries_by_parent",
    sum(screen_drug$drug_significant & screen_drug$parent_class_significant) == 337L &&
      sum(screen_drug$drug_significant & !screen_drug$parent_class_significant) == 54L,
    "337 in significant parent classes; 54 outside",
    sprintf("%d in; %d outside",
            sum(screen_drug$drug_significant & screen_drug$parent_class_significant),
            sum(screen_drug$drug_significant & !screen_drug$parent_class_significant))
  )
  add_stage4_check(
    "legacy_family_audit_only",
    sum(screen_drug$parent_class_significant) == 658L &&
      sum(screen_drug$conditional_drug_significant_legacy) == 350L,
    "658 selected-parent drugs; 350 legacy flags",
    sprintf("%d selected-parent drugs; %d legacy flags",
            sum(screen_drug$parent_class_significant),
            sum(screen_drug$conditional_drug_significant_legacy)),
    "Legacy fields are explicitly named and excluded from characterisation and publication counts."
  )
  add_stage4_check(
    "models_complete",
    all(screen_class$converged) && all(screen_drug$converged) && nrow(model_failures) == 0L,
    "1194 completed models; 0 failures",
    sprintf("%d completed; %d failures",
            sum(screen_class$converged) + sum(screen_drug$converged), nrow(model_failures))
  )
  add_stage4_check(
    "characterisation_scope",
    nrow(char_class) == 125L && nrow(char_drug) == 391L &&
      setequal(char_drug$bnf_drug_code,
               screen_drug$bnf_drug_code[screen_drug$drug_significant]),
    "125 class discoveries and exactly 391 all-drug discoveries characterised",
    sprintf("%d classes; %d drugs", nrow(char_class), nrow(char_drug))
  )
  add_stage4_check(
    "meaningful_counts",
    sum(results_class$meaningful) == 30L && sum(results_drug$meaningful) == 88L &&
      sum(results_drug$meaningful & results_drug$parent_class_significant) == 87L &&
      sum(results_drug$meaningful & !results_drug$parent_class_significant) == 1L,
    "30 classes; 88 exploratory drugs (87 in significant parents, 1 outside)",
    sprintf("%d classes; %d drugs (%d in, %d outside)",
            sum(results_class$meaningful), sum(results_drug$meaningful),
            sum(results_drug$meaningful & results_drug$parent_class_significant),
            sum(results_drug$meaningful & !results_drug$parent_class_significant))
  )

  ref_class <- fread(file.path(stage3_snapshot_dir, "screen_class.csv"))
  ref_class[, bnf_class_code := pad_class_code(bnf_class_code)]
  new_class <- as.data.table(copy(screen_class))
  new_class[, bnf_class_code := pad_class_code(bnf_class_code)]
  class_compare <- merge(
    new_class[, .(bnf_class_code, p_new = p_value, q_new = class_q_bh,
                  sig_new = class_significant, route_new = route)],
    ref_class[, .(bnf_class_code, p_old = p_value, q_old = p_adj,
                  sig_old = significant, route_old = route)],
    by = "bnf_class_code"
  )
  add_stage4_check(
    "class_inference_matches_stage3",
    nrow(class_compare) == 220L &&
      max(abs(class_compare$p_new - class_compare$p_old), na.rm = TRUE) < 1e-10 &&
      max(abs(class_compare$q_new - class_compare$q_old), na.rm = TRUE) < 1e-10 &&
      all(class_compare$sig_new == class_compare$sig_old) &&
      all(class_compare$route_new == class_compare$route_old),
    "same raw p, class BH q, flags and routes as Stage 3",
    sprintf("%d rows; max p delta %.3g; max q delta %.3g",
            nrow(class_compare), max(abs(class_compare$p_new - class_compare$p_old)),
            max(abs(class_compare$q_new - class_compare$q_old)))
  )

  ref_drug <- fread(
    file.path(stage3_snapshot_dir, "screen_drug.csv"),
    colClasses = list(character = "bnf_drug_code")
  )
  new_drug <- as.data.table(copy(screen_drug))
  drug_compare <- merge(
    new_drug[, .(bnf_drug_code, p_new = p_value, q_new = drug_all_q_bh,
                 sig_new = drug_significant, parent_new = parent_class_significant,
                 legacy_q_new = conditional_drug_q_legacy,
                 legacy_sig_new = conditional_drug_significant_legacy, route_new = route)],
    ref_drug[, .(bnf_drug_code, p_old = p_value, q_old = p_adj_all,
                 sig_old = sig_all, parent_old = parent_class_sig,
                 legacy_q_old = p_adj, legacy_sig_old = significant, route_old = route)],
    by = "bnf_drug_code"
  )
  legacy_q_equal <- with(drug_compare,
    all((is.na(legacy_q_new) & is.na(legacy_q_old)) |
          (!is.na(legacy_q_new) & !is.na(legacy_q_old) &
             abs(legacy_q_new - legacy_q_old) < 1e-10)))
  add_stage4_check(
    "all_drug_inference_matches_stage3",
    nrow(drug_compare) == 974L &&
      max(abs(drug_compare$p_new - drug_compare$p_old), na.rm = TRUE) < 1e-10 &&
      max(abs(drug_compare$q_new - drug_compare$q_old), na.rm = TRUE) < 1e-10 &&
      all(drug_compare$sig_new == drug_compare$sig_old) &&
      all(drug_compare$parent_new == drug_compare$parent_old) &&
      all(drug_compare$route_new == drug_compare$route_old) &&
      legacy_q_equal && all(drug_compare$legacy_sig_new == drug_compare$legacy_sig_old),
    "Stage 4 q/flags equal Stage 3 all-drug fields; legacy audit equals old selected-parent fields",
    sprintf("%d rows; max p delta %.3g; max all-drug q delta %.3g",
            nrow(drug_compare), max(abs(drug_compare$p_new - drug_compare$p_old)),
            max(abs(drug_compare$q_new - drug_compare$q_old)))
  )

  ref_results_drug <- fread(
    file.path(stage3_snapshot_dir, "results_drug.csv"),
    colClasses = list(character = "bnf_drug_code")
  )[sig_all == TRUE]
  stage4_drug_membership_ok <- setequal(results_drug$bnf_drug_code,
                                        ref_results_drug$bnf_drug_code)
  result_compare <- merge(
    as.data.table(results_drug)[, .(bnf_drug_code, peak_new = peak_trough_ratio,
                                    lo_new = ptr_lci, hi_new = ptr_uci,
                                    strength_new = stl_seasonal_strength,
                                    meaningful_new = meaningful)],
    ref_results_drug[, .(bnf_drug_code, peak_old = peak_trough_ratio,
                         lo_old = ptr_lci, hi_old = ptr_uci,
                         strength_old = stl_seasonal_strength,
                         meaningful_old = meaningful)],
    by = "bnf_drug_code"
  )
  max_characterisation_delta <- max(
    abs(result_compare$peak_new - result_compare$peak_old),
    abs(result_compare$lo_new - result_compare$lo_old),
    abs(result_compare$hi_new - result_compare$hi_old),
    abs(result_compare$strength_new - result_compare$strength_old),
    na.rm = TRUE
  )
  add_stage4_check(
    "all_drug_characterisation_matches_stage3",
    stage4_drug_membership_ok && nrow(result_compare) == 391L &&
      max_characterisation_delta < 1e-10 &&
      all(result_compare$meaningful_new == result_compare$meaningful_old),
    "Stage 3 all-drug subset reproduced exactly under the same frozen environment",
    sprintf("membership equal=%s; %d rows; max reported-value delta %.3g",
            stage4_drug_membership_ok, nrow(result_compare), max_characterisation_delta)
  )

  machine_output_files <- c(
    "screen_class.csv", "screen_drug.csv", "characterisation_class.csv",
    "characterisation_drug.csv", "results_class.csv", "results_drug.csv"
  )
  forbidden_names <- c(
    "p_adj", "significant", "p_adj_all", "sig_all", "parent_class_sig",
    "p_adj_primary", "significant_primary", "class_p_adj"
  )
  ambiguous_columns <- rbindlist(lapply(machine_output_files, function(relative_path) {
    nm <- names(fread(file.path(out_dir, relative_path), nrows = 0L))
    bad_names <- intersect(nm, forbidden_names)
    if (length(bad_names)) {
      data.table(file = relative_path, column = bad_names)
    } else {
      data.table(file = character(), column = character())
    }
  }), fill = TRUE)
  add_stage4_check(
    "no_ambiguous_machine_columns", nrow(ambiguous_columns) == 0L,
    "no unqualified/legacy ambiguous column names",
    if (nrow(ambiguous_columns)) paste(ambiguous_columns$file, ambiguous_columns$column,
                                       collapse = ";") else "none"
  )

  table2_names <- names(fread(file.path(tab_dir, "table2_meaningful_drugs.csv"), nrows = 0L))
  appendix_drug_names <- names(fread(file.path(tab_dir, "appendixA2_all_drugs.csv"), nrows = 0L))
  publication_scope_ok <-
    "Drug q (BH; all eligible drugs)" %in% table2_names &&
    "Drug q (BH; all eligible drugs)" %in% appendix_drug_names &&
    !any(grepl("within class|all-drug scan|^Adjusted p$", c(table2_names, appendix_drug_names),
               ignore.case = TRUE)) &&
    nrow(main_table_drugs) == 88L && nrow(appendix_table_drugs) == 391L &&
    accounting$Drugs[2] == 391L
  add_stage4_check(
    "publication_tables_use_all_drug_family", publication_scope_ok,
    "all drug tables use only the 974-drug BH family and explicit labels",
    sprintf("table2 rows=%d; appendix rows=%d; accounting discoveries=%d",
            nrow(main_table_drugs), nrow(appendix_table_drugs), accounting$Drugs[2])
  )

  metadata_ok <- nrow(inference_scope_metadata) == 3L &&
    inference_scope_metadata[level == "class", analysis_role] == "primary_inferential" &&
    inference_scope_metadata[level == "drug", analysis_role] == "secondary_exploratory" &&
    grepl("whole hierarchy", inference_scope_metadata[level == "hierarchy", interpretation])
  add_stage4_check(
    "inference_scope_metadata", metadata_ok,
    "primary class family, secondary exploratory drug family, no whole-hierarchy claim",
    paste(inference_scope_metadata$level, inference_scope_metadata$analysis_role,
          sep = "=", collapse = ";")
  )

  stage4_comparison <- data.table(
    metric = c(
      "class_family_size", "class_discoveries", "drug_family_size", "drug_discoveries",
      "drug_discoveries_parent_significant", "drug_discoveries_parent_not_significant",
      "legacy_conditional_discoveries", "characterised_drugs", "meaningful_classes",
      "meaningful_drugs"
    ),
    stage3_value = c(220, 125, 974, 391, 337, 54, 350, 404, 30, 88),
    stage4_value = c(
      nrow(screen_class), sum(screen_class$class_significant), nrow(screen_drug),
      sum(screen_drug$drug_significant),
      sum(screen_drug$drug_significant & screen_drug$parent_class_significant),
      sum(screen_drug$drug_significant & !screen_drug$parent_class_significant),
      sum(screen_drug$conditional_drug_significant_legacy), nrow(char_drug),
      sum(results_class$meaningful), sum(results_drug$meaningful)
    ),
    interpretation = c(
      rep("unchanged", 7),
      "13 legacy-conditional-only drugs are no longer characterised/reported",
      "unchanged", "unchanged among all-drug discoveries"
    )
  )
  atomic_fwrite(stage4_comparison, file.path(stage4_dir, "stage3_stage4_comparison.csv"))

  stage4_qc_summary <- rbindlist(stage4_checks, use.names = TRUE, fill = TRUE)
  atomic_fwrite(stage4_qc_summary, file.path(stage4_dir, "stage4_qc_summary.csv"))
  if (any(!stage4_qc_summary$pass)) {
    stop("Stage 4 completion gate failed: ",
         paste(stage4_qc_summary[pass == FALSE, check_id], collapse = ", "),
         ". See ", file.path(stage4_dir, "stage4_qc_summary.csv"), ".")
  }

  snapshot_files <- c(
    "screen_class.csv", "screen_drug.csv", "screen_class_route_summary.csv",
    "screen_drug_route_summary.csv", "model_failures.csv",
    "characterisation_class.csv", "characterisation_drug.csv",
    "results_class.csv", "results_drug.csv", "appendix_route_proportions.csv",
    "inference_scope_metadata.csv",
    file.path("results", "tables", "table1_meaningful_classes.csv"),
    file.path("results", "tables", "table2_meaningful_drugs.csv"),
    file.path("results", "tables", "table_accounting.csv"),
    file.path("results", "tables", "appendixA1_all_classes.csv"),
    file.path("results", "tables", "appendixA2_all_drugs.csv"),
    file.path("results", "tables", "appendix_route_proportions.csv")
  )
  for (relative_path in snapshot_files) {
    source_path <- file.path(out_dir, relative_path)
    target_path <- file.path(snapshot_dir, relative_path)
    dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
    temporary_path <- paste0(target_path, ".tmp")
    if (!file.copy(source_path, temporary_path, overwrite = TRUE) ||
        !file.rename(temporary_path, target_path)) {
      stop("Could not seal Stage 4 snapshot file: ", relative_path)
    }
  }
  stage4_manifest <- data.table(
    relative_path = snapshot_files,
    bytes = as.numeric(file.info(file.path(snapshot_dir, snapshot_files))$size),
    sha256 = vapply(file.path(snapshot_dir, snapshot_files), sha256_file, character(1))
  )
  atomic_fwrite(stage4_manifest, file.path(stage4_dir, "stage4_snapshot_manifest.csv"))
  stage4_completion <- data.table(
    stage = "stage4", status = "PASS",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    checks_passed = sum(stage4_qc_summary$pass), checks_total = nrow(stage4_qc_summary),
    class_discoveries = sum(screen_class$class_significant),
    drug_discoveries = sum(screen_drug$drug_significant),
    meaningful_classes = sum(results_class$meaningful),
    meaningful_drugs = sum(results_drug$meaningful),
    snapshot_files = nrow(stage4_manifest),
    analysis_script_sha256 = sha256_file(script_path),
    renv_lock_sha256 = sha256_file(lock_path)
  )
  atomic_fwrite(stage4_completion, file.path(stage4_dir, "stage4_completion.csv"))
  message(
    "Stage 4 complete: one all-eligible-drug BH family is authoritative; ",
    nrow(stage4_qc_summary), " checks passed; stopping before all-series plotting."
  )
  return(invisible(TRUE))
}

## ---- STAGE 5.1: secular-trend specification sensitivity -------------------
if (run_stage == "stage5_trend") {
  trend_dir <- file.path(out_dir, "qc", "stage5", "trend")
  trend_snapshot_dir <- file.path(trend_dir, "trend_snapshot")
  stage4_dir <- file.path(out_dir, "qc", "stage4")
  stage4_snapshot_dir <- file.path(stage4_dir, "stage4_snapshot")
  dir.create(trend_snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  trend_checks <- list()
  add_trend_check <- function(check_id, pass, expected, observed, details = "") {
    trend_checks[[length(trend_checks) + 1L]] <<- data.table(
      check_id = check_id, pass = isTRUE(pass), expected = as.character(expected),
      observed = as.character(observed), details = as.character(details)
    )
  }

  # The sensitivity starts from the sealed Stage 4 authority. It adds trend
  # bases only to a local copy so the Stage 2 covariate artefact stays unchanged.
  stage4_manifest_path <- file.path(stage4_dir, "stage4_snapshot_manifest.csv")
  if (!file.exists(stage4_manifest_path)) {
    stop("Stage 5.1 requires a completed Stage 4 snapshot.")
  }
  stage4_manifest <- fread(stage4_manifest_path)
  stage4_paths <- file.path(stage4_snapshot_dir, stage4_manifest$relative_path)
  stage4_hash_now <- vapply(stage4_paths, sha256_file, character(1))
  stage4_intact <- all(file.exists(stage4_paths)) &&
    all(as.numeric(file.info(stage4_paths)$size) == stage4_manifest$bytes) &&
    identical(unname(stage4_hash_now), stage4_manifest$sha256)
  add_trend_check(
    "stage4_authority_unchanged", stage4_intact,
    paste(nrow(stage4_manifest), "unchanged Stage 4 files"),
    paste(sum(stage4_hash_now == stage4_manifest$sha256), "matching hashes")
  )
  if (!stage4_intact) stop("The sealed Stage 4 authority failed verification.")

  covar_trend <- copy(covar)
  covar_trend[, trend_linear := t]
  trend4 <- splines::ns(covar_trend$t, df = 4)
  covar_trend[, `:=`(
    trend4_1 = trend4[, 1], trend4_2 = trend4[, 2],
    trend4_3 = trend4[, 3], trend4_4 = trend4[, 4]
  )]

  f_linear_full <- items ~ trend_linear + sin12 + cos12 + sin6 + cos6 + offset(off)
  f_linear_red  <- items ~ trend_linear + offset(off)
  f_spline4_full <- items ~ trend4_1 + trend4_2 + trend4_3 + trend4_4 +
    sin12 + cos12 + sin6 + cos6 + offset(off)
  f_spline4_red <- items ~ trend4_1 + trend4_2 + trend4_3 + trend4_4 + offset(off)

  standardise_primary_class <- screen_class |>
    transmute(
      bnf_class_code, bnf_class_name,
      analysis_specification = "trend_spline3",
      trend_specification = "ns_df3",
      multiplicity_family = "all_eligible_classes",
      p_value, seasonality_q_bh = class_q_bh,
      seasonality_detected = class_significant,
      distribution, route, disp_ratio, disp_p, lb_p, nw_lag, hac_capped,
      theta, b_sin12, b_cos12, b_sin6, b_cos6, converged, note
    )
  standardise_primary_drug <- screen_drug |>
    transmute(
      bnf_class_code, bnf_class_name, bnf_drug_code, bnf_drug_name,
      analysis_specification = "trend_spline3",
      trend_specification = "ns_df3",
      multiplicity_family = "all_eligible_drugs",
      p_value, seasonality_q_bh = drug_all_q_bh,
      seasonality_detected = drug_significant,
      parent_class_q_bh = class_q_bh,
      parent_class_detected = parent_class_significant,
      distribution, route, disp_ratio, disp_p, lb_p, nw_lag, hac_capped,
      theta, b_sin12, b_cos12, b_sin6, b_cos6, converged, note
    )

  fit_trend_screen <- function(monthly, id_cols, f_full, f_red,
                               analysis_specification, trend_specification,
                               multiplicity_family) {
    fitted <- monthly |>
      select(all_of(c(id_cols, "year_month", "items"))) |>
      group_by(across(all_of(id_cols))) |>
      group_modify(~ fit_test_series(
        .x, covar_trend, f_full = f_full, f_red = f_red,
        offset_col = "offset_log_patient_days"
      )) |>
      ungroup()
    fitted |>
      mutate(
        analysis_specification = analysis_specification,
        trend_specification = trend_specification,
        multiplicity_family = multiplicity_family,
        seasonality_q_bh = p.adjust(p_value, method = "BH"),
        seasonality_detected = !is.na(seasonality_q_bh) &
          seasonality_q_bh < fdr_alpha
      ) |>
      arrange(seasonality_q_bh)
  }

  class_linear <- fit_trend_screen(
    class_monthly_elig,
    c("bnf_class_code", "bnf_class_name"),
    f_linear_full, f_linear_red, "trend_linear", "linear",
    "all_eligible_classes"
  )
  class_spline4 <- fit_trend_screen(
    class_monthly_elig,
    c("bnf_class_code", "bnf_class_name"),
    f_spline4_full, f_spline4_red, "trend_spline4", "ns_df4",
    "all_eligible_classes"
  )
  trend_screen_class <- bind_rows(
    standardise_primary_class, class_linear, class_spline4
  ) |>
    arrange(analysis_specification, seasonality_q_bh, bnf_class_code)

  drug_linear <- fit_trend_screen(
    drug_monthly_elig,
    c("bnf_class_code", "bnf_drug_code", "bnf_drug_name"),
    f_linear_full, f_linear_red, "trend_linear", "linear",
    "all_eligible_drugs"
  ) |>
    left_join(
      class_linear |>
        select(bnf_class_code, bnf_class_name,
               parent_class_q_bh = seasonality_q_bh,
               parent_class_detected = seasonality_detected),
      by = "bnf_class_code"
    )
  drug_spline4 <- fit_trend_screen(
    drug_monthly_elig,
    c("bnf_class_code", "bnf_drug_code", "bnf_drug_name"),
    f_spline4_full, f_spline4_red, "trend_spline4", "ns_df4",
    "all_eligible_drugs"
  ) |>
    left_join(
      class_spline4 |>
        select(bnf_class_code, bnf_class_name,
               parent_class_q_bh = seasonality_q_bh,
               parent_class_detected = seasonality_detected),
      by = "bnf_class_code"
    )
  trend_screen_drug <- bind_rows(
    standardise_primary_drug, drug_linear, drug_spline4
  ) |>
    arrange(analysis_specification, seasonality_q_bh, bnf_drug_code)

  characterise_trend_level <- function(scr, monthly, id_col, f_full) {
    detected <- scr |> filter(seasonality_detected)
    if (!nrow(detected)) return(tibble())
    calculated <- detected |>
      group_by(across(all_of(id_col))) |>
      group_modify(function(row, key) {
        series <- monthly |> semi_join(key, by = id_col)
        d <- merge(covar_trend, series[, c("year_month", "items")], by = "year_month")
        d <- d[order(d$t), ]
        d$off <- d$offset_log_patient_days
        fit <- if (identical(row$distribution[1], "negbin")) {
          .nb_fit(f_full, d)
        } else {
          glm(f_full, poisson, data = d)
        }
        if (is.null(fit)) {
          return(tibble(
            peak_trough_ratio_raw = NA_real_, ptr_lci_raw = NA_real_,
            ptr_uci_raw = NA_real_, peak_month = NA_character_,
            trough_month = NA_character_, stl_seasonal_strength_raw = NA_real_
          ))
        }
        seasonal <- .seasonal_curve(
          row$b_sin12[1], row$b_cos12[1], row$b_sin6[1], row$b_cos6[1]
        )
        ci <- .ptr_ci(fit, row$route[1])
        strength <- .stl_strength(series, covar_trend)
        tibble(
          peak_trough_ratio_raw = exp(max(seasonal) - min(seasonal)),
          ptr_lci_raw = ci[1], ptr_uci_raw = ci[2],
          peak_month = month.abb[which.max(seasonal)],
          trough_month = month.abb[which.min(seasonal)],
          stl_seasonal_strength_raw = unname(strength["seasonal"])
        )
      }) |>
      ungroup()
    calculated |>
      left_join(scr, by = id_col) |>
      mutate(
        peak_trough_ratio = round(peak_trough_ratio_raw, 3),
        ptr_lci = round(ptr_lci_raw, 3),
        ptr_uci = round(ptr_uci_raw, 3),
        stl_seasonal_strength = round(stl_seasonal_strength_raw, 3),
        meaningful = !is.na(ptr_lci) & !is.na(stl_seasonal_strength) &
          ptr_lci >= meaningful_threshold &
          stl_seasonal_strength >= stl_strength_threshold
      ) |>
      select(
        all_of(id_col), any_of(c("bnf_class_code", "bnf_class_name",
                                 "bnf_drug_name")),
        analysis_specification, trend_specification,
        multiplicity_family, seasonality_q_bh, seasonality_detected,
        any_of(c("parent_class_q_bh", "parent_class_detected")),
        peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
        stl_seasonal_strength, meaningful, distribution, route
      ) |>
      distinct()
  }

  primary_char_class <- results_class |>
    transmute(
      bnf_class_code, bnf_class_name,
      analysis_specification = "trend_spline3", trend_specification = "ns_df3",
      multiplicity_family = "all_eligible_classes",
      seasonality_q_bh = class_q_bh, seasonality_detected = class_significant,
      peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
      stl_seasonal_strength, meaningful, distribution, route
    )
  primary_char_drug <- results_drug |>
    left_join(
      screen_drug |> select(bnf_drug_code, bnf_class_code),
      by = "bnf_drug_code"
    ) |>
    transmute(
      bnf_class_code, bnf_class_name, bnf_drug_code, bnf_drug_name,
      analysis_specification = "trend_spline3", trend_specification = "ns_df3",
      multiplicity_family = "all_eligible_drugs",
      seasonality_q_bh = drug_all_q_bh, seasonality_detected = drug_significant,
      parent_class_q_bh = class_q_bh,
      parent_class_detected = parent_class_significant,
      peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
      stl_seasonal_strength, meaningful, distribution, route
    )
  trend_char_class <- bind_rows(
    primary_char_class,
    characterise_trend_level(class_linear, class_monthly_elig,
                              "bnf_class_code", f_linear_full),
    characterise_trend_level(class_spline4, class_monthly_elig,
                              "bnf_class_code", f_spline4_full)
  ) |>
    arrange(analysis_specification, desc(meaningful), desc(peak_trough_ratio))
  trend_char_drug <- bind_rows(
    primary_char_drug,
    characterise_trend_level(drug_linear, drug_monthly_elig,
                              "bnf_drug_code", f_linear_full),
    characterise_trend_level(drug_spline4, drug_monthly_elig,
                              "bnf_drug_code", f_spline4_full)
  ) |>
    arrange(analysis_specification, desc(meaningful), desc(peak_trough_ratio))

  summarise_trend <- function(scr, chr, level) {
    detected_summary <- scr |>
      group_by(analysis_specification, trend_specification) |>
      summarise(
        family_size = n(), model_failures = sum(!converged),
        discoveries = sum(seasonality_detected), .groups = "drop"
      )
    meaningful_summary <- chr |>
      group_by(analysis_specification, trend_specification) |>
      summarise(
        characterised = n(), meaningful = sum(meaningful), .groups = "drop"
      )
    detected_summary |>
      left_join(meaningful_summary,
                by = c("analysis_specification", "trend_specification")) |>
      mutate(level = level, .before = 1)
  }
  trend_summary <- as.data.table(bind_rows(
    summarise_trend(trend_screen_class, trend_char_class, "class"),
    summarise_trend(trend_screen_drug, trend_char_drug, "drug")
  ))

  compare_trend_sets <- function(scr, chr, id_col, level) {
    primary_discovery <- as.character(scr |>
      filter(analysis_specification == "trend_spline3", seasonality_detected) |>
      pull(all_of(id_col)))
    primary_meaningful_codes <- as.character(chr |>
      filter(analysis_specification == "trend_spline3", meaningful) |>
      pull(all_of(id_col)))
    bind_rows(lapply(c("trend_linear", "trend_spline4"), function(spec) {
      alternative_discovery <- as.character(scr |>
        filter(analysis_specification == spec, seasonality_detected) |>
        pull(all_of(id_col)))
      alternative_meaningful_codes <- as.character(chr |>
        filter(analysis_specification == spec, meaningful) |>
        pull(all_of(id_col)))
      tibble(
        level = level, comparison = paste0(spec, "_vs_trend_spline3"),
        primary_discoveries = length(primary_discovery),
        alternative_discoveries = length(alternative_discovery),
        discovery_overlap = length(base::intersect(primary_discovery, alternative_discovery)),
        discoveries_added = length(base::setdiff(alternative_discovery, primary_discovery)),
        discoveries_lost = length(base::setdiff(primary_discovery, alternative_discovery)),
        primary_meaningful = length(primary_meaningful_codes),
        alternative_meaningful = length(alternative_meaningful_codes),
        meaningful_overlap = length(base::intersect(
          primary_meaningful_codes, alternative_meaningful_codes
        )),
        meaningful_added = length(base::setdiff(
          alternative_meaningful_codes, primary_meaningful_codes
        )),
        meaningful_lost = length(base::setdiff(
          primary_meaningful_codes, alternative_meaningful_codes
        ))
      )
    }))
  }
  trend_overlap_summary <- bind_rows(
    compare_trend_sets(trend_screen_class, trend_char_class,
                       "bnf_class_code", "class"),
    compare_trend_sets(trend_screen_drug, trend_char_drug,
                       "bnf_drug_code", "drug")
  )

  make_trend_changes <- function(scr, chr, id_col, name_col, level) {
    characterised <- chr |>
      select(all_of(id_col), analysis_specification,
             peak_trough_ratio, ptr_lci, ptr_uci, peak_month,
             stl_seasonal_strength, meaningful)
    full <- scr |>
      left_join(characterised, by = c(id_col, "analysis_specification")) |>
      mutate(meaningful = coalesce(meaningful, FALSE))
    primary <- full |>
      filter(analysis_specification == "trend_spline3") |>
      select(
        all_of(id_col), primary_q_bh = seasonality_q_bh,
        primary_detected = seasonality_detected,
        primary_peak_trough_ratio = peak_trough_ratio,
        primary_ptr_lci = ptr_lci, primary_peak_month = peak_month,
        primary_seasonal_strength = stl_seasonal_strength,
        primary_meaningful = meaningful
      )
    full |>
      filter(analysis_specification != "trend_spline3") |>
      left_join(primary, by = id_col) |>
      mutate(
        discovery_changed = seasonality_detected != primary_detected,
        meaningfulness_changed = meaningful != primary_meaningful,
        level = level,
        code = as.character(.data[[id_col]]),
        series_name = .data[[name_col]]
      ) |>
      filter(discovery_changed | meaningfulness_changed) |>
      select(
        level, code, series_name, analysis_specification,
        seasonality_q_bh, primary_q_bh, seasonality_detected,
        primary_detected, discovery_changed,
        peak_trough_ratio, primary_peak_trough_ratio, ptr_lci,
        primary_ptr_lci, peak_month, primary_peak_month,
        stl_seasonal_strength, primary_seasonal_strength,
        meaningful, primary_meaningful, meaningfulness_changed
      )
  }
  trend_classification_changes <- bind_rows(
    make_trend_changes(trend_screen_class, trend_char_class,
                       "bnf_class_code", "bnf_class_name", "class"),
    make_trend_changes(trend_screen_drug, trend_char_drug,
                       "bnf_drug_code", "bnf_drug_name", "drug")
  ) |>
    arrange(level, analysis_specification, code)

  trend_q_boundary <- bind_rows(
    trend_screen_class |>
      mutate(level = "class", code = as.character(bnf_class_code),
             series_name = bnf_class_name),
    trend_screen_drug |>
      mutate(level = "drug", code = as.character(bnf_drug_code),
             series_name = bnf_drug_name)
  ) |>
    group_by(level, analysis_specification) |>
    slice_min(abs(seasonality_q_bh - fdr_alpha), n = 20, with_ties = FALSE) |>
    ungroup() |>
    select(level, analysis_specification, trend_specification, code,
           series_name, p_value, seasonality_q_bh, seasonality_detected)

  output_files <- c(
    "trend_screen_class.csv", "trend_screen_drug.csv",
    "trend_characterisation_class.csv", "trend_characterisation_drug.csv",
    "trend_summary.csv", "trend_overlap_summary.csv",
    "trend_classification_changes.csv", "trend_q_boundary.csv"
  )
  atomic_fwrite(as.data.table(trend_screen_class), file.path(trend_dir, output_files[1]))
  atomic_fwrite(as.data.table(trend_screen_drug), file.path(trend_dir, output_files[2]))
  atomic_fwrite(as.data.table(trend_char_class), file.path(trend_dir, output_files[3]))
  atomic_fwrite(as.data.table(trend_char_drug), file.path(trend_dir, output_files[4]))
  atomic_fwrite(as.data.table(trend_summary), file.path(trend_dir, output_files[5]))
  atomic_fwrite(as.data.table(trend_overlap_summary), file.path(trend_dir, output_files[6]))
  atomic_fwrite(as.data.table(trend_classification_changes), file.path(trend_dir, output_files[7]))
  atomic_fwrite(as.data.table(trend_q_boundary), file.path(trend_dir, output_files[8]))

  add_trend_check(
    "complete_families",
    all(trend_summary[level == "class", family_size] == 220L) &&
      all(trend_summary[level == "drug", family_size] == 974L),
    "three 220-class and three 974-drug families",
    paste(trend_summary$level, trend_summary$analysis_specification,
          trend_summary$family_size, sep = "=", collapse = ";")
  )
  add_trend_check(
    "all_models_complete", all(trend_summary$model_failures == 0L),
    "0 failed models", sum(trend_summary$model_failures)
  )
  add_trend_check(
    "valid_q_values",
    all(is.finite(trend_screen_class$seasonality_q_bh)) &&
      all(trend_screen_class$seasonality_q_bh >= 0 &
            trend_screen_class$seasonality_q_bh <= 1) &&
      all(is.finite(trend_screen_drug$seasonality_q_bh)) &&
      all(trend_screen_drug$seasonality_q_bh >= 0 &
            trend_screen_drug$seasonality_q_bh <= 1),
    "all q-values finite and in [0,1]", "checked"
  )
  add_trend_check(
    "all_discoveries_characterised",
    all(trend_summary$discoveries == trend_summary$characterised) &&
      !anyNA(trend_char_class$ptr_lci) && !anyNA(trend_char_drug$ptr_lci),
    "every discovery has amplitude CI and characterisation",
    paste(trend_summary$level, trend_summary$analysis_specification,
          paste0(trend_summary$characterised, "/", trend_summary$discoveries),
          sep = "=", collapse = ";")
  )
  primary_summary <- trend_summary[trend_summary$analysis_specification == "trend_spline3", ]
  primary_class_summary <- primary_summary[level == "class"]
  primary_drug_summary <- primary_summary[level == "drug"]
  add_trend_check(
    "primary_counts_preserved",
    primary_class_summary$discoveries == 125L &&
      primary_drug_summary$discoveries == 391L &&
      primary_class_summary$meaningful == 30L &&
      primary_drug_summary$meaningful == 88L,
    "primary 125/391 discoveries and 30/88 meaningful",
    paste(primary_summary$level, primary_summary$discoveries,
          primary_summary$meaningful, sep = "=", collapse = ";")
  )
  add_trend_check(
    "meaningful_subset_of_discoveries",
    all(trend_char_class$seasonality_detected[trend_char_class$meaningful]) &&
      all(trend_char_drug$seasonality_detected[trend_char_drug$meaningful]),
    "all meaningful series are discoveries in the same specification", "checked"
  )
  add_trend_check(
    "comparison_outputs_present",
    all(file.exists(file.path(trend_dir, output_files))) &&
      nrow(trend_overlap_summary) == 4L && nrow(trend_q_boundary) == 120L,
    "8 result files; 4 overlap rows; 120 boundary rows",
    sprintf("%d files; %d overlap rows; %d boundary rows",
            sum(file.exists(file.path(trend_dir, output_files))),
            nrow(trend_overlap_summary), nrow(trend_q_boundary))
  )
  add_trend_check(
    "set_comparisons_reconcile",
    all(trend_overlap_summary$discovery_overlap +
          trend_overlap_summary$discoveries_added ==
          trend_overlap_summary$alternative_discoveries) &&
      all(trend_overlap_summary$discovery_overlap +
            trend_overlap_summary$discoveries_lost ==
            trend_overlap_summary$primary_discoveries) &&
      all(trend_overlap_summary$meaningful_overlap +
            trend_overlap_summary$meaningful_added ==
            trend_overlap_summary$alternative_meaningful) &&
      all(trend_overlap_summary$meaningful_overlap +
            trend_overlap_summary$meaningful_lost ==
            trend_overlap_summary$primary_meaningful),
    "all overlap/add/loss counts reconcile to both compared sets", "checked"
  )

  trend_qc_summary <- rbindlist(trend_checks, use.names = TRUE, fill = TRUE)
  atomic_fwrite(trend_qc_summary, file.path(trend_dir, "stage5_trend_qc_summary.csv"))
  if (any(!trend_qc_summary$pass)) {
    stop("Stage 5.1 completion gate failed: ",
         paste(trend_qc_summary[pass == FALSE, check_id], collapse = ", "),
         ". See ", file.path(trend_dir, "stage5_trend_qc_summary.csv"), ".")
  }

  for (relative_path in output_files) {
    source_path <- file.path(trend_dir, relative_path)
    target_path <- file.path(trend_snapshot_dir, relative_path)
    temporary_path <- paste0(target_path, ".tmp")
    if (!file.copy(source_path, temporary_path, overwrite = TRUE) ||
        !file.rename(temporary_path, target_path)) {
      stop("Could not seal Stage 5.1 snapshot file: ", relative_path)
    }
  }
  trend_manifest <- data.table(
    relative_path = output_files,
    bytes = as.numeric(file.info(file.path(trend_snapshot_dir, output_files))$size),
    sha256 = vapply(file.path(trend_snapshot_dir, output_files),
                    sha256_file, character(1))
  )
  atomic_fwrite(trend_manifest, file.path(trend_dir, "stage5_trend_snapshot_manifest.csv"))
  trend_completion <- data.table(
    stage = "stage5_trend", status = "PASS",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    checks_passed = sum(trend_qc_summary$pass), checks_total = nrow(trend_qc_summary),
    classification_change_rows = nrow(trend_classification_changes),
    snapshot_files = nrow(trend_manifest),
    analysis_script_sha256 = sha256_file(script_path),
    renv_lock_sha256 = sha256_file(lock_path)
  )
  atomic_fwrite(trend_completion, file.path(trend_dir, "stage5_trend_completion.csv"))
  message(
    "Stage 5.1 complete: trend sensitivity fitted for all class/drug families; ",
    nrow(trend_qc_summary), " structural checks passed; stopping for result review."
  )
  return(invisible(TRUE))
}

## ---- STAGE 5.2: uniform Poisson-HAC sensitivity ---------------------------
if (run_stage == "stage5_hac") {
  hac_dir <- file.path(out_dir, "qc", "stage5", "uniform_hac")
  hac_snapshot_dir <- file.path(hac_dir, "uniform_hac_snapshot")
  stage4_dir <- file.path(out_dir, "qc", "stage4")
  stage4_snapshot_dir <- file.path(stage4_dir, "stage4_snapshot")
  trend_dir <- file.path(out_dir, "qc", "stage5", "trend")
  trend_snapshot_dir <- file.path(trend_dir, "trend_snapshot")
  dir.create(hac_snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  hac_checks <- list()
  add_hac_check <- function(check_id, pass, expected, observed, details = "") {
    hac_checks[[length(hac_checks) + 1L]] <<- data.table(
      check_id = check_id, pass = isTRUE(pass), expected = as.character(expected),
      observed = as.character(observed), details = as.character(details)
    )
  }
  verify_snapshot <- function(manifest_path, snapshot_path) {
    if (!file.exists(manifest_path)) return(list(ok = FALSE, n = 0L, matching = 0L))
    manifest <- fread(manifest_path)
    paths <- file.path(snapshot_path, manifest$relative_path)
    if (!all(file.exists(paths))) {
      return(list(ok = FALSE, n = nrow(manifest), matching = sum(file.exists(paths))))
    }
    hashes <- vapply(paths, sha256_file, character(1))
    size_col <- if ("bytes" %in% names(manifest)) "bytes" else "size_bytes"
    matching <- hashes == manifest$sha256 &
      as.numeric(file.info(paths)$size) == manifest[[size_col]]
    list(ok = all(matching), n = nrow(manifest), matching = sum(matching))
  }
  stage4_verification <- verify_snapshot(
    file.path(stage4_dir, "stage4_snapshot_manifest.csv"), stage4_snapshot_dir
  )
  trend_verification <- verify_snapshot(
    file.path(trend_dir, "stage5_trend_snapshot_manifest.csv"), trend_snapshot_dir
  )
  add_hac_check(
    "prior_authorities_unchanged",
    stage4_verification$ok && trend_verification$ok,
    "17 Stage 4 and 8 Stage 5.1 files unchanged",
    sprintf("Stage 4 %d/%d; Stage 5.1 %d/%d",
            stage4_verification$matching, stage4_verification$n,
            trend_verification$matching, trend_verification$n)
  )
  if (!stage4_verification$ok || !trend_verification$ok) {
    stop("Stage 5.2 requires intact Stage 4 and Stage 5.1 snapshots.")
  }

  fit_uniform_hac <- function(series, covar) {
    res <- data.frame(
      p_value = NA_real_, distribution = "poisson",
      route = "HAC-Wald-Uniform", disp_ratio = NA_real_, disp_p = NA_real_,
      lb_p = NA_real_, nw_lag = NA_real_, hac_capped = NA,
      b_sin12 = NA_real_, b_cos12 = NA_real_, b_sin6 = NA_real_,
      b_cos6 = NA_real_, converged = FALSE, note = NA_character_,
      stringsAsFactors = FALSE
    )
    tryCatch({
      d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
      d <- d[order(d$t), ]
      if (nrow(d) != nrow(covar) || anyNA(d$items)) {
        stop("series does not cover all 48 months")
      }
      d$off <- d$offset_log_patient_days
      fp <- glm(.f_full, family = poisson, data = d)
      if (!isTRUE(fp$converged)) stop("Poisson full model did not converge")
      pear <- residuals(fp, type = "pearson")
      res$disp_ratio <- sum(pear^2) / df.residual(fp)
      res$disp_p <- tryCatch(
        AER::dispersiontest(fp)$p.value, error = function(e) NA_real_
      )
      res$lb_p <- tryCatch(
        Box.test(pear, lag = harmonic_period_months[1],
                 type = "Ljung-Box")$p.value,
        error = function(e) NA_real_
      )
      hac <- .do_hac(fp)
      coefficients <- coef(fp)
      res$p_value <- hac$p
      res$nw_lag <- hac$bw
      res$hac_capped <- hac$capped
      res$b_sin12 <- unname(coefficients["sin12"])
      res$b_cos12 <- unname(coefficients["cos12"])
      res$b_sin6 <- unname(coefficients["sin6"])
      res$b_cos6 <- unname(coefficients["cos6"])
      res$converged <- TRUE
    }, error = function(e) res$note <<- conditionMessage(e))
    res
  }

  fit_uniform_level <- function(monthly, id_cols, multiplicity_family) {
    monthly |>
      select(all_of(c(id_cols, "year_month", "items"))) |>
      group_by(across(all_of(id_cols))) |>
      group_modify(~ fit_uniform_hac(.x, covar)) |>
      ungroup() |>
      mutate(
        analysis_specification = "uniform_poisson_hac",
        inference_specification = "uniform_poisson_hac",
        multiplicity_family = multiplicity_family,
        seasonality_q_bh = p.adjust(p_value, method = "BH"),
        seasonality_detected = !is.na(seasonality_q_bh) &
          seasonality_q_bh < fdr_alpha
      ) |>
      arrange(seasonality_q_bh)
  }

  uniform_class <- fit_uniform_level(
    class_monthly_elig, c("bnf_class_code", "bnf_class_name"),
    "all_eligible_classes"
  )
  uniform_drug <- fit_uniform_level(
    drug_monthly_elig,
    c("bnf_class_code", "bnf_drug_code", "bnf_drug_name"),
    "all_eligible_drugs"
  ) |>
    left_join(
      uniform_class |>
        select(bnf_class_code, bnf_class_name,
               parent_class_q_bh = seasonality_q_bh,
               parent_class_detected = seasonality_detected),
      by = "bnf_class_code"
    )

  primary_hac_class <- screen_class |>
    transmute(
      bnf_class_code, bnf_class_name,
      analysis_specification = "primary_routed",
      inference_specification = "diagnostic_routed",
      multiplicity_family = "all_eligible_classes",
      p_value, seasonality_q_bh = class_q_bh,
      seasonality_detected = class_significant,
      distribution, route, disp_ratio, disp_p, lb_p, nw_lag, hac_capped,
      b_sin12, b_cos12, b_sin6, b_cos6, converged, note
    )
  primary_hac_drug <- screen_drug |>
    transmute(
      bnf_class_code, bnf_class_name, bnf_drug_code, bnf_drug_name,
      analysis_specification = "primary_routed",
      inference_specification = "diagnostic_routed",
      multiplicity_family = "all_eligible_drugs",
      p_value, seasonality_q_bh = drug_all_q_bh,
      seasonality_detected = drug_significant,
      parent_class_q_bh = class_q_bh,
      parent_class_detected = parent_class_significant,
      distribution, route, disp_ratio, disp_p, lb_p, nw_lag, hac_capped,
      b_sin12, b_cos12, b_sin6, b_cos6, converged, note
    )
  hac_screen_class <- bind_rows(primary_hac_class, uniform_class) |>
    arrange(analysis_specification, seasonality_q_bh, bnf_class_code)
  hac_screen_drug <- bind_rows(primary_hac_drug, uniform_drug) |>
    arrange(analysis_specification, seasonality_q_bh, bnf_drug_code)

  characterise_uniform_level <- function(scr, monthly, id_col) {
    detected <- scr |> filter(seasonality_detected)
    if (!nrow(detected)) return(tibble())
    calculated <- detected |>
      group_by(across(all_of(id_col))) |>
      group_modify(function(row, key) {
        series <- monthly |> semi_join(key, by = id_col)
        d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
        d <- d[order(d$t), ]
        d$off <- d$offset_log_patient_days
        fit <- glm(.f_full, poisson, data = d)
        seasonal <- .seasonal_curve(
          row$b_sin12[1], row$b_cos12[1], row$b_sin6[1], row$b_cos6[1]
        )
        ci <- .ptr_ci(fit, "HAC-Wald-Uniform")
        strength <- .stl_strength(series, covar)
        tibble(
          peak_trough_ratio_raw = exp(max(seasonal) - min(seasonal)),
          ptr_lci_raw = ci[1], ptr_uci_raw = ci[2],
          peak_month = month.abb[which.max(seasonal)],
          trough_month = month.abb[which.min(seasonal)],
          stl_seasonal_strength_raw = unname(strength["seasonal"])
        )
      }) |>
      ungroup()
    calculated |>
      left_join(scr, by = id_col) |>
      mutate(
        peak_trough_ratio = round(peak_trough_ratio_raw, 3),
        ptr_lci = round(ptr_lci_raw, 3), ptr_uci = round(ptr_uci_raw, 3),
        stl_seasonal_strength = round(stl_seasonal_strength_raw, 3),
        meaningful = !is.na(ptr_lci) & !is.na(stl_seasonal_strength) &
          ptr_lci >= meaningful_threshold &
          stl_seasonal_strength >= stl_strength_threshold
      ) |>
      select(
        all_of(id_col), any_of(c("bnf_class_code", "bnf_class_name",
                                 "bnf_drug_name")),
        analysis_specification, inference_specification,
        multiplicity_family, seasonality_q_bh, seasonality_detected,
        any_of(c("parent_class_q_bh", "parent_class_detected")),
        peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
        stl_seasonal_strength, meaningful, distribution, route
      ) |>
      distinct()
  }

  primary_hac_char_class <- results_class |>
    transmute(
      bnf_class_code, bnf_class_name,
      analysis_specification = "primary_routed",
      inference_specification = "diagnostic_routed",
      multiplicity_family = "all_eligible_classes",
      seasonality_q_bh = class_q_bh, seasonality_detected = class_significant,
      peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
      stl_seasonal_strength, meaningful, distribution, route
    )
  primary_hac_char_drug <- results_drug |>
    left_join(
      screen_drug |> select(bnf_drug_code, bnf_class_code),
      by = "bnf_drug_code"
    ) |>
    transmute(
      bnf_class_code, bnf_class_name, bnf_drug_code, bnf_drug_name,
      analysis_specification = "primary_routed",
      inference_specification = "diagnostic_routed",
      multiplicity_family = "all_eligible_drugs",
      seasonality_q_bh = drug_all_q_bh, seasonality_detected = drug_significant,
      parent_class_q_bh = class_q_bh,
      parent_class_detected = parent_class_significant,
      peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
      stl_seasonal_strength, meaningful, distribution, route
    )
  hac_char_class <- bind_rows(
    primary_hac_char_class,
    characterise_uniform_level(uniform_class, class_monthly_elig,
                               "bnf_class_code")
  ) |>
    arrange(analysis_specification, desc(meaningful), desc(peak_trough_ratio))
  hac_char_drug <- bind_rows(
    primary_hac_char_drug,
    characterise_uniform_level(uniform_drug, drug_monthly_elig,
                               "bnf_drug_code")
  ) |>
    arrange(analysis_specification, desc(meaningful), desc(peak_trough_ratio))

  summarise_hac <- function(scr, chr, level) {
    screen_summary <- scr |>
      group_by(analysis_specification, inference_specification) |>
      summarise(
        family_size = n(), model_failures = sum(!converged),
        discoveries = sum(seasonality_detected), .groups = "drop"
      )
    char_summary <- chr |>
      group_by(analysis_specification, inference_specification) |>
      summarise(characterised = n(), meaningful = sum(meaningful),
                .groups = "drop")
    screen_summary |>
      left_join(char_summary,
                by = c("analysis_specification", "inference_specification")) |>
      mutate(level = level, .before = 1)
  }
  hac_summary <- as.data.table(bind_rows(
    summarise_hac(hac_screen_class, hac_char_class, "class"),
    summarise_hac(hac_screen_drug, hac_char_drug, "drug")
  ))

  compare_hac_sets <- function(scr, chr, id_col, level) {
    primary_discovery_codes <- as.character(scr |>
      filter(analysis_specification == "primary_routed", seasonality_detected) |>
      pull(all_of(id_col)))
    uniform_discovery_codes <- as.character(scr |>
      filter(analysis_specification == "uniform_poisson_hac", seasonality_detected) |>
      pull(all_of(id_col)))
    primary_meaningful_codes <- as.character(chr |>
      filter(analysis_specification == "primary_routed", meaningful) |>
      pull(all_of(id_col)))
    uniform_meaningful_codes <- as.character(chr |>
      filter(analysis_specification == "uniform_poisson_hac", meaningful) |>
      pull(all_of(id_col)))
    tibble(
      level = level, comparison = "uniform_poisson_hac_vs_primary_routed",
      primary_discoveries = length(primary_discovery_codes),
      alternative_discoveries = length(uniform_discovery_codes),
      discovery_overlap = length(base::intersect(
        primary_discovery_codes, uniform_discovery_codes
      )),
      discoveries_added = length(base::setdiff(
        uniform_discovery_codes, primary_discovery_codes
      )),
      discoveries_lost = length(base::setdiff(
        primary_discovery_codes, uniform_discovery_codes
      )),
      primary_meaningful = length(primary_meaningful_codes),
      alternative_meaningful = length(uniform_meaningful_codes),
      meaningful_overlap = length(base::intersect(
        primary_meaningful_codes, uniform_meaningful_codes
      )),
      meaningful_added = length(base::setdiff(
        uniform_meaningful_codes, primary_meaningful_codes
      )),
      meaningful_lost = length(base::setdiff(
        primary_meaningful_codes, uniform_meaningful_codes
      ))
    )
  }
  hac_overlap_summary <- bind_rows(
    compare_hac_sets(hac_screen_class, hac_char_class,
                     "bnf_class_code", "class"),
    compare_hac_sets(hac_screen_drug, hac_char_drug,
                     "bnf_drug_code", "drug")
  )

  make_hac_changes <- function(scr, chr, id_col, name_col, level) {
    characterised <- chr |>
      select(all_of(id_col), analysis_specification,
             peak_trough_ratio, ptr_lci, ptr_uci, peak_month,
             stl_seasonal_strength, meaningful)
    full <- scr |>
      left_join(characterised, by = c(id_col, "analysis_specification")) |>
      mutate(meaningful = coalesce(meaningful, FALSE))
    primary <- full |>
      filter(analysis_specification == "primary_routed") |>
      select(
        all_of(id_col), primary_q_bh = seasonality_q_bh,
        primary_detected = seasonality_detected,
        primary_peak_trough_ratio = peak_trough_ratio,
        primary_ptr_lci = ptr_lci, primary_peak_month = peak_month,
        primary_seasonal_strength = stl_seasonal_strength,
        primary_meaningful = meaningful, primary_route = route
      )
    full |>
      filter(analysis_specification == "uniform_poisson_hac") |>
      left_join(primary, by = id_col) |>
      mutate(
        discovery_changed = seasonality_detected != primary_detected,
        meaningfulness_changed = meaningful != primary_meaningful,
        level = level, code = as.character(.data[[id_col]]),
        series_name = .data[[name_col]]
      ) |>
      filter(discovery_changed | meaningfulness_changed) |>
      select(
        level, code, series_name, seasonality_q_bh, primary_q_bh,
        seasonality_detected, primary_detected, discovery_changed,
        route, primary_route, peak_trough_ratio, primary_peak_trough_ratio,
        ptr_lci, primary_ptr_lci, peak_month, primary_peak_month,
        stl_seasonal_strength, primary_seasonal_strength,
        meaningful, primary_meaningful, meaningfulness_changed
      )
  }
  hac_classification_changes <- bind_rows(
    make_hac_changes(hac_screen_class, hac_char_class,
                     "bnf_class_code", "bnf_class_name", "class"),
    make_hac_changes(hac_screen_drug, hac_char_drug,
                     "bnf_drug_code", "bnf_drug_name", "drug")
  ) |>
    arrange(level, code)

  hac_q_boundary <- bind_rows(
    hac_screen_class |>
      mutate(level = "class", code = as.character(bnf_class_code),
             series_name = bnf_class_name),
    hac_screen_drug |>
      mutate(level = "drug", code = as.character(bnf_drug_code),
             series_name = bnf_drug_name)
  ) |>
    group_by(level, analysis_specification) |>
    slice_min(abs(seasonality_q_bh - fdr_alpha), n = 20, with_ties = FALSE) |>
    ungroup() |>
    select(level, analysis_specification, inference_specification, code,
           series_name, p_value, seasonality_q_bh, seasonality_detected, route)

  hac_bandwidth_summary <- bind_rows(
    uniform_class |> mutate(level = "class"),
    uniform_drug |> mutate(level = "drug")
  ) |>
    group_by(level) |>
    summarise(
      n = n(), n_capped = sum(hac_capped),
      min_lag = min(nw_lag), median_lag = median(nw_lag),
      max_lag = max(nw_lag), .groups = "drop"
    )

  output_files <- c(
    "uniform_hac_screen_class.csv", "uniform_hac_screen_drug.csv",
    "uniform_hac_characterisation_class.csv",
    "uniform_hac_characterisation_drug.csv", "uniform_hac_summary.csv",
    "uniform_hac_overlap_summary.csv", "uniform_hac_classification_changes.csv",
    "uniform_hac_q_boundary.csv", "uniform_hac_bandwidth_summary.csv"
  )
  atomic_fwrite(as.data.table(hac_screen_class), file.path(hac_dir, output_files[1]))
  atomic_fwrite(as.data.table(hac_screen_drug), file.path(hac_dir, output_files[2]))
  atomic_fwrite(as.data.table(hac_char_class), file.path(hac_dir, output_files[3]))
  atomic_fwrite(as.data.table(hac_char_drug), file.path(hac_dir, output_files[4]))
  atomic_fwrite(as.data.table(hac_summary), file.path(hac_dir, output_files[5]))
  atomic_fwrite(as.data.table(hac_overlap_summary), file.path(hac_dir, output_files[6]))
  atomic_fwrite(as.data.table(hac_classification_changes),
                file.path(hac_dir, output_files[7]))
  atomic_fwrite(as.data.table(hac_q_boundary), file.path(hac_dir, output_files[8]))
  atomic_fwrite(as.data.table(hac_bandwidth_summary),
                file.path(hac_dir, output_files[9]))

  add_hac_check(
    "complete_families",
    all(hac_summary[level == "class", family_size] == 220L) &&
      all(hac_summary[level == "drug", family_size] == 974L),
    "two 220-class and two 974-drug families",
    paste(hac_summary$level, hac_summary$analysis_specification,
          hac_summary$family_size, sep = "=", collapse = ";")
  )
  add_hac_check(
    "all_models_complete", all(hac_summary$model_failures == 0L),
    "0 failed models", sum(hac_summary$model_failures)
  )
  add_hac_check(
    "uniform_route_applied",
    all(uniform_class$route == "HAC-Wald-Uniform") &&
      all(uniform_drug$route == "HAC-Wald-Uniform") &&
      all(uniform_class$distribution == "poisson") &&
      all(uniform_drug$distribution == "poisson"),
    "all 1194 sensitivity models use Poisson-HAC", "checked"
  )
  add_hac_check(
    "valid_uniform_hac_settings",
    all(is.finite(c(uniform_class$nw_lag, uniform_drug$nw_lag))) &&
      all(c(uniform_class$nw_lag, uniform_drug$nw_lag) >= 0) &&
      all(c(uniform_class$nw_lag, uniform_drug$nw_lag) <= 12),
    "all bandwidths finite and capped at 12",
    sprintf("range %.3f to %.3f; capped=%d",
            min(c(uniform_class$nw_lag, uniform_drug$nw_lag)),
            max(c(uniform_class$nw_lag, uniform_drug$nw_lag)),
            sum(c(uniform_class$hac_capped, uniform_drug$hac_capped)))
  )
  add_hac_check(
    "valid_q_values",
    all(is.finite(uniform_class$seasonality_q_bh)) &&
      all(uniform_class$seasonality_q_bh >= 0 & uniform_class$seasonality_q_bh <= 1) &&
      all(is.finite(uniform_drug$seasonality_q_bh)) &&
      all(uniform_drug$seasonality_q_bh >= 0 & uniform_drug$seasonality_q_bh <= 1),
    "all uniform-HAC q-values finite and in [0,1]", "checked"
  )
  add_hac_check(
    "all_discoveries_characterised",
    all(hac_summary$discoveries == hac_summary$characterised) &&
      !anyNA(hac_char_class$ptr_lci) && !anyNA(hac_char_drug$ptr_lci),
    "every discovery has amplitude CI and characterisation",
    paste(hac_summary$level, hac_summary$analysis_specification,
          paste0(hac_summary$characterised, "/", hac_summary$discoveries),
          sep = "=", collapse = ";")
  )
  primary_hac_summary <- hac_summary[analysis_specification == "primary_routed"]
  add_hac_check(
    "primary_counts_preserved",
    primary_hac_summary[level == "class", discoveries] == 125L &&
      primary_hac_summary[level == "drug", discoveries] == 391L &&
      primary_hac_summary[level == "class", meaningful] == 30L &&
      primary_hac_summary[level == "drug", meaningful] == 88L,
    "primary 125/391 discoveries and 30/88 meaningful",
    paste(primary_hac_summary$level, primary_hac_summary$discoveries,
          primary_hac_summary$meaningful, sep = "=", collapse = ";")
  )
  add_hac_check(
    "set_comparisons_reconcile",
    all(hac_overlap_summary$discovery_overlap +
          hac_overlap_summary$discoveries_added ==
          hac_overlap_summary$alternative_discoveries) &&
      all(hac_overlap_summary$discovery_overlap +
            hac_overlap_summary$discoveries_lost ==
            hac_overlap_summary$primary_discoveries) &&
      all(hac_overlap_summary$meaningful_overlap +
            hac_overlap_summary$meaningful_added ==
            hac_overlap_summary$alternative_meaningful) &&
      all(hac_overlap_summary$meaningful_overlap +
            hac_overlap_summary$meaningful_lost ==
            hac_overlap_summary$primary_meaningful),
    "all overlap/add/loss counts reconcile to both compared sets", "checked"
  )
  add_hac_check(
    "comparison_outputs_present",
    all(file.exists(file.path(hac_dir, output_files))) &&
      nrow(hac_overlap_summary) == 2L && nrow(hac_q_boundary) == 80L,
    "9 result files; 2 overlap rows; 80 boundary rows",
    sprintf("%d files; %d overlap rows; %d boundary rows",
            sum(file.exists(file.path(hac_dir, output_files))),
            nrow(hac_overlap_summary), nrow(hac_q_boundary))
  )

  hac_qc_summary <- rbindlist(hac_checks, use.names = TRUE, fill = TRUE)
  atomic_fwrite(hac_qc_summary, file.path(hac_dir, "stage5_hac_qc_summary.csv"))
  if (any(!hac_qc_summary$pass)) {
    stop("Stage 5.2 completion gate failed: ",
         paste(hac_qc_summary[pass == FALSE, check_id], collapse = ", "),
         ". See ", file.path(hac_dir, "stage5_hac_qc_summary.csv"), ".")
  }

  for (relative_path in output_files) {
    source_path <- file.path(hac_dir, relative_path)
    target_path <- file.path(hac_snapshot_dir, relative_path)
    temporary_path <- paste0(target_path, ".tmp")
    if (!file.copy(source_path, temporary_path, overwrite = TRUE) ||
        !file.rename(temporary_path, target_path)) {
      stop("Could not seal Stage 5.2 snapshot file: ", relative_path)
    }
  }
  hac_manifest <- data.table(
    relative_path = output_files,
    bytes = as.numeric(file.info(file.path(hac_snapshot_dir, output_files))$size),
    sha256 = vapply(file.path(hac_snapshot_dir, output_files),
                    sha256_file, character(1))
  )
  atomic_fwrite(hac_manifest, file.path(hac_dir, "stage5_hac_snapshot_manifest.csv"))
  hac_completion <- data.table(
    stage = "stage5_hac", status = "PASS",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    checks_passed = sum(hac_qc_summary$pass), checks_total = nrow(hac_qc_summary),
    classification_change_rows = nrow(hac_classification_changes),
    snapshot_files = nrow(hac_manifest),
    analysis_script_sha256 = sha256_file(script_path),
    renv_lock_sha256 = sha256_file(lock_path)
  )
  atomic_fwrite(hac_completion, file.path(hac_dir, "stage5_hac_completion.csv"))
  message(
    "Stage 5.2 complete: uniform Poisson-HAC fitted for all class/drug series; ",
    nrow(hac_qc_summary), " checks passed; stopping for result review."
  )
  return(invisible(TRUE))
}

## ---- STAGE 5.3: working-day offset sensitivity ----------------------------
if (run_stage == "stage5_working_days") {
  working_dir <- file.path(out_dir, "qc", "stage5", "working_days")
  working_snapshot_dir <- file.path(working_dir, "working_days_snapshot")
  dir.create(working_snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  working_checks <- list()
  add_working_check <- function(check_id, pass, expected, observed, details = "") {
    working_checks[[length(working_checks) + 1L]] <<- data.table(
      check_id = check_id, pass = isTRUE(pass), expected = as.character(expected),
      observed = as.character(observed), details = as.character(details)
    )
  }
  verify_prior_snapshot <- function(manifest_path, snapshot_path) {
    if (!file.exists(manifest_path)) return(c(ok = FALSE, n = 0, matching = 0))
    manifest <- fread(manifest_path)
    paths <- file.path(snapshot_path, manifest$relative_path)
    if (!all(file.exists(paths))) {
      return(c(ok = FALSE, n = nrow(manifest), matching = sum(file.exists(paths))))
    }
    hashes <- vapply(paths, sha256_file, character(1))
    size_col <- if ("bytes" %in% names(manifest)) "bytes" else "size_bytes"
    matching <- hashes == manifest$sha256 &
      as.numeric(file.info(paths)$size) == manifest[[size_col]]
    c(ok = all(matching), n = nrow(manifest), matching = sum(matching))
  }
  prior_specs <- list(
    stage4 = c(file.path(out_dir, "qc", "stage4", "stage4_snapshot_manifest.csv"),
               file.path(out_dir, "qc", "stage4", "stage4_snapshot")),
    trend = c(file.path(out_dir, "qc", "stage5", "trend",
                        "stage5_trend_snapshot_manifest.csv"),
              file.path(out_dir, "qc", "stage5", "trend", "trend_snapshot")),
    hac = c(file.path(out_dir, "qc", "stage5", "uniform_hac",
                      "stage5_hac_snapshot_manifest.csv"),
            file.path(out_dir, "qc", "stage5", "uniform_hac",
                      "uniform_hac_snapshot"))
  )
  prior_verification <- lapply(prior_specs, function(x) verify_prior_snapshot(x[1], x[2]))
  prior_ok <- all(vapply(prior_verification, function(x) as.logical(x["ok"]), logical(1)))
  add_working_check(
    "prior_authorities_unchanged", prior_ok,
    "17 Stage 4, 8 trend and 9 uniform-HAC files unchanged",
    paste(names(prior_verification), vapply(prior_verification, function(x) {
      sprintf("%d/%d", as.integer(x["matching"]), as.integer(x["n"]))
    }, character(1)), collapse = ";")
  )
  if (!prior_ok) stop("Stage 5.3 requires intact prior-stage snapshots.")

  fit_working_level <- function(monthly, id_cols, multiplicity_family) {
    monthly |>
      select(all_of(c(id_cols, "year_month", "items"))) |>
      group_by(across(all_of(id_cols))) |>
      group_modify(~ fit_test_series(
        .x, covar, f_full = .f_full, f_red = .f_red,
        offset_col = "offset_log_patient_working_days"
      )) |>
      ungroup() |>
      mutate(
        analysis_specification = "working_day_offset",
        offset_specification = "working_days",
        multiplicity_family = multiplicity_family,
        seasonality_q_bh = p.adjust(p_value, method = "BH"),
        seasonality_detected = !is.na(seasonality_q_bh) &
          seasonality_q_bh < fdr_alpha
      ) |>
      arrange(seasonality_q_bh)
  }

  working_class <- fit_working_level(
    class_monthly_elig, c("bnf_class_code", "bnf_class_name"),
    "all_eligible_classes"
  )
  working_drug <- fit_working_level(
    drug_monthly_elig,
    c("bnf_class_code", "bnf_drug_code", "bnf_drug_name"),
    "all_eligible_drugs"
  ) |>
    left_join(
      working_class |>
        select(bnf_class_code, bnf_class_name,
               parent_class_q_bh = seasonality_q_bh,
               parent_class_detected = seasonality_detected),
      by = "bnf_class_code"
    )

  primary_working_class <- screen_class |>
    transmute(
      bnf_class_code, bnf_class_name,
      analysis_specification = "primary_calendar_day_offset",
      offset_specification = "calendar_days",
      multiplicity_family = "all_eligible_classes",
      p_value, seasonality_q_bh = class_q_bh,
      seasonality_detected = class_significant,
      distribution, route, disp_ratio, disp_p, lb_p, nw_lag, hac_capped,
      theta, b_sin12, b_cos12, b_sin6, b_cos6, converged, note
    )
  primary_working_drug <- screen_drug |>
    transmute(
      bnf_class_code, bnf_class_name, bnf_drug_code, bnf_drug_name,
      analysis_specification = "primary_calendar_day_offset",
      offset_specification = "calendar_days",
      multiplicity_family = "all_eligible_drugs",
      p_value, seasonality_q_bh = drug_all_q_bh,
      seasonality_detected = drug_significant,
      parent_class_q_bh = class_q_bh,
      parent_class_detected = parent_class_significant,
      distribution, route, disp_ratio, disp_p, lb_p, nw_lag, hac_capped,
      theta, b_sin12, b_cos12, b_sin6, b_cos6, converged, note
    )
  working_screen_class <- bind_rows(primary_working_class, working_class) |>
    arrange(analysis_specification, seasonality_q_bh, bnf_class_code)
  working_screen_drug <- bind_rows(primary_working_drug, working_drug) |>
    arrange(analysis_specification, seasonality_q_bh, bnf_drug_code)

  characterise_working_level <- function(scr, monthly, id_col) {
    detected <- scr |> filter(seasonality_detected)
    if (!nrow(detected)) return(tibble())
    calculated <- detected |>
      group_by(across(all_of(id_col))) |>
      group_modify(function(row, key) {
        series <- monthly |> semi_join(key, by = id_col)
        d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
        d <- d[order(d$t), ]
        d$off <- d$offset_log_patient_working_days
        fit <- if (identical(row$distribution[1], "negbin")) {
          .nb_fit(.f_full, d)
        } else {
          glm(.f_full, poisson, data = d)
        }
        if (is.null(fit)) {
          return(tibble(
            peak_trough_ratio_raw = NA_real_, ptr_lci_raw = NA_real_,
            ptr_uci_raw = NA_real_, peak_month = NA_character_,
            trough_month = NA_character_, stl_seasonal_strength_raw = NA_real_
          ))
        }
        seasonal <- .seasonal_curve(
          row$b_sin12[1], row$b_cos12[1], row$b_sin6[1], row$b_cos6[1]
        )
        ci <- .ptr_ci(fit, row$route[1])
        strength <- .stl_strength(
          series, covar, offset_col = "offset_log_patient_working_days"
        )
        tibble(
          peak_trough_ratio_raw = exp(max(seasonal) - min(seasonal)),
          ptr_lci_raw = ci[1], ptr_uci_raw = ci[2],
          peak_month = month.abb[which.max(seasonal)],
          trough_month = month.abb[which.min(seasonal)],
          stl_seasonal_strength_raw = unname(strength["seasonal"])
        )
      }) |>
      ungroup()
    calculated |>
      left_join(scr, by = id_col) |>
      mutate(
        peak_trough_ratio = round(peak_trough_ratio_raw, 3),
        ptr_lci = round(ptr_lci_raw, 3), ptr_uci = round(ptr_uci_raw, 3),
        stl_seasonal_strength = round(stl_seasonal_strength_raw, 3),
        meaningful = !is.na(ptr_lci) & !is.na(stl_seasonal_strength) &
          ptr_lci >= meaningful_threshold &
          stl_seasonal_strength >= stl_strength_threshold
      ) |>
      select(
        all_of(id_col), any_of(c("bnf_class_code", "bnf_class_name",
                                 "bnf_drug_name")),
        analysis_specification, offset_specification,
        multiplicity_family, seasonality_q_bh, seasonality_detected,
        any_of(c("parent_class_q_bh", "parent_class_detected")),
        peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
        stl_seasonal_strength, meaningful, distribution, route
      ) |>
      distinct()
  }

  primary_working_char_class <- results_class |>
    transmute(
      bnf_class_code, bnf_class_name,
      analysis_specification = "primary_calendar_day_offset",
      offset_specification = "calendar_days",
      multiplicity_family = "all_eligible_classes",
      seasonality_q_bh = class_q_bh, seasonality_detected = class_significant,
      peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
      stl_seasonal_strength, meaningful, distribution, route
    )
  primary_working_char_drug <- results_drug |>
    left_join(
      screen_drug |> select(bnf_drug_code, bnf_class_code),
      by = "bnf_drug_code"
    ) |>
    transmute(
      bnf_class_code, bnf_class_name, bnf_drug_code, bnf_drug_name,
      analysis_specification = "primary_calendar_day_offset",
      offset_specification = "calendar_days",
      multiplicity_family = "all_eligible_drugs",
      seasonality_q_bh = drug_all_q_bh, seasonality_detected = drug_significant,
      parent_class_q_bh = class_q_bh,
      parent_class_detected = parent_class_significant,
      peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
      stl_seasonal_strength, meaningful, distribution, route
    )
  working_char_class <- bind_rows(
    primary_working_char_class,
    characterise_working_level(working_class, class_monthly_elig,
                               "bnf_class_code")
  ) |>
    arrange(analysis_specification, desc(meaningful), desc(peak_trough_ratio))
  working_char_drug <- bind_rows(
    primary_working_char_drug,
    characterise_working_level(working_drug, drug_monthly_elig,
                               "bnf_drug_code")
  ) |>
    arrange(analysis_specification, desc(meaningful), desc(peak_trough_ratio))

  summarise_working <- function(scr, chr, level) {
    scr |>
      group_by(analysis_specification, offset_specification) |>
      summarise(
        family_size = n(), model_failures = sum(!converged),
        discoveries = sum(seasonality_detected), .groups = "drop"
      ) |>
      left_join(
        chr |>
          group_by(analysis_specification, offset_specification) |>
          summarise(characterised = n(), meaningful = sum(meaningful),
                    .groups = "drop"),
        by = c("analysis_specification", "offset_specification")
      ) |>
      mutate(level = level, .before = 1)
  }
  working_summary <- as.data.table(bind_rows(
    summarise_working(working_screen_class, working_char_class, "class"),
    summarise_working(working_screen_drug, working_char_drug, "drug")
  ))

  compare_working_sets <- function(scr, chr, id_col, level) {
    primary_discovery_codes <- as.character(scr |>
      filter(analysis_specification == "primary_calendar_day_offset",
             seasonality_detected) |>
      pull(all_of(id_col)))
    working_discovery_codes <- as.character(scr |>
      filter(analysis_specification == "working_day_offset",
             seasonality_detected) |>
      pull(all_of(id_col)))
    primary_meaningful_codes <- as.character(chr |>
      filter(analysis_specification == "primary_calendar_day_offset", meaningful) |>
      pull(all_of(id_col)))
    working_meaningful_codes <- as.character(chr |>
      filter(analysis_specification == "working_day_offset", meaningful) |>
      pull(all_of(id_col)))
    tibble(
      level = level, comparison = "working_days_vs_calendar_days",
      primary_discoveries = length(primary_discovery_codes),
      alternative_discoveries = length(working_discovery_codes),
      discovery_overlap = length(base::intersect(
        primary_discovery_codes, working_discovery_codes
      )),
      discoveries_added = length(base::setdiff(
        working_discovery_codes, primary_discovery_codes
      )),
      discoveries_lost = length(base::setdiff(
        primary_discovery_codes, working_discovery_codes
      )),
      primary_meaningful = length(primary_meaningful_codes),
      alternative_meaningful = length(working_meaningful_codes),
      meaningful_overlap = length(base::intersect(
        primary_meaningful_codes, working_meaningful_codes
      )),
      meaningful_added = length(base::setdiff(
        working_meaningful_codes, primary_meaningful_codes
      )),
      meaningful_lost = length(base::setdiff(
        primary_meaningful_codes, working_meaningful_codes
      ))
    )
  }
  working_overlap_summary <- bind_rows(
    compare_working_sets(working_screen_class, working_char_class,
                         "bnf_class_code", "class"),
    compare_working_sets(working_screen_drug, working_char_drug,
                         "bnf_drug_code", "drug")
  )

  make_working_changes <- function(scr, chr, id_col, name_col, level) {
    characterised <- chr |>
      select(all_of(id_col), analysis_specification,
             peak_trough_ratio, ptr_lci, ptr_uci, peak_month,
             stl_seasonal_strength, meaningful)
    full <- scr |>
      left_join(characterised, by = c(id_col, "analysis_specification")) |>
      mutate(meaningful = coalesce(meaningful, FALSE))
    primary <- full |>
      filter(analysis_specification == "primary_calendar_day_offset") |>
      select(
        all_of(id_col), primary_q_bh = seasonality_q_bh,
        primary_detected = seasonality_detected,
        primary_peak_trough_ratio = peak_trough_ratio,
        primary_ptr_lci = ptr_lci, primary_peak_month = peak_month,
        primary_seasonal_strength = stl_seasonal_strength,
        primary_meaningful = meaningful, primary_route = route
      )
    full |>
      filter(analysis_specification == "working_day_offset") |>
      left_join(primary, by = id_col) |>
      mutate(
        discovery_changed = seasonality_detected != primary_detected,
        meaningfulness_changed = meaningful != primary_meaningful,
        level = level, code = as.character(.data[[id_col]]),
        series_name = .data[[name_col]]
      ) |>
      filter(discovery_changed | meaningfulness_changed) |>
      select(
        level, code, series_name, seasonality_q_bh, primary_q_bh,
        seasonality_detected, primary_detected, discovery_changed,
        route, primary_route, peak_trough_ratio, primary_peak_trough_ratio,
        ptr_lci, primary_ptr_lci, peak_month, primary_peak_month,
        stl_seasonal_strength, primary_seasonal_strength,
        meaningful, primary_meaningful, meaningfulness_changed
      )
  }
  working_classification_changes <- bind_rows(
    make_working_changes(working_screen_class, working_char_class,
                         "bnf_class_code", "bnf_class_name", "class"),
    make_working_changes(working_screen_drug, working_char_drug,
                         "bnf_drug_code", "bnf_drug_name", "drug")
  ) |>
    arrange(level, code)

  working_q_boundary <- bind_rows(
    working_screen_class |>
      mutate(level = "class", code = as.character(bnf_class_code),
             series_name = bnf_class_name),
    working_screen_drug |>
      mutate(level = "drug", code = as.character(bnf_drug_code),
             series_name = bnf_drug_name)
  ) |>
    group_by(level, analysis_specification) |>
    slice_min(abs(seasonality_q_bh - fdr_alpha), n = 20, with_ties = FALSE) |>
    ungroup() |>
    select(level, analysis_specification, offset_specification, code,
           series_name, p_value, seasonality_q_bh, seasonality_detected, route)

  working_route_summary <- bind_rows(
    working_screen_class |> mutate(level = "class"),
    working_screen_drug |> mutate(level = "drug")
  ) |>
    count(level, analysis_specification, offset_specification,
          distribution, route, name = "n") |>
    group_by(level, analysis_specification) |>
    mutate(pct = round(100 * n / sum(n), 1)) |>
    ungroup()
  working_offset_monthly <- covar |>
    transmute(
      year_month, month_date, list_size, days_in_month, working_days,
      calendar_to_working_exposure_ratio = days_in_month / working_days,
      offset_log_difference = offset_log_patient_working_days -
        offset_log_patient_days
    )

  output_files <- c(
    "working_day_screen_class.csv", "working_day_screen_drug.csv",
    "working_day_characterisation_class.csv",
    "working_day_characterisation_drug.csv", "working_day_summary.csv",
    "working_day_overlap_summary.csv", "working_day_classification_changes.csv",
    "working_day_q_boundary.csv", "working_day_route_summary.csv",
    "working_day_offset_monthly.csv"
  )
  atomic_fwrite(as.data.table(working_screen_class),
                file.path(working_dir, output_files[1]))
  atomic_fwrite(as.data.table(working_screen_drug),
                file.path(working_dir, output_files[2]))
  atomic_fwrite(as.data.table(working_char_class),
                file.path(working_dir, output_files[3]))
  atomic_fwrite(as.data.table(working_char_drug),
                file.path(working_dir, output_files[4]))
  atomic_fwrite(as.data.table(working_summary),
                file.path(working_dir, output_files[5]))
  atomic_fwrite(as.data.table(working_overlap_summary),
                file.path(working_dir, output_files[6]))
  atomic_fwrite(as.data.table(working_classification_changes),
                file.path(working_dir, output_files[7]))
  atomic_fwrite(as.data.table(working_q_boundary),
                file.path(working_dir, output_files[8]))
  atomic_fwrite(as.data.table(working_route_summary),
                file.path(working_dir, output_files[9]))
  atomic_fwrite(as.data.table(working_offset_monthly),
                file.path(working_dir, output_files[10]))

  add_working_check(
    "complete_families",
    all(working_summary[level == "class", family_size] == 220L) &&
      all(working_summary[level == "drug", family_size] == 974L),
    "two 220-class and two 974-drug families",
    paste(working_summary$level, working_summary$analysis_specification,
          working_summary$family_size, sep = "=", collapse = ";")
  )
  add_working_check(
    "all_models_complete", all(working_summary$model_failures == 0L),
    "0 failed models", sum(working_summary$model_failures)
  )
  add_working_check(
    "working_offset_valid",
    nrow(working_offset_monthly) == 48L &&
      all(working_offset_monthly$working_days >= 18L &
            working_offset_monthly$working_days <= 23L) &&
      all(is.finite(working_offset_monthly$offset_log_difference)) &&
      all(working_offset_monthly$offset_log_difference != 0),
    "48 months; 18-23 working days; finite offset differs from calendar offset",
    sprintf("%d months; range %d-%d",
            nrow(working_offset_monthly), min(working_offset_monthly$working_days),
            max(working_offset_monthly$working_days))
  )
  add_working_check(
    "valid_q_values",
    all(is.finite(working_class$seasonality_q_bh)) &&
      all(working_class$seasonality_q_bh >= 0 & working_class$seasonality_q_bh <= 1) &&
      all(is.finite(working_drug$seasonality_q_bh)) &&
      all(working_drug$seasonality_q_bh >= 0 & working_drug$seasonality_q_bh <= 1),
    "all working-day q-values finite and in [0,1]", "checked"
  )
  add_working_check(
    "all_discoveries_characterised",
    all(working_summary$discoveries == working_summary$characterised) &&
      !anyNA(working_char_class$ptr_lci) && !anyNA(working_char_drug$ptr_lci),
    "every discovery has amplitude CI and characterisation",
    paste(working_summary$level, working_summary$analysis_specification,
          paste0(working_summary$characterised, "/", working_summary$discoveries),
          sep = "=", collapse = ";")
  )
  primary_working_summary <- working_summary[
    analysis_specification == "primary_calendar_day_offset"
  ]
  add_working_check(
    "primary_counts_preserved",
    primary_working_summary[level == "class", discoveries] == 125L &&
      primary_working_summary[level == "drug", discoveries] == 391L &&
      primary_working_summary[level == "class", meaningful] == 30L &&
      primary_working_summary[level == "drug", meaningful] == 88L,
    "primary 125/391 discoveries and 30/88 meaningful",
    paste(primary_working_summary$level, primary_working_summary$discoveries,
          primary_working_summary$meaningful, sep = "=", collapse = ";")
  )
  add_working_check(
    "set_comparisons_reconcile",
    all(working_overlap_summary$discovery_overlap +
          working_overlap_summary$discoveries_added ==
          working_overlap_summary$alternative_discoveries) &&
      all(working_overlap_summary$discovery_overlap +
            working_overlap_summary$discoveries_lost ==
            working_overlap_summary$primary_discoveries) &&
      all(working_overlap_summary$meaningful_overlap +
            working_overlap_summary$meaningful_added ==
            working_overlap_summary$alternative_meaningful) &&
      all(working_overlap_summary$meaningful_overlap +
            working_overlap_summary$meaningful_lost ==
            working_overlap_summary$primary_meaningful),
    "all overlap/add/loss counts reconcile to both compared sets", "checked"
  )
  add_working_check(
    "comparison_outputs_present",
    all(file.exists(file.path(working_dir, output_files))) &&
      nrow(working_overlap_summary) == 2L && nrow(working_q_boundary) == 80L,
    "10 result files; 2 overlap rows; 80 boundary rows",
    sprintf("%d files; %d overlap rows; %d boundary rows",
            sum(file.exists(file.path(working_dir, output_files))),
            nrow(working_overlap_summary), nrow(working_q_boundary))
  )

  working_qc_summary <- rbindlist(working_checks, use.names = TRUE, fill = TRUE)
  atomic_fwrite(working_qc_summary,
                file.path(working_dir, "stage5_working_days_qc_summary.csv"))
  if (any(!working_qc_summary$pass)) {
    stop("Stage 5.3 completion gate failed: ",
         paste(working_qc_summary[pass == FALSE, check_id], collapse = ", "),
         ". See ",
         file.path(working_dir, "stage5_working_days_qc_summary.csv"), ".")
  }

  for (relative_path in output_files) {
    source_path <- file.path(working_dir, relative_path)
    target_path <- file.path(working_snapshot_dir, relative_path)
    temporary_path <- paste0(target_path, ".tmp")
    if (!file.copy(source_path, temporary_path, overwrite = TRUE) ||
        !file.rename(temporary_path, target_path)) {
      stop("Could not seal Stage 5.3 snapshot file: ", relative_path)
    }
  }
  working_manifest <- data.table(
    relative_path = output_files,
    bytes = as.numeric(file.info(file.path(working_snapshot_dir, output_files))$size),
    sha256 = vapply(file.path(working_snapshot_dir, output_files),
                    sha256_file, character(1))
  )
  atomic_fwrite(
    working_manifest,
    file.path(working_dir, "stage5_working_days_snapshot_manifest.csv")
  )
  working_completion <- data.table(
    stage = "stage5_working_days", status = "PASS",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    checks_passed = sum(working_qc_summary$pass),
    checks_total = nrow(working_qc_summary),
    classification_change_rows = nrow(working_classification_changes),
    snapshot_files = nrow(working_manifest),
    analysis_script_sha256 = sha256_file(script_path),
    renv_lock_sha256 = sha256_file(lock_path)
  )
  atomic_fwrite(
    working_completion,
    file.path(working_dir, "stage5_working_days_completion.csv")
  )
  message(
    "Stage 5.3 complete: working-day offset fitted for all class/drug series; ",
    nrow(working_qc_summary), " checks passed; stopping for result review."
  )
  return(invisible(TRUE))
}

## ---- STAGE 5.4: STL-strength threshold sensitivity -----------------------
if (run_stage == "stage5_threshold") {
  threshold_dir <- file.path(out_dir, "qc", "stage5", "threshold")
  threshold_snapshot_dir <- file.path(threshold_dir, "threshold_snapshot")
  dir.create(threshold_snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  threshold_checks <- list()
  add_threshold_check <- function(check_id, pass, expected, observed, details = "") {
    threshold_checks[[length(threshold_checks) + 1L]] <<- data.table(
      check_id = check_id, pass = isTRUE(pass), expected = as.character(expected),
      observed = as.character(observed), details = as.character(details)
    )
  }
  verify_threshold_prior <- function(manifest_path, snapshot_path) {
    if (!file.exists(manifest_path)) {
      return(list(ok = FALSE, n = 0L, matching = 0L))
    }
    manifest <- fread(manifest_path)
    paths <- file.path(snapshot_path, manifest$relative_path)
    if (!all(file.exists(paths))) {
      return(list(ok = FALSE, n = nrow(manifest), matching = sum(file.exists(paths))))
    }
    hashes <- vapply(paths, sha256_file, character(1))
    size_col <- if ("bytes" %in% names(manifest)) "bytes" else "size_bytes"
    matching <- hashes == manifest$sha256 &
      as.numeric(file.info(paths)$size) == manifest[[size_col]]
    list(ok = all(matching), n = nrow(manifest), matching = sum(matching))
  }

  prior_specs <- list(
    stage4 = list(
      file.path(out_dir, "qc", "stage4", "stage4_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage4", "stage4_snapshot")
    ),
    trend = list(
      file.path(out_dir, "qc", "stage5", "trend", "stage5_trend_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage5", "trend", "trend_snapshot")
    ),
    hac = list(
      file.path(out_dir, "qc", "stage5", "uniform_hac", "stage5_hac_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage5", "uniform_hac", "uniform_hac_snapshot")
    ),
    working_days = list(
      file.path(out_dir, "qc", "stage5", "working_days",
                "stage5_working_days_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage5", "working_days", "working_days_snapshot")
    )
  )
  prior_verification <- lapply(
    prior_specs, function(x) verify_threshold_prior(x[[1]], x[[2]])
  )
  prior_ok <- all(vapply(prior_verification, function(x) x$ok, logical(1)))
  add_threshold_check(
    "prior_authorities_unchanged", prior_ok,
    "17 Stage 4, 8 trend, 9 uniform-HAC and 10 working-day files unchanged",
    paste(names(prior_verification), vapply(prior_verification, function(x) {
      sprintf("%d/%d", x$matching, x$n)
    }, character(1)), collapse = "; ")
  )
  if (!prior_ok) {
    stop("Stage 5.4 requires intact Stage 4 and Stage 5.1--5.3 snapshots.")
  }

  # Stage 5.4 is deliberately post-processing only. It reads the sealed Stage 4
  # characterisations and changes the descriptive STL-strength threshold; the
  # discovery families, amplitude rule and primary 0.50 threshold stay fixed.
  stage4_snapshot_dir <- prior_specs$stage4[[2]]
  stage4_char_class <- fread(
    file.path(stage4_snapshot_dir, "characterisation_class.csv"),
    colClasses = list(character = "bnf_class_code")
  )
  stage4_char_drug <- fread(
    file.path(stage4_snapshot_dir, "characterisation_drug.csv"),
    colClasses = list(character = c("bnf_drug_code", "bnf_class_code"))
  )
  stage4_results_class <- fread(
    file.path(stage4_snapshot_dir, "results_class.csv"),
    colClasses = list(character = "bnf_class_code")
  )
  stage4_results_drug <- fread(
    file.path(stage4_snapshot_dir, "results_drug.csv"),
    colClasses = list(character = "bnf_drug_code")
  )

  threshold_grid <- c(0.40, 0.50, 0.60)
  class_source <- stage4_char_class[, .(
    series_code = bnf_class_code, series_name = bnf_class_name,
    parent_class_significant = NA,
    peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
    stl_seasonal_strength, stl_trend_strength,
    seasonality_q_bh = class_q_bh,
    seasonality_detected = class_significant
  )]
  drug_source <- stage4_char_drug[, .(
    series_code = bnf_drug_code, series_name = bnf_drug_name,
    parent_class_significant,
    peak_trough_ratio, ptr_lci, ptr_uci, peak_month, trough_month,
    stl_seasonal_strength, stl_trend_strength,
    seasonality_q_bh = drug_all_q_bh,
    seasonality_detected = drug_significant
  )]

  expand_thresholds <- function(source, level_name) {
    rbindlist(lapply(threshold_grid, function(cut) {
      data.table(
        level = level_name,
        series_code = source$series_code,
        series_name = source$series_name,
        inference_scope = if (level_name == "class") {
          "primary_inferential"
        } else {
          "secondary_exploratory"
        },
        classification_role = "descriptive_threshold_sensitivity",
        amplitude_lci_threshold = meaningful_threshold,
        stl_strength_threshold = cut,
        primary_stl_strength_threshold = stl_strength_threshold,
        peak_trough_ratio = source$peak_trough_ratio,
        ptr_lci = source$ptr_lci,
        ptr_uci = source$ptr_uci,
        peak_month = source$peak_month,
        trough_month = source$trough_month,
        stl_seasonal_strength = source$stl_seasonal_strength,
        stl_trend_strength = source$stl_trend_strength,
        seasonality_q_bh = source$seasonality_q_bh,
        seasonality_detected = source$seasonality_detected,
        parent_class_significant = as.logical(source$parent_class_significant),
        amplitude_qualified = source$ptr_lci >= meaningful_threshold,
        meaningful = source$ptr_lci >= meaningful_threshold &
          source$stl_seasonal_strength >= cut
      )
    }), use.names = TRUE)
  }

  threshold_class <- expand_thresholds(class_source, "class")
  threshold_drug <- expand_thresholds(drug_source, "drug")

  threshold_class_summary <- threshold_class[, .(
    n_discoveries = uniqueN(series_code),
    n_amplitude_qualified = sum(amplitude_qualified),
    n_meaningful = sum(meaningful),
    n_meaningful_in_significant_parent = NA_integer_,
    n_meaningful_outside_significant_parent = NA_integer_
  ), by = .(level, stl_strength_threshold)]
  threshold_drug_summary <- threshold_drug[, .(
    n_discoveries = uniqueN(series_code),
    n_amplitude_qualified = sum(amplitude_qualified),
    n_meaningful = sum(meaningful),
    n_meaningful_in_significant_parent = sum(meaningful & parent_class_significant),
    n_meaningful_outside_significant_parent = sum(meaningful & !parent_class_significant)
  ), by = .(level, stl_strength_threshold)]
  threshold_summary <- rbindlist(
    list(threshold_class_summary, threshold_drug_summary), use.names = TRUE
  )
  setorder(threshold_summary, level, stl_strength_threshold)

  make_threshold_changes <- function(panel) {
    reference <- panel[abs(stl_strength_threshold - 0.50) < 1e-12, .(
      series_code, series_name, parent_class_significant,
      ptr_lci, stl_seasonal_strength,
      reference_meaningful = meaningful
    )]
    rbindlist(lapply(c(0.40, 0.60), function(cut) {
      alternative <- panel[abs(stl_strength_threshold - cut) < 1e-12, .(
        series_code, alternative_meaningful = meaningful
      )]
      comparison <- merge(reference, alternative, by = "series_code", all = TRUE)
      comparison[reference_meaningful != alternative_meaningful, .(
        level = unique(panel$level), series_code, series_name,
        comparison = sprintf("%.2f_vs_0.50", cut),
        reference_stl_threshold = 0.50,
        alternative_stl_threshold = cut,
        reference_meaningful, alternative_meaningful,
        direction = ifelse(alternative_meaningful, "added", "lost"),
        ptr_lci, stl_seasonal_strength, parent_class_significant
      )]
    }), use.names = TRUE)
  }
  threshold_classification_changes <- rbindlist(
    list(make_threshold_changes(threshold_class),
         make_threshold_changes(threshold_drug)),
    use.names = TRUE
  )
  setorder(threshold_classification_changes, level, alternative_stl_threshold,
           direction, series_code)

  make_threshold_boundary <- function(panel, n = 20L) {
    boundary <- panel[
      abs(stl_strength_threshold - 0.50) < 1e-12 & amplitude_qualified
    ][order(abs(stl_seasonal_strength - 0.50), series_code)]
    boundary <- head(boundary, n)
    boundary[, .(
      level, series_code, series_name, parent_class_significant,
      peak_trough_ratio, ptr_lci, ptr_uci,
      stl_seasonal_strength,
      distance_from_primary_threshold = abs(stl_seasonal_strength - 0.50),
      meaningful_at_0_40 = ptr_lci >= meaningful_threshold &
        stl_seasonal_strength >= 0.40,
      meaningful_at_0_50 = ptr_lci >= meaningful_threshold &
        stl_seasonal_strength >= 0.50,
      meaningful_at_0_60 = ptr_lci >= meaningful_threshold &
        stl_seasonal_strength >= 0.60
    )]
  }
  threshold_boundary <- rbindlist(
    list(make_threshold_boundary(threshold_class),
         make_threshold_boundary(threshold_drug)),
    use.names = TRUE
  )

  add_threshold_check(
    "sealed_source_scope",
    nrow(stage4_char_class) == 125L && nrow(stage4_char_drug) == 391L &&
      !anyDuplicated(stage4_char_class$bnf_class_code) &&
      !anyDuplicated(stage4_char_drug$bnf_drug_code) &&
      all(stage4_char_class$class_significant) &&
      all(stage4_char_drug$drug_significant),
    "125 class and 391 drug discoveries from the sealed Stage 4 authority",
    sprintf("%d classes; %d drugs", nrow(stage4_char_class), nrow(stage4_char_drug))
  )
  add_threshold_check(
    "current_primary_alignment",
    setequal(as.character(results_class$bnf_class_code),
             stage4_results_class$bnf_class_code) &&
      setequal(as.character(results_drug$bnf_drug_code),
               stage4_results_drug$bnf_drug_code) &&
      setequal(as.character(results_class$bnf_class_code[results_class$meaningful]),
               stage4_results_class$bnf_class_code[stage4_results_class$meaningful]) &&
      setequal(as.character(results_drug$bnf_drug_code[results_drug$meaningful]),
               stage4_results_drug$bnf_drug_code[stage4_results_drug$meaningful]),
    "current primary discovery and meaningful sets match sealed Stage 4",
    sprintf("%d/%d current discoveries; %d/%d current meaningful",
            nrow(results_class), nrow(results_drug),
            sum(results_class$meaningful), sum(results_drug$meaningful))
  )
  add_threshold_check(
    "threshold_grid_complete",
    nrow(threshold_class) == 125L * length(threshold_grid) &&
      nrow(threshold_drug) == 391L * length(threshold_grid) &&
      setequal(unique(threshold_class$stl_strength_threshold), threshold_grid) &&
      setequal(unique(threshold_drug$stl_strength_threshold), threshold_grid),
    "0.40, 0.50 and 0.60 for every sealed discovery",
    sprintf("%d class rows; %d drug rows", nrow(threshold_class), nrow(threshold_drug))
  )
  add_threshold_check(
    "continuous_values_complete",
    all(is.finite(threshold_class$ptr_lci)) &&
      all(is.finite(threshold_drug$ptr_lci)) &&
      all(threshold_class$stl_seasonal_strength >= 0 &
            threshold_class$stl_seasonal_strength <= 1) &&
      all(threshold_drug$stl_seasonal_strength >= 0 &
            threshold_drug$stl_seasonal_strength <= 1) &&
      all(threshold_class$seasonality_detected) &&
      all(threshold_drug$seasonality_detected),
    "finite amplitude limits, seasonal strengths in [0,1], discoveries only",
    "checked"
  )

  primary_class_codes <- threshold_class[
    abs(stl_strength_threshold - 0.50) < 1e-12 & meaningful, series_code
  ]
  primary_drug_codes <- threshold_drug[
    abs(stl_strength_threshold - 0.50) < 1e-12 & meaningful, series_code
  ]
  add_threshold_check(
    "primary_threshold_reproduced",
    length(primary_class_codes) == 30L && length(primary_drug_codes) == 88L &&
      setequal(primary_class_codes,
               stage4_results_class$bnf_class_code[stage4_results_class$meaningful]) &&
      setequal(primary_drug_codes,
               stage4_results_drug$bnf_drug_code[stage4_results_drug$meaningful]),
    "0.50 exactly reproduces 30 primary classes and 88 exploratory drugs",
    sprintf("%d classes; %d drugs", length(primary_class_codes), length(primary_drug_codes))
  )
  nested_sets <- function(panel) {
    s40 <- panel[abs(stl_strength_threshold - 0.40) < 1e-12 & meaningful, series_code]
    s50 <- panel[abs(stl_strength_threshold - 0.50) < 1e-12 & meaningful, series_code]
    s60 <- panel[abs(stl_strength_threshold - 0.60) < 1e-12 & meaningful, series_code]
    all(s60 %in% s50) && all(s50 %in% s40)
  }
  add_threshold_check(
    "meaningful_sets_nested",
    nested_sets(threshold_class) && nested_sets(threshold_drug),
    "0.60 set is within 0.50, which is within 0.40, at both levels",
    "checked"
  )

  detail_counts <- rbindlist(list(
    threshold_class[, .(
      n_discoveries = uniqueN(series_code),
      n_amplitude_qualified = sum(amplitude_qualified),
      n_meaningful = sum(meaningful)
    ), by = .(level, stl_strength_threshold)],
    threshold_drug[, .(
      n_discoveries = uniqueN(series_code),
      n_amplitude_qualified = sum(amplitude_qualified),
      n_meaningful = sum(meaningful)
    ), by = .(level, stl_strength_threshold)]
  ))
  summary_comparison <- merge(
    threshold_summary[, .(level, stl_strength_threshold,
                          n_discoveries, n_amplitude_qualified, n_meaningful)],
    detail_counts,
    by = c("level", "stl_strength_threshold"), suffixes = c("_summary", "_detail")
  )
  add_threshold_check(
    "summary_reconciles",
    nrow(summary_comparison) == 6L &&
      all(summary_comparison$n_discoveries_summary ==
            summary_comparison$n_discoveries_detail) &&
      all(summary_comparison$n_amplitude_qualified_summary ==
            summary_comparison$n_amplitude_qualified_detail) &&
      all(summary_comparison$n_meaningful_summary ==
            summary_comparison$n_meaningful_detail),
    "six summary rows reconcile to detailed flags",
    paste(threshold_summary$level,
          sprintf("%.2f=%d", threshold_summary$stl_strength_threshold,
                  threshold_summary$n_meaningful), collapse = "; ")
  )

  expected_change_count <- function(panel) {
    reference <- panel[abs(stl_strength_threshold - 0.50) < 1e-12,
                       .(series_code, reference = meaningful)]
    sum(vapply(c(0.40, 0.60), function(cut) {
      alternative <- panel[abs(stl_strength_threshold - cut) < 1e-12,
                           .(series_code, alternative = meaningful)]
      comparison <- merge(reference, alternative, by = "series_code")
      sum(comparison$reference != comparison$alternative)
    }, integer(1)))
  }
  expected_changes <- expected_change_count(threshold_class) +
    expected_change_count(threshold_drug)
  add_threshold_check(
    "classification_changes_reconcile",
    nrow(threshold_classification_changes) == expected_changes &&
      !anyDuplicated(threshold_classification_changes[,
        .(level, comparison, series_code)]) &&
      all((threshold_classification_changes$direction == "added") ==
            threshold_classification_changes$alternative_meaningful),
    "one unique row for every identity change relative to 0.50",
    sprintf("%d rows", nrow(threshold_classification_changes))
  )
  add_threshold_check(
    "boundary_rows_complete",
    nrow(threshold_boundary[level == "class"]) == 20L &&
      nrow(threshold_boundary[level == "drug"]) == 20L &&
      all(threshold_boundary$ptr_lci >= meaningful_threshold),
    "20 amplitude-qualified boundary rows per level",
    sprintf("%d classes; %d drugs",
            nrow(threshold_boundary[level == "class"]),
            nrow(threshold_boundary[level == "drug"]))
  )

  output_files <- c(
    "threshold_characterisation_class.csv",
    "threshold_characterisation_drug.csv",
    "threshold_summary.csv",
    "threshold_classification_changes.csv",
    "threshold_boundary.csv"
  )
  atomic_fwrite(threshold_class, file.path(threshold_dir, output_files[1]))
  atomic_fwrite(threshold_drug, file.path(threshold_dir, output_files[2]))
  atomic_fwrite(threshold_summary, file.path(threshold_dir, output_files[3]))
  atomic_fwrite(threshold_classification_changes,
                file.path(threshold_dir, output_files[4]))
  atomic_fwrite(threshold_boundary, file.path(threshold_dir, output_files[5]))

  add_threshold_check(
    "outputs_written",
    all(file.exists(file.path(threshold_dir, output_files))) &&
      all(file.info(file.path(threshold_dir, output_files))$size > 0),
    paste(length(output_files), "non-empty analytical files"),
    paste(sum(file.exists(file.path(threshold_dir, output_files))), "present")
  )
  threshold_qc_summary <- rbindlist(threshold_checks)
  atomic_fwrite(
    threshold_qc_summary,
    file.path(threshold_dir, "stage5_threshold_qc_summary.csv")
  )
  if (!all(threshold_qc_summary$pass)) {
    failed <- threshold_qc_summary[!pass, check_id]
    stop("Stage 5.4 completion gate failed: ", paste(failed, collapse = ", "),
         ". See ", file.path(threshold_dir, "stage5_threshold_qc_summary.csv"), ".")
  }

  for (relative_path in output_files) {
    source_path <- file.path(threshold_dir, relative_path)
    target_path <- file.path(threshold_snapshot_dir, relative_path)
    temporary_path <- paste0(target_path, ".tmp")
    if (!file.copy(source_path, temporary_path, overwrite = TRUE) ||
        !file.rename(temporary_path, target_path)) {
      stop("Could not seal Stage 5.4 snapshot file: ", relative_path)
    }
  }
  threshold_manifest <- data.table(
    relative_path = output_files,
    bytes = as.numeric(file.info(file.path(threshold_snapshot_dir, output_files))$size),
    sha256 = vapply(file.path(threshold_snapshot_dir, output_files),
                    sha256_file, character(1))
  )
  atomic_fwrite(
    threshold_manifest,
    file.path(threshold_dir, "stage5_threshold_snapshot_manifest.csv")
  )
  threshold_completion <- data.table(
    stage = "stage5_threshold", status = "PASS",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    checks_passed = sum(threshold_qc_summary$pass),
    checks_total = nrow(threshold_qc_summary),
    classification_change_rows = nrow(threshold_classification_changes),
    snapshot_files = nrow(threshold_manifest),
    analysis_script_sha256 = sha256_file(script_path),
    renv_lock_sha256 = sha256_file(lock_path)
  )
  atomic_fwrite(
    threshold_completion,
    file.path(threshold_dir, "stage5_threshold_completion.csv")
  )
  message(
    "Stage 5.4 complete: STL-strength thresholds summarised from the sealed ",
    "Stage 4 authority; ", nrow(threshold_qc_summary),
    " checks passed; stopping before Stage 6 diagnostics."
  )
  return(invisible(TRUE))
}

## ---- STAGE 6: targeted diagnostics, cohort flow and missingness -----------
if (run_stage == "stage6") {
  stage6_dir <- file.path(out_dir, "qc", "stage6")
  stage6_snapshot_dir <- file.path(stage6_dir, "stage6_snapshot")
  dir.create(stage6_snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  stage6_checks <- list()
  add_stage6_check <- function(check_id, pass, expected, observed, details = "") {
    stage6_checks[[length(stage6_checks) + 1L]] <<- data.table(
      check_id = check_id, pass = isTRUE(pass), expected = as.character(expected),
      observed = as.character(observed), details = as.character(details)
    )
  }
  verify_stage6_prior <- function(manifest_path, snapshot_path) {
    if (!file.exists(manifest_path)) {
      return(list(ok = FALSE, n = 0L, matching = 0L))
    }
    manifest <- fread(manifest_path)
    paths <- file.path(snapshot_path, manifest$relative_path)
    if (!all(file.exists(paths))) {
      return(list(ok = FALSE, n = nrow(manifest), matching = sum(file.exists(paths))))
    }
    hashes <- vapply(paths, sha256_file, character(1))
    size_col <- if ("bytes" %in% names(manifest)) "bytes" else "size_bytes"
    matching <- hashes == manifest$sha256 &
      as.numeric(file.info(paths)$size) == manifest[[size_col]]
    list(ok = all(matching), n = nrow(manifest), matching = sum(matching))
  }
  prior_specs <- list(
    stage4 = list(
      file.path(out_dir, "qc", "stage4", "stage4_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage4", "stage4_snapshot")
    ),
    trend = list(
      file.path(out_dir, "qc", "stage5", "trend", "stage5_trend_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage5", "trend", "trend_snapshot")
    ),
    hac = list(
      file.path(out_dir, "qc", "stage5", "uniform_hac", "stage5_hac_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage5", "uniform_hac", "uniform_hac_snapshot")
    ),
    working_days = list(
      file.path(out_dir, "qc", "stage5", "working_days",
                "stage5_working_days_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage5", "working_days", "working_days_snapshot")
    ),
    threshold = list(
      file.path(out_dir, "qc", "stage5", "threshold",
                "stage5_threshold_snapshot_manifest.csv"),
      file.path(out_dir, "qc", "stage5", "threshold", "threshold_snapshot")
    )
  )
  prior_verification <- lapply(
    prior_specs, function(x) verify_stage6_prior(x[[1]], x[[2]])
  )
  prior_ok <- all(vapply(prior_verification, function(x) x$ok, logical(1)))
  add_stage6_check(
    "prior_authorities_unchanged", prior_ok,
    "17 Stage 4 and 8/9/10/5 Stage 5 files unchanged",
    paste(names(prior_verification), vapply(prior_verification, function(x) {
      sprintf("%d/%d", x$matching, x$n)
    }, character(1)), collapse = "; ")
  )
  if (!prior_ok) stop("Stage 6 requires intact Stage 4 and Stage 5 snapshots.")

  stage4_results_class <- fread(
    file.path(prior_specs$stage4[[2]], "results_class.csv"),
    colClasses = list(character = "bnf_class_code")
  )
  stage4_results_drug <- fread(
    file.path(prior_specs$stage4[[2]], "results_drug.csv"),
    colClasses = list(character = "bnf_drug_code")
  )
  primary_alignment <-
    setequal(as.character(results_class$bnf_class_code),
             stage4_results_class$bnf_class_code) &&
    setequal(as.character(results_drug$bnf_drug_code),
             stage4_results_drug$bnf_drug_code) &&
    setequal(as.character(results_class$bnf_class_code[results_class$meaningful]),
             stage4_results_class$bnf_class_code[stage4_results_class$meaningful]) &&
    setequal(as.character(results_drug$bnf_drug_code[results_drug$meaningful]),
             stage4_results_drug$bnf_drug_code[stage4_results_drug$meaningful])
  add_stage6_check(
    "primary_authority_alignment", primary_alignment,
    "current primary discovery and meaningful sets match sealed Stage 4",
    sprintf("%d/%d discoveries; %d/%d meaningful",
            nrow(results_class), nrow(results_drug),
            sum(results_class$meaningful), sum(results_drug$meaningful))
  )
  if (!primary_alignment) stop("Current primary results do not match Stage 4 authority.")

  fit_primary_mean <- function(series, screen_row) {
    tryCatch({
      d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
      d <- d[order(d$t), ]
      if (nrow(d) != nrow(covar) || anyNA(d$items)) {
        stop("series does not form one complete 48-month model input")
      }
      d$off <- d$offset_log_patient_days
      distribution <- as.character(screen_row$distribution[[1]])
      fit <- if (identical(distribution, "negbin")) .nb_fit(.f_full, d) else NULL
      if (identical(distribution, "negbin") && is.null(fit)) {
        stop("negative-binomial diagnostic refit failed")
      }
      if (is.null(fit)) fit <- glm(.f_full, family = poisson, data = d)
      if (!isTRUE(fit$converged)) stop("diagnostic mean-model refit did not converge")
      list(ok = TRUE, data = d, fit = fit, note = NA_character_)
    }, error = function(e) {
      list(ok = FALSE, data = NULL, fit = NULL, note = conditionMessage(e))
    })
  }

  diagnose_primary_series <- function(series, screen_row) {
    fitted_object <- fit_primary_mean(series, screen_row)
    if (!fitted_object$ok) {
      return(data.table(
        diagnostic_refit_converged = FALSE,
        diagnostic_refit_note = fitted_object$note,
        n_model_months = NA_integer_, pearson_residual_rmse = NA_real_,
        max_abs_pearson_residual = NA_real_, max_abs_pearson_month = NA_integer_,
        max_abs_residual_acf = NA_real_, max_abs_residual_acf_lag = NA_integer_,
        max_cooks_distance = NA_real_, max_cooks_month = NA_integer_,
        max_segment_mean_shift_sd = NA_real_, segment_shift_after_month = NA_integer_,
        observed_fitted_log_rate_r2 = NA_real_, fitted_to_observed_total_ratio = NA_real_,
        harmonic_peak_trough_ratio = NA_real_,
        diagnostic_stl_seasonal_strength = NA_real_,
        diagnostic_stl_trend_strength = NA_real_
      ))
    }
    d <- fitted_object$data
    fit <- fitted_object$fit
    pearson <- as.numeric(residuals(fit, type = "pearson"))
    fitted_count <- as.numeric(fitted(fit))
    abs_pearson <- abs(pearson)
    acf_values <- as.numeric(stats::acf(
      pearson, lag.max = harmonic_period_months[1], plot = FALSE,
      demean = TRUE, na.action = na.pass
    )$acf)[-1]
    cooks <- suppressWarnings(as.numeric(cooks.distance(fit)))
    candidate_splits <- 12L:(length(pearson) - 12L)
    residual_sd <- sd(pearson)
    split_shift <- if (is.finite(residual_sd) && residual_sd > 0) {
      vapply(candidate_splits, function(split) {
        abs(mean(pearson[seq_len(split)]) -
              mean(pearson[(split + 1L):length(pearson)])) / residual_sd
      }, numeric(1))
    } else {
      rep(0, length(candidate_splits))
    }
    observed_rate <- d$items / d$list_size
    fitted_rate <- fitted_count / d$list_size
    rate_cor <- suppressWarnings(cor(log1p(observed_rate), log1p(fitted_rate)))
    seasonal_curve <- .seasonal_curve(
      screen_row$b_sin12[[1]], screen_row$b_cos12[[1]],
      screen_row$b_sin6[[1]], screen_row$b_cos6[[1]]
    )
    stl_values <- .stl_strength(series, covar)
    data.table(
      diagnostic_refit_converged = TRUE,
      diagnostic_refit_note = NA_character_,
      n_model_months = nrow(d),
      pearson_residual_rmse = sqrt(mean(pearson^2)),
      max_abs_pearson_residual = max(abs_pearson),
      max_abs_pearson_month = d$year_month[which.max(abs_pearson)],
      max_abs_residual_acf = max(abs(acf_values)),
      max_abs_residual_acf_lag = which.max(abs(acf_values)),
      max_cooks_distance = max(cooks, na.rm = TRUE),
      max_cooks_month = d$year_month[which.max(cooks)],
      max_segment_mean_shift_sd = max(split_shift),
      segment_shift_after_month = d$year_month[candidate_splits[which.max(split_shift)]],
      observed_fitted_log_rate_r2 = if (is.finite(rate_cor)) rate_cor^2 else NA_real_,
      fitted_to_observed_total_ratio = sum(fitted_count) / sum(d$items),
      harmonic_peak_trough_ratio = exp(diff(range(seasonal_curve))),
      diagnostic_stl_seasonal_strength = unname(stl_values["seasonal"]),
      diagnostic_stl_trend_strength = unname(stl_values["trend"])
    )
  }

  diagnose_level <- function(monthly, screen, id_col, name_col, level_name) {
    rbindlist(lapply(seq_len(nrow(screen)), function(i) {
      screen_row <- screen[i, , drop = FALSE]
      code <- as.character(screen_row[[id_col]][[1]])
      series <- monthly[as.character(monthly[[id_col]]) == code,
                        c("year_month", "items")]
      metrics <- diagnose_primary_series(series, screen_row)
      data.table(
        level = level_name,
        series_code = code,
        series_name = as.character(screen_row[[name_col]][[1]]),
        inference_scope = if (level_name == "class") {
          "primary_inferential"
        } else {
          "secondary_exploratory"
        },
        seasonality_q_bh = if (level_name == "class") {
          as.numeric(screen_row$class_q_bh[[1]])
        } else {
          as.numeric(screen_row$drug_all_q_bh[[1]])
        },
        seasonality_detected = if (level_name == "class") {
          as.logical(screen_row$class_significant[[1]])
        } else {
          as.logical(screen_row$drug_significant[[1]])
        },
        parent_class_significant = if (level_name == "drug") {
          as.logical(screen_row$parent_class_significant[[1]])
        } else {
          NA
        },
        distribution = as.character(screen_row$distribution[[1]]),
        route = as.character(screen_row$route[[1]]),
        disp_ratio = as.numeric(screen_row$disp_ratio[[1]]),
        disp_p = as.numeric(screen_row$disp_p[[1]]),
        lb_p = as.numeric(screen_row$lb_p[[1]]),
        nw_lag = as.numeric(screen_row$nw_lag[[1]]),
        hac_capped = as.logical(screen_row$hac_capped[[1]]),
        theta = as.numeric(screen_row$theta[[1]]),
        original_model_converged = as.logical(screen_row$converged[[1]]),
        original_model_note = as.character(screen_row$note[[1]])
      )[, cbind(.SD, metrics)]
    }), use.names = TRUE)
  }

  diagnostic_class <- diagnose_level(
    class_monthly_elig, screen_class,
    "bnf_class_code", "bnf_class_name", "class"
  )
  diagnostic_drug <- diagnose_level(
    drug_monthly_elig, screen_drug,
    "bnf_drug_code", "bnf_drug_name", "drug"
  )
  diagnostic_inventory <- rbindlist(
    list(diagnostic_class, diagnostic_drug), use.names = TRUE
  )
  setorder(diagnostic_inventory, level, seasonality_q_bh, series_code)

  add_stage6_check(
    "diagnostic_inventory_complete",
    nrow(diagnostic_class) == 220L && nrow(diagnostic_drug) == 974L &&
      !anyDuplicated(diagnostic_inventory[, .(level, series_code)]) &&
      all(diagnostic_inventory$n_model_months == 48L) &&
      all(diagnostic_inventory$diagnostic_refit_converged),
    "220 class and 974 drug diagnostic refits, each using 48 months",
    sprintf("%d classes; %d drugs; %d failed refits",
            nrow(diagnostic_class), nrow(diagnostic_drug),
            sum(!diagnostic_inventory$diagnostic_refit_converged))
  )
  diagnostic_values_finite <- diagnostic_inventory[
    diagnostic_refit_converged == TRUE,
    all(is.finite(pearson_residual_rmse)) &&
      all(is.finite(max_abs_pearson_residual)) &&
      all(is.finite(max_abs_residual_acf)) &&
      all(is.finite(max_cooks_distance)) &&
      all(is.finite(max_segment_mean_shift_sd)) &&
      all(is.finite(harmonic_peak_trough_ratio)) &&
      all(is.finite(diagnostic_stl_seasonal_strength)) &&
      all(is.finite(diagnostic_stl_trend_strength))
  ]
  add_stage6_check(
    "diagnostic_metrics_valid", diagnostic_values_finite,
    "finite residual, influence, ACF, amplitude and STL metrics for every refit",
    "checked"
  )

  # Selection is explicit and reproducible. Automated metrics are retained for
  # all series; detailed month-level residuals and ACFs are limited to the
  # prespecified scientific and diagnostic review groups below.
  selection_parts <- list()
  add_selection <- function(inventory, codes, reason) {
    rows <- inventory[series_code %in% as.character(codes),
                      .(level, series_code, series_name)]
    if (nrow(rows)) {
      rows[, selection_reason := reason]
      selection_parts[[length(selection_parts) + 1L]] <<- rows
    }
  }
  top_metric_codes <- function(inventory, metric, n = 5L) {
    values <- inventory[is.finite(get(metric))][order(-get(metric))]
    head(values$series_code, n)
  }

  add_selection(
    diagnostic_class,
    results_class$bnf_class_code[results_class$meaningful],
    "all primary meaningful classes"
  )
  add_selection(
    diagnostic_class,
    head(diagnostic_class[order(abs(seasonality_q_bh - fdr_alpha))]$series_code, 10L),
    "10 classes closest to the class BH boundary"
  )
  class_strength_boundary <- as.data.table(char_class)[
    ptr_lci >= meaningful_threshold
  ][order(abs(stl_seasonal_strength - stl_strength_threshold))]
  add_selection(
    diagnostic_class, head(class_strength_boundary$bnf_class_code, 10L),
    "10 amplitude-qualified classes closest to the STL boundary"
  )
  for (metric in c("max_abs_pearson_residual", "max_abs_residual_acf",
                   "max_cooks_distance", "max_segment_mean_shift_sd")) {
    add_selection(
      diagnostic_class, top_metric_codes(diagnostic_class, metric),
      paste("top five classes by", metric)
    )
  }
  add_selection(
    diagnostic_class,
    diagnostic_class[harmonic_peak_trough_ratio >= 10, series_code],
    "extreme harmonic amplitude (peak-to-trough ratio at least 10)"
  )
  add_selection(
    diagnostic_class,
    head(diagnostic_class[
      diagnostic_stl_trend_strength >= 0.80 &
        diagnostic_stl_seasonal_strength < stl_strength_threshold
    ][order(-diagnostic_stl_trend_strength)]$series_code, 5L),
    "five strongest trend-dominated class profiles"
  )
  add_selection(
    diagnostic_class,
    diagnostic_class[
      !original_model_converged |
        (!is.na(original_model_note) & nzchar(original_model_note)) |
        grepl("NBfail", route, fixed = TRUE),
      series_code
    ],
    "model failure, warning note or negative-binomial fallback"
  )

  add_selection(
    diagnostic_drug,
    head(as.data.table(results_drug)[meaningful == TRUE][order(-peak_trough_ratio)]$bnf_drug_code,
         10L),
    "10 highest-amplitude meaningful exploratory drugs"
  )
  add_selection(
    diagnostic_drug,
    head(diagnostic_drug[order(abs(seasonality_q_bh - fdr_alpha))]$series_code, 10L),
    "10 drugs closest to the all-drug BH boundary"
  )
  drug_strength_boundary <- as.data.table(char_drug)[
    ptr_lci >= meaningful_threshold
  ][order(abs(stl_seasonal_strength - stl_strength_threshold))]
  add_selection(
    diagnostic_drug, head(drug_strength_boundary$bnf_drug_code, 20L),
    "20 amplitude-qualified drugs closest to the STL boundary"
  )
  threshold_changes <- fread(
    file.path(prior_specs$threshold[[2]], "threshold_classification_changes.csv"),
    colClasses = list(character = "series_code")
  )
  add_selection(
    diagnostic_drug,
    threshold_changes[level == "drug", series_code],
    "drug meaningfulness changed in the 0.40/0.60 STL sensitivity"
  )
  working_changes <- fread(
    file.path(prior_specs$working_days[[2]], "working_day_classification_changes.csv"),
    colClasses = list(character = "code")
  )
  add_selection(
    diagnostic_drug,
    working_changes[
      level == "drug" & primary_meaningful & !meaningful, code
    ],
    "primary meaningful drug lost under the working-day sensitivity"
  )
  for (metric in c("max_abs_pearson_residual", "max_abs_residual_acf",
                   "max_cooks_distance", "max_segment_mean_shift_sd")) {
    add_selection(
      diagnostic_drug, top_metric_codes(diagnostic_drug, metric),
      paste("top five drugs by", metric)
    )
  }
  add_selection(
    diagnostic_drug,
    diagnostic_drug[harmonic_peak_trough_ratio >= 10, series_code],
    "extreme harmonic amplitude (peak-to-trough ratio at least 10)"
  )
  add_selection(
    diagnostic_drug,
    head(diagnostic_drug[
      diagnostic_stl_trend_strength >= 0.80 &
        diagnostic_stl_seasonal_strength < stl_strength_threshold
    ][order(-diagnostic_stl_trend_strength)]$series_code, 5L),
    "five strongest trend-dominated drug profiles"
  )
  add_selection(
    diagnostic_drug,
    diagnostic_drug[
      !original_model_converged |
        (!is.na(original_model_note) & nzchar(original_model_note)) |
        grepl("NBfail", route, fixed = TRUE),
      series_code
    ],
    "model failure, warning note or negative-binomial fallback"
  )

  diagnostic_selection <- rbindlist(selection_parts, use.names = TRUE)[, .(
    selection_reason = paste(sort(unique(selection_reason)), collapse = "; ")
  ), by = .(level, series_code, series_name)]
  diagnostic_selection <- merge(
    diagnostic_selection,
    diagnostic_inventory,
    by = c("level", "series_code", "series_name"), all.x = TRUE
  )
  setorder(diagnostic_selection, level, series_code)

  extract_targeted_details <- function(selection) {
    time_parts <- list()
    acf_parts <- list()
    for (i in seq_len(nrow(selection))) {
      level_name <- selection$level[i]
      code <- selection$series_code[i]
      if (level_name == "class") {
        monthly <- class_monthly_elig[
          as.character(class_monthly_elig$bnf_class_code) == code,
          c("year_month", "items")
        ]
        screen_row <- screen_class[
          as.character(screen_class$bnf_class_code) == code, , drop = FALSE
        ]
      } else {
        monthly <- drug_monthly_elig[
          as.character(drug_monthly_elig$bnf_drug_code) == code,
          c("year_month", "items")
        ]
        screen_row <- screen_drug[
          as.character(screen_drug$bnf_drug_code) == code, , drop = FALSE
        ]
      }
      fitted_object <- fit_primary_mean(monthly, screen_row)
      if (!fitted_object$ok) next
      d <- fitted_object$data
      fit <- fitted_object$fit
      pearson <- as.numeric(residuals(fit, type = "pearson"))
      cooks <- suppressWarnings(as.numeric(cooks.distance(fit)))
      fitted_count <- as.numeric(fitted(fit))
      time_parts[[length(time_parts) + 1L]] <- data.table(
        level = level_name, series_code = code,
        series_name = selection$series_name[i],
        selection_reason = selection$selection_reason[i],
        year_month = d$year_month, month_date = d$month_date,
        observed_items = d$items, fitted_items = fitted_count,
        observed_items_per_1000_patients = d$items / d$list_size * 1000,
        fitted_items_per_1000_patients = fitted_count / d$list_size * 1000,
        pearson_residual = pearson, cooks_distance = cooks
      )
      acf_values <- as.numeric(stats::acf(
        pearson, lag.max = harmonic_period_months[1], plot = FALSE,
        demean = TRUE, na.action = na.pass
      )$acf)[-1]
      acf_parts[[length(acf_parts) + 1L]] <- data.table(
        level = level_name, series_code = code,
        series_name = selection$series_name[i],
        selection_reason = selection$selection_reason[i],
        lag_months = seq_along(acf_values), residual_acf = acf_values
      )
    }
    list(
      timeseries = rbindlist(time_parts, use.names = TRUE),
      acf = rbindlist(acf_parts, use.names = TRUE)
    )
  }
  targeted_details <- extract_targeted_details(diagnostic_selection)
  targeted_residual_timeseries <- targeted_details$timeseries
  targeted_residual_acf <- targeted_details$acf

  add_stage6_check(
    "targeted_selection_complete",
    all(results_class$bnf_class_code[results_class$meaningful] %in%
          diagnostic_selection[level == "class", series_code]) &&
      nrow(diagnostic_selection[level == "class"]) < 220L &&
      nrow(diagnostic_selection[level == "drug"]) < 974L,
    "all 30 meaningful classes plus targeted boundaries/flags; not all drugs",
    sprintf("%d classes; %d drugs selected",
            nrow(diagnostic_selection[level == "class"]),
            nrow(diagnostic_selection[level == "drug"]))
  )
  selected_n <- nrow(diagnostic_selection)
  add_stage6_check(
    "targeted_details_complete",
    nrow(targeted_residual_timeseries) == selected_n * 48L &&
      nrow(targeted_residual_acf) == selected_n * 12L &&
      !anyNA(targeted_residual_timeseries$pearson_residual) &&
      !anyNA(targeted_residual_acf$residual_acf),
    "48 monthly and 12 residual-ACF rows for every selected series",
    sprintf("%d monthly rows; %d ACF rows for %d series",
            nrow(targeted_residual_timeseries), nrow(targeted_residual_acf), selected_n)
  )

  cook_review_thresholds <- diagnostic_inventory[, .(
    cook_distance_review_threshold = as.numeric(
      quantile(max_cooks_distance, 0.99, na.rm = TRUE, names = FALSE)
    )
  ), by = level]
  diagnostic_review_flags <- merge(
    diagnostic_selection, cook_review_thresholds, by = "level", all.x = TRUE
  )[, .(
    level, series_code, series_name, inference_scope, selection_reason,
    seasonality_q_bh, seasonality_detected,
    upper_tail_influence = max_cooks_distance >= cook_distance_review_threshold,
    cook_distance_review_threshold,
    very_high_residual_acf = max_abs_residual_acf >= 0.80,
    large_segment_mean_shift = max_segment_mean_shift_sd >= 2,
    extreme_harmonic_amplitude = harmonic_peak_trough_ratio >= 10,
    diagnostic_refit_failed = !diagnostic_refit_converged,
    max_abs_pearson_residual, max_abs_pearson_month,
    max_abs_residual_acf, max_abs_residual_acf_lag,
    max_cooks_distance, max_cooks_month,
    max_segment_mean_shift_sd, segment_shift_after_month,
    observed_fitted_log_rate_r2, harmonic_peak_trough_ratio,
    diagnostic_stl_seasonal_strength, diagnostic_stl_trend_strength
  )]
  diagnostic_review_flags[, any_exception :=
    upper_tail_influence | very_high_residual_acf | large_segment_mean_shift |
      extreme_harmonic_amplitude | diagnostic_refit_failed]
  setorder(diagnostic_review_flags, -any_exception, level, series_code)

  make_cohort_flow <- function(eligibility, screen, results, id_col, level_name) {
    tab <- as.data.table(copy(eligibility))
    tab[, series_code := as.character(get(id_col))]
    observed_items <- sum(tab$total_items)
    screen_codes <- as.character(screen[[id_col]])
    complete_codes <- as.character(screen[[id_col]][screen$converged])
    result_codes <- as.character(results[[id_col]])
    meaningful_codes <- as.character(results[[id_col]][results$meaningful])
    detected_codes <- if (level_name == "class") {
      as.character(screen[[id_col]][screen$class_significant])
    } else {
      as.character(screen[[id_col]][screen$drug_significant])
    }
    row_for <- function(step_order, step, rule_type, codes, notes) {
      codes <- unique(as.character(codes))
      included_items <- sum(tab$total_items[tab$series_code %in% codes])
      data.table(
        level = level_name,
        inference_scope = if (level_name == "class") {
          "primary_inferential"
        } else {
          "secondary_exploratory"
        },
        step_order, step, rule_type,
        n_series = length(codes),
        n_not_in_step_from_observed = nrow(tab) - length(codes),
        included_items,
        included_item_share = included_items / observed_items,
        notes
      )
    }
    rbindlist(list(
      row_for(1L, "observed_after_chapter_restriction", "starting_universe",
              tab$series_code,
              "Series accounting begins after restriction to BNF chapters 01-14."),
      row_for(2L, "available_in_all_48_months", "marginal_eligibility_rule",
              tab$series_code[tab$rule_every_month],
              "Marginal rule; exclusions overlap with the annual-volume rule."),
      row_for(3L, "met_annual_volume_rule_in_all_4_years", "marginal_eligibility_rule",
              tab$series_code[tab$rule_min_volume],
              "Marginal rule; at least 1,000 items in each calendar year."),
      row_for(4L, "eligible_both_rules", "combined_eligibility",
              tab$series_code[tab$eligible],
              "Intersection of complete monthly coverage and annual volume."),
      row_for(5L, "model_attempted", "downstream_flow", screen_codes,
              "Every eligible series entered its complete testing family."),
      row_for(6L, "model_completed", "downstream_flow", complete_codes,
              "Convergence/failure status retained in screening outputs."),
      row_for(7L, "included_in_complete_bh_family", "downstream_flow", complete_codes,
              if (level_name == "class") {
                "BH across all eligible classes."
              } else {
                "BH across all eligible drugs, independently of parent class."
              }),
      row_for(8L, "detected_at_family_fdr_5_percent", "analysis_result", detected_codes,
              "FDR claim applies within this complete family only."),
      row_for(9L, "meaningful_descriptive_classification", "analysis_result",
              meaningful_codes,
              "Detected series also meeting the amplitude and STL rules.")
    ), use.names = TRUE)
  }
  cohort_flow <- rbindlist(list(
    make_cohort_flow(elig_class, screen_class, results_class,
                     "bnf_class_code", "class"),
    make_cohort_flow(elig_drug, screen_drug, results_drug,
                     "bnf_drug_code", "drug")
  ), use.names = TRUE)

  make_exclusion_overlap <- function(eligibility, level_name) {
    as.data.table(copy(eligibility))[, .(
      n_series = .N,
      total_items = sum(total_items)
    ), by = .(rule_every_month, rule_min_volume)][, `:=`(
      level = level_name,
      status = fifelse(
        rule_every_month & rule_min_volume, "eligible",
        fifelse(!rule_every_month & !rule_min_volume,
                "failed both rules",
                fifelse(!rule_every_month, "failed monthly coverage only",
                        "failed annual volume only"))
      )
    )]
  }
  exclusion_overlap_summary <- rbindlist(list(
    make_exclusion_overlap(elig_class, "class"),
    make_exclusion_overlap(elig_drug, "drug")
  ), use.names = TRUE)
  exclusion_overlap_summary[, item_share_within_level :=
    total_items / sum(total_items), by = level]
  setcolorder(exclusion_overlap_summary,
              c("level", "status", "rule_every_month", "rule_min_volume",
                "n_series", "total_items", "item_share_within_level"))
  setorder(exclusion_overlap_summary, level, -rule_every_month, -rule_min_volume)

  epd_main <- epd_file_qc[analytical_role == "main_2022_2025"]
  source_scope_flow <- data.table(
    step_order = 1:3,
    step = c("raw_EPD_source", "retained_BNF_chapters_01_14",
             "removed_outside_BNF_chapters_01_14"),
    n_records = c(
      sum(epd_main$raw_row_count),
      sum(epd_main$retained_chapter_01_14_rows),
      sum(epd_main$removed_chapter_rows)
    ),
    n_items = c(
      sum(epd_main$raw_total_items),
      sum(epd_main$retained_chapter_01_14_items),
      sum(epd_main$removed_chapter_items)
    )
  )
  source_scope_flow[, item_share_of_raw := n_items / n_items[step_order == 1L]]
  source_scope_flow[, notes := c(
    "Forty-eight monthly EPD archives in the primary 2022-2025 window.",
    "Starting universe for class/drug series accounting.",
    "Outside the declared medicines scope; not treated as missing data."
  )]

  add_stage6_check(
    "cohort_flow_reconciles",
    identical(cohort_flow[level == "class", n_series],
              c(344L, 236L, 220L, 220L, 220L, 220L, 220L, 125L, 30L)) &&
      identical(cohort_flow[level == "drug", n_series],
                c(2155L, 1194L, 975L, 974L, 974L, 974L, 974L, 391L, 88L)) &&
      sum(exclusion_overlap_summary[level == "class", n_series]) == 344L &&
      sum(exclusion_overlap_summary[level == "drug", n_series]) == 2155L,
    "344/220/125/30 classes and 2155/974/391/88 drugs reconcile",
    paste(cohort_flow[step_order %in% c(1L, 4L, 8L, 9L),
                      paste(level, step, n_series, sep = "=")], collapse = "; ")
  )
  add_stage6_check(
    "source_scope_reconciles",
    source_scope_flow$n_records[1] ==
      source_scope_flow$n_records[2] + source_scope_flow$n_records[3] &&
      source_scope_flow$n_items[1] ==
        source_scope_flow$n_items[2] + source_scope_flow$n_items[3] &&
      source_scope_flow$n_items[2] == 4641000754,
    "raw EPD records/items partition exactly into retained and out-of-scope",
    sprintf("%s raw items; %s retained; %s removed",
            format(source_scope_flow$n_items[1], scientific = FALSE),
            format(source_scope_flow$n_items[2], scientific = FALSE),
            format(source_scope_flow$n_items[3], scientific = FALSE))
  )

  list_size_main <- list_size_source_qc[analytical_role == "main_2022_2025"]
  invalid_epd_items <- epd_main[, sum(
    missing_items_rows + nonnumeric_items_rows + negative_items_rows +
      noninteger_items_rows
  )]
  missing_epd_keys <- epd_main[, sum(
    month_mismatch_rows + missing_practice_rows + missing_presentation_rows
  )]
  invalid_list_size <- list_size_main[, sum(
    missing_patient_counts + nonnumeric_patient_counts +
      negative_patient_counts + noninteger_patient_counts
  )]
  missing_model_inputs_class <- class_monthly_elig |>
    group_by(bnf_class_code) |>
    summarise(bad = n_distinct(year_month) != 48L || anyNA(items), .groups = "drop") |>
    summarise(n = sum(bad)) |>
    pull(n)
  missing_model_inputs_drug <- drug_monthly_elig |>
    group_by(bnf_drug_code) |>
    summarise(bad = n_distinct(year_month) != 48L || anyNA(items), .groups = "drop") |>
    summarise(n = sum(bad)) |>
    pull(n)

  missingness_summary <- rbindlist(list(
    data.table(
      stage = "raw_EPD_binding", variable_or_check = "required month/practice/presentation keys",
      n_missing_or_invalid = missing_epd_keys, denominator = sum(epd_main$raw_row_count),
      handling = "Required analytical keys; zero required for release.",
      affects_analysis = missing_epd_keys > 0
    ),
    data.table(
      stage = "raw_EPD_binding", variable_or_check = "ITEMS missing/nonnumeric/negative/noninteger",
      n_missing_or_invalid = invalid_epd_items, denominator = sum(epd_main$raw_row_count),
      handling = "Invalid item values would stop the pipeline.",
      affects_analysis = invalid_epd_items > 0
    ),
    data.table(
      stage = "raw_EPD_binding", variable_or_check = "SNOMED code missing",
      n_missing_or_invalid = sum(epd_main$missing_snomed_rows),
      denominator = sum(epd_main$raw_row_count),
      handling = "Retained: SNOMED is neither an analytical key nor an analysis variable.",
      affects_analysis = FALSE
    ),
    data.table(
      stage = "descriptor_lookup", variable_or_check = "May-2025 BNF class descriptor unmatched",
      n_missing_or_invalid = sum(epd_main$unmatched_lookup_raw_rows),
      denominator = sum(epd_main$retained_chapter_01_14_rows),
      handling = paste0("Retained by code; ", sum(epd_main$unmatched_lookup_items),
                        " items have missing lookup-derived descriptors."),
      affects_analysis = FALSE
    ),
    data.table(
      stage = "list_size_binding", variable_or_check = "practice code missing",
      n_missing_or_invalid = sum(list_size_main$missing_practice_codes),
      denominator = sum(list_size_main$source_rows),
      handling = "Required source key; zero required for release.",
      affects_analysis = sum(list_size_main$missing_practice_codes) > 0
    ),
    data.table(
      stage = "list_size_binding", variable_or_check = "patient count missing/invalid",
      n_missing_or_invalid = invalid_list_size,
      denominator = sum(list_size_main$source_rows),
      handling = "Invalid denominator values would stop the pipeline.",
      affects_analysis = invalid_list_size > 0
    ),
    data.table(
      stage = "monthly_aggregate", variable_or_check = "duplicate class/drug month keys",
      n_missing_or_invalid = sum(epd_main$drug_month_duplicate_rows) +
        sum(epd_main$class_month_duplicate_rows),
      denominator = 48L,
      handling = "National aggregates must be unique by series and month.",
      affects_analysis = FALSE
    ),
    data.table(
      stage = "monthly_covariates", variable_or_check = "missing/nonfinite list size or primary offset",
      n_missing_or_invalid = sum(!is.finite(covar$list_size)) +
        sum(!is.finite(covar$offset_log_patient_days)),
      denominator = nrow(covar) * 2L,
      handling = "Every primary month requires both values.",
      affects_analysis = FALSE
    ),
    data.table(
      stage = "model_input", variable_or_check = "eligible series lacking 48 complete item months",
      n_missing_or_invalid = missing_model_inputs_class + missing_model_inputs_drug,
      denominator = nrow(screen_class) + nrow(screen_drug),
      handling = "Eligibility and fit functions independently require 48 months.",
      affects_analysis = FALSE
    ),
    data.table(
      stage = "model_output", variable_or_check = "missing p/q value or non-converged model",
      n_missing_or_invalid = sum(!screen_class$converged | is.na(screen_class$p_value) |
                                   is.na(screen_class$class_q_bh)) +
        sum(!screen_drug$converged | is.na(screen_drug$p_value) |
              is.na(screen_drug$drug_all_q_bh)),
      denominator = nrow(screen_class) + nrow(screen_drug),
      handling = "Retained explicitly in screening and model-failure outputs.",
      affects_analysis = FALSE
    ),
    data.table(
      stage = "characterisation", variable_or_check = "missing amplitude CI or STL strength",
      n_missing_or_invalid = sum(!is.finite(char_class$ptr_lci) |
                                   !is.finite(char_class$ptr_uci) |
                                   !is.finite(char_class$stl_seasonal_strength)) +
        sum(!is.finite(char_drug$ptr_lci) |
              !is.finite(char_drug$ptr_uci) |
              !is.finite(char_drug$stl_seasonal_strength)),
      denominator = nrow(char_class) + nrow(char_drug),
      handling = "A failed uncertainty or STL calculation remains visible and cannot be meaningful.",
      affects_analysis = FALSE
    )
  ), use.names = TRUE)

  warning_summary <- data.table(
    category = c(
      "captured_runtime_warning_messages", "class_model_failures",
      "drug_model_failures", "class_negative_binomial_fallbacks",
      "drug_negative_binomial_fallbacks", "class_HAC_bandwidth_capped",
      "drug_HAC_bandwidth_capped", "failed_amplitude_confidence_intervals",
      "missing_STL_strength_values"
    ),
    n = c(
      length(unique(warnings_seen)), sum(!screen_class$converged),
      sum(!screen_drug$converged), sum(grepl("NBfail", screen_class$route, fixed = TRUE)),
      sum(grepl("NBfail", screen_drug$route, fixed = TRUE)),
      sum(screen_class$hac_capped %in% TRUE), sum(screen_drug$hac_capped %in% TRUE),
      sum(!is.finite(char_class$ptr_lci) | !is.finite(char_class$ptr_uci)) +
        sum(!is.finite(char_drug$ptr_lci) | !is.finite(char_drug$ptr_uci)),
      sum(!is.finite(char_class$stl_seasonal_strength)) +
        sum(!is.finite(char_drug$stl_seasonal_strength))
    ),
    interpretation = c(
      "Unique warnings reaching the analysis-level warning handler before Stage 6 completion.",
      "Primary class model failures.", "Exploratory drug model failures.",
      "Overdispersed class series whose NB route failed and used HAC.",
      "Overdispersed drug series whose NB route failed and used HAC.",
      "Class HAC bandwidth reached the declared cap of 12.",
      "Drug HAC bandwidth reached the declared cap of 12.",
      "Detected series without a finite coefficient-draw amplitude interval.",
      "Detected series without a finite STL seasonal strength."
    )
  )
  add_stage6_check(
    "missingness_and_warnings_accounted",
    all(missingness_summary[affects_analysis == TRUE, n_missing_or_invalid] == 0) &&
      warning_summary[category == "class_model_failures", n] == 0L &&
      warning_summary[category == "drug_model_failures", n] == 0L &&
      warning_summary[category == "failed_amplitude_confidence_intervals", n] == 0L &&
      warning_summary[category == "missing_STL_strength_values", n] == 0L &&
      warning_summary[category == "drug_negative_binomial_fallbacks", n] == 2L,
    "no missing required analytical values/failures; two documented drug NB fallbacks",
    paste(warning_summary$category, warning_summary$n, sep = "=", collapse = "; ")
  )

  output_files <- c(
    "model_diagnostic_inventory.csv", "diagnostic_selection.csv",
    "targeted_residual_timeseries.csv", "targeted_residual_acf.csv",
    "diagnostic_review_flags.csv", "cohort_flow.csv",
    "exclusion_overlap_summary.csv", "source_scope_flow.csv",
    "missingness_summary.csv", "warning_summary.csv"
  )
  atomic_fwrite(diagnostic_inventory, file.path(stage6_dir, output_files[1]))
  atomic_fwrite(diagnostic_selection, file.path(stage6_dir, output_files[2]))
  atomic_fwrite(targeted_residual_timeseries, file.path(stage6_dir, output_files[3]))
  atomic_fwrite(targeted_residual_acf, file.path(stage6_dir, output_files[4]))
  atomic_fwrite(diagnostic_review_flags, file.path(stage6_dir, output_files[5]))
  atomic_fwrite(cohort_flow, file.path(stage6_dir, output_files[6]))
  atomic_fwrite(exclusion_overlap_summary, file.path(stage6_dir, output_files[7]))
  atomic_fwrite(source_scope_flow, file.path(stage6_dir, output_files[8]))
  atomic_fwrite(missingness_summary, file.path(stage6_dir, output_files[9]))
  atomic_fwrite(warning_summary, file.path(stage6_dir, output_files[10]))
  add_stage6_check(
    "outputs_written",
    all(file.exists(file.path(stage6_dir, output_files))) &&
      all(file.info(file.path(stage6_dir, output_files))$size > 0),
    paste(length(output_files), "non-empty Stage 6 analytical files"),
    paste(sum(file.exists(file.path(stage6_dir, output_files))), "present")
  )

  stage6_qc_summary <- rbindlist(stage6_checks)
  atomic_fwrite(stage6_qc_summary, file.path(stage6_dir, "stage6_qc_summary.csv"))
  if (!all(stage6_qc_summary$pass)) {
    failed <- stage6_qc_summary[!pass, check_id]
    stop("Stage 6 completion gate failed: ", paste(failed, collapse = ", "),
         ". See ", file.path(stage6_dir, "stage6_qc_summary.csv"), ".")
  }

  for (relative_path in output_files) {
    source_path <- file.path(stage6_dir, relative_path)
    target_path <- file.path(stage6_snapshot_dir, relative_path)
    temporary_path <- paste0(target_path, ".tmp")
    if (!file.copy(source_path, temporary_path, overwrite = TRUE) ||
        !file.rename(temporary_path, target_path)) {
      stop("Could not seal Stage 6 snapshot file: ", relative_path)
    }
  }
  stage6_manifest <- data.table(
    relative_path = output_files,
    bytes = as.numeric(file.info(file.path(stage6_snapshot_dir, output_files))$size),
    sha256 = vapply(file.path(stage6_snapshot_dir, output_files),
                    sha256_file, character(1))
  )
  atomic_fwrite(stage6_manifest, file.path(stage6_dir, "stage6_snapshot_manifest.csv"))
  stage6_completion <- data.table(
    stage = "stage6", status = "PASS",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    checks_passed = sum(stage6_qc_summary$pass),
    checks_total = nrow(stage6_qc_summary),
    diagnostic_inventory_rows = nrow(diagnostic_inventory),
    targeted_series = nrow(diagnostic_selection),
    diagnostic_exception_rows = sum(diagnostic_review_flags$any_exception),
    snapshot_files = nrow(stage6_manifest),
    analysis_script_sha256 = sha256_file(script_path),
    renv_lock_sha256 = sha256_file(lock_path)
  )
  atomic_fwrite(stage6_completion, file.path(stage6_dir, "stage6_completion.csv"))
  message(
    "Stage 6 complete: all-series diagnostic inventory, targeted residual review, ",
    "cohort flow and missingness/warning accounting passed ",
    nrow(stage6_qc_summary), " checks; stopping before Stage 7 reporting crosswalk."
  )
  return(invisible(TRUE))
}

## ---- shared plotting theme --------------------------------------------------

theme_pub <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# observed + fitted rate (per 1000 patients) for one class or drug series
.series_fit <- function(mon, id_col, id_val) {
  ser <- mon |> filter(.data[[id_col]] == id_val) |> select(year_month, items)
  d   <- covar |> left_join(ser, by = "year_month") |> arrange(t)
  d$off <- d$offset_log_patient_days
  fit <- suppressWarnings(glm(.f_full, poisson, data = d))
  tibble(month_date = d$month_date,
         observed = d$items / d$list_size * 1000,
         fitted   = as.numeric(fitted(fit)) / d$list_size * 1000)
}

## ---- MAIN TEXT figures: exemplars ------------------------------------------

exemplar_codes <- results_class |>
  filter(meaningful, peak_trough_ratio < 10) |>
  slice_head(n = 6) |> pull(bnf_class_code)

if (length(exemplar_codes)) {
  ex_names <- setNames(char_class$bnf_class_name[match(exemplar_codes, char_class$bnf_class_code)],
                       exemplar_codes)
  lvl <- results_class |> filter(bnf_class_code %in% exemplar_codes) |>
    arrange(desc(peak_trough_ratio)) |> pull(bnf_class_name)
  
  fit_df <- bind_rows(lapply(exemplar_codes, function(cd)
    .series_fit(class_monthly_elig, "bnf_class_code", cd) |> mutate(class = ex_names[as.character(cd)])))
  fit_df$class <- factor(fit_df$class, levels = lvl)
  
  p1 <- ggplot(fit_df, aes(month_date)) +
    geom_point(aes(y = observed), size = 1, colour = "grey45") +
    geom_line(aes(y = fitted), colour = "#B2182B", linewidth = 0.8) +
    facet_wrap(~ class, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Prescription items per 1000 registered patients",
         title = "Observed monthly prescribing with fitted trend and seasonal model") +
    theme_pub
  ggsave(file.path(fig_dir, "Figure2_exemplar_observed_fitted.png"), p1, width = 9, height = 8, dpi = 300)
  
  shape_df <- char_class |> filter(bnf_class_code %in% exemplar_codes) |>
    group_by(bnf_class_name) |>
    reframe(month = 1:12,
            factor = { s <- .seasonal_curve(b_sin12, b_cos12, b_sin6, b_cos6, 1:12); exp(s)/mean(exp(s)) })
  shape_df$bnf_class_name <- factor(shape_df$bnf_class_name, levels = lvl)
  p2 <- ggplot(shape_df, aes(month, factor)) +
    geom_hline(yintercept = 1, colour = "grey70", linetype = 2) +
    geom_line(colour = "#2166AC", linewidth = 0.9) +
    facet_wrap(~ bnf_class_name, ncol = 2) +
    scale_x_continuous(breaks = c(1,4,7,10), labels = month.abb[c(1,4,7,10)]) +
    labs(x = NULL, y = "Seasonal factor (relative to annual mean)",
         title = "Fitted seasonal shape by calendar month") +
    theme_pub
  ggsave(file.path(fig_dir, "Figure3_exemplar_seasonal_shapes.png"), p2, width = 9, height = 8, dpi = 300)
}

## ---- MAIN TEXT figure: amplitude x peak-month landscape --------------------


land <- results_class |>
  filter(meaningful) |>
  left_join(chap_class, by = "bnf_class_code") |>
  mutate(peak_month = factor(peak_month, levels = month.abb),
         chapter = bnf_chapter_name)
# visible ceiling: classes above y_cap are drawn as triangles pinned at the top
# with their true ratio in the label, so the bulk of classes spread out below.
y_cap <- 3
land <- land |>
  mutate(above  = peak_trough_ratio > y_cap,
         y_plot = pmin(peak_trough_ratio, y_cap),
         lab    = ifelse(above,
                         sprintf("%s (%.0f\u00d7)", bnf_class_name, peak_trough_ratio),
                         bnf_class_name))
if (nrow(land)) {
  pal <- grDevices::hcl.colors(max(3, length(unique(land$chapter))), "Dark 3")
  set.seed(plot_jitter_seed)
  jit <- position_jitter(width = 0.18, height = 0, seed = plot_jitter_seed)
  p3 <- ggplot(land, aes(peak_month, y_plot, colour = chapter)) +
    geom_hline(yintercept = 1, colour = "grey80") +
    geom_point(aes(shape = above), size = 2.8, alpha = 0.9, position = jit) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_y_log10(breaks = c(1, 1.25, 1.5, 2, y_cap),
                  labels = c("1", "1.25", "1.5", "2", paste0("\u2265", y_cap)),
                  expand = expansion(mult = c(0.03, 0.14))) +
    scale_x_discrete(drop = FALSE) +
    scale_colour_manual(values = pal, name = "BNF chapter") +
    labs(x = "Month of peak prescribing", y = "Peak-to-trough ratio (log scale)",
         title = "Seasonality landscape: amplitude and timing of meaningfully seasonal drug classes",
         caption = paste0("\u25b2 = class exceeds the axis ceiling (", y_cap,
                          "\u00d7); actual peak-to-trough ratio shown in the label.")) +
    theme_pub + theme(legend.position = "right", legend.text = element_text(size = 8),
                      plot.caption = element_text(size = 8, hjust = 0))
  if (has_repel)
    p3 <- p3 + ggrepel::geom_text_repel(aes(label = lab), size = 2.4,
                                        max.overlaps = Inf, seed = plot_jitter_seed, position = jit, show.legend = FALSE, force = 2, point.padding = NA, min.segment.length = 0, box.padding = 0.5)
  ggsave(file.path(fig_dir, "Figure1_seasonality_landscape.png"), p3, width = 11, height = 8, dpi = 300)
}

## ---- APPENDIX figure: genuine vs excluded (why exclusions were made) --------

pick <- function(nm) results_class$bnf_class_code[match(nm, results_class$bnf_class_name)]
compare_codes <- na.omit(c(
  kept    = pick("Tetracyclines"), kept2 = pick("Antihistamines"), kept3 = pick("Penicillins"),
  break1  = pick("Amino acids and nutritional agents"),
  break2  = pick("Digestive aids"),
  trend1  = pick("Filaricides")))
if (length(compare_codes) >= 2) {
  cmp_df <- bind_rows(lapply(compare_codes, function(cd) {
    nm <- results_class$bnf_class_name[results_class$bnf_class_code == cd][1]
    keep <- results_class$meaningful[results_class$bnf_class_code == cd][1]
    .series_fit(class_monthly_elig, "bnf_class_code", cd) |>
      mutate(class = sprintf("%s (%s)", nm, ifelse(keep, "retained", "excluded")))
  }))
  p4 <- ggplot(cmp_df, aes(month_date)) +
    geom_point(aes(y = observed), size = 0.9, colour = "grey45") +
    geom_line(aes(y = fitted), colour = "#B2182B", linewidth = 0.7) +
    facet_wrap(~ class, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Prescription items per 1000 registered patients",
         title = "Genuine seasonality versus excluded structural-break and trend-dominated series") +
    theme_pub
  ggsave(file.path(fig_dir, "FigureS2_retained_vs_nonmeaningful.png"), p4, width = 9, height = 8, dpi = 300)
}

## ---- APPENDIX figures: ALL series, observed + fitted, paginated by chapter --
# Every eligible series is plotted. Facet-strip colour flags analytic status:
#   dark red + white text  = meaningful seasonality
#   pale red + dark text    = statistically significant but not meaningful
#   default grey            = eligible but not significant
# Strips are recoloured on the plot's grob, matched to each panel by its label.

# read a facet strip's panel label
.strip_label <- function(sg) {
  gt <- sg$grobs[[1]]
  for (ch in gt$children) if (!is.null(ch$children))
    for (tx in ch$children) if (!is.null(tx$label)) return(tx$label)
  NA_character_
}
# recolour one strip: background fill + text colour
.recolour_strip <- function(sg, fill, textcol) {
  gt <- sg$grobs[[1]]
  for (i in seq_along(gt$children)) {
    ch <- gt$children[[i]]
    if (inherits(ch, "rect")) {
      if (is.null(ch$gp)) ch$gp <- grid::gpar(); ch$gp$fill <- fill; ch$gp$col <- NA
    } else if (!is.null(ch$children)) {
      for (j in seq_along(ch$children))
        if (inherits(ch$children[[j]], "text")) {
          if (is.null(ch$children[[j]]$gp)) ch$children[[j]]$gp <- grid::gpar()
          ch$children[[j]]$gp$col <- textcol
        }
    }
    gt$children[[i]] <- ch
  }
  sg$grobs[[1]] <- gt; sg
}
# apply a tier lookup (named list: label -> c(fill, textcol)) to all strips
.style_strips <- function(p, tier_map) {
  g <- ggplotGrob(p)
  for (si in which(grepl("^strip", g$layout$name))) {
    lab <- .strip_label(g$grobs[[si]])
    if (!is.na(lab) && lab %in% names(tier_map))
      g$grobs[[si]] <- .recolour_strip(g$grobs[[si]], tier_map[[lab]][1], tier_map[[lab]][2])
  }
  g
}

# Draw a partially filled appendix page without enlarging its panels. The plot
# is scaled into the occupied cells of a fixed 4 x 4 page grid; unused cells
# remain blank. This keeps a one-panel page at the same panel size as a full
# 16-panel page.
.draw_fixed_4x4_page <- function(plot_grob, page_title, n_panels, first_page,
                                 y_label) {
  n_panels <- as.integer(n_panels)
  if (length(n_panels) != 1L || is.na(n_panels) ||
      n_panels < 1L || n_panels > 16L) {
    stop("Fixed appendix pages require between 1 and 16 panels.")
  }
  if (!isTRUE(first_page)) grid::grid.newpage()

  rows_used <- ceiling(n_panels / 4L)
  cols_used <- if (rows_used == 1L) n_panels else 4L
  plot_height <- 0.90 * rows_used / 4

  grid::grid.text(
    page_title,
    x = grid::unit(0.035, "npc"), y = grid::unit(0.985, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontface = "bold", fontsize = 9)
  )
  grid::grid.text(
    y_label,
    x = grid::unit(0.012, "npc"),
    y = grid::unit(0.93 - plot_height / 2, "npc"),
    rot = 90, gp = grid::gpar(fontsize = 7)
  )

  plot_vp <- grid::viewport(
    x = grid::unit(0.035, "npc"), y = grid::unit(0.93, "npc"),
    width = grid::unit(0.955 * cols_used / 4, "npc"),
    height = grid::unit(plot_height, "npc"),
    just = c("left", "top")
  )
  grid::pushViewport(plot_vp)
  grid::grid.draw(plot_grob)
  grid::popViewport()
  invisible(NULL)
}

# tier colour lookups keyed by series NAME (meaningful takes precedence over
# significant). results_* contain only significant/characterised series; the
# `meaningful` column marks the top tier. Stage 4/full runs use the complete
# all-eligible-drug BH family as the sole drug-level significance definition.
FILL_MEAN <- "#8B0000"; TXT_MEAN <- "white"
FILL_SIG  <- "#E8A0A0"; TXT_SIG  <- "grey15"
tier_class <- {
  m <- setNames(vector("list", nrow(results_class)), results_class$bnf_class_name)
  for (i in seq_len(nrow(results_class)))
    m[[i]] <- if (isTRUE(results_class$meaningful[i])) c(FILL_MEAN, TXT_MEAN) else c(FILL_SIG, TXT_SIG)
  m
}
tier_drug <- {
  sig <- if (stage3_legacy_inference) {
    (results_drug$significant_primary %in% TRUE) | (results_drug$sig_all %in% TRUE)
  } else {
    results_drug$drug_significant %in% TRUE
  }
  m <- setNames(vector("list", nrow(results_drug)), results_drug$bnf_drug_name)
  for (i in seq_len(nrow(results_drug)))
    m[[i]] <- if (isTRUE(results_drug$meaningful[i])) c(FILL_MEAN, TXT_MEAN)
  else if (sig[i])                        c(FILL_SIG, TXT_SIG)
  else                                    NULL   # eligible-only: default grey
  m[!vapply(m, is.null, logical(1))]
}

# Same tiering, blue palette, used only for the seasonal-SHAPE appendix figures
# (kept separate from tier_class/tier_drug so the observed+fitted appendix
# stays red and the shape appendix is blue).
FILL_MEAN_SHAPE <- "#08306B"; TXT_MEAN_SHAPE <- "white"   # dark blue
FILL_SIG_SHAPE  <- "#9ECAE1"; TXT_SIG_SHAPE  <- "grey15"   # pale blue
tier_class_shape <- {
  m <- setNames(vector("list", nrow(results_class)), results_class$bnf_class_name)
  for (i in seq_len(nrow(results_class)))
    m[[i]] <- if (isTRUE(results_class$meaningful[i])) c(FILL_MEAN_SHAPE, TXT_MEAN_SHAPE)
  else c(FILL_SIG_SHAPE, TXT_SIG_SHAPE)
  m
}
tier_drug_shape <- {
  sig <- if (stage3_legacy_inference) {
    (results_drug$significant_primary %in% TRUE) | (results_drug$sig_all %in% TRUE)
  } else {
    results_drug$drug_significant %in% TRUE
  }
  m <- setNames(vector("list", nrow(results_drug)), results_drug$bnf_drug_name)
  for (i in seq_len(nrow(results_drug)))
    m[[i]] <- if (isTRUE(results_drug$meaningful[i])) c(FILL_MEAN_SHAPE, TXT_MEAN_SHAPE)
  else if (sig[i])                        c(FILL_SIG_SHAPE, TXT_SIG_SHAPE)
  else                                    NULL
  m[!vapply(m, is.null, logical(1))]
}

render_all_series <- function(mon, id_col, name_col, chap_lookup, out_pdf,
                              tier_map = list(), per_page = 16L, ncol = 4L) {
  if (per_page != 16L || ncol != 4L) {
    stop("Paginated appendix PDFs use a fixed 4 x 4 layout (16 panels per page).")
  }
  ids <- chap_lookup |> arrange(bnf_chapter_code) |>
    left_join(distinct(mon, .data[[id_col]], .data[[name_col]]), by = id_col)
  pdf(out_pdf, width = 11.7, height = 8.3)   # A4 landscape
  on.exit(dev.off())
  first_page <- TRUE
  for (ch in unique(ids$bnf_chapter_code)) {
    chn <- ids |> filter(bnf_chapter_code == ch)
    pages <- split(seq_len(nrow(chn)), ceiling(seq_len(nrow(chn)) / per_page))
    for (pi in seq_along(pages)) {
      rows <- chn[pages[[pi]], ]
      # names are not guaranteed unique (e.g. two drug codes sharing a display
      # name) - facet levels must be unique, so disambiguate any duplicate on
      # THIS page by appending its code, and build the tier lookup for the page
      # using the same disambiguated label so strip colouring still matches.
      nm <- rows[[name_col]]
      dup <- nm %in% nm[duplicated(nm)]
      panel_label <- ifelse(dup, sprintf("%s [%s]", nm, rows[[id_col]]), nm)
      rows$panel_label <- panel_label
      page_tier_map <- setNames(tier_map[nm], panel_label)
      page_tier_map <- page_tier_map[!vapply(page_tier_map, is.null, logical(1))]
      
      dd <- bind_rows(lapply(seq_len(nrow(rows)), function(i) {
        .series_fit(mon, id_col, rows[[id_col]][i]) |>
          mutate(panel = rows$panel_label[i])
      }))
      dd$panel <- factor(dd$panel, levels = rows$panel_label)
      ptitle <- sprintf("Chapter %02d - %s  (page %d of %d)",
                        as.integer(ch), rows$bnf_chapter_name[1], pi, length(pages))
      p <- ggplot(dd, aes(month_date)) +
        geom_point(aes(y = observed), size = 0.35, colour = "grey55") +
        geom_line(aes(y = fitted), colour = "#B2182B", linewidth = 0.4) +
        facet_wrap(~ panel, scales = "free_y", ncol = ncol) +
        labs(x = NULL, y = NULL, title = NULL) +
        theme_bw(base_size = 7) +
        theme(strip.background = element_rect(fill = "grey92", colour = NA),
              strip.text = element_text(size = 6, face = "bold"), panel.grid.minor = element_blank(),
              plot.margin = margin(2, 2, 2, 2))
      # Recolour flagged strips, then draw into the occupied cells of the fixed
      # 4 x 4 page. Unused cells remain blank rather than being stretched.
      .draw_fixed_4x4_page(
        .style_strips(p, page_tier_map), ptitle, nrow(rows), first_page,
        "Items per 1000 patients"
      )
      first_page <- FALSE
    }
  }
}

render_all_series(class_monthly_elig, "bnf_class_code", "bnf_class_name",
                  chap_class, file.path(fig_dir, "S3 appendix_all_classes_by_chapter.pdf"),
                  tier_map = tier_class, per_page = 16L, ncol = 4L)
render_all_series(drug_monthly_elig, "bnf_drug_code", "bnf_drug_name",
                  chap_drug, file.path(fig_dir, "S4. appendix_all_drugs_by_chapter.pdf"),
                  tier_map = tier_drug, per_page = 16L, ncol = 4L)

## ---- APPENDIX figures: ALL fitted seasonal SHAPES, paginated by chapter ----
# Same tiered strip colouring, but showing the fitted seasonal curve (the
# multiplicative factor by calendar month, as in figure2_seasonal_shape.png)
# for every eligible series rather than observed+fitted counts. Fit is Poisson
# throughout regardless of a series' inferential route - this is a visual shape
# summary, not the inferential model, exactly as for the observed+fitted
# appendix. A flat/noisy curve for a non-significant series is itself
# informative (it shows there is nothing there), so all eligible series are
# included, not just significant ones.
FILL_MEAN_SHAPE <- "#08306B"; TXT_MEAN_SHAPE <- "white"   # dark blue
FILL_SIG_SHAPE  <- "#9ECAE1"; TXT_SIG_SHAPE  <- "grey15"   # pale blue

.series_shape <- function(mon, id_col, id_val) {
  ser <- mon |> filter(.data[[id_col]] == id_val) |> select(year_month, items)
  d   <- covar |> left_join(ser, by = "year_month") |> arrange(t)
  d$off <- d$offset_log_patient_days
  fit <- suppressWarnings(glm(.f_full, poisson, data = d))
  b <- coef(fit)
  s <- .seasonal_curve(b[["sin12"]], b[["cos12"]], b[["sin6"]], b[["cos6"]], 1:12)
  f <- exp(s)
  tibble(month = 1:12, factor = f / mean(f))
}

# shared per-page label/tier assembly (duplicate names disambiguated by code)
.page_labels <- function(rows, name_col, id_col, tier_map) {
  nm  <- rows[[name_col]]
  dup <- nm %in% nm[duplicated(nm)]
  panel_label <- ifelse(dup, sprintf("%s [%s]", nm, rows[[id_col]]), nm)
  tier <- setNames(tier_map[nm], panel_label)
  list(labels = panel_label, tier = tier[!vapply(tier, is.null, logical(1))])
}

render_all_shapes <- function(mon, id_col, name_col, chap_lookup, out_pdf,
                              tier_map = list(), per_page = 16L, ncol = 4L) {
  if (per_page != 16L || ncol != 4L) {
    stop("Paginated appendix PDFs use a fixed 4 x 4 layout (16 panels per page).")
  }
  ids <- chap_lookup |> arrange(bnf_chapter_code) |>
    left_join(distinct(mon, .data[[id_col]], .data[[name_col]]), by = id_col)
  pdf(out_pdf, width = 11.7, height = 8.3)   # A4 landscape
  on.exit(dev.off())
  first_page <- TRUE
  for (ch in unique(ids$bnf_chapter_code)) {
    chn <- ids |> filter(bnf_chapter_code == ch)
    pages <- split(seq_len(nrow(chn)), ceiling(seq_len(nrow(chn)) / per_page))
    for (pi in seq_along(pages)) {
      rows <- chn[pages[[pi]], ]
      pl <- .page_labels(rows, name_col, id_col, tier_map)
      rows$panel_label <- pl$labels
      
      dd <- bind_rows(lapply(seq_len(nrow(rows)), function(i) {
        .series_shape(mon, id_col, rows[[id_col]][i]) |>
          mutate(panel = rows$panel_label[i])
      }))
      dd$panel <- factor(dd$panel, levels = rows$panel_label)
      ptitle <- sprintf("Chapter %02d - %s  (page %d of %d) - fitted seasonal shape",
                        as.integer(ch), rows$bnf_chapter_name[1], pi, length(pages))
      p <- ggplot(dd, aes(month, factor)) +
        geom_hline(yintercept = 1, colour = "grey75", linewidth = 0.3) +
        geom_line(colour = "#2166AC", linewidth = 0.5) +
        facet_wrap(~ panel, ncol = ncol) +
        scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c("Jan", "Apr", "Jul", "Oct")) +
        labs(x = NULL, y = NULL, title = NULL) +
        theme_bw(base_size = 7) +
        theme(strip.background = element_rect(fill = "grey92", colour = NA),
              strip.text = element_text(size = 6, face = "bold"), panel.grid.minor = element_blank(),
              plot.margin = margin(2, 2, 2, 2))
      .draw_fixed_4x4_page(
        .style_strips(p, pl$tier), ptitle, nrow(rows), first_page,
        "Seasonal factor (relative to annual mean)"
      )
      first_page <- FALSE
    }
  }
}

render_all_shapes(class_monthly_elig, "bnf_class_code", "bnf_class_name", chap_class,
                  file.path(fig_dir, "S5. appendix_all_classes_seasonal_shape_by_chapter.pdf"),
                  tier_map = tier_class_shape, per_page = 16L, ncol = 4L)
render_all_shapes(drug_monthly_elig, "bnf_drug_code", "bnf_drug_name", chap_drug,
                  file.path(fig_dir, "S6. appendix_all_drugs_seasonal_shape_by_chapter.pdf"),
                  tier_map = tier_drug_shape, per_page = 16L, ncol = 4L)

## ---- report -----------------------------------------------------------------

cat(sprintf(paste0(
  "Publication outputs written to %s\n",
  "  Tables (tables/): table1_meaningful_classes, table_accounting,\n",
  "    appendixA1_all_classes, appendixA2_all_drugs, coverage, exclusions, crosswalk, routes\n",
  "  Figures (figures/): Figure1_seasonality_landscape, Figure2_exemplar_observed_fitted,\n",
  "    Figure3_exemplar_seasonal_shapes, FigureS2_retained_vs_nonmeaningful,\n",
  "    S3 appendix_all_classes_by_chapter.pdf, S4. appendix_all_drugs_by_chapter.pdf\n",
  "  Main-text meaningful classes: %d | eligible classes plotted: %d | eligible drugs plotted: %d\n"),
  res_dir, nrow(main_table_classes),
  n_distinct(class_monthly_elig$bnf_class_code), n_distinct(drug_monthly_elig$bnf_drug_code)))


### 13. Publication tables as Word documents (one .docx per table) -----------
# Reads the CSVs already written by Section 12 (no re-derivation), and stages
# a JSON manifest per table for a Node/docx-js script to render as a formatted,
# captioned Word document. Row shading mirrors the appendix-figure tiering:
# dark red = meaningful, pale red = significant but not meaningful, white =
# neither (kept for tables that carry a "Meaningful seasonality" column).

suppressMessages({library(jsonlite)})

word_dir   <- file.path(res_dir, "word")
stage_dir  <- file.path(res_dir, "word_stage")
for (d in c(word_dir, stage_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FILL_MEAN_HEX <- "8B0000"; TXT_MEAN_HEX <- "FFFFFF"
FILL_SIG_HEX  <- "E8A0A0"; TXT_SIG_HEX  <- "262626"

## table manifest: id, source csv, output filename, caption, notes, orientation
tbl_spec <- list(
  list(id = "table1", file = "table1_meaningful_classes.csv",
       out = "Table1_meaningful_classes.docx",
       caption = "Table 1. Drug classes with meaningful seasonal variation in community prescribing, England, January 2022 to December 2025.",
       notes = "Meaningful seasonality was defined as a peak-to-trough ratio whose lower 95% confidence bound was at least 1.10 and an STL seasonal strength of at least 0.50 (see Methods). Seasonal strength summarises how much of the non-trend variation in the series is explained by a recurring within-year cycle after flexible detrending (1.00 = fully explained). Adjusted p-values were derived from Benjamini-Hochberg correction across all screened drug classes.",
       orientation = "landscape", shade_col = NULL),
  list(id = "table2", file = "table2_meaningful_drugs.csv",
       out = "Table2_meaningful_drugs.docx",
       caption = "Table 2. Individual drugs with meaningful seasonal variation in community prescribing, England, January 2022 to December 2025.",
       notes = "Meaningful seasonality was defined as for Table 1. Drug q-values were derived from Benjamini-Hochberg correction across all eligible drugs. Parent-class significance is shown as descriptive context and was not used to select drugs for testing or reporting. Drug-level findings are secondary and exploratory.",
       orientation = "landscape", shade_col = NULL),
  list(id = "table3", file = "table_accounting.csv",
       out = "Table3_analytic_accounting.docx",
       caption = "Table 3. Analytic accounting of drug classes and individual drugs, from eligibility to meaningful seasonality.",
       notes = "Eligible series were prescribed in every month of the study window and reached at least 1,000 items dispensed nationally in each calendar year. Benjamini-Hochberg correction was applied separately across all eligible classes (primary analysis) and all eligible drugs (secondary exploratory analysis). False discovery rate control is claimed within each complete family, not across the hierarchy as a whole. Trivial-amplitude series are significant series with a peak-to-trough ratio under 1.05.",
       orientation = "portrait", shade_col = NULL),
  list(id = "s1", file = "appendixA1_all_classes.csv",
       out = "TableS1_all_classes.docx",
       caption = "Supplementary Table S1. Seasonality characteristics for drug classes with detected seasonality.",
       notes = "The table contains classes meeting the class-family 5% FDR threshold. Shading: dark red rows show meaningful seasonality (as defined for Table 1); pale red rows are statistically significant but did not meet the amplitude or seasonal-strength threshold. HAC = Newey-West heteroscedasticity- and autocorrelation-consistent standard errors; LRT = likelihood ratio test; STL = seasonal-trend decomposition using Loess.",
       orientation = "landscape", shade_col = "Meaningful seasonality",
       sig_col = NULL),
  list(id = "s2", file = "appendixA2_all_drugs.csv",
       out = "TableS2_all_drugs.docx",
       caption = "Supplementary Table S2. Seasonality characteristics for individual drugs with detected seasonality.",
       notes = "The table contains drugs meeting the complete drug-family 5% FDR threshold. Shading as for Table S1. Drug-level findings are secondary and exploratory. Drug q-values were derived from Benjamini-Hochberg correction across all eligible drugs, regardless of parent-class significance; parent-class significance is descriptive context only.",
       orientation = "landscape", shade_col = "Meaningful seasonality",
       sig_col = NULL),
  list(id = "s3", file = "appendix_coverage_summary.csv",
       out = "TableS3_coverage_summary.docx",
       caption = "Supplementary Table S3. Coverage of eligible drug classes and individual drugs.",
       notes = "Excluded item share is the proportion of total prescribing volume accounted for by series that failed the eligibility criteria (see Table 3 notes).",
       orientation = "portrait", shade_col = NULL),
  list(id = "s4", file = "appendix_exclusions_class.csv",
       out = "TableS4_excluded_classes.docx",
       caption = "Supplementary Table S4. Drug classes excluded from analysis, with reasons.",
       notes = "'Not every month' indicates the class was not prescribed in every month of the study window; '<1000 items/year' indicates at least one calendar year fell below the national volume threshold.",
       orientation = "landscape", shade_col = NULL),
  list(id = "s5", file = "appendix_exclusions_drug.csv",
       out = "TableS5_excluded_drugs.docx",
       caption = "Supplementary Table S5. Individual drugs excluded from analysis, with reasons.",
       notes = "Exclusion reasons as for Table S4, applied at the individual-drug (chemical substance) level.",
       orientation = "landscape", shade_col = NULL),
  list(id = "s6", file = "appendix_recode_crosswalk.csv",
       out = "TableS6_recode_crosswalk.docx",
       caption = "Supplementary Table S6. BNF code reconciliation crosswalk for products reclassified during the study period.",
       notes = "Products whose BNF code changed mid-window (identified by matching chemical-substance names across non-overlapping time periods and confirmed by inspection) were reconciled to a single series under the code shown in 'to_code'; 'from_code' gives the retired code.",
       orientation = "portrait", shade_col = NULL),
  list(id = "s7", file = "appendix_route_proportions.csv",
       out = "TableS7_route_proportions.docx",
       caption = "Supplementary Table S7. Distribution and inferential method used for each analysed series, by class and drug level.",
       notes = "Route was selected per series from residual diagnostics (autocorrelation first): Poisson with Newey-West HAC standard errors and a joint Wald test where residual autocorrelation was present; negative binomial with a likelihood ratio test where overdispersion was present without autocorrelation; Poisson with a likelihood ratio test otherwise.",
       orientation = "portrait", shade_col = NULL)
)

## sanitise text for Word/LibreOffice rendering: some Unicode symbols (e.g. the
## en-dash and >= sign baked into earlier CSVs) do not reliably render through
## the docx-js -> LibreOffice PDF path, so swap them for plain ASCII regardless
## of which upstream step introduced them.
.ascii_safe <- function(x) {
  x |>
    gsub("\u2265", ">=", x = _, useBytes = TRUE, fixed = TRUE) |>
    gsub("\u2264", "<=", x = _, useBytes = TRUE, fixed = TRUE) |>
    gsub("\u2013", "-",  x = _, useBytes = TRUE, fixed = TRUE) |>
    gsub("\u2014", "-",  x = _, useBytes = TRUE, fixed = TRUE) |>
    gsub("\u00d7", "x",  x = _, useBytes = TRUE, fixed = TRUE) |>
    gsub("\u2018", "'",  x = _, useBytes = TRUE, fixed = TRUE) |>
    gsub("\u2019", "'",  x = _, useBytes = TRUE, fixed = TRUE) |>
    gsub("\u201c", '"',  x = _, useBytes = TRUE, fixed = TRUE) |>
    gsub("\u201d", '"',  x = _, useBytes = TRUE, fixed = TRUE)
}

## build one JSON manifest entry per table that has a source file present
manifest <- list()
for (spec in tbl_spec) {
  f <- file.path(tab_dir, spec$file)
  if (!file.exists(f)) { message("skip (not found): ", spec$file); next }
  # read as character throughout: these CSVs already carry publication-formatted
  # values (dotted codes, formatted p-values etc.) - character preserves them
  # exactly rather than risking any column being re-coerced to numeric on read.
  dt <- fread(f, colClasses = "character")
  if (nrow(dt) == 0) { message("skip (empty): ", spec$file); next }
  
  # appendixA1/A2 (results_class/results_drug) contain only statistically
  # significant series by construction, so shading is simply: meaningful ->
  # dark red, otherwise -> pale red (significant but not meaningful). Tables
  # with no shade_col (e.g. Table 1/2, already pre-filtered to meaningful) are
  # left unshaded.
  row_shade <- rep("none", nrow(dt))
  if (!is.null(spec$shade_col) && spec$shade_col %in% names(dt)) {
    row_shade <- ifelse(dt[[spec$shade_col]] == "Yes", "meaningful", "significant")
  }
  
  out <- list(
    id = jsonlite::unbox(spec$id), out_file = jsonlite::unbox(spec$out),
    caption = jsonlite::unbox(.ascii_safe(spec$caption)), notes = jsonlite::unbox(.ascii_safe(spec$notes)),
    orientation = jsonlite::unbox(spec$orientation),
    headers = I(.ascii_safe(names(dt))),
    rows = unname(lapply(seq_len(nrow(dt)), function(i) .ascii_safe(as.character(unlist(dt[i, ]))))),
    row_shade = I(row_shade),
    fill_mean = jsonlite::unbox(FILL_MEAN_HEX), txt_mean = jsonlite::unbox(TXT_MEAN_HEX),
    fill_sig  = jsonlite::unbox(FILL_SIG_HEX),  txt_sig  = jsonlite::unbox(TXT_SIG_HEX))
  write_json(out, file.path(stage_dir, paste0(spec$id, ".json")),
             auto_unbox = FALSE, null = "null")
  manifest[[length(manifest) + 1]] <- list(id = jsonlite::unbox(spec$id),
                                           out_file = jsonlite::unbox(spec$out))
}
write_json(manifest, file.path(stage_dir, "manifest.json"), auto_unbox = FALSE)

cat(sprintf("Staged %d table(s) for Word export to %s\n", length(manifest), stage_dir))
cat("Run the Node renderer next: node render_word_tables.js '", stage_dir, "' '", word_dir, "'\n", sep = "")

### 14. Window-justification figure: where does the pandemic disruption sit? -----
# Descriptive only (no model is refitted). Reads the 2021 EPD and 2021 list-size
# data, aggregates 2021 to national monthly totals by BNF class, and plots the
# raw monthly prescribing rate across 2021-2025 for a set of representative
# classes, with 2021 (the excluded year) shaded. This shows empirically that the
# most severe disruption is confined to 2021 and that a regular annual pattern is
# re-established from 2022, supporting the choice of analytic window. It does not
# alter any object produced by Sections 2-13.
#
# Requires Sections 3-11 objects: covar, class_monthly, class_monthly_elig,
# results_class. Reads 2021 raw data from the same drives as Section 2. Caches
# the (small) 2021 class aggregates so the 1 GB monthly files are read only once.

stopifnot(exists("out_dir"), exists("covar"), exists("results_class"),
          exists("class_monthly_elig"))
if (!exists("data_dir")) data_dir <- out_dir
if (!exists("class_monthly"))
  class_monthly <- readRDS(file.path(data_dir, "class_monthly.rds"))

# --- editable: representative seasonal classes (matched by exact BNF class name).
#     A high-volume non-significant class is appended automatically as a flat
#     control, so leave room for it (aim for <= 5 named here).
exemplar_names <- c(
  "Penicillins",                 # winter antibiotic (COVID-suppressed 2020-21)
  "Macrolides",                  # winter antibiotic
  "Antihistamines",              # summer allergy (opposite-phase cycle)
  "Cough suppressants",          # winter respiratory
  "Sunscreening preparations"    # summer dermatological
)

## 14a. Aggregate the 2021 EPD to national monthly totals per BNF class ---------
# Cached per month (all chapter 1-14 classes, so the cache is reusable even if
# the exemplar set is changed later); the heavy CSV is read only when uncached.

cache_dir_2021 <- file.path(data_dir, "cache_2021")
dir.create(cache_dir_2021, showWarnings = FALSE, recursive = TRUE)

epd_2021 <- epd_2021_files

for (f in epd_2021) {
  ym      <- str_extract(basename(f), "\\d{6}")
  cache_f <- file.path(cache_dir_2021, paste0("data_byclass_2021_", ym, ".csv"))
  if (file.exists(cache_f)) { message("Skipping 2021 ", ym, " (cached)"); next }
  message("Aggregating 2021 EPD ", ym)
  
  # read only the presentation code and item count (schema-robust)
  hdr      <- names(read_epd_archive(f, nrows = 0))
  code_col <- if ("BNF_CODE" %in% hdr) "BNF_CODE"
  else if ("BNF_PRESENTATION_CODE" %in% hdr) "BNF_PRESENTATION_CODE"
  else stop("no recognised presentation-code column in ", basename(f))
  d <- read_epd_archive(f, select = c(code_col, "ITEMS"))
  setnames(d, c(code_col, "ITEMS"), c("bnf_code", "items"))
  
  d[, bnf_chapter_code := substr(bnf_code, 1, 2)]
  d <- d[!is.na(bnf_chapter_code) & bnf_chapter_code %in% valid_chapters]
  d[, bnf_class_code := substr(bnf_code, 1, 6)]
  agg <- d[, .(items = sum(items)), by = bnf_class_code]
  agg[, year_month := as.integer(ym)]
  
  fwrite(agg[, .(year_month, bnf_class_code, items)], cache_f)
  rm(d, agg); gc()                                                # free before next file
}

cache_files <- list.files(cache_dir_2021,
                          pattern = "^data_byclass_2021_\\d{6}\\.csv$", full.names = TRUE)
class_2021 <- rbindlist(lapply(cache_files, function(f) fread(
  f, colClasses = list(character = c(
    "bnf_chapter_code", "bnf_section_code", "bnf_class_code"
  ))
)), use.names = TRUE)
n_months_2021 <- uniqueN(class_2021$year_month)
if (n_months_2021 < 12L)
  warning(sprintf("Only %d of 12 months of 2021 EPD are present.", n_months_2021))

## 14b. National monthly list size for 2021 (quarterly -> monthly by LOCF) ------

ls_2021_files <- list_size_2021_files

ls2021_totals <- rbindlist(lapply(ls_2021_files, function(f) {
  ym <- str_extract(basename(f), "\\d{6}")
  d  <- fread(f, select = "NUMBER_OF_PATIENTS")
  data.table(year_month = as.integer(ym), list_size = sum(d$NUMBER_OF_PATIENTS, na.rm = TRUE))
}))

# complete 2021 monthly grid, carry forward, back-fill any leading gap
ls2021 <- data.table(date = seq(window_start, window_end, by = "month")) |>
  mutate(year_month = as.integer(format(date, "%Y%m"))) |>
  left_join(ls2021_totals, by = "year_month") |>
  arrange(date) |>
  tidyr::fill(list_size, .direction = "downup") |>
  transmute(year_month, list_size)

# combine with the 2022-2025 denominator already used in the analysis (covar)
ls_all <- bind_rows(ls2021, covar |> transmute(year_month, list_size)) |>
  arrange(year_month)

## 14c. Choose exemplar classes (named seasonal + one flat control) -------------

ex <- results_class |>
  distinct(bnf_class_code, bnf_class_name) |>
  filter(bnf_class_name %in% exemplar_names)
missing_names <- setdiff(exemplar_names, ex$bnf_class_name)
if (length(missing_names))
  warning("exemplar class name(s) not found and skipped: ",
          paste(missing_names, collapse = ", "),
          " - check spelling against results_class$bnf_class_name.")

# flat control: highest-volume eligible class that was NOT significant
flat <- class_monthly_elig |>
  filter(!bnf_class_code %in% results_class$bnf_class_code) |>
  group_by(bnf_class_code, bnf_class_name) |>
  summarise(items = sum(items), .groups = "drop") |>
  arrange(desc(items)) |>
  slice_head(n = 1) |>
  select(bnf_class_code, bnf_class_name)

exemplars <- bind_rows(ex, flat) |> distinct(bnf_class_code, bnf_class_name)
if (nrow(exemplars) < 1L) stop("No exemplar classes resolved.")

# facet order: named seasonal classes by amplitude (largest first), control last
ex_ord <- ex |>
  left_join(results_class |> select(bnf_class_code, peak_trough_ratio), by = "bnf_class_code") |>
  arrange(desc(peak_trough_ratio))
facet_levels <- unique(c(ex_ord$bnf_class_name, flat$bnf_class_name))

## 14d. Build the 2021-2025 rate series and plot --------------------------------

items_all <- bind_rows(
  class_2021    |> transmute(year_month, bnf_class_code = as.character(bnf_class_code), items),
  class_monthly |> transmute(year_month, bnf_class_code = as.character(bnf_class_code), items))

plot_df <- exemplars |>
  inner_join(items_all, by = "bnf_class_code") |>
  left_join(ls_all, by = "year_month") |>
  mutate(date = as.Date(paste0(year_month, "01"), "%Y%m%d"),
         rate = items / list_size * 1000,
         class = factor(bnf_class_name, levels = facet_levels)) |>
  arrange(class, date)

res_dir <- file.path(data_dir, "results")
fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!exists("theme_pub")) {
  theme_pub <- theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold", size = 8),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
}

p_win <- ggplot(plot_df, aes(date, rate)) +
  annotate("rect", xmin = as.Date("2020-12-16"), xmax = as.Date("2021-12-16"),
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_vline(xintercept = study_start, linetype = 2,
             colour = "grey40", linewidth = 0.4) +
  geom_line(colour = "#B2182B", linewidth = 0.5) +
  geom_point(size = 0.7, colour = "grey30") +
  facet_wrap(~ class, scales = "free_y", ncol = 2) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = "Prescription items per 1000 registered patients",
       title = "Monthly prescribing, 2021-2025, for representative drug classes",
       subtitle = "Shaded year (2021) is excluded from the analytic window; dashed line marks the January 2022 start",
       caption = paste0("Raw monthly rates (no model fitted). The most severe pandemic disruption is confined to 2021; ",
                        "a regular annual pattern is re-established from 2022 onward.")) +
  theme_pub +
  theme(plot.caption = element_text(size = 8, hjust = 0))

ggsave(file.path(fig_dir, "figure_window_justification_2021_2025.png"),
       p_win, width = 10, height = 8, dpi = 300)
fwrite(plot_df |> select(bnf_class_code, bnf_class_name, year_month, items, list_size, rate),
       file.path(res_dir, "tables", "window_justification_rates.csv"))

cat(sprintf(paste0(
  "Section 14 complete (window-justification figure).\n",
  "  2021 EPD months aggregated: %d/12\n",
  "  Exemplar classes plotted:   %s\n",
  "  Figure: %s\n"),
  n_months_2021,
  paste(facet_levels, collapse = "; "),
  file.path(fig_dir, "figure_window_justification_2021_2025.png")))


### 14e. Systematic 2021-anomaly scan across all eligible drug classes ----------
# Diagnostic only. For every eligible class, quantify how far 2021 departs from
# the pattern established over 2022-2025, to locate where pandemic residue sits.
#   level: 2021 mean log-rate vs a LINEAR back-extrapolation of the 2022-2025
#          trend (so a merely growing/shrinking class is not flagged); reported
#          as % deviation and in residual-SD units.
#   shape: correlation of 2021's within-year (mean-centred) profile with the mean
#          2022-2025 profile, judged against the 2022-2025 year-to-year baseline.
# Reuses the 14a all-class 2021 cache and the 14b denominator - no new data.

stopifnot(exists("ls_all"), exists("res_dir"))          # from 14b / 14d
if (!exists("class_2021")) {
  cf <- list.files(cache_dir_2021, pattern = "^data_byclass_2021_\\d{6}\\.csv$", full.names = TRUE)
  if (!length(cf)) stop("No 2021 class cache found - run Section 14a first.")
  class_2021 <- rbindlist(lapply(cf, function(f) fread(
    f, colClasses = list(character = c(
      "bnf_chapter_code", "bnf_section_code", "bnf_class_code"
    ))
  )), use.names = TRUE)
}
if (!exists("class_monthly"))
  class_monthly <- readRDS(file.path(data_dir, "class_monthly.rds"))
if (!exists("class_monthly_elig"))
  class_monthly_elig <- readRDS(file.path(data_dir, "class_monthly_eligible.rds"))

cm  <- as.data.table(class_monthly)
cme <- as.data.table(class_monthly_elig)
elig_codes <- unique(cme$bnf_class_code)
meta_lu <- unique(cm[, .(bnf_class_code, bnf_class_name, bnf_chapter_code, bnf_chapter_name)])

# complete 2021-2025 monthly panel (eligible classes); 2021 gaps -> 0 items
grid21 <- CJ(bnf_class_code = elig_codes, year_month = 202101:202112)
it21 <- merge(grid21, class_2021[, .(bnf_class_code, year_month, items)],
              by = c("bnf_class_code", "year_month"), all.x = TRUE)
it21[is.na(items), items := 0]
itrest <- cm[bnf_class_code %in% elig_codes, .(bnf_class_code, year_month, items)]
panel <- rbind(it21, itrest)
panel <- merge(panel, as.data.table(ls_all)[, .(year_month, list_size)], by = "year_month")
panel[, `:=`(year = year_month %/% 100L, month = year_month %% 100L,
             lr = log(items + 0.5) - log(list_size),
             rate = items / list_size * 1000)]
setorder(panel, bnf_class_code, year_month)

idx <- function(ym) (ym %/% 100L - 2022L) * 12L + (ym %% 100L)   # 2022-01 -> 1

scan_one <- function(d) {
  d <- d[order(year_month)]
  rest <- d[year >= 2022]; y21 <- d[year == 2021]
  if (nrow(rest) < 24 || nrow(y21) < 6) return(NULL)
  rt <- idx(rest$year_month); yt <- idx(y21$year_month)
  fit <- lm(lr ~ rt, data = data.frame(lr = rest$lr, rt = rt))
  pred21 <- predict(fit, newdata = data.frame(rt = yt))
  rsd <- sd(residuals(fit))
  lev_log <- mean(y21$lr) - mean(pred21)
  prof <- function(sub){ v <- sub$lr[order(sub$month)]; if (length(v) != 12) return(rep(NA, 12)); v - mean(v) }
  ry <- 2022:2025
  P   <- sapply(ry, function(y) prof(d[year == y]))
  ref <- rowMeans(P)
  scor <- suppressWarnings(stats::cor(prof(y21), ref))
  bcor <- mean(sapply(seq_along(ry), function(i) suppressWarnings(stats::cor(P[, i], ref))))
  data.table(level_dev_pct = 100 * (exp(lev_log) - 1),
             anomaly_sd    = if (rsd > 0) lev_log / rsd else NA_real_,
             shape_cor_21  = scor, base_shape_cor = bcor,
             shape_gap     = bcor - scor)
}

scan <- panel[, scan_one(.SD), by = bnf_class_code]
scan <- merge(scan, meta_lu, by = "bnf_class_code")
scan[, direction   := fifelse(level_dev_pct < 0, "2021 lower", "2021 higher")]
scan[, flag_level  := !is.na(anomaly_sd) & abs(anomaly_sd) >= 2]
scan[, flag_shape  := !is.na(shape_gap)  & shape_gap >= 0.30]
scan[, flagged     := flag_level | flag_shape]
setorder(scan, level_dev_pct)
fwrite(scan, file.path(res_dir, "tables", "window_2021_anomaly_scan.csv"))

## summary ---------------------------------------------------------------------
cat(sprintf("\n2021 anomaly scan: %d eligible classes\n", nrow(scan)))
cat(sprintf("  flagged (level >=2 SD off trend OR shape gap >=0.30): %d\n", sum(scan$flagged)))
cat(sprintf("    of which 2021 lower: %d | 2021 higher: %d\n",
            sum(scan$flagged & scan$level_dev_pct < 0),
            sum(scan$flagged & scan$level_dev_pct >= 0)))
cat("  flagged by BNF chapter:\n")
print(scan[flagged == TRUE, .N, by = bnf_chapter_name][order(-N)])
cat("\n  Most-below-trend classes (2021 vs 2022-2025):\n")
print(scan[order(level_dev_pct)][1:12,
                                 .(bnf_class_name, level_dev_pct = round(level_dev_pct, 0),
                                   anomaly_sd = round(anomaly_sd, 1), shape_gap = round(shape_gap, 2))])
cat("\n  Most shape-disrupted classes:\n")
print(scan[order(-shape_gap)][1:12,
                              .(bnf_class_name, shape_gap = round(shape_gap, 2),
                                level_dev_pct = round(level_dev_pct, 0))])

## plot the 9 classes whose 2021 level departs most (either direction) ----------
sel <- scan[order(-abs(level_dev_pct))][1:9]
lab <- setNames(sprintf("%s (2021 %+.0f%%)", sel$bnf_class_name, sel$level_dev_pct),
                sel$bnf_class_code)
pd  <- panel[bnf_class_code %in% sel$bnf_class_code]
pd[, date  := as.Date(paste0(year_month, "01"), "%Y%m%d")]
pd[, class := factor(lab[as.character(bnf_class_code)],
                     levels = lab[as.character(sel$bnf_class_code)])]

p_scan <- ggplot(pd, aes(date, rate)) +
  annotate("rect", xmin = as.Date("2020-12-16"), xmax = as.Date("2021-12-16"),
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_vline(xintercept = study_start, linetype = 2, colour = "grey40", linewidth = 0.4) +
  geom_line(colour = "#B2182B", linewidth = 0.5) + geom_point(size = 0.6, colour = "grey30") +
  facet_wrap(~ class, scales = "free_y", ncol = 3) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = "Prescription items per 1000 registered patients",
       title = "Eligible classes whose 2021 prescribing departs most from the 2022-2025 pattern",
       subtitle = "Shaded year (2021) is excluded; percentage = 2021 level vs back-extrapolated 2022-2025 trend") +
  theme_pub + theme(strip.text = element_text(size = 7))
ggsave(file.path(res_dir, "figures", "figure_window_2021_anomaly_top9.png"),
       p_scan, width = 11, height = 8, dpi = 300)

cat(sprintf("\nSaved: window_2021_anomaly_scan.csv and figure_window_2021_anomaly_top9.png\n"))

### 14f. Plot every eligible class, 2021-2025, for visual window selection ------
# Diagnostic. One paginated PDF (A4 landscape) with every eligible class's raw
# monthly rate across 2021-2025, 2021 shaded, faceted and paginated by BNF
# chapter. Facet strips are recoloured by the 14e flags so disrupted classes are
# findable at a glance:
#   red   = flagged on LEVEL (2021 >= 2 residual-SD off the 2022-2025 trend)
#   blue  = flagged on SHAPE (within-year profile gap >= 0.30) only
#   grey  = not flagged
# Each strip shows the class name plus its signed 2021 level deviation.
# Requires 14e objects: panel, scan. Reuses .style_strips / .recolour_strip /
# .strip_label from Section 12 if present (redefined here if not).

stopifnot(exists("panel"), exists("scan"), exists("res_dir"))
suppressMessages(library(ggplot2))

# --- strip-recolouring helpers (reuse Section 12's if already defined) --------
if (!exists(".strip_label")) {
  .strip_label <- function(sg) {
    gt <- sg$grobs[[1]]
    for (ch in gt$children) if (!is.null(ch$children))
      for (tx in ch$children) if (!is.null(tx$label)) return(tx$label)
    NA_character_
  }
}
if (!exists(".recolour_strip")) {
  .recolour_strip <- function(sg, fill, textcol) {
    gt <- sg$grobs[[1]]
    for (i in seq_along(gt$children)) {
      ch <- gt$children[[i]]
      if (inherits(ch, "rect")) {
        if (is.null(ch$gp)) ch$gp <- grid::gpar(); ch$gp$fill <- fill; ch$gp$col <- NA
      } else if (!is.null(ch$children)) {
        for (j in seq_along(ch$children))
          if (inherits(ch$children[[j]], "text")) {
            if (is.null(ch$children[[j]]$gp)) ch$children[[j]]$gp <- grid::gpar()
            ch$children[[j]]$gp$col <- textcol
          }
      }
      gt$children[[i]] <- ch
    }
    sg$grobs[[1]] <- gt; sg
  }
}
if (!exists(".style_strips")) {
  .style_strips <- function(p, tier_map) {
    g <- ggplotGrob(p)
    for (si in which(grepl("^strip", g$layout$name))) {
      lab <- .strip_label(g$grobs[[si]])
      if (!is.na(lab) && lab %in% names(tier_map))
        g$grobs[[si]] <- .recolour_strip(g$grobs[[si]], tier_map[[lab]][1], tier_map[[lab]][2])
    }
    g
  }
}

# tier colours
FILL_LEVEL <- "#B2182B"; TXT_LEVEL <- "white"    # level-disrupted (strongest signal)
FILL_SHAPE <- "#2166AC"; TXT_SHAPE <- "white"    # shape-disrupted only
FILL_NONE  <- "grey92";  TXT_NONE  <- "grey15"   # not flagged

pd <- copy(panel)
pd[, date := as.Date(paste0(year_month, "01"), "%Y%m%d")]
pd[, rate := items / list_size * 1000]

# panel label per class: name + signed 2021 level deviation (disambiguate any
# duplicate names by appending the code, so facet levels stay unique per page)
sc <- as.data.table(scan)[, .(bnf_class_code, bnf_class_name, bnf_chapter_code,
                              bnf_chapter_name, level_dev_pct, flag_level, flag_shape)]
sc[, base_label := sprintf("%s (2021 %+.0f%%)", bnf_class_name, round(level_dev_pct))]
sc[, tier := fifelse(flag_level, "level", fifelse(flag_shape, "shape", "none"))]

order_key <- function(dt) dt[order(bnf_chapter_code, -abs(level_dev_pct))]  # worst first within chapter
sc <- order_key(sc)

out_pdf <- file.path(res_dir, "figures", "S1. figure_window_all_classes_2021_2025.pdf")
per_page <- 16L; ncol <- 4L

if (!exists("theme_pub")) {
  theme_pub <- theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          panel.grid.minor = element_blank())
}

pdf(out_pdf, width = 11.7, height = 8.3)   # A4 landscape
first_page <- TRUE
for (ch in unique(sc$bnf_chapter_code)) {
  chn <- sc[bnf_chapter_code == ch]
  pages <- split(seq_len(nrow(chn)), ceiling(seq_len(nrow(chn)) / per_page))
  for (pi in seq_along(pages)) {
    rows <- chn[pages[[pi]]]
    # unique facet labels on this page
    lab <- rows$base_label
    dup <- lab %in% lab[duplicated(lab)]
    lab <- ifelse(dup, sprintf("%s [%s]", lab, rows$bnf_class_code), lab)
    rows[, panel_label := lab]
    
    tier_map <- setNames(
      lapply(rows$tier, function(t)
        switch(t, level = c(FILL_LEVEL, TXT_LEVEL),
               shape = c(FILL_SHAPE, TXT_SHAPE),
               c(FILL_NONE,  TXT_NONE))),
      rows$panel_label)
    
    dd <- merge(pd[bnf_class_code %in% rows$bnf_class_code],
                rows[, .(bnf_class_code, panel_label)], by = "bnf_class_code")
    dd[, panel := factor(panel_label, levels = rows$panel_label)]
    
    ptitle <- sprintf("Chapter %02d - %s  (page %d of %d)  |  red = 2021 level-disrupted, blue = shape-disrupted",
                      as.integer(ch), rows$bnf_chapter_name[1], pi, length(pages))
    p <- ggplot(dd, aes(date, rate)) +
      annotate("rect", xmin = as.Date("2020-12-16"), xmax = as.Date("2021-12-16"),
               ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
      geom_vline(xintercept = study_start, linetype = 2, colour = "grey45", linewidth = 0.3) +
      geom_line(colour = "grey20", linewidth = 0.35) +
      geom_point(size = 0.3, colour = "grey35") +
      facet_wrap(~ panel, scales = "free_y", ncol = ncol) +
      scale_x_date(date_breaks = "1 year", date_labels = "'%y") +
      labs(x = NULL, y = NULL, title = NULL) +
      theme_bw(base_size = 7) +
      theme(strip.background = element_rect(fill = "grey92", colour = NA),
            strip.text = element_text(size = 6, face = "bold"),
            panel.grid.minor = element_blank(),
            plot.margin = margin(2, 2, 2, 2))
    
    .draw_fixed_4x4_page(
      .style_strips(p, tier_map), ptitle, nrow(rows), first_page,
      "Items per 1000 patients"
    )
    first_page <- FALSE
  }
}
dev.off()

cat(sprintf(paste0(
  "Section 14f complete.\n",
  "  All %d eligible classes plotted 2021-2025 -> %s\n",
  "  Strip colour: red = level-disrupted (%d), blue = shape-only (%d), grey = stable (%d)\n"),
  nrow(sc), out_pdf,
  sum(sc$tier == "level"), sum(sc$tier == "shape"), sum(sc$tier == "none")))

}, warning = function(w) {
  warning_text <- conditionMessage(w)
  warnings_seen <<- c(warnings_seen, warning_text)
  message("WARNING: ", warning_text)
  invokeRestart("muffleWarning")
})

invisible(TRUE)
}

.run_analysis()
