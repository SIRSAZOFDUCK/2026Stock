# Module 01: frozen-source import, aggregation, recoding, and quality control.
run_import_qc <- function(state) {
  # Evaluate in the coordinator context so downstream modules receive these objects.
  evalq({
## 01.1 Read the BNF lookup and monthly prescribing sources -------------------

### BNF class lookup

# May 2025 BNF hierarchy snapshot from the NHSBSA open-data release.

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


### English Prescribing Dataset
# Per-month aggregates and QC files are resumable checkpoints, avoiding repeated
# expansion of the 60 large frozen archives.

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

# Resolve all checkpoint and QC paths for one prescribing month.
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

  }, envir = state)

  if (isTRUE(state$partial_epd_run)) {
    message("Requested EPD checkpoint month(s) completed; stopping before cross-month Stage 2 checks.")
    return(invisible(FALSE))
  }

  evalq({

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


## 01.2 Read and interpolate registered-patient denominators ------------------

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


## 01.3 Assemble and validate monthly analysis panels -------------------------

# The coordinator defines the output directory and complete month sequence.

fails <- character(0)   # validation failures collected, reported together

# Convert accepted source month formats to the common integer YYYYMM key.
normalise_ym <- function(x, fname) {
  if (is.numeric(x)) return(as.integer(x))
  x <- trimws(as.character(x))
  if (all(grepl("^\\d{6}$", x))) return(as.integer(x))
  if (all(grepl("^\\d{4}-\\d{2}$", x)))
    return(as.integer(paste0(substr(x, 1, 4), substr(x, 6, 7))))
  stop(sprintf("%s: unrecognised year_month format (e.g. %s).",
               fname, paste(head(unique(x), 3), collapse = ", ")))
}

### Bind the monthly class and drug checkpoints

# Validate each checkpoint before combining months.
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


## 01.4 Build the shared time, holiday, and offset covariates -----------------

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

## 01.5 Record source flow and input-level QC ---------------------------------

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

## 01.6 Save validated panels before code reconciliation ----------------------

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

## 01.7 Reconcile BNF codes that changed during 2022–2025 ---------------------

# Normalise substance names so minor spelling drift maps to one candidate.
normalise_name <- function(x) {
  x |>
    str_to_lower() |>
    str_remove("\\s*\\([^()]*\\)\\s*$") |>   # drop a single trailing "(qualifier)"
    str_replace_all("[[:punct:]]", " ") |>
    str_replace_all("\\s+", " ") |>
    str_trim()
}

# Resolve the code and label columns for one reconciliation level.
.level_cols <- function(level) {
  switch(level,
         drug  = list(code = "bnf_drug_code",  name = "bnf_drug_name"),
         class = list(code = "bnf_class_code", name = "bnf_class_name"),
         stop("level must be 'drug' or 'class'"))
}

# Find same-name code histories that do not overlap in time.
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

# Apply an approved crosswalk and rebuild unique monthly series.
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

# Identify non-overlapping code histories that may describe one substance.
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

# Apply the prespecified, manually reviewed beclometasone code reconciliation.
cand_drug
cand_class

cand_drug |> filter(heals_to_full, combined_items >= 4000)   # keeps the only candidate we need

# Map the retired code to the retained current code.
cand_drug |> slice(1) |> pull(codes)     # the two full codes for the line below
drug_monthly |> filter(bnf_drug_code %in% c("0301011AB","0302000AA")) |>
  group_by(bnf_drug_code) |>
  summarise(first = min(year_month), last = max(year_month), items = sum(items), .groups = "drop")

# Rebuild the drug series after applying the approved crosswalk.
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

# Re-derive class totals from the reconciled drug panel.
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

## 01.8 Apply the final import and source-QC gate ------------------------------

# Represent one result in the final import/QC gate.
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


  }, envir = state)

  if (identical(state$run_stage, "stage2")) {
    evalq({
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
    }, envir = state)
    return(invisible(FALSE))
  }

  invisible(TRUE)
}
