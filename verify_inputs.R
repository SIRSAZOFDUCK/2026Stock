#!/usr/bin/env Rscript

# Verify that independently acquired source files match the frozen study
# manifest before starting the expensive import and modelling workflow.

suppressPackageStartupMessages(library(dplyr))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

manifest_path <- file.path(repo_dir, "reproducibility", "input_manifest.csv")
if (!file.exists(manifest_path)) stop("Missing source manifest: ", manifest_path)

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("source_id", "file_size_bytes", "sha256")
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop("Source manifest lacks required columns: ", paste(missing_columns, collapse = ", "))
}

bundled_only <- "--bundled-only" %in% commandArgs(trailingOnly = TRUE)
if (bundled_only) {
  manifest <- manifest %>%
    filter(grepl("^(BNF_REFERENCE|BANK_HOLIDAYS)", source_id)) # retain bundled references
}

epd_dir <- Sys.getenv(
  "STOCK2026_EPD_DIR",
  unset = file.path(repo_dir, "data", "epd")
)
list_size_dir <- Sys.getenv(
  "STOCK2026_LIST_SIZE_DIR",
  unset = file.path(repo_dir, "data", "list_size")
)
bnf_path <- Sys.getenv(
  "STOCK2026_BNF_PATH",
  unset = file.path(repo_dir, "bnf_code_current_202505_version_88.csv")
)
calendar_path <- file.path(
  repo_dir, "reproducibility", "source_snapshots",
  "bank-holidays_2026-08-05.json"
)

resolve_path <- function(source_id) {
  if (grepl("^EPD_[0-9]{6}$", source_id)) {
    return(file.path(epd_dir, paste0("EPD_SNOMED_", sub("^EPD_", "", source_id), ".ZIP")))
  }
  if (grepl("^LIST_SIZE_[0-9]{6}$", source_id)) {
    return(file.path(
      list_size_dir,
      paste0("gp-reg-pat-prac-all_", sub("^LIST_SIZE_", "", source_id), ".csv")
    ))
  }
  if (grepl("^BNF_REFERENCE", source_id)) return(bnf_path)
  if (grepl("^BANK_HOLIDAYS", source_id)) return(calendar_path)
  stop("No release path rule for source_id: ", source_id)
}

sha256_file <- function(path) {
  shasum <- Sys.which("shasum")
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(shasum)) {
    output <- system2(shasum, c("-a", "256", shQuote(path)), stdout = TRUE)
  } else if (nzchar(sha256sum)) {
    output <- system2(sha256sum, shQuote(path), stdout = TRUE)
  } else {
    stop("Install either shasum or sha256sum to verify source hashes.")
  }
  sub("[[:space:]].*$", "", output[1])
}

paths <- vapply(manifest$source_id, resolve_path, character(1))
exists <- file.exists(paths)
observed_size <- rep(NA_real_, length(paths))
observed_hash <- rep(NA_character_, length(paths))

if (any(exists)) {
  observed_size[exists] <- unname(file.info(paths[exists])$size)
}

size_matches <- exists & observed_size == as.numeric(manifest$file_size_bytes)
for (index in which(size_matches)) {
  message(sprintf("Hashing %d/%d: %s", index, nrow(manifest), basename(paths[index])))
  observed_hash[index] <- sha256_file(paths[index])
}
hash_matches <- size_matches & observed_hash == manifest$sha256

verification <- tibble(
  source_id = manifest$source_id,
  path = paths,
  exists = exists,
  expected_size_bytes = as.numeric(manifest$file_size_bytes),
  observed_size_bytes = observed_size,
  size_matches = size_matches,
  expected_sha256 = manifest$sha256,
  observed_sha256 = observed_hash,
  sha256_matches = hash_matches
)

failed <- verification %>%
  filter(!sha256_matches)                # retain missing or altered sources
if (nrow(failed)) {
  failed %>%
    select(                               # show fields needed to diagnose failures
      source_id, path, exists, expected_size_bytes,
      observed_size_bytes, size_matches, sha256_matches
    ) %>%
    print(n = Inf)
  stop(
    nrow(failed), " of ", nrow(verification),
    " required sources are missing or differ from the frozen manifest."
  )
}

cat(
  "PASS:", nrow(verification),
  if (bundled_only) "bundled reference files" else "required source files",
  "match the frozen sizes and SHA-256 hashes.\n"
)
