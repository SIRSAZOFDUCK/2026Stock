#### SEASONALITY IN PRIMARY CARE PRESCRIBING IN ENGLAND PRIMARY CARE
#### Nadine Stock, Gillian Carr, Islam Omar, Saran Shantikumar, July 2026

### 1. Set up----------

# Clear environment
rm(list=ls())

# Required CRAN packages
required_packages <- c(
  "data.table",
  "dplyr",
  "stringr",
  "tidyr",
  "conflicted",
  "purrr",
  "lubridate",
  "readr",
  "tibble",
  "MASS",
  "AER",
  "sandwich",
  "lmtest",
  "car",
  "ggplot2",
  "ggrepel",
  "jsonlite"
)

# Identify packages that are not currently installed
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

# Install any missing packages
if (length(missing_packages) > 0L) {
  message(
    "Installing missing package(s): ",
    paste(missing_packages, collapse = ", ")
  )
  
  install.packages(
    missing_packages,
    dependencies = TRUE
  )
}

# Load all packages
invisible(
  lapply(
    required_packages,
    function(package) {
      suppressPackageStartupMessages(
        library(package, character.only = TRUE)
      )
    }
  )
)

message("All required packages are installed and loaded.")

# Resolve known function-name conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("first", "dplyr")

# Specify output directory
out_dir <- "OUTPUT DIRECTORY"

# Set working directory
setwd("WORKING DIRECTORY")


### 2. Get data -------------

## BNF drug class lookup table

# Retrived from https://opendata.nhsbsa.net/dataset/bnf-code-information-current-year (May 2025 version)

bnf_ref <- fread("bnf_code_current_202505_version_88.csv", colClasses = "character") %>% # Read in and force character to avoid dropping leading zeros
  select(BNF_SECTION,
         BNF_SECTION_CODE,
         BNF_PARAGRAPH,
         BNF_PARAGRAPH_CODE) %>% # Keep required columns
  rename(BNF_SECTION_NAME = BNF_SECTION,
         BNF_CLASS_NAME = BNF_PARAGRAPH,
         BNF_CLASS_CODE = BNF_PARAGRAPH_CODE) %>% # Rename
  distinct() %>% # Remove duplicate rows
  rename_with(tolower) # lower case column headers


## EPD data

# Set working directory
setwd("EPD DIRECTORY")


# Fetch names of EPD files
epd_files <- list.files(pattern = "\\.csv$")

for (i in seq_along(epd_files)) {
  
  setwd("EPD DIRECTORY")
  
  ym <- str_extract(epd_files[i], "\\d{6}")  # Extract year-month of dataset
  
  # Skip if already processed and output file is in output directory
  if (file.exists(file.path(out_dir, paste0("data_bydrug_", ym, ".csv"))) &&
      file.exists(file.path(out_dir, paste0("data_byclass_", ym, ".csv")))) {
    message("Skipping ", ym, " (already done)")
    next
  }
  
  message("Processing ", ym, " (", i, " of ", length(epd_files), ")")
  
  # Read in one EPD file
  data <- fread(epd_files[i])
  
  # Standardise drug-name column across schema versions (pre-/post-202503)
  if ("CHEMICAL_SUBSTANCE_BNF_DESCR" %in% names(data)) {
    setnames(data, "CHEMICAL_SUBSTANCE_BNF_DESCR", "BNF_DRUG_NAME")
  } else if ("BNF_CHEMICAL_SUBSTANCE" %in% names(data)) {
    setnames(data, "BNF_CHEMICAL_SUBSTANCE", "BNF_DRUG_NAME")
  } else {
    stop("No recognised drug-name column in ", epd_files[i])
  }
  
  # Standardise presentation-code column across schema versions (pre-/post-202503)
  if ("BNF_CODE" %in% names(data)) {
    # already canonical, nothing to do
  } else if ("BNF_PRESENTATION_CODE" %in% names(data)) {
    setnames(data, "BNF_PRESENTATION_CODE", "BNF_CODE")
  } else {
    stop("No recognised presentation-code column in ", epd_files[i])
  }
  
  # Filter and aggregate to drug level
  data1 <- data %>% select(YEAR_MONTH,
                           PRACTICE_CODE,
                           BNF_DRUG_NAME,
                           BNF_CODE,
                           BNF_CHAPTER_PLUS_CODE,
                           ITEMS) %>%
    mutate(
      BNF_CLASS_CODE   = substr(BNF_CODE, 1, 6),
      BNF_DRUG_CODE    = substr(BNF_CODE, 1, 9),
      BNF_CHAPTER_CODE = as.numeric(substr(BNF_CHAPTER_PLUS_CODE, 1, 2)),
      BNF_CHAPTER_NAME = substring(BNF_CHAPTER_PLUS_CODE, 5)
    ) %>%
    select(-c(BNF_CODE, BNF_CHAPTER_PLUS_CODE)) %>%
    filter(BNF_CHAPTER_CODE < 15) %>%
    mutate(BNF_CHAPTER_CODE = sprintf("%02d", BNF_CHAPTER_CODE)) %>%
    select(YEAR_MONTH, PRACTICE_CODE, BNF_CHAPTER_CODE, BNF_CHAPTER_NAME,
           BNF_CLASS_CODE, BNF_DRUG_CODE, BNF_DRUG_NAME, ITEMS) %>%
    rename_with(tolower) %>%
    left_join(bnf_ref, by = "bnf_class_code") %>%
    select(year_month, practice_code, bnf_chapter_code, bnf_chapter_name,
           bnf_section_code, bnf_section_name, bnf_class_code, bnf_class_name,
           bnf_drug_code, bnf_drug_name, items) %>%
    group_by(across(-items)) %>%
    summarise(items = sum(items), .groups = "drop") %>%
    group_by(year_month, bnf_chapter_code, bnf_chapter_name,
             bnf_section_code, bnf_section_name,
             bnf_class_code, bnf_class_name,
             bnf_drug_code, bnf_drug_name) %>%
    summarise(items = sum(items), .groups = "drop")
  
    # Save drug-level
  write.csv(data1, file = file.path(out_dir, paste0("data_bydrug_", ym, ".csv")), row.names = FALSE)
  
  # Make and save class-level version
  data2 <- data1 %>%
    group_by(year_month, bnf_chapter_code, bnf_chapter_name,
             bnf_section_code, bnf_section_name,
             bnf_class_code, bnf_class_name) %>%
    summarise(items = sum(items), .groups = "drop")
  
  write.csv(data2, file = file.path(out_dir, paste0("data_byclass_", ym, ".csv")), row.names = FALSE)
  
  
  # Remove unneeded data objects
  rm(data, data1, data2)
  gc()  # reclaim memory before the next 17M-row file
}


## List size data

# Set working directory
setwd("LIST SIZE DIRECTORY")

# Find all file names. List size data total for each GP practice, from https://digital.nhs.uk/data-and-information/publications/statistical/patients-registered-at-a-gp-practice
files <- list.files(pattern = "^gp-reg-pat-prac-all_\\d{6}\\.csv$", full.names = TRUE)

# Create empty list size table
list.size <- data.table()  # will grow one row of total list size per file

# Loop to process list size data to get total number of patients in England per month
for (f in files) {
  ym    <- str_extract(basename(f), "\\d{6}")            # yearmonth from filename
  d     <- fread(f, select = c("CODE", "NUMBER_OF_PATIENTS"))
  total <- sum(d$NUMBER_OF_PATIENTS, na.rm = TRUE)        # aggregate across all GP practices
  list.size <- rbind(list.size,
                       data.table(year_month = ym, list_size = total)) # Add to list size table
}

# Create a dataframe with monthly list size data, where the list size is carried forward for a three-month period

# Parse year_month to a first-of-month Date
list.size <- list.size %>%
  mutate(date = as.Date(paste0(year_month, "01"), format = "%Y%m%d"))

# Complete monthly grid from the earliest month to 2025-12
full <- tibble(date = seq(min(list.size$date), as.Date("2025-12-01"), by = "month"))

# Join and carry the last value forward into the empty months
list_size <- full %>%
  left_join(list.size %>% select(date, list_size), by = "date") %>%
  arrange(date) %>%
  fill(list_size, .direction = "down") %>%        # Last observation carried forward
  mutate(year_month = format(date, "%Y%m")) %>%
  select(year_month, list_size)

# Save
setwd(out_dir)
write.csv(list_size,"listsize.csv", row.names = FALSE)

# Remove initial list size data objects
rm(list = ls(pattern = "^ls\\d{6}$"))


### 3. Process data -----------

# Set working directory
data_dir <- out_dir
setwd(data_dir)

# Define timeframe expected
expected_ym <- as.integer(format(seq(as.Date("2022-01-01"),
                                     as.Date("2025-12-01"), by = "month"), "%Y%m"))

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
  cc  <- if (prefix == "data_bydrug") list(character = "bnf_drug_code") else NULL
  
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
    dt[, `:=`(bnf_chapter_code = as.integer(bnf_chapter_code),
              bnf_class_code   = as.integer(bnf_class_code),
              items            = as.integer(items))]
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

if (!all(class_monthly$bnf_chapter_code %in% 1:14) ||
    !all(drug_monthly$bnf_chapter_code %in% 1:14)) {
  fails <- c(fails, "chapters outside 1-14 present - import filter has slipped.")
}

## drug code prefix must reproduce the class code 
n_prefix_bad <- drug_monthly[as.integer(substr(bnf_drug_code, 1, 6)) != bnf_class_code, .N]
if (n_prefix_bad > 0)
  fails <- c(fails, sprintf("drug_monthly: %d rows where substr(drug code, 1, 6) != class code.",
                            n_prefix_bad))

## Validate the denominator

listsize <- fread(file.path(data_dir, "listsize.csv"),
                  select = c("year_month", "list_size"))

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


covar <- data.table(year_month = expected_ym)
covar[, month_date := as.Date(paste0(year_month, "01"), "%Y%m%d")]
covar[, t := seq_len(.N)]
covar[, days_in_month := as.integer(lubridate::days_in_month(month_date))]

# Fourier terms: annual (12-month) and first harmonic (6-month) periods
covar[, `:=`(sin12 = sin(2 * pi * t / 12), cos12 = cos(2 * pi * t / 12),
             sin6  = sin(4 * pi * t / 12), cos6  = cos(4 * pi * t / 12))]

# Spline basis stored as fixed columns so every series shares identical knots
S <- splines::ns(covar$t, df = 3)
covar[, `:=`(trend1 = S[, 1], trend2 = S[, 2], trend3 = S[, 3])]

# Denominator and combined offset (sum of logs; no large product formed)
covar <- merge(covar, listsize, by = "year_month", sort = TRUE)
covar[, offset_log_patient_days := log(list_size) + log(days_in_month)]

## Validate the frame 
stopifnot(
  nrow(covar) == 48,
  !anyNA(covar),
  sum(covar$days_in_month) == 1461,                      # incl. 29 Feb 2024
  covar[year_month == 202402, days_in_month] == 29,
  all(diff(covar$t) == 1)
)

## Save derived objects and quality control outcomes  ----------------------------------

saveRDS(class_monthly, file.path(data_dir, "class_monthly.rds"))
saveRDS(drug_monthly,  file.path(data_dir, "drug_monthly.rds"))
saveRDS(covar,         file.path(data_dir, "covariate_frame.rds"))
fwrite(qc_month,       file.path(data_dir, "input_qc_by_month.csv"))

cat(sprintf(paste0(
  "Complete.\n",
  "  Months:            %d (%d - %d)\n",
  "  Classes observed:  %d\n",
  "  Drugs observed:    %d\n",
  "  Total items:       %s\n",
  "  Listsize:          %s - %s (%d distinct values)\n",
  "  Saved: class_monthly.rds, drug_monthly.rds, covariate_frame.rds, input_qc_by_month.csv\n"),
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
  x %>%
    str_to_lower() %>%
    str_remove("\\s*\\([^()]*\\)\\s*$") %>%   # drop a single trailing "(qualifier)"
    str_replace_all("[[:punct:]]", " ") %>%
    str_replace_all("\\s+", " ") %>%
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
  x <- x %>% rename(.code = all_of(cols$code), .name = all_of(cols$name))
  
  code_activity <- x %>%
    group_by(.code) %>%
    summarise(name           = .name[which.max(year_month)],
              months_present = n_distinct(year_month),
              total_items    = sum(items),
              active         = list(sort(unique(year_month))),
              .groups = "drop") %>%
    mutate(nname = normalise_name(name))
  
  code_activity %>%
    group_by(nname) %>%
    filter(n() > 1) %>%
    summarise(
      n_codes            = n(),
      codes              = paste(.code, collapse = " + "),
      example_name       = dplyr::first(name),
      combined_items     = sum(total_items),
      best_single_months = max(months_present),
      union_months       = n_distinct(unlist(active)),
      max_simultaneous   = max(as.integer(table(unlist(active)))),
      .groups = "drop") %>%
    mutate(level         = level,
           clean_recode  = max_simultaneous == 1L,
           recovers      = union_months > best_single_months,
           heals_to_full = union_months == n_months_req) %>%
    filter(recovers) %>%
    arrange(desc(combined_items))
}

apply_recode_crosswalk <- function(x, crosswalk, level = c("drug", "class")) {
  level <- match.arg(level)
  stopifnot(all(c("from_code", "to_code") %in% names(crosswalk)))
  
  if (level == "class") {
    class_desc <- x %>%
      distinct(bnf_class_code, bnf_chapter_code, bnf_chapter_name,
               bnf_section_code, bnf_section_name, bnf_class_name)
    out <- x %>%
      left_join(crosswalk, by = c("bnf_class_code" = "from_code")) %>%
      mutate(bnf_class_code = coalesce(to_code, bnf_class_code)) %>%
      select(year_month, bnf_class_code, items) %>%
      group_by(year_month, bnf_class_code) %>%
      summarise(items = sum(items), .groups = "drop") %>%
      left_join(class_desc, by = "bnf_class_code")
    if (anyNA(out$bnf_class_name))
      warning("apply_recode_crosswalk: a canonical class code has no descriptor - check crosswalk targets.")
    return(out %>%
             select(year_month, bnf_chapter_code, bnf_chapter_name, bnf_section_code,
                    bnf_section_name, bnf_class_code, bnf_class_name, items))
  }
  
  ## drug level: remap the 9-char code, re-derive class from its prefix
  class_desc <- x %>%
    distinct(bnf_class_code, bnf_chapter_code, bnf_chapter_name,
             bnf_section_code, bnf_section_name, bnf_class_name)
  out <- x %>%
    left_join(crosswalk, by = c("bnf_drug_code" = "from_code")) %>%
    mutate(bnf_drug_code  = coalesce(to_code, bnf_drug_code),
           bnf_class_code = as.integer(substr(bnf_drug_code, 1, 6))) %>%
    group_by(bnf_drug_code) %>%
    mutate(bnf_drug_name = bnf_drug_name[which.max(year_month)]) %>%
    ungroup() %>%
    group_by(year_month, bnf_class_code, bnf_drug_code, bnf_drug_name) %>%
    summarise(items = sum(items), .groups = "drop") %>%
    left_join(class_desc, by = "bnf_class_code")
  if (anyNA(out$bnf_class_name))
    warning("apply_recode_crosswalk: a canonical class code has no descriptor - check crosswalk targets.")
  out %>%
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

cand_drug %>% filter(heals_to_full, combined_items >= 4000)   # keeps the only candidate we need

# Select which of the codes we want to keep for this
cand_drug %>% slice(1) %>% pull(codes)     # the two full codes for the line below
drug_monthly %>% filter(bnf_drug_code %in% c("0301011AB","0302000AA")) %>%
  group_by(bnf_drug_code) %>%
  summarise(first = min(year_month), last = max(year_month), items = sum(items), .groups = "drop")

# Unify codes so same code is used for the above candidate drug
xwalk_drug <- tibble::tribble(
  ~from_code,   ~to_code,
  "0301011AB",  "0302000AA"   # Trimbow: reconciled to the current (2025) corticosteroid code
)

drug_monthly <- apply_recode_crosswalk(drug_monthly, xwalk_drug, "drug")

# Rederive class_monthly
class_monthly <- drug_monthly %>%
  group_by(year_month, bnf_chapter_code, bnf_chapter_name, bnf_section_code,
           bnf_section_name, bnf_class_code, bnf_class_name) %>%
  summarise(items = sum(items), .groups = "drop")

### 6. Eligibility and coverage

# Eligible if prescribed in every month AND at least 1000 items dispensed nationally each calendar year
if (!exists("data_dir")) stop("data_dir not set - run the setup/import section first.")
if (!exists("expected_ym")) {
  expected_ym <- as.integer(format(seq(as.Date("2022-01-01"),
                                       as.Date("2025-12-01"), by = "month"), "%Y%m"))
}
if (!exists("class_monthly")) class_monthly <- readRDS(file.path(data_dir, "class_monthly.rds"))
if (!exists("drug_monthly"))  drug_monthly  <- readRDS(file.path(data_dir, "drug_monthly.rds"))

n_months_req <- length(expected_ym)                    # 48
n_years_req  <- n_distinct(expected_ym %/% 100L)       # 4
min_items_yr <- 1000L

## Eligibility flags per drug class

elig_class <- class_monthly %>%
  group_by(bnf_class_code, bnf_class_name) %>%
  summarise(n_months = n_distinct(year_month), total_items = sum(items), .groups = "drop") %>%
  left_join(
    class_monthly %>%
      group_by(bnf_class_code, year = year_month %/% 100L) %>%
      summarise(items_yr = sum(items), .groups = "drop") %>%
      group_by(bnf_class_code) %>%
      summarise(n_years = n(), min_year_items = min(items_yr), .groups = "drop"),
    by = "bnf_class_code") %>%
  mutate(rule_every_month = n_months == n_months_req,
         rule_min_volume  = n_years == n_years_req & min_year_items >= min_items_yr,
         eligible         = rule_every_month & rule_min_volume)


# Eligibility flags per drug )note that drug names for the same code may have changed, so aggregate by code, and use more recent name)

drug_meta <- drug_monthly %>%
  group_by(bnf_drug_code) %>%
  summarise(bnf_class_code = bnf_class_code[which.max(year_month)],
            bnf_class_name = bnf_class_name[which.max(year_month)],
            bnf_drug_name  = bnf_drug_name[which.max(year_month)],
            .groups = "drop")

elig_drug <- drug_meta %>%
  left_join(
    drug_monthly %>%
      group_by(bnf_drug_code) %>%
      summarise(n_months = n_distinct(year_month), total_items = sum(items), .groups = "drop"),
    by = "bnf_drug_code") %>%
  left_join(
    drug_monthly %>%
      group_by(bnf_drug_code, year = year_month %/% 100L) %>%
      summarise(items_yr = sum(items), .groups = "drop") %>%
      group_by(bnf_drug_code) %>%
      summarise(n_years = n(), min_year_items = min(items_yr), .groups = "drop"),
    by = "bnf_drug_code") %>%
  mutate(rule_every_month = n_months == n_months_req,
         rule_min_volume  = n_years == n_years_req & min_year_items >= min_items_yr,
         eligible         = rule_every_month & rule_min_volume)

## one row per series (guards against grouping anomalies)
stopifnot(
  nrow(elig_class) == n_distinct(class_monthly$bnf_class_code),
  nrow(elig_drug)  == n_distinct(drug_monthly$bnf_drug_code)
)

## List exclusions and summarise coverage

coverage <- bind_rows(
  elig_class %>% summarise(level = "class", total_series = n(), eligible = sum(eligible),
                          excluded = sum(!eligible),
                          excluded_item_share = sum(total_items[!eligible]) / sum(total_items)),
  elig_drug  %>% summarise(level = "drug",  total_series = n(), eligible = sum(eligible),
                          excluded = sum(!eligible),
                          excluded_item_share = sum(total_items[!eligible]) / sum(total_items))
)

add_reason <- function(tab) {
  tab %>%
    filter(!eligible) %>%
    mutate(reason = case_when(
      !rule_every_month & !rule_min_volume ~ "not every month; <1000 items/year",
      !rule_every_month                    ~ "not every month",
      TRUE                                 ~ "<1000 items/year")) %>%
    arrange(desc(total_items))
}
excl_class <- add_reason(elig_class)
excl_drug  <- add_reason(elig_drug)

# Summarise reasons for exclusion
reason_breakdown <- function(excl, elig, level_name) {
  lvls <- c("<1000 items/year", "not every month", "not every month; <1000 items/year")
  n_total <- nrow(elig)
  
  tab <- excl %>%
    count(reason, name = "n") %>%
    right_join(tibble(reason = lvls), by = "reason") %>%   # keep empty categories, fix order
    mutate(n = coalesce(n, 0L),
           pct_of_excluded = 100 * n / sum(n),
           pct_of_all      = 100 * n / n_total) %>%
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

class_monthly_elig <- class_monthly %>%
  semi_join(filter(elig_class, eligible), by = "bnf_class_code")

drug_monthly_elig <- drug_monthly %>%
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
  print(excl_drug %>% slice_head(n = 20) %>%
          select(bnf_drug_code, bnf_drug_name, total_items, reason))
}
if (nrow(excl_class)) {
  cat("  Largest excluded classes:\n")
  print(excl_class %>% slice_head(n = 20) %>%
          select(bnf_class_code, bnf_class_name, total_items, reason))
}

### 7. Model fitting Poisson GLM (fixed spline + Fourier + patient-days offset) -------
# Core engine: one row per series (p-value, distribution, route, diagnostics, harmonic coefficients)

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
                            alpha_disp = 0.05, alpha_lb = 0.05, lb_lag = 12L) {
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
    d$off <- d$offset_log_patient_days
    
    fp <- glm(.f_full, family = poisson, data = d)
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
      fnb_full <- .nb_fit(.f_full, d)
      fnb_red  <- .nb_fit(.f_red,  d)
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
      fr <- glm(.f_red, family = poisson, data = d)
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

screen_class <- class_monthly_elig %>%
  select(bnf_class_code, bnf_class_name, year_month, items) %>%
  group_by(bnf_class_code, bnf_class_name) %>%
  group_modify(~ fit_test_series(.x, covar)) %>%
  ungroup()

# BH across the class family; a class is significant at FDR 5%
screen_class <- screen_class %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         significant = !is.na(p_adj) & p_adj < 0.05) %>%
  arrange(p_adj)

# proportions the methods commit to reporting
route_summary_class <- screen_class %>%
  count(distribution, route, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1))

n_fail_class <- sum(!screen_class$converged)

saveRDS(screen_class,        file.path(data_dir, "screen_class.rds"))
fwrite(screen_class,         file.path(data_dir, "screen_class.csv"))
fwrite(route_summary_class,  file.path(data_dir, "screen_class_route_summary.csv"))

cat(sprintf("Class screen: %d classes | %d significant (BH 5%%) | %d non-converged\n",
            nrow(screen_class), sum(screen_class$significant), n_fail_class))
print(route_summary_class)

# Output key: p_adj = BH-adjusted seasonality p (<0.05 significant); distribution/route = chosen model and test; disp_*/lb_p = overdispersion/autocorrelation diagnostics; sin/cos = harmonic coefficients

### 9. Run fit over every eligible drug, applying BH at 5% FDR, reports route/distribution properties
# Fitted for all drugs; a drug is flagged only if its parent class is significantly seasonal

stopifnot(exists("fit_test_series"), exists("covar"),
          exists("drug_monthly_elig"), exists("screen_class"))

sig_classes <- screen_class %>% filter(significant) %>% pull(bnf_class_code)

# fit every eligible drug
screen_drug <- drug_monthly_elig %>%
  select(bnf_class_code, bnf_drug_code, bnf_drug_name, year_month, items) %>%
  group_by(bnf_class_code, bnf_drug_code, bnf_drug_name) %>%
  group_modify(~ fit_test_series(.x, covar)) %>%
  ungroup() %>%
  mutate(parent_class_sig = bnf_class_code %in% sig_classes)

# PRIMARY: pooled BH within the significant-class group only 
sig_family <- screen_drug %>%
  filter(parent_class_sig) %>%
  mutate(p_adj = p.adjust(p_value, "BH")) %>%
  select(bnf_drug_code, p_adj)

screen_drug <- screen_drug %>%
  left_join(sig_family, by = "bnf_drug_code") %>%
  mutate(significant = parent_class_sig & !is.na(p_adj) & p_adj < 0.05,
         # EXPLORATORY: BH across all eligible drugs
         p_adj_all   = p.adjust(p_value, "BH"),
         sig_all     = !is.na(p_adj_all) & p_adj_all < 0.05) %>%
  left_join(screen_class %>% select(bnf_class_code, bnf_class_name, class_p_adj = p_adj),
            by = "bnf_class_code") %>%
  arrange(!parent_class_sig, p_adj, p_adj_all)   # primary-family first, by significance

route_summary_drug <- screen_drug %>%
  count(distribution, route, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1))

n_fail_drug <- sum(!screen_drug$converged)

saveRDS(screen_drug,       file.path(data_dir, "screen_drug.rds"))
fwrite(screen_drug,        file.path(data_dir, "screen_drug.csv"))
fwrite(route_summary_drug, file.path(data_dir, "screen_drug_route_summary.csv"))

cat(sprintf(paste0(
  "Drug analysis: %d eligible drugs fitted | %d non-converged\n",
  "  PRIMARY (in significant classes): %d of %d significant (pooled BH 5%%)\n",
  "  EXPLORATORY (all drugs):          %d of %d significant (BH 5%%)\n",
  "    of which in NON-significant classes: %d\n"),
  nrow(screen_drug), n_fail_drug,
  sum(screen_drug$significant), sum(screen_drug$parent_class_sig),
  sum(screen_drug$sig_all), nrow(screen_drug),
  sum(screen_drug$sig_all & !screen_drug$parent_class_sig)))
print(route_summary_drug)

### 10. Characterise significant seasonality ------
# Peak:trough ratio (95% CI), peak/trough months, modality (one or two annual cycles)
stopifnot(exists(".nb_fit"), exists(".f_full"), exists("covar"),
          exists("screen_class"), exists("screen_drug"),
          exists("class_monthly_elig"), exists("drug_monthly_elig"))

.f_1h <- items ~ trend1 + trend2 + trend3 + sin12 + cos12 + offset(off)
.harm <- c("sin12", "cos12", "sin6", "cos6")

# 12-month harmonic basis (sin12, cos12, sin6, cos6), reused by the point curve and the bootstrap
.harm_basis <- cbind(sin(2*pi*(1:12)/12), cos(2*pi*(1:12)/12),
                     sin(4*pi*(1:12)/12), cos(4*pi*(1:12)/12))

# Seasonal curve over calendar months 1-12 from the four harmonic coefficients (t = 1 is January)
.seasonal_curve <- function(b_sin12, b_cos12, b_sin6, b_cos6, m = 1:12) {
  b_sin12 * sin(2*pi*m/12) + b_cos12 * cos(2*pi*m/12) +
    b_sin6  * sin(4*pi*m/12) + b_cos6  * cos(4*pi*m/12)
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

# 95% CI for the peak:trough ratio by parametric bootstrap of the harmonic coefficients (route-appropriate covariance)
.ptr_ci <- function(f_full, route, seed = 1L, B = 2000L) {
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

# Seasonal reproducibility: mean pairwise correlation of each year's detrended 12-month profile (~1 if the cycle recurs, ~0 for a one-off break)
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

# STL seasonal/trend strength (Wang-Smith-Hyndman): variance share of each component under a flexible loess trend, so trend-dominated series score low on seasonality
.stl_strength <- function(series, covar) {
  tryCatch({
    d <- merge(covar, series[, c("year_month", "items")], by = "year_month")
    d <- d[order(d$t), ]
    x <- ts(log(d$items) - d$offset_log_patient_days, frequency = 12)  # log rate/patient-day
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
  
  # Refit 1- and 2-harmonic models under the chosen distribution for AIC modality; reuse the 2-harmonic fit for the CI
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
  sig <- scr %>% filter({{keep_flag}})
  if (nrow(sig) == 0) return(tibble())
  sig %>%
    group_by(across(all_of(id_cols))) %>%
    group_modify(function(row, key) {
      ser <- monthly %>% semi_join(key, by = id_cols)
      characterise_one(row, ser, covar)
    }) %>%
    ungroup() %>%
    left_join(scr, by = id_cols) %>%
    arrange(desc(peak_trough_ratio))
}

char_class <- characterise_level(screen_class, class_monthly_elig,
                                 c("bnf_class_code"), significant)
char_drug  <- characterise_level(screen_drug,  drug_monthly_elig,
                                 c("bnf_drug_code"), significant | sig_all)

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
print(char_class %>% slice_head(n = 10) %>%
        mutate(ptr = sprintf("%.2f (%.2f-%.2f)", peak_trough_ratio, ptr_lci, ptr_uci)) %>%
        select(bnf_class_name, ptr, peak_month, modality, p_adj))


### 11. Reporting --------

suppressMessages(library(ggplot2))
stopifnot(exists("char_class"), exists("char_drug"), exists(".seasonal_curve"),
          exists(".f_full"), exists(".nb_fit"), exists("covar"),
          exists("class_monthly_elig"))

# "Meaningful" = appreciable, precisely-estimated amplitude (lower CI of peak:trough) AND a genuine recurring cycle surviving flexible detrending (STL); the 0.50 STL cut sits mid-gap in the real data (robust over 0.40-0.60)
meaningful_threshold    <- 1.10   # lower CI of peak:trough ratio
stl_strength_threshold  <- 0.50   # STL seasonal strength (mid-gap on real data)

fig_dir <- file.path(data_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

## Results table - classes -----------------------------------------------

results_class <- char_class %>%
  transmute(bnf_class_code, bnf_class_name,
            peak_trough_ratio = round(peak_trough_ratio, 3),
            ptr_lci = round(ptr_lci, 3), ptr_uci = round(ptr_uci, 3),
            peak_month, trough_month, modality, n_peaks,
            seasonal_reproducibility = round(seasonal_reproducibility, 3),
            stl_seasonal_strength = round(stl_seasonal_strength, 3),
            stl_trend_strength    = round(stl_trend_strength, 3),
            distribution, route, hac_capped, p_adj,
            meaningful = ptr_lci >= meaningful_threshold &
              stl_seasonal_strength >= stl_strength_threshold) %>%
  arrange(desc(meaningful), desc(peak_trough_ratio))

## Results table - drugs -------------------------------------------------

results_drug <- char_drug %>%
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
              stl_seasonal_strength >= stl_strength_threshold) %>%
  arrange(desc(meaningful), desc(peak_trough_ratio))

## Figures ---------------------------------------------------------------

# observed + fitted monthly rate (per 1000 registered patients) for one class
.fit_frame <- function(code) {
  ser  <- class_monthly_elig %>% filter(bnf_class_code == code) %>% select(year_month, items)
  d    <- covar %>% left_join(ser, by = "year_month") %>% arrange(t)
  d$off <- d$offset_log_patient_days
  dist <- char_class$distribution[char_class$bnf_class_code == code][1]
  fit  <- if (identical(dist, "negbin")) .nb_fit(.f_full, d) else NULL
  if (is.null(fit)) fit <- glm(.f_full, poisson, data = d)
  nm <- char_class$bnf_class_name[char_class$bnf_class_code == code][1]
  d %>% transmute(class = nm, month_date,
                 observed = items / list_size * 1000,
                 fitted   = as.numeric(fitted(fit)) / list_size * 1000)
}

# Exemplars: strongest meaningful classes, excluding the extreme vaccine/antiviral series that would flatten the shared axes
exemplar_class_codes <- results_class %>%
  filter(meaningful, peak_trough_ratio < 10) %>%
  slice_head(n = 6) %>% pull(bnf_class_code)

if (length(exemplar_class_codes)) {
  fit_df <- bind_rows(lapply(exemplar_class_codes, .fit_frame)) %>%
    arrange(class)
  fit_df$class <- factor(fit_df$class,
                         levels = results_class %>% filter(bnf_class_code %in% exemplar_class_codes) %>%
                           arrange(desc(peak_trough_ratio)) %>% pull(bnf_class_name))
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
  shape_df <- char_class %>%
    filter(bnf_class_code %in% exemplar_class_codes) %>%
    group_by(bnf_class_name) %>%
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
  screen_class %>% count(distribution, route, name = "n") %>% mutate(level = "class"),
  screen_drug  %>% count(distribution, route, name = "n") %>% mutate(level = "drug")) %>%
  group_by(level) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup() %>%
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
print(results_class %>% filter(meaningful) %>%
        mutate(ptr = sprintf("%.2f (%.2f-%.2f)", peak_trough_ratio, ptr_lci, ptr_uci)) %>%
        select(bnf_class_name, ptr, peak_month, modality, stl_seasonal_strength) %>%
        slice_head(n = 20))


### 12. Publication outputs (tables and figures) ------------------------------
# Publication-ready tables and figures into <data_dir>/results; appendix plots every eligible series, paginated by BNF chapter (requires Sections 4-11 objects)

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

# Render integer class codes zero-padded and DOTTED ("03.01.01") so they survive re-reading in Excel/fread as non-numeric; drug codes contain letters and are already safe
bnf_class_dotted <- function(code) {
  s <- sprintf("%06d", as.integer(code))
  sub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1.\\2.\\3", s)
}

# BNF chapter lookup from the eligible monthly frames
chap_class <- class_monthly_elig %>%
  distinct(bnf_class_code, bnf_chapter_code, bnf_chapter_name)
chap_drug <- drug_monthly_elig %>%
  distinct(bnf_drug_code, bnf_chapter_code, bnf_chapter_name)

## ---- MAIN TEXT: Table 1, meaningful-seasonal classes ------------------------

main_table_classes <- results_class %>%
  filter(meaningful) %>%
  transmute(
    `BNF class code`                = bnf_class_dotted(bnf_class_code),
    `BNF class`                     = bnf_class_name,
    `Peak-to-trough ratio (95% CI)` = fmt_ci(peak_trough_ratio, ptr_lci, ptr_uci),
    `Peak month`                    = peak_month,
    `Trough month`                  = trough_month,
    `Seasonal shape`                = relabel_modality(modality),
    `Seasonal strength`             = sprintf("%.2f", stl_seasonal_strength),
    `Adjusted p`                    = fmt_p(p_adj))
fwrite(main_table_classes, file.path(tab_dir, "table1_meaningful_classes.csv"))

## ---- MAIN TEXT: Table 2, meaningful-seasonal drugs -------------------------
# Drug-level parallel to Table 1; adjusted p is the within-class (confirmatory) value where available, else the all-drug (exploratory) value
main_table_drugs <- results_drug %>%
  filter(meaningful) %>%
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
fwrite(main_table_drugs, file.path(tab_dir, "table2_meaningful_drugs.csv"))

## ---- MAIN TEXT: analytic accounting (funnel) -------------------------------

sig_class_n <- sum(screen_class$significant, na.rm = TRUE)
sig_drug_n  <- if (exists("screen_drug")) sum(screen_drug$significant, na.rm = TRUE) else NA
accounting <- tibble::tibble(
  Stage = c("Eligible series",
            "Statistically significant (BH FDR 5%)",
            "Significant but trivial amplitude (ratio < 1.05)",
            "Meaningful seasonality (amplitude CI \u2265 threshold and seasonal strength \u2265 threshold)"),
  Classes = c(
    if (exists("elig_class")) sum(elig_class$eligible) else NA_integer_,
    sig_class_n,
    sum(results_class$peak_trough_ratio < 1.05, na.rm = TRUE),
    sum(results_class$meaningful, na.rm = TRUE)),
  Drugs = c(
    if (exists("elig_drug")) sum(elig_drug$eligible) else NA_integer_,
    sig_drug_n,
    sum(results_drug$peak_trough_ratio < 1.05, na.rm = TRUE),
    sum(results_drug$meaningful, na.rm = TRUE)))
fwrite(accounting, file.path(tab_dir, "table_accounting.csv"))

## ---- APPENDIX: full class and drug results ---------------------------------

appendix_table_classes <- results_class %>%
  left_join(chap_class, by = "bnf_class_code") %>%
  transmute(
    `BNF chapter`                   = bnf_chapter_name,
    `BNF class code`                = bnf_class_dotted(bnf_class_code),
    `BNF class`                     = bnf_class_name,
    `Peak-to-trough ratio`          = sprintf("%.3f", peak_trough_ratio),
    `95% CI lower`                  = sprintf("%.3f", ptr_lci),
    `95% CI upper`                  = sprintf("%.3f", ptr_uci),
    `Peak month`                    = peak_month,
    `Trough month`                  = trough_month,
    `Seasonal shape`                = relabel_modality(modality),
    `Seasonal strength (STL)`       = sprintf("%.3f", stl_seasonal_strength),
    `Trend strength (STL)`          = sprintf("%.3f", stl_trend_strength),
    `Cross-year reproducibility`    = sprintf("%.3f", seasonal_reproducibility),
    `Distribution`                  = tools::toTitleCase(distribution),
    `Inference method`              = relabel_route(route),
    `HAC bandwidth capped`          = ifelse(is.na(hac_capped), "", ifelse(hac_capped, "Yes", "No")),
    `Adjusted p`                    = fmt_p(p_adj),
    `Meaningful seasonality`        = ifelse(meaningful, "Yes", "No"))
fwrite(appendix_table_classes, file.path(tab_dir, "appendixA1_all_classes.csv"))

appendix_table_drugs <- results_drug %>%
  left_join(chap_drug, by = "bnf_drug_code") %>%
  transmute(
    `BNF chapter`                   = bnf_chapter_name,
    `BNF class code`                = sub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1.\\2.\\3",
                                          substr(as.character(bnf_drug_code), 1, 6)),
    `BNF class`                     = bnf_class_name,
    `BNF drug code`                 = as.character(bnf_drug_code),
    `Drug (chemical substance)`     = bnf_drug_name,
    `Peak-to-trough ratio`          = sprintf("%.3f", peak_trough_ratio),
    `95% CI lower`                  = sprintf("%.3f", ptr_lci),
    `95% CI upper`                  = sprintf("%.3f", ptr_uci),
    `Peak month`                    = peak_month,
    `Trough month`                  = trough_month,
    `Seasonal shape`                = relabel_modality(modality),
    `Seasonal strength (STL)`       = sprintf("%.3f", stl_seasonal_strength),
    `Trend strength (STL)`          = sprintf("%.3f", stl_trend_strength),
    `Distribution`                  = tools::toTitleCase(distribution),
    `Inference method`              = relabel_route(route),
    `In significant class`          = ifelse(parent_class_sig, "Yes", "No"),
    `Significant (within class)`    = ifelse(significant_primary, "Yes", "No"),
    `Adjusted p (within class)`     = fmt_p(p_adj_primary),
    `Significant (all-drug scan)`   = ifelse(sig_all, "Yes", "No"),
    `Adjusted p (all-drug scan)`    = fmt_p(p_adj_all),
    `Meaningful seasonality`        = ifelse(meaningful, "Yes", "No"))
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
  fwrite(excl_class %>% mutate(bnf_class_code = bnf_class_dotted(bnf_class_code)),
         file.path(tab_dir, "appendix_exclusions_class.csv"))
if (exists("excl_drug") && nrow(excl_drug))
  fwrite(excl_drug %>% mutate(bnf_drug_code = as.character(bnf_drug_code)),
         file.path(tab_dir, "appendix_exclusions_drug.csv"))
if (exists("xwalk_drug")) fwrite(xwalk_drug, file.path(tab_dir, "appendix_recode_crosswalk.csv"))
if (exists("appendix_routes")) fwrite(appendix_routes, file.path(tab_dir, "appendix_route_proportions.csv"))

## ---- shared plotting theme --------------------------------------------------

theme_pub <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# observed + fitted rate (per 1000 patients) for one class or drug series
.series_fit <- function(mon, id_col, id_val) {
  ser <- mon %>% filter(.data[[id_col]] == id_val) %>% select(year_month, items)
  d   <- covar %>% left_join(ser, by = "year_month") %>% arrange(t)
  d$off <- d$offset_log_patient_days
  fit <- suppressWarnings(glm(.f_full, poisson, data = d))
  tibble(month_date = d$month_date,
         observed = d$items / d$list_size * 1000,
         fitted   = as.numeric(fitted(fit)) / d$list_size * 1000)
}

## ---- MAIN TEXT figures: exemplars ------------------------------------------

exemplar_codes <- results_class %>%
  filter(meaningful, peak_trough_ratio < 10) %>%
  slice_head(n = 6) %>% pull(bnf_class_code)

if (length(exemplar_codes)) {
  ex_names <- setNames(char_class$bnf_class_name[match(exemplar_codes, char_class$bnf_class_code)],
                       exemplar_codes)
  lvl <- results_class %>% filter(bnf_class_code %in% exemplar_codes) %>%
    arrange(desc(peak_trough_ratio)) %>% pull(bnf_class_name)
  
  fit_df <- bind_rows(lapply(exemplar_codes, function(cd)
    .series_fit(class_monthly_elig, "bnf_class_code", cd) %>% mutate(class = ex_names[as.character(cd)])))
  fit_df$class <- factor(fit_df$class, levels = lvl)
  
  p1 <- ggplot(fit_df, aes(month_date)) +
    geom_point(aes(y = observed), size = 1, colour = "grey45") +
    geom_line(aes(y = fitted), colour = "#B2182B", linewidth = 0.8) +
    facet_wrap(~ class, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Prescription items per 1000 registered patients",
         title = "Observed monthly prescribing with fitted trend and seasonal model") +
    theme_pub
  ggsave(file.path(fig_dir, "figure1_exemplar_observed_fitted.png"), p1, width = 9, height = 8, dpi = 300)
  
  shape_df <- char_class %>% filter(bnf_class_code %in% exemplar_codes) %>%
    group_by(bnf_class_name) %>%
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
  ggsave(file.path(fig_dir, "figure2_seasonal_shape.png"), p2, width = 9, height = 8, dpi = 300)
}

## ---- MAIN TEXT figure: amplitude x peak-month landscape --------------------


land <- results_class %>%
  filter(meaningful) %>%
  left_join(chap_class, by = "bnf_class_code") %>%
  mutate(peak_month = factor(peak_month, levels = month.abb),
         chapter = bnf_chapter_name)
# Visible ceiling: classes above y_cap are drawn as triangles pinned at the top, with the true ratio in the label
y_cap <- 3
land <- land %>%
  mutate(above  = peak_trough_ratio > y_cap,
         y_plot = pmin(peak_trough_ratio, y_cap),
         lab    = ifelse(above,
                         sprintf("%s (%.0f\u00d7)", bnf_class_name, peak_trough_ratio),
                         bnf_class_name))
if (nrow(land)) {
  pal <- grDevices::hcl.colors(max(3, length(unique(land$chapter))), "Dark 3")
  set.seed(1)
  jit <- position_jitter(width = 0.18, height = 0, seed = 1)
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
                                        max.overlaps = Inf, seed = 1, position = jit, show.legend = FALSE, force = 2, point.padding = NA, min.segment.length = 0, box.padding = 0.5)
  ggsave(file.path(fig_dir, "figure3_seasonality_landscape.png"), p3, width = 11, height = 8, dpi = 300)
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
    .series_fit(class_monthly_elig, "bnf_class_code", cd) %>%
      mutate(class = sprintf("%s (%s)", nm, ifelse(keep, "retained", "excluded")))
  }))
  p4 <- ggplot(cmp_df, aes(month_date)) +
    geom_point(aes(y = observed), size = 0.9, colour = "grey45") +
    geom_line(aes(y = fitted), colour = "#B2182B", linewidth = 0.7) +
    facet_wrap(~ class, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Prescription items per 1000 registered patients",
         title = "Genuine seasonality versus excluded structural-break and trend-dominated series") +
    theme_pub
  ggsave(file.path(fig_dir, "figure_appendix_genuine_vs_excluded.png"), p4, width = 9, height = 8, dpi = 300)
}

## ---- APPENDIX figures: ALL series, observed + fitted, paginated by chapter --
# Every eligible series plotted; strip colour = analytic tier (dark red meaningful, pale red significant, grey neither), recoloured on the grob by panel label

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

# Tier lookups keyed by series name (meaningful outranks significant); for drugs "significant" = within-class flag OR all-drug scan
FILL_MEAN <- "#8B0000"; TXT_MEAN <- "white"
FILL_SIG  <- "#E8A0A0"; TXT_SIG  <- "grey15"
tier_class <- {
  m <- setNames(vector("list", nrow(results_class)), results_class$bnf_class_name)
  for (i in seq_len(nrow(results_class)))
    m[[i]] <- if (isTRUE(results_class$meaningful[i])) c(FILL_MEAN, TXT_MEAN) else c(FILL_SIG, TXT_SIG)
  m
}
tier_drug <- {
  sig <- (results_drug$significant_primary %in% TRUE) | (results_drug$sig_all %in% TRUE)
  m <- setNames(vector("list", nrow(results_drug)), results_drug$bnf_drug_name)
  for (i in seq_len(nrow(results_drug)))
    m[[i]] <- if (isTRUE(results_drug$meaningful[i])) c(FILL_MEAN, TXT_MEAN)
  else if (sig[i])                        c(FILL_SIG, TXT_SIG)
  else                                    NULL   # eligible-only: default grey
  m[!vapply(m, is.null, logical(1))]
}

# Same tiering in a blue palette, used only for the seasonal-shape appendix figures
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
  sig <- (results_drug$significant_primary %in% TRUE) | (results_drug$sig_all %in% TRUE)
  m <- setNames(vector("list", nrow(results_drug)), results_drug$bnf_drug_name)
  for (i in seq_len(nrow(results_drug)))
    m[[i]] <- if (isTRUE(results_drug$meaningful[i])) c(FILL_MEAN_SHAPE, TXT_MEAN_SHAPE)
  else if (sig[i])                        c(FILL_SIG_SHAPE, TXT_SIG_SHAPE)
  else                                    NULL
  m[!vapply(m, is.null, logical(1))]
}

render_all_series <- function(mon, id_col, name_col, chap_lookup, out_pdf,
                              tier_map = list(), per_page = 20, ncol = 4) {
  ids <- chap_lookup %>% arrange(bnf_chapter_code) %>%
    left_join(distinct(mon, .data[[id_col]], .data[[name_col]]), by = id_col)
  pdf(out_pdf, width = 11.7, height = 8.3)   # A4 landscape
  on.exit(dev.off())
  first_page <- TRUE
  for (ch in unique(ids$bnf_chapter_code)) {
    chn <- ids %>% filter(bnf_chapter_code == ch)
    pages <- split(seq_len(nrow(chn)), ceiling(seq_len(nrow(chn)) / per_page))
    for (pi in seq_along(pages)) {
      rows <- chn[pages[[pi]], ]
      # Disambiguate any duplicate panel name on this page by appending its code, and key the page tier lookup by the same label
      nm <- rows[[name_col]]
      dup <- nm %in% nm[duplicated(nm)]
      panel_label <- ifelse(dup, sprintf("%s [%s]", nm, rows[[id_col]]), nm)
      rows$panel_label <- panel_label
      page_tier_map <- setNames(tier_map[nm], panel_label)
      page_tier_map <- page_tier_map[!vapply(page_tier_map, is.null, logical(1))]
      
      dd <- bind_rows(lapply(seq_len(nrow(rows)), function(i) {
        .series_fit(mon, id_col, rows[[id_col]][i]) %>%
          mutate(panel = rows$panel_label[i])
      }))
      dd$panel <- factor(dd$panel, levels = rows$panel_label)
      ptitle <- sprintf("Chapter %02d - %s  (page %d of %d)",
                        ch, rows$bnf_chapter_name[1], pi, length(pages))
      p <- ggplot(dd, aes(month_date)) +
        geom_point(aes(y = observed), size = 0.35, colour = "grey55") +
        geom_line(aes(y = fitted), colour = "#B2182B", linewidth = 0.4) +
        facet_wrap(~ panel, scales = "free_y", ncol = ncol) +
        labs(x = NULL, y = "Items per 1000 patients", title = ptitle) +
        theme_bw(base_size = 7) +
        theme(strip.background = element_rect(fill = "grey92", colour = NA),
              strip.text = element_text(size = 6, face = "bold"), panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = 9))
      # recolour flagged strips on the grob, then draw (print() does not draw a gtable)
      if (!first_page) grid::grid.newpage()
      grid::grid.draw(.style_strips(p, page_tier_map))
      first_page <- FALSE
    }
  }
}

render_all_series(class_monthly_elig, "bnf_class_code", "bnf_class_name",
                  chap_class, file.path(fig_dir, "appendix_all_classes_by_chapter.pdf"),
                  tier_map = tier_class, per_page = 20, ncol = 4)
render_all_series(drug_monthly_elig, "bnf_drug_code", "bnf_drug_name",
                  chap_drug, file.path(fig_dir, "appendix_all_drugs_by_chapter.pdf"),
                  tier_map = tier_drug, per_page = 20, ncol = 4)

## ---- APPENDIX figures: ALL fitted seasonal SHAPES, paginated by chapter ----
# Fitted seasonal shape (Poisson throughout) for every eligible series; a flat curve for a null series is itself informative
FILL_MEAN_SHAPE <- "#08306B"; TXT_MEAN_SHAPE <- "white"   # dark blue
FILL_SIG_SHAPE  <- "#9ECAE1"; TXT_SIG_SHAPE  <- "grey15"   # pale blue

.series_shape <- function(mon, id_col, id_val) {
  ser <- mon %>% filter(.data[[id_col]] == id_val) %>% select(year_month, items)
  d   <- covar %>% left_join(ser, by = "year_month") %>% arrange(t)
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
                              tier_map = list(), per_page = 20, ncol = 4) {
  ids <- chap_lookup %>% arrange(bnf_chapter_code) %>%
    left_join(distinct(mon, .data[[id_col]], .data[[name_col]]), by = id_col)
  pdf(out_pdf, width = 11.7, height = 8.3)   # A4 landscape
  on.exit(dev.off())
  first_page <- TRUE
  for (ch in unique(ids$bnf_chapter_code)) {
    chn <- ids %>% filter(bnf_chapter_code == ch)
    pages <- split(seq_len(nrow(chn)), ceiling(seq_len(nrow(chn)) / per_page))
    for (pi in seq_along(pages)) {
      rows <- chn[pages[[pi]], ]
      pl <- .page_labels(rows, name_col, id_col, tier_map)
      rows$panel_label <- pl$labels
      
      dd <- bind_rows(lapply(seq_len(nrow(rows)), function(i) {
        .series_shape(mon, id_col, rows[[id_col]][i]) %>%
          mutate(panel = rows$panel_label[i])
      }))
      dd$panel <- factor(dd$panel, levels = rows$panel_label)
      ptitle <- sprintf("Chapter %02d - %s  (page %d of %d) - fitted seasonal shape",
                        ch, rows$bnf_chapter_name[1], pi, length(pages))
      p <- ggplot(dd, aes(month, factor)) +
        geom_hline(yintercept = 1, colour = "grey75", linewidth = 0.3) +
        geom_line(colour = "#2166AC", linewidth = 0.5) +
        facet_wrap(~ panel, ncol = ncol) +
        scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c("Jan", "Apr", "Jul", "Oct")) +
        labs(x = NULL, y = "Seasonal factor (relative to annual mean)", title = ptitle) +
        theme_bw(base_size = 7) +
        theme(strip.background = element_rect(fill = "grey92", colour = NA),
              strip.text = element_text(size = 6, face = "bold"), panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = 9))
      if (!first_page) grid::grid.newpage()
      grid::grid.draw(.style_strips(p, pl$tier))
      first_page <- FALSE
    }
  }
}

render_all_shapes(class_monthly_elig, "bnf_class_code", "bnf_class_name", chap_class,
                  file.path(fig_dir, "appendix_all_classes_seasonal_shape_by_chapter.pdf"),
                  tier_map = tier_class_shape, per_page = 20, ncol = 4)
render_all_shapes(drug_monthly_elig, "bnf_drug_code", "bnf_drug_name", chap_drug,
                  file.path(fig_dir, "appendix_all_drugs_seasonal_shape_by_chapter.pdf"),
                  tier_map = tier_drug_shape, per_page = 20, ncol = 4)

## ---- report -----------------------------------------------------------------

cat(sprintf(paste0(
  "Publication outputs written to %s\n",
  "  Tables (tables/): table1_meaningful_classes, table_accounting,\n",
  "    appendixA1_all_classes, appendixA2_all_drugs, coverage, exclusions, crosswalk, routes\n",
  "  Figures (figures/): figure1_exemplar_observed_fitted, figure2_seasonal_shape,\n",
  "    figure3_seasonality_landscape, figure_appendix_genuine_vs_excluded,\n",
  "    appendix_all_classes_by_chapter.pdf, appendix_all_drugs_by_chapter.pdf\n",
  "  Main-text meaningful classes: %d | eligible classes plotted: %d | eligible drugs plotted: %d\n"),
  res_dir, nrow(main_table_classes),
  n_distinct(class_monthly_elig$bnf_class_code), n_distinct(drug_monthly_elig$bnf_drug_code)))


### 13. Publication tables as Word documents (one .docx per table) -----------
# Stages a JSON manifest per Section-12 CSV for a Node/docx-js renderer; row shading mirrors the appendix tiering

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
       notes = "Meaningful seasonality defined as for Table 1. 'In significant class' indicates whether the drug's parent BNF class was itself significantly seasonal after correction (the primary, conditional analysis); adjusted p is the within-class value where available, otherwise the value from the exploratory all-drug analysis. Drug-level findings are a descriptive drill-down and are not an independent confirmatory analysis.",
       orientation = "landscape", shade_col = NULL),
  list(id = "table3", file = "table_accounting.csv",
       out = "Table3_analytic_accounting.docx",
       caption = "Table 3. Analytic accounting of drug classes and individual drugs, from eligibility to meaningful seasonality.",
       notes = "Eligible series were prescribed in every month of the study window and reached at least 1,000 items dispensed nationally in each calendar year. Significant series passed the class- or drug-level seasonality test at a 5% false discovery rate. Trivial-amplitude series are significant series with a peak-to-trough ratio under 1.05.",
       orientation = "portrait", shade_col = NULL),
  list(id = "s1", file = "appendixA1_all_classes.csv",
       out = "TableS1_all_classes.docx",
       caption = "Supplementary Table S1. Seasonality characteristics for all eligible drug classes.",
       notes = "Shading: dark red rows show meaningful seasonality (as defined for Table 1); pale red rows are statistically significant (5% FDR) but did not meet the amplitude or seasonal-strength threshold; unshaded rows were eligible but not statistically significant. HAC = Newey-West heteroscedasticity- and autocorrelation-consistent standard errors; LRT = likelihood ratio test; STL = seasonal-trend decomposition using Loess.",
       orientation = "landscape", shade_col = "Meaningful seasonality",
       sig_col = NULL),
  list(id = "s2", file = "appendixA2_all_drugs.csv",
       out = "TableS2_all_drugs.docx",
       caption = "Supplementary Table S2. Seasonality characteristics for all eligible individual drugs.",
       notes = "Shading as for Table S1. 'Significant (within class)' is the primary, conditional test (drugs in classes found significant at class level, pooled Benjamini-Hochberg correction); 'Significant (all-drug scan)' is an exploratory analysis correcting across every eligible drug regardless of parent-class significance, included because a class can be non-significant on average while containing an individually seasonal drug.",
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

## Swap Unicode symbols (en-dash, >=, etc.) for plain ASCII so they render through the docx-js -> LibreOffice path
.ascii_safe <- function(x) {
  x %>%
    gsub("\u2265", ">=", x = ., useBytes = TRUE, fixed = TRUE) %>%
    gsub("\u2264", "<=", x = ., useBytes = TRUE, fixed = TRUE) %>%
    gsub("\u2013", "-",  x = ., useBytes = TRUE, fixed = TRUE) %>%
    gsub("\u2014", "-",  x = ., useBytes = TRUE, fixed = TRUE) %>%
    gsub("\u00d7", "x",  x = ., useBytes = TRUE, fixed = TRUE) %>%
    gsub("\u2018", "'",  x = ., useBytes = TRUE, fixed = TRUE) %>%
    gsub("\u2019", "'",  x = ., useBytes = TRUE, fixed = TRUE) %>%
    gsub("\u201c", '"',  x = ., useBytes = TRUE, fixed = TRUE) %>%
    gsub("\u201d", '"',  x = ., useBytes = TRUE, fixed = TRUE)
}

## build one JSON manifest entry per table that has a source file present
manifest <- list()
for (spec in tbl_spec) {
  f <- file.path(tab_dir, spec$file)
  if (!file.exists(f)) { message("skip (not found): ", spec$file); next }
  # Read as character so publication-formatted values (dotted codes, formatted p) are preserved exactly
  dt <- fread(f, colClasses = "character")
  if (nrow(dt) == 0) { message("skip (empty): ", spec$file); next }
  
  # Shade rows: meaningful -> dark red, else significant -> pale red; tables without a shade column stay unshaded
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
# Descriptive only: aggregates 2021 EPD by class and plots raw 2021-2025 rates (2021 shaded) to show disruption is confined to 2021 (requires Sections 3-11 objects)

stopifnot(exists("out_dir"), exists("covar"), exists("results_class"),
          exists("class_monthly_elig"))
if (!exists("data_dir")) data_dir <- out_dir
if (!exists("class_monthly"))
  class_monthly <- readRDS(file.path(data_dir, "class_monthly.rds"))

# Editable: representative seasonal classes (exact BNF names); a high-volume non-significant class is appended automatically as a flat control
exemplar_names <- c(
  "Penicillins",                 # winter antibiotic (COVID-suppressed 2020-21)
  "Macrolides",                  # winter antibiotic
  "Antihistamines",              # summer allergy (opposite-phase cycle)
  "Cough suppressants",          # winter respiratory
  "Sunscreening preparations"    # summer dermatological
)

epd_dir <- "EPD DIRECTORY"
ls_dir  <- "LIST SIZE DIRECTORY"
if (!dir.exists(epd_dir)) stop("EPD directory not found: ", epd_dir,
                               " - mount the drive or adjust the path.")
if (!dir.exists(ls_dir))  stop("List-size directory not found: ", ls_dir,
                               " - mount the drive or adjust the path.")

## 14a. Aggregate the 2021 EPD to national monthly totals per BNF class ---------
# Cached per month (all chapters 1-14), so the 1 GB CSVs are read only once

cache_dir_2021 <- file.path(data_dir, "cache_2021")
dir.create(cache_dir_2021, showWarnings = FALSE, recursive = TRUE)

epd_all  <- list.files(epd_dir, pattern = "\\.csv$", full.names = TRUE)
epd_2021 <- epd_all[as.integer(str_extract(basename(epd_all), "\\d{6}")) %in% 202101:202112]
if (!length(epd_2021))
  stop("No 2021 EPD monthly files found in ", epd_dir,
       " - download the Jan-Dec 2021 SNOMED EPD files first.")

for (f in epd_2021) {
  ym      <- str_extract(basename(f), "\\d{6}")
  cache_f <- file.path(cache_dir_2021, paste0("data_byclass_2021_", ym, ".csv"))
  if (file.exists(cache_f)) { message("Skipping 2021 ", ym, " (cached)"); next }
  message("Aggregating 2021 EPD ", ym)
  
  # read only the presentation code and item count (schema-robust)
  hdr      <- names(fread(f, nrows = 0))
  code_col <- if ("BNF_CODE" %in% hdr) "BNF_CODE"
  else if ("BNF_PRESENTATION_CODE" %in% hdr) "BNF_PRESENTATION_CODE"
  else stop("no recognised presentation-code column in ", basename(f))
  d <- fread(f, select = c(code_col, "ITEMS"))
  setnames(d, c(code_col, "ITEMS"), c("bnf_code", "items"))
  
  d[, bnf_chapter_code := as.integer(substr(bnf_code, 1, 2))]
  d <- d[!is.na(bnf_chapter_code) & bnf_chapter_code <= 14]      # chapters 1-14
  d[, bnf_class_code := as.integer(substr(bnf_code, 1, 6))]       # paragraph = class
  agg <- d[, .(items = sum(items)), by = bnf_class_code]
  agg[, year_month := as.integer(ym)]
  
  fwrite(agg[, .(year_month, bnf_class_code, items)], cache_f)
  rm(d, agg); gc()                                                # free before next file
}

cache_files <- list.files(cache_dir_2021,
                          pattern = "^data_byclass_2021_\\d{6}\\.csv$", full.names = TRUE)
class_2021  <- rbindlist(lapply(cache_files, fread), use.names = TRUE)
n_months_2021 <- uniqueN(class_2021$year_month)
if (n_months_2021 < 12L)
  warning(sprintf("Only %d of 12 months of 2021 EPD are present.", n_months_2021))

## 14b. National monthly list size for 2021 (quarterly -> monthly by LOCF) ------

ls_files <- list.files(ls_dir, pattern = "^gp-reg-pat-prac-all_\\d{6}\\.csv$",
                       full.names = TRUE)
ls_2021_files <- ls_files[as.integer(str_extract(basename(ls_files), "\\d{6}")) %in% 202101:202112]
if (!length(ls_2021_files))
  stop("No 2021 list-size files found in ", ls_dir,
       " (expected gp-reg-pat-prac-all_YYYYMM.csv). Adjust the pattern if the ",
       "2021 filenames differ.")

ls2021_totals <- rbindlist(lapply(ls_2021_files, function(f) {
  ym <- str_extract(basename(f), "\\d{6}")
  d  <- fread(f, select = "NUMBER_OF_PATIENTS")
  data.table(year_month = as.integer(ym), list_size = sum(d$NUMBER_OF_PATIENTS, na.rm = TRUE))
}))

# complete 2021 monthly grid, carry forward, back-fill any leading gap
ls2021 <- data.table(date = seq(as.Date("2021-01-01"), as.Date("2021-12-01"), by = "month")) %>%
  mutate(year_month = as.integer(format(date, "%Y%m"))) %>%
  left_join(ls2021_totals, by = "year_month") %>%
  arrange(date) %>%
  tidyr::fill(list_size, .direction = "downup") %>%
  transmute(year_month, list_size)

# combine with the 2022-2025 denominator already used in the analysis (covar)
ls_all <- bind_rows(ls2021, covar %>% transmute(year_month, list_size)) %>%
  arrange(year_month)

## 14c. Choose exemplar classes (named seasonal + one flat control) -------------

ex <- results_class %>%
  distinct(bnf_class_code, bnf_class_name) %>%
  filter(bnf_class_name %in% exemplar_names)
missing_names <- setdiff(exemplar_names, ex$bnf_class_name)
if (length(missing_names))
  warning("exemplar class name(s) not found and skipped: ",
          paste(missing_names, collapse = ", "),
          " - check spelling against results_class$bnf_class_name.")

# flat control: highest-volume eligible class that was NOT significant
flat <- class_monthly_elig %>%
  filter(!bnf_class_code %in% results_class$bnf_class_code) %>%
  group_by(bnf_class_code, bnf_class_name) %>%
  summarise(items = sum(items), .groups = "drop") %>%
  arrange(desc(items)) %>%
  slice_head(n = 1) %>%
  select(bnf_class_code, bnf_class_name)

exemplars <- bind_rows(ex, flat) %>% distinct(bnf_class_code, bnf_class_name)
if (nrow(exemplars) < 1L) stop("No exemplar classes resolved.")

# facet order: named seasonal classes by amplitude (largest first), control last
ex_ord <- ex %>%
  left_join(results_class %>% select(bnf_class_code, peak_trough_ratio), by = "bnf_class_code") %>%
  arrange(desc(peak_trough_ratio))
facet_levels <- unique(c(ex_ord$bnf_class_name, flat$bnf_class_name))

## 14d. Build the 2021-2025 rate series and plot --------------------------------

items_all <- bind_rows(
  class_2021    %>% transmute(year_month, bnf_class_code = as.integer(bnf_class_code), items),
  class_monthly %>% transmute(year_month, bnf_class_code = as.integer(bnf_class_code), items))

plot_df <- exemplars %>%
  inner_join(items_all, by = "bnf_class_code") %>%
  left_join(ls_all, by = "year_month") %>%
  mutate(date = as.Date(paste0(year_month, "01"), "%Y%m%d"),
         rate = items / list_size * 1000,
         class = factor(bnf_class_name, levels = facet_levels)) %>%
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
  geom_vline(xintercept = as.Date("2022-01-01"), linetype = 2,
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
fwrite(plot_df %>% select(bnf_class_code, bnf_class_name, year_month, items, list_size, rate),
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
# Diagnostic: per class, quantify 2021's departure from the 2022-2025 pattern in level (vs back-extrapolated trend) and shape (within-year profile correlation)

stopifnot(exists("ls_all"), exists("res_dir"))          # from 14b / 14d
if (!exists("class_2021")) {
  cf <- list.files(cache_dir_2021, pattern = "^data_byclass_2021_\\d{6}\\.csv$", full.names = TRUE)
  if (!length(cf)) stop("No 2021 class cache found - run Section 14a first.")
  class_2021 <- rbindlist(lapply(cf, fread), use.names = TRUE)
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
  geom_vline(xintercept = as.Date("2022-01-01"), linetype = 2, colour = "grey40", linewidth = 0.4) +
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
# Diagnostic: paginated PDF of every eligible class's raw 2021-2025 rate (2021 shaded); strips recoloured by the 14e flags (red level, blue shape, grey none)

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

# Panel label per class: name + signed 2021 level deviation (duplicates disambiguated by code)
sc <- as.data.table(scan)[, .(bnf_class_code, bnf_class_name, bnf_chapter_code,
                              bnf_chapter_name, level_dev_pct, flag_level, flag_shape)]
sc[, base_label := sprintf("%s (2021 %+.0f%%)", bnf_class_name, round(level_dev_pct))]
sc[, tier := fifelse(flag_level, "level", fifelse(flag_shape, "shape", "none"))]

order_key <- function(dt) dt[order(bnf_chapter_code, -abs(level_dev_pct))]  # worst first within chapter
sc <- order_key(sc)

out_pdf <- file.path(res_dir, "figures", "figure_window_all_classes_2021_2025.pdf")
per_page <- 20L; ncol <- 4L

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
                      ch, rows$bnf_chapter_name[1], pi, length(pages))
    p <- ggplot(dd, aes(date, rate)) +
      annotate("rect", xmin = as.Date("2020-12-16"), xmax = as.Date("2021-12-16"),
               ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
      geom_vline(xintercept = as.Date("2022-01-01"), linetype = 2, colour = "grey45", linewidth = 0.3) +
      geom_line(colour = "grey20", linewidth = 0.35) +
      geom_point(size = 0.3, colour = "grey35") +
      facet_wrap(~ panel, scales = "free_y", ncol = ncol) +
      scale_x_date(date_breaks = "1 year", date_labels = "'%y") +
      labs(x = NULL, y = "Items per 1000 patients", title = ptitle) +
      theme_bw(base_size = 7) +
      theme(strip.background = element_rect(fill = "grey92", colour = NA),
            strip.text = element_text(size = 6, face = "bold"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 8))
    
    if (!first_page) grid::grid.newpage()
    grid::grid.draw(.style_strips(p, tier_map))
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
