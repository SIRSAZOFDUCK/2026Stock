# Prespecified sensitivity analyses, each operating on the shared run context.

run_trend_sensitivity <- function(state) {
  evalq({
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
  }, envir = state)
  invisible(TRUE)
}


run_hac_sensitivity <- function(state) {
  evalq({
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
  }, envir = state)
  invisible(TRUE)
}


run_working_days_sensitivity <- function(state) {
  evalq({
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
  }, envir = state)
  invisible(TRUE)
}


run_threshold_sensitivity <- function(state) {
  evalq({
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
  }, envir = state)
  invisible(TRUE)
}
