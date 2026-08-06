# Module 04: model diagnostics, cohort flow, missingness, and warnings.
## 04.1 Verify upstream snapshots and primary-result alignment ----------------
run_diagnostics <- function(state) {
  # Diagnostics run against the exact sealed primary and sensitivity outputs.
  evalq({
  stage6_dir <- file.path(out_dir, "qc", "stage6")
  stage6_snapshot_dir <- file.path(stage6_dir, "stage6_snapshot")
  dir.create(stage6_snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  # Collect all diagnostic-release gates in a single table.
  stage6_checks <- list()
  add_stage6_check <- function(check_id, pass, expected, observed, details = "") {
    stage6_checks[[length(stage6_checks) + 1L]] <<- data.table(
      check_id = check_id, pass = isTRUE(pass), expected = as.character(expected),
      observed = as.character(observed), details = as.character(details)
    )
  }
  # Verify each upstream file by both size and SHA-256 hash.
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

  # Confirm that the current primary results reproduce the sealed reference sets.
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

  ## 04.2 Refit every eligible series and calculate diagnostic metrics ---------
  # Recreate the mean model using the distribution selected in the primary screen.
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

  # Return one complete diagnostic row, including explicit refit failures.
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

  # Apply the same diagnostic calculation to each class or drug series.
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

  ## 04.3 Select a focused set for detailed residual review -------------------
  # Automated metrics cover all series; detailed residuals use prespecified groups.
  selection_parts <- list()
  # Add existing requested codes to the detailed-review selection.
  add_selection <- function(inventory, codes, reason) {
    rows <- inventory[series_code %in% as.character(codes),
                      .(level, series_code, series_name)]
    if (nrow(rows)) {
      rows[, selection_reason := reason]
      selection_parts[[length(selection_parts) + 1L]] <<- rows
    }
  }
  # Select the largest finite values of one diagnostic metric.
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

  ## 04.4 Extract month-level residuals and autocorrelation --------------------
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

  ## 04.5 Reconcile the sequential cohort flow -------------------------------
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
    # Construct one ordered cohort-flow row from a unique set of codes.
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

  # Count overlaps between the two eligibility rules without double-counting.
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
  ## 04.6 Account for missing, invalid, and unmatched analytical values --------
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

  ## 04.7 Write diagnostics and seal the stage after all gates pass ------------
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
  }, envir = state)
  invisible(TRUE)
}
