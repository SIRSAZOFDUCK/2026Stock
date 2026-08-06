# Primary eligibility, modelling, characterisation, and tabular outputs.
run_primary_analysis <- function(state) {
  evalq({
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

  }, envir = state)
  invisible(TRUE)
}
