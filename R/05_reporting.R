# Module 05: publication figures and the prespecified 2021 window assessment.
run_reporting <- function(state) {
  # Reuse primary objects from the coordinator context; no model definitions change here.
  evalq({
## 05.1 Shared plotting theme and fit-frame helper ----------------------------

theme_pub <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# Return observed and fitted rates per 1,000 registered patients for one series.
.series_fit <- function(mon, id_col, id_val) {
  ser <- mon |> filter(.data[[id_col]] == id_val) |> select(year_month, items)
  d   <- covar |> left_join(ser, by = "year_month") |> arrange(t)
  d$off <- d$offset_log_patient_days
  fit <- suppressWarnings(glm(.f_full, poisson, data = d))
  tibble(month_date = d$month_date,
         observed = d$items / d$list_size * 1000,
         fitted   = as.numeric(fitted(fit)) / d$list_size * 1000)
}

## 05.2 Main-text exemplar and seasonality-landscape figures -----------------

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

## 05.3 Appendix comparison of retained and excluded discoveries --------------

# Resolve named class exemplars to their stable BNF codes.
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

## 05.4 Paginated observed-and-fitted appendix figures ------------------------
# Plot every eligible series. Strip colours distinguish meaningful, significant,
# and eligible-only series, matched to panels by their display labels.

# Read a facet strip's panel label.
.strip_label <- function(sg) {
  gt <- sg$grobs[[1]]
  for (ch in gt$children) if (!is.null(ch$children))
    for (tx in ch$children) if (!is.null(tx$label)) return(tx$label)
  NA_character_
}
# Set one facet strip's background and text colours.
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
# Apply a named label-to-colour lookup across all facet strips.
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

# Tier lookups are keyed by display name; meaningful status takes precedence.
# Drug significance always comes from the complete eligible-drug BH family.
FILL_MEAN <- "#8B0000"; TXT_MEAN <- "white"
FILL_SIG  <- "#E8A0A0"; TXT_SIG  <- "grey15"
tier_class <- {
  m <- setNames(vector("list", nrow(results_class)), results_class$bnf_class_name)
  for (i in seq_len(nrow(results_class)))
    m[[i]] <- if (isTRUE(results_class$meaningful[i])) c(FILL_MEAN, TXT_MEAN) else c(FILL_SIG, TXT_SIG)
  m
}
tier_drug <- {
  sig <- results_drug$drug_significant %in% TRUE
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
  sig <- results_drug$drug_significant %in% TRUE
  m <- setNames(vector("list", nrow(results_drug)), results_drug$bnf_drug_name)
  for (i in seq_len(nrow(results_drug)))
    m[[i]] <- if (isTRUE(results_drug$meaningful[i])) c(FILL_MEAN_SHAPE, TXT_MEAN_SHAPE)
  else if (sig[i])                        c(FILL_SIG_SHAPE, TXT_SIG_SHAPE)
  else                                    NULL
  m[!vapply(m, is.null, logical(1))]
}

# Render every eligible series on fixed 4 × 4 chapter-specific PDF pages.
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

## 05.5 Paginated fitted-seasonal-shape appendix figures ----------------------
# Show the fitted monthly factor for every eligible series. Poisson fits provide
# a common visual summary only; inferential routes remain those selected earlier.
FILL_MEAN_SHAPE <- "#08306B"; TXT_MEAN_SHAPE <- "white"   # dark blue
FILL_SIG_SHAPE  <- "#9ECAE1"; TXT_SIG_SHAPE  <- "grey15"   # pale blue

# Convert one fitted harmonic curve to factors relative to its annual mean.
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

# Build unique per-page labels and matching status colours.
.page_labels <- function(rows, name_col, id_col, tier_map) {
  nm  <- rows[[name_col]]
  dup <- nm %in% nm[duplicated(nm)]
  panel_label <- ifelse(dup, sprintf("%s [%s]", nm, rows[[id_col]]), nm)
  tier <- setNames(tier_map[nm], panel_label)
  list(labels = panel_label, tier = tier[!vapply(tier, is.null, logical(1))])
}

# Render all seasonal-shape panels with the same fixed pagination.
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

### Report created files and headline counts

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
## 05.6 Assess the exclusion of the pandemic-disrupted 2021 window ------------
# This descriptive extension combines cached 2021 class totals with the
# 2022–2025 analysis panel. It does not refit or alter the primary models.

stopifnot(exists("out_dir"), exists("covar"), exists("results_class"),
          exists("class_monthly_elig"))
if (!exists("data_dir")) data_dir <- out_dir
if (!exists("class_monthly"))
  class_monthly <- readRDS(file.path(data_dir, "class_monthly.rds"))

# Named seasonal exemplars are matched exactly; one high-volume flat control is added.
exemplar_names <- c(
  "Penicillins",                 # winter antibiotic (COVID-suppressed 2020-21)
  "Macrolides",                  # winter antibiotic
  "Antihistamines",              # summer allergy (opposite-phase cycle)
  "Cough suppressants",          # winter respiratory
  "Sunscreening preparations"    # summer dermatological
)

### 05.6a Aggregate 2021 prescribing to monthly class totals
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

### 05.6b Convert quarterly 2021 list size to monthly denominators

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

### 05.6c Select named seasonal exemplars and one flat control

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

### 05.6d Build and plot the 2021–2025 rate series

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
  "Window-justification figure complete.\n",
  "  2021 EPD months aggregated: %d/12\n",
  "  Exemplar classes plotted:   %s\n",
  "  Figure: %s\n"),
  n_months_2021,
  paste(facet_levels, collapse = "; "),
  file.path(fig_dir, "figure_window_justification_2021_2025.png")))


### 05.6e Scan 2021 anomalies across all eligible classes
# Compare 2021 with the 2022–2025 level trend and mean within-year shape.
# Reuse the 05.6a class cache and 05.6b denominator without reading new data.

stopifnot(exists("ls_all"), exists("res_dir"))          # built in 05.6b / 05.6d
if (!exists("class_2021")) {
  cf <- list.files(cache_dir_2021, pattern = "^data_byclass_2021_\\d{6}\\.csv$", full.names = TRUE)
  if (!length(cf)) stop("No 2021 class cache found; run subsection 05.6a first.")
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

# Map YYYYMM to a continuous month index beginning at January 2022.
idx <- function(ym) (ym %/% 100L - 2022L) * 12L + (ym %% 100L)

# Compare 2021 level and within-year shape with the 2022–2025 reference pattern.
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

#### Summarise level and shape flags
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

#### Plot the nine largest absolute 2021 level departures
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

### 05.6f Plot every eligible class across 2021–2025
# Plot every eligible class by chapter with 2021 shaded. Strips are red for
# level disruption, blue for shape-only disruption, and grey otherwise.

stopifnot(exists("panel"), exists("scan"), exists("res_dir"))
suppressMessages(library(ggplot2))

# Define strip helpers only if the appendix setup did not already create them.
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

# Colour level disruptions red, shape-only disruptions blue, and others grey.
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

# Put the largest absolute departures first within each chapter.
order_key <- function(dt) dt[order(bnf_chapter_code, -abs(level_dev_pct))]
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
  "All-class window plots complete.\n",
  "  All %d eligible classes plotted 2021-2025 -> %s\n",
  "  Strip colour: red = level-disrupted (%d), blue = shape-only (%d), grey = stable (%d)\n"),
  nrow(sc), out_pdf,
  sum(sc$tier == "level"), sum(sc$tier == "shape"), sum(sc$tier == "none")))
  }, envir = state)
  invisible(TRUE)
}
