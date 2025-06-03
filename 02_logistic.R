
# setup ------------------------------------------------------------------------

## libraries -------------------------------------------------------------------

library(broom.mixed)
library(gtsummary)
library(patchwork)
library(tidytable)
library(collapse)
library(stringr)
library(ggplot2)
library(scales)
library(arrow)
library(here)
library(lme4)
library(gt)

## helpers and functions ---------------------------------------------------

today  = format(Sys.Date(), "%y%m%d")
files  = list.files(here("clean"), pattern = "^analysis_\\d{6}\\.csv$", full.names = T)
dates  = str_extract(files, "\\d{6}")
latest = files[which.max(dates)]

## load data -------------------------------------------------------------------

df_main_90 =
  fread(latest) |>
  fmutate(
    cancer_cat_4 = case_when(
      cancer_cat %in% c("MCL", "Follicular Lymphoma") ~ "Other Lymphomas",
      TRUE                                            ~ cancer_cat
    ),
    cancer_cat_3 = case_when(
      cancer_cat %in% c("DLBCL", "Other Lymphomas") ~ "Lymphoma",
      TRUE                                          ~ cancer_cat
    ),
    cancer_cat_2 = case_when(
      cancer_cat == "Multiple Myeloma" ~ "Multiple Myeloma",
      TRUE                             ~ "Leukemia/Lymphoma"
    )
  )

df_main_365 = fsubset(df_main_90, admin_date <= as.Date("2024-04-01"))

rm(files, dates, latest)
gc()

# model fitting ----------------------------------------------------------------

## variables -------------------------------------------------------------------

exposure_vars = c(
  "e_hours_precise",
  "e_hours_sunrise"
)

outcome_vars = c(
  "o_os_90",
  "o_efs_90",
  "o_icu_01",
  "o_imv_01",
  "o_vasoactive_01",
  "o_infx_any_01",
  "o_crs34_01",
  "o_icans34_01",
  "tocilizumab_01",
  "anakinra_01"
)

outcome_365 = c(
  "o_os_365",
  "o_efs_365"
)

covariates = c(
  "age",
  "female_01",
  "nhw_01",
  "ecog",
  "vw",
  "saps2",
  "cancer_cat",
  "product_cat",
  "cell_dose_std",
  "conditioning_flucy_01",
  "outpt_01",
  "los_pre_cart",
  "year_cat",
  "season"
)

## fit models ------------------------------------------------------------------

for (e in exposure_vars) {
  for (o in c(outcome_vars, outcome_365)) {

    df = if (o %in% outcome_365) df_main_365 else df_main_90
    f  = paste0( o, " ~ ", e, " + ", paste(covariates, collapse = " + "), " + (1 | hospital)")
    m  = glmer(formula = as.formula(f), data = df, family = binomial)
    os = gsub("o_|_01|", "", o)
    es = gsub("e_hours_", "", e)
    mn = paste0("m_", os, "_", es)
    assign(mn, m, envir = .GlobalEnv)
  }
}

## ARR for 90-day OS -----------------------------------------------------------

### template -------------------------------------------------------------------

tpl                 = model.frame(m_os_90_precise)[1, ]
tpl$e_hours_precise = 8
tpl$hospital        = NA
hours               = seq(8, 15, by = 1)
crit                = qnorm(0.975)

newdata =
  bind_rows(
    lapply(hours, function(h) {
      row = tpl
      row$e_hours_precise = h
      row
    })
  ) |>
  fmutate(Hour = hours)

### make predictions -----------------------------------------------------------

pred_link =
  predict(
    m_os_90_precise,
    newdata = newdata,
    type    = "link",
    se.fit  = T,
    re.form = NA
  )

preds =
  newdata |>
  fmutate(
    fit_link  = pred_link$fit,
    se_link   = pred_link$se.fit,
    prob      = plogis(fit_link),
    prob_lo   = plogis(fit_link - crit*se_link),
    prob_hi   = plogis(fit_link + crit*se_link)
  )

### arr vs baseline (0800) -----------------------------------------------------

baseline = preds |> fsubset(Hour == 8) |> pull(prob)
preds    = preds |> fmutate(ARR = prob-baseline)
mean_arr = fmean(preds$ARR)

### bootstrap CI ---------------------------------------------------------------

set.seed(2025)

arr_boot =
  replicate(1000, {
    sample_rows = sample(nrow(preds), replace = T)
    fmean(preds$ARR[sample_rows])
  })

ci_arr = quantile(arr_boot, c(0.025, 0.975), na.rm = T)

arr_os90 =
  tidytable(
    arr_mean = mean_arr,
    ci_lo    = ci_arr[1],
    ci_hi    = ci_arr[2]
  )

arr_os90

# marginal effects plots -------------------------------------------------------

## set up grids ----------------------------------------------------------------

hour_grids =
  list(
    e_hours_precise = seq(8, 18, by = 0.1),
    e_hours_sunrise = seq(0, 11, by = 0.1)
  )

model_names      = ls(pattern = "^m_")
pred_list        = vector("list", length(model_names))
names(pred_list) = model_names

## make predictions for each model ---------------------------------------------

### extract model details ------------------------------------------------------

for (i in seq_along(model_names)) {

  model_name    = model_names[i]
  model         = get(model_name)
  exposure_var  = grep("e_hours_", names(model.frame(model)), value = T)
  stopifnot(length(exposure_var) == 1)
  outcome_var   = all.vars(formula(model))[1]
  fixed_effects = attr(terms(model), "term.labels")
  covariates    = setdiff(fixed_effects, exposure_var)

  ref =
    model.frame(model) |>
    summarise(across(
      all_of(covariates),
      ~ if (is.numeric(.x)) fmedian(.x) else factor(names(which.max(table(.x))), levels = levels(.x))
    ))

  time_seq = hour_grids[[exposure_var]]

  newdat =
    ref[rep(1, length(time_seq)), ] |>
    mutate(!!exposure_var := time_seq)

### make predictions -----------------------------------------------------------

  pred = predict(model, newdata = newdat, type = "link", se.fit = T, re.form = NA)
  crit = qnorm(0.975)

### exposure aOR p-values ------------------------------------------------------

  pval =
    broom.mixed::tidy(model, effects = "fixed") |>
    fsubset(term == exposure_var) |>
    pull(p.value)

### data frame of predictions --------------------------------------------------

  pred_df =
    tidytable(
      exposure = exposure_var,
      outcome  = outcome_var,
      time     = time_seq,
      prob     = plogis(pred$fit),
      ci_lo    = plogis(pred$fit - crit * pred$se.fit),
      ci_hi    = plogis(pred$fit + crit * pred$se.fit),
      p_value  = pval
    )

  pred_list[[i]] = pred_df
}

all_preds = rowbind(pred_list)


## make marginal effects plots -------------------------------------------------

### outcome levels -------------------------------------------------------------

surv_levels = c(
  "o_os_90",
  "o_os_365",
  "o_efs_90",
  "o_efs_365"
)

comp_levels = c(
  "o_crs34_01",
  "o_icans34_01",
  "o_infx_any_01",
  "o_icu_01",
  "o_imv_01",
  "o_vasoactive_01"
)

### outcome labels -------------------------------------------------------------

surv_labels = c(
  "90-Day Overall Survival",
  "365-Day Overall Survival",
  "90-Day Event-Free Survival",
  "365-Day Event-Free Survival"
)

comp_labels = c(
  "Severe CRS",
  "Severe ICANS",
  "Infection",
  "ICU Admission",
  "Mechanical Ventilation",
  "Vasopressor Use"
)

### separate survival and complication outcomes --------------------------------

preds_surv =
  fsubset(all_preds, str_detect(outcome, "s_")) |>
  fmutate(outcome  = factor(outcome,  levels = surv_levels, labels = surv_labels)) |>
  fmutate(exposure = factor(exposure, labels = c("Hour of Day", "Hours Since Sunrise")))

preds_comp =
  fsubset(all_preds, outcome %in% comp_levels) |>
  fmutate(outcome  = factor(outcome,  levels = comp_levels, labels = comp_labels)) |>
  fmutate(exposure = factor(exposure, labels = c("Hour of Day", "Hours Since Sunrise")))

### extract p-values for plots -------------------------------------------------

pval_surv =
  distinct(preds_surv, exposure, outcome, p_value) |>
  fmutate(lab = paste0("p = ", signif(p_value, 2))) |>
  fmutate(x   = -Inf) |>
  fmutate(y   =  Inf)

pval_comp =
  distinct(preds_comp, exposure, outcome, p_value) |>
  fmutate(lab = paste0("p = ", signif(p_value, 2))) |>
  fmutate(x   = -Inf) |>
  fmutate(y   =  Inf)

### plot survival outcomes -----------------------------------------------------

me_surv =
  ggplot(preds_surv, aes(x = time, y = prob)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2) +
  geom_line() +
  geom_text(
    data = pval_surv, aes(x = x, y = y, label = lab),
    hjust = -0.1,
    vjust = 1.1,
    size  = 3
  ) +
  facet_grid(rows = vars(outcome), cols = vars(exposure), scales = "free_x") +
  labs(x = "Time", y = "Adjusted Probability of Outcome") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

me_surv

ggsave(
  here("figures", paste0("marginal_effects_surv_", today, ".pdf")),
  height = 11,
  width  = 08,
  units  = "in",
  dpi    = 600
)

### plot complications ---------------------------------------------------------

me_comp =
  ggplot(preds_comp, aes(x = time, y = prob)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2) +
  geom_line() +
  geom_text(
    data = pval_comp, aes(x = x, y = y, label = lab),
    hjust = -0.1,
    vjust = 1.1,
    size  = 3
  ) +
  facet_grid(rows = vars(outcome), cols = vars(exposure), scales = "free_x") +
  labs(x = "Time", y = "Adjusted Probability of Outcome") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

me_comp

ggsave(
  here("figures", paste0("marginal_effects_comp_", today, ".pdf")),
  height = 11,
  width  = 08,
  units  = "in",
  dpi    = 600
)

# plot 90-day OS individually for Figure 1 -------------------------------------

p_value =
  fsubset(preds_surv, outcome == "90-Day Overall Survival" & exposure == "Hour of Day") |>
  pull(p_value) |>
  funique()

figure_1b =
  fsubset(preds_surv, outcome == "90-Day Overall Survival" & exposure == "Hour of Day") |>
  ggplot(aes(x = time, y = prob)) +
  geom_line(size = 1.25, color = "black") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "black", alpha = 0.15) +
  annotate(
    "text",
    x      = 8,
    y      = 0.22,
    label  = paste0("p = ", signif(p_value, 2)),
    hjust  = 0,
    vjust  = 1,
    fontface = "bold"
  ) +
  scale_x_continuous(limits = c(7, 18), breaks = seq(8, 18, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(xlim = c(8, 17), ylim = c(0.2, 1)) +
  labs(
    x     = "Infusion Time of Day (Hours)",
    y     = "Adjusted 90-Day Overall Survival",
    title = "B. Adjusted 90-Day Overall Survival"
  ) +
  theme_bw(base_size = 14)

figure_1b

ggsave(
  here("figures", paste0("01b_marginal_surv90_tod_", today, ".pdf")),
  width  = 8,
  height = 8,
  units  = "in",
  dpi    = 600
)

rm(figure_1b, p_value, pval_surv, pval_comp)
gc()

# save odds ratios and e-values for tables/figures -----------------------------

## e-value function ------------------------------------------------------------

library(EValue)

compute_evalue_row = function(or_est, ci_lo, ci_hi) {

  or_num = as.numeric(or_est)
  lo_num = as.numeric(ci_lo)
  hi_num = as.numeric(ci_hi)

  if (any(is.na(c(or_num, lo_num, hi_num))) || lo_num == hi_num) {
    return(c(e_value = NA_real_, e_value_lo = NA_real_, e_value_hi = NA_real_))
  }

  tryCatch({
    ev_mat = EValue::evalues.OR(est = or_num, lo = lo_num, hi = hi_num, rare = F)
    c(
      e_value    = ev_mat["E-values", "point"],
      e_value_lo = ev_mat["E-values", "lower"],
      e_value_hi = ev_mat["E-values", "upper"]
    )
  }, error = function(e) {
    c(e_value = NA_real_, e_value_lo = NA_real_, e_value_hi = NA_real_)
  })
}

## list of model coefficients --------------------------------------------------

results = list()

for (i in seq_along(model_names)) {

  model_name    = model_names[i]
  model         = get(model_name)
  exposure_var  = grep("e_hours_", names(model.frame(model)), value = TRUE)
  stopifnot(length(exposure_var) == 1)
  outcome_var   = all.vars(formula(model))[1]
  fixed_effects = attr(terms(model), "term.labels")
  covariates    = setdiff(fixed_effects, exposure_var)
  tidy_model    = broom.mixed::tidy(model, effects = "fixed", exponentiate = T, conf.int = T)
  exposure_row  = fsubset(tidy_model, term == exposure_var)

  if (nrow(exposure_row) != 1) {
    warning("Could not find unique exposure row in model: ", model_name)
    next
  }

## compile ORs and E-values ----------------------------------------------------

  aor     = exposure_row$estimate
  ci_lo   = exposure_row$conf.low
  ci_hi   = exposure_row$conf.high
  p_value = exposure_row$p.value

  ev = compute_evalue_row(aor, ci_lo, ci_hi)

  results[[i]] = tidytable(
    model_name,
    outcome    = outcome_var,
    exposure   = exposure_var,
    aor        = aor,
    ci_lo      = ci_lo,
    ci_hi      = ci_hi,
    p_value    = p_value,
    e_value    = ev["e_value"],
    e_value_lo = ev["e_value_lo"],
    e_value_hi = ev["e_value_hi"]
  )
}

## table of results ------------------------------------------------------------

final_results =
  bind_rows(results) |>
  select(outcome, exposure, aor, ci_lo, ci_hi, p_value, e_value) |>
  mutate(across(
    .cols = where(is.numeric),
    .fns  = ~signif(.x, 3)
  )) |>
  fmutate(e_value = if_else(p_value <= 0.05, e_value, NA_real_))

fwrite(final_results, here("output", paste0("logistic_aor_", today, ".csv")))

# dot-whisker plot of adjusted odds ratios -------------------------------------

## outcome labels --------------------------------------------------------------

outcome_labels = c(
  o_os_90         = "Overall Survival at 90 Days",
  o_os_365        = "Overall Survival at 365 Days",
  o_efs_90        = "Event-Free Survival at 90 Days",
  o_efs_365       = "Event-Free Survival at 365 Days",
  o_icu_01        = "ICU Admission",
  o_imv_01        = "Mechanical Ventilation",
  o_vasoactive_01 = "Vasopressor Use",
  o_infx_any_01   = "Any Infection",
  o_crs34_01      = "Grade 3/4 CRS",
  o_icans34_01    = "Grade 3/4 ICANS",
  tocilizumab_01  = "Tocilizumab",
  anakinra_01     = "Anakinra"
)

ordered_outcomes = c(
  "o_os_90",
  "o_os_365",
  "o_efs_90",
  "o_efs_365",
  "o_icu_01",
  "o_imv_01",        # Vent
  "o_vasoactive_01", # Vaso
  "o_infx_any_01",   # Infection
  "o_crs34_01",      # CRS
  "o_icans34_01",    # ICANS
  "tocilizumab_01",  # Tocilizumab
  "anakinra_01"      # Anakinra
)

## prepare data for plotting ---------------------------------------------------

final_results =
  fsubset(final_results, outcome %in% ordered_outcomes) |>
  fmutate(
    category = case_when(
      outcome %in% surv_levels                        ~ "Survival Outcomes",
      outcome %in% c("tocilizumab_01", "anakinra_01") ~ "Medications",
      TRUE                                            ~ "Complications"
    ),
    category = factor(category, levels = c("Survival Outcomes", "Complications", "Medications")),
    outcome  = factor(outcome,  levels = rev(ordered_outcomes)) |>
      forcats::fct_relabel(~outcome_labels[.x]),
    exposure = case_when(
      exposure == "e_hours_precise" ~ "Each Hour of Day",
      TRUE                          ~ "Each Hour Since Sunrise"
    )
  ) |>
  tidytable::slice(n = 1, .by = c(outcome, exposure))

# plot only adjusted -----------------------------------------------------------

## just hour of day ------------------------------------------------------------

plot_clock =
  fsubset(final_results, exposure == "Each Hour of Day") |>
  ggplot(aes(x = aor, y = outcome, xmin = ci_lo, xmax = ci_hi)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(x = "Adjusted Odds Ratio (Per Hour)", y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.placement    = "outside"   # move row‐strips to the left
  ) +
  facet_grid(
    rows   = vars(category),
    cols   = vars(exposure),
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  )

plot_clock

ggsave(
  here("figures", paste0("dw_clock_", today, ".pdf")),
  width  = 11,
  height = 8,
  units  = "in",
  dpi    = 600
)

## just hour since sunrise -----------------------------------------------------

plot_sunrise =
  fsubset(final_results, exposure == "Each Hour Since Sunrise") |>
  ggplot(aes(x = aor, y = outcome, xmin = ci_lo, xmax = ci_hi)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(x = "Adjusted Odds Ratio (Per Hour)", y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.placement    = "outside"   # move row‐strips to the left
  ) +
  facet_grid(
    rows   = vars(category),
    cols   = vars(exposure),
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  )

plot_sunrise

ggsave(
  here("figures", paste0("dw_sunrise_", today, ".pdf")),
  width  = 11,
  height = 8,
  units  = "in",
  dpi    = 600
)

## both combined ---------------------------------------------------------------

plot_adj =
  final_results |>
  ggplot(aes(x = aor, y = outcome, xmin = ci_lo, xmax = ci_hi, group = exposure, color = exposure)) +
  geom_point(size       = 3.0, position = position_dodge(width = 0.5)) +
  geom_errorbarh(height = 0.2, position = position_dodge(width = 0.5)) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(x = "Adjusted Odds Ratio (Per Hour)", y = NULL, color = "Exposure") +
  scale_color_manual(
    values = c(
      "Each Hour Since Sunrise" = "gray70",
      "Each Hour of Day"        = "black"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.placement     = "outside"   # move row‐strips to the left
  ) +
  facet_grid(
    rows   = vars(category),
    cols   = vars(exposure),
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  )

plot_adj

ggsave(
  here("figures", paste0("dw_outcomes_both_", today, ".pdf")),
  width  = 11,
  height = 8,
  units  = "in",
  dpi    = 600
)
