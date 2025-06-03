

# Subgroup Analysis by Cancer / CAR T-cell target


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
    ),
    target = if_else(cancer_cat_2 == "Multiple Myeloma", "BCMA", "CD-19")
  )

df_main_365 = fsubset(df_main_90, admin_date <= as.Date("2024-04-01"))

rm(files, dates, latest)
gc()

# interaction model ------------------------------------------------------------

## more libraries --------------------------------------------------------------

library(margins)
library(sandwich)

## variable vectors ------------------------------------------------------------

exposures = c(
  "e_hours_precise",
  "e_hours_sunrise"
)

outcomes = c(
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
  "nhw_01",
  "ecog",
  "vw",
  "saps2",
  "female_01",
  # "cell_dose_std",
  # "conditioning_flucy_01",
  "outpt_01",
  "los_pre_cart",
  "year_cat",
  "season",
  "hospital"
)

## function to fit logit model + cluster SEs + extract interaction term --------

fit_marginal_logit = function(o, e, data) {

  ### fit model
  rhs          = paste(c(e, "target", paste0(e, ":target"), covariates), collapse = " + ")
  f            = as.formula(paste(o, "~", rhs))
  model        = glm(f, data = data, family = binomial(link = "logit"))
  vcov_robust  = sandwich::vcovCL(model, cluster = ~hospital)
  term_pattern = paste0(e, ":targetCD-19")
  tidy_all     = broom::tidy(model, conf.int = T, vcov = vcov_robust)

  tidy_interact =
    fsubset(tidy_all, term == term_pattern) |>
    mutate(outcome = o, exposure = e)

  ### marginal predicted probabilities
  me     = margins(model, variables = e, at = list(target = c("CD-19", "BCMA")))
  me_sum = summary(me)

  return(list(
    interaction_term = tidy_interact,
    marginal_summary = me_sum
  ))
}

## run function over all permutations ------------------------------------------

model_grid = tidyr::crossing(o = outcomes, e = exposures)

interaction_table =
  model_grid |>
  mutate(res = map2(o, e, ~fit_marginal_logit(.x, .y, df_main_90)))

interaction_term =
  interaction_table |>
  mutate(out = map(res, ~ .x$interaction_term)) |>
  unnest(out) |>
  select(outcome, exposure, term, estimate, p.value)

marginal_summary =
  interaction_table |>
  mutate(margins = map(res, ~ .x$marginal_summary)) |>
  select(o, e, margins) |>
  unnest(margins)

fsubset(marginal_summary, o == "o_os_90")

# plot marginal effects by target -------------------------------------------------

library(rlang)

plot_predicted_probs_by_ca = function(outcome, exposure, data, covariates) {
  # 1. Build formula with interaction
  rhs = paste(c(exposure, "target", paste0(exposure, ":target"), covariates), collapse = " + ")
  f   = as.formula(paste(outcome, "~", rhs))

  # 2. Fit model
  model = glm(f, data = data, family = binomial(link = "logit"))

  # 3. Create prediction grid
  exposure_seq = seq(min(data[[exposure]], na.rm = TRUE),
                     max(data[[exposure]], na.rm = TRUE), by = 0.1)

  # Typical values for covariates
  cov_summary = data |>
    summarise(across(all_of(covariates), \(x) {
      if (is.numeric(x)) median(x, na.rm = TRUE)
      else {
        ux = na.omit(unique(x))
        ux[which.max(tabulate(match(x, ux)))]
      }
    }))

  # Expand grid by exposure and target
  newdata = expand_grid(
    !!sym(exposure) := exposure_seq,
    target = c("CD-19", "BCMA")
  )

  # Add covariates to the grid (replicate values)
  n_rows = nrow(newdata)
  newdata = bind_cols(newdata, cov_summary[rep(1, n_rows), ])

  # 4. Predict: get fit + SE on link scale
  pred_link = predict(model, newdata = newdata, type = "link", se.fit = TRUE)

  newdata = newdata |>
    mutate(
      fit_link = pred_link$fit,
      se_link  = pred_link$se.fit,
      ci_lo    = plogis(fit_link - 1.96*se_link),
      ci_hi    = plogis(fit_link + 1.96*se_link),
      prob     = plogis(fit_link)
    )

  # 5. Plot with ribbon
  ggplot(newdata, aes_string(x = exposure, y = "prob", color = "target", fill = "target")) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw(base_size = 14) +
    facet_wrap(~target) +
    theme(legend.position = "NULL")
}

p_value =
  interaction_term |>
  fsubset(outcome == "o_os_90" & exposure == "e_hours_precise") |>
  pull(p.value) |>
  signif(3)

t_marginal_effects =
  plot_predicted_probs_by_ca(
    outcome    = "o_os_90",
    exposure   = "e_hours_precise",
    data       = df_main_90,
    covariates = covariates
  ) +
  xlab("Hour of Day") +
  ylab("Adjusted 90-Day Overall Survival") +
  scale_x_continuous(breaks = seq(8, 18, 1)) +
  scale_fill_manual(values  = c("BCMA" = "#2BAE66", "CD-19" = "#EEA47F")) +
  scale_color_manual(values = c("BCMA" = "#2BAE66", "CD-19" = "#EEA47F")) +
  coord_cartesian(xlim = c(8, 17)) +
  annotate(
    "text",
    x     = 8.0,
    y     = 0.4,
    label = paste0("p (interaction) = ", p_value),
    hjust = 0,
    vjust = 1,
    size  = 5
  )

t_marginal_effects

ggsave(
  here("figures", paste0("me_target_", today, ".pdf")),
  height = 8,
  width  = 11,
  units  = "in",
  dpi    = 600
)

# dot-whisker plot by target ---------------------------------------------------

## function to run independent models by target --------------------------------

fit_glmer_by_ = function(l, e, o, c, data) {

  df_target = data |> filter(target == l)

  rhs    = paste(c(e, c, "(1 | hospital)"), collapse = " + ")
  f      = as.formula(paste(o, "~", rhs))
  m      = glmer(f, data = df_target, family = binomial)
  tidy_m = broom.mixed::tidy(m, effects = "fixed", exponentiate = T, conf.int = T)
  out_row = fsubset(tidy_m, term == e)

  tidytable(
    outcome  = o,
    exposure = e,
    target   = l,
    aor      = out_row$estimate,
    ci_lo    = out_row$conf.low,
    ci_hi    = out_row$conf.high,
    p_value  = out_row$p.value

  )
}

## fit models for main binary outcomes -----------------------------------------

c  = c(
  "age",
  "nhw_01",
  "ecog",
  "vw",
  "saps2",
  "female_01",
  # "cell_dose_std",
  # "conditioning_flucy_01",
  "outpt_01",
  "los_pre_cart",
  "year_cat",
  "season"
)

grid = expand.grid(
  target           = c("CD-19", "BCMA"),
  exposure         = exposures,
  outcome          = outcomes,
  stringsAsFactors = FALSE
)

results_stratified =
  pmap_dfr(
    grid,
    function(target, exposure, outcome) {
      tryCatch(
        fit_glmer_by_(target, exposure, outcome, covariates, df_main_90),
        error = function(e) tidytable(
          outcome  = outcome,
          exposure = exposure,
          target   = target,
          aor      = NA_real_,
          ci_lo    = NA_real_,
          ci_hi    = NA_real_
        )
      )
    }
  )

## fit models for 365-day binary outcomes --------------------------------------

grid_365 = expand.grid(
  target           = c("CD-19", "BCMA"),
  exposure         = exposures,
  outcome          = c("o_os_365", "o_efs_365"),
  stringsAsFactors = FALSE
)

results_365 =
  pmap_dfr(
    grid_365,
    function(target, exposure, outcome) {
      tryCatch(
        fit_glmer_by_(target, exposure, outcome, covariates, df_main_90),
        error = function(e) tidytable(
          outcome  = outcome,
          exposure = exposure,
          target   = target,
          aor      = NA_real_,
          ci_lo    = NA_real_,
          ci_hi    = NA_real_
        )
      )
    }
  )

surv_labels = c(
  "90-Day OS",
  "365-Day OS",
  "90-Day EFS",
  "365-Day EFS"
)

results_stratified =
  rowbind(results_stratified, results_365) |>
  mutate(
    outcome = recode(
      outcome,
      "o_os_90"         = "90-Day OS",
      "o_os_365"        = "365-Day OS",
      "o_efs_90"        = "90-Day EFS",
      "o_efs_365"       = "365-Day EFS",
      "o_crs34_01"      = "Severe CRS",
      "o_icans34_01"    = "Severe ICANS",
      "o_infx_any_01"   = "Infection",
      "o_icu_01"        = "ICU Admission",
      "o_imv_01"        = "Mechanical Ventilation",
      "o_vasoactive_01" = "Vasopressor Use",
      "anakinra_01"     = "Anakinra",
      "tocilizumab_01"  = "Tocilizumab"
    ),
    outcome  = factor(outcome, levels = rev(unique(outcome))),
    exposure = recode(
      exposure,
      "e_hours_precise" = "Each Hour of Day",
      "e_hours_sunrise" = "Each Hour Since Sunrise"
    ),
    category = case_when(
      outcome %in% surv_labels                  ~ "Survival Outcomes",
      outcome %in% c("Tocilizumab", "Anakinra") ~ "Medications",
      TRUE                                      ~ "Complications"
    ),
    category = factor(category, levels = c("Survival Outcomes", "Complications", "Medications"))
  )

plot_target =
  ggplot(results_stratified, aes(x = aor, y = outcome, color = target)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi),
    height   = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_x_log10() +
  scale_color_manual(values = c("BCMA" = "#2BAE66", "CD-19" = "#EEA47F")) +
  labs(
    x     = "Adjusted Odds Ratio (per Hour)",
    y     = NULL,
    color = "CAR T-cell Target"
  ) +
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

plot_target

ggsave(
  here("figures", paste0("dw_target_", today, ".pdf")),
  height = 8,
  width  = 11,
  units  = "in",
  dpi    = 600
)

# Kaplan-Meier by target -------------------------------------------------------

## additional libraries --------------------------------------------------------

library(survminer)
library(survival)
library(cmprsk)
library(coxme)

## most recent survival data ---------------------------------------------------

files  = list.files(here("clean"), pattern = "^df_survival_\\d{6}\\.parquet$", full.names = TRUE)
dates  = str_extract(files, "\\d{6}")
latest = files[which.max(dates)]
df_s   = read_parquet(latest)
message("Loaded file: ", latest)

## make separate sex datasets --------------------------------------------------

target = select(df_main_90, mrn, target)
df_s   = join(target, df_s, how = "left", multiple = F)
bcma   = fsubset(df_s, target == "BCMA")
cd19   = fsubset(df_s, target == "CD-19")

## create survival objects for KM analysis -------------------------------------

b_surv = Surv(time = bcma$overall_time, event = bcma$overall_status)
b_fit  = survfit(b_surv ~ e_late_f15_01, data = bcma)
c_surv = Surv(time = cd19$overall_time, event = cd19$overall_status)
c_fit  = survfit(c_surv ~ e_late_f15_01, data = cd19)

## plot side-by-side -----------------------------------------------------------

### female ---------------------------------------------------------------------

s_b =
  ggsurvplot(
    b_fit,
    data         = bcma,
    risk.table   = TRUE,
    pval         = TRUE,
    conf.int     = FALSE,
    palette      = c("#F1A226", "#800072"),
    legend.title = "Time of Infusion",
    legend.labs  = c("Before 15:00", "After 15:00"),
    xlim         = c(0, 720),
    break.x.by   = 180,
    title        = "Overall Survival (BCMA)"
  )

### male -----------------------------------------------------------------------

s_c =
  ggsurvplot(
    c_fit,
    data         = cd19,
    risk.table   = TRUE,
    pval         = TRUE,
    conf.int     = FALSE,
    palette      = c("#F1A226", "#800072"),
    legend.title = "Time of Infusion",
    legend.labs  = c("Before 15:00", "After 15:00"),
    xlim         = c(0, 720),
    break.x.by   = 180,
    title        = "Overall Survival (CD19)"
  )

### pair and save --------------------------------------------------------------

surv_pair =
  arrange_ggsurvplots(
    list(s_b, s_c),
    ncol              = 2,
    nrow              = 1,
    risk.table.height = 0.25,
    common.legend     = TRUE
  )

ggsave(
  filename = here("figures", paste0("km_target_", today, ".pdf")),
  plot     = surv_pair[[1]],
  width    = 11,
  height   = 8,
  units    = "in",
  dpi      = 600
)

# RMST by target ------------------------------------------------------------------

## load most recent survival data ----------------------------------------------

files  = list.files(here("clean"), pattern = "^df_survival_\\d{6}\\.parquet$", full.names = T)
dates  = str_extract(files, "\\d{6}")
latest = files[which.max(dates)]

df =
  read_parquet(latest) |>
  join(df_main_90, how = "left", multiple = F) |>
  select(
    target,
    all_of(exposures),
    all_of(covariates),
    ends_with("_time"),
    ends_with("_status")
  )

## time/status vars ------------------------------------------------------------

rmst_sets = list(
  overall = c(time = "overall_time", status = "overall_status"),
  efs     = c(time = "efs_time",     status = "efs_status")
)

## rmst interaction function ---------------------------------------------------

rmst_interaction = function(e, time_v, status_v, c, data) {
  data$phi = pseudo::pseudomean(data[[time_v]], data[[status_v]], tmax = 730)

  f = as.formula(paste0("phi ~ ", e, " * target + ", paste(c, collapse = " + "), " + (1|hospital)"))
  m = lmer(f, data = data)

  terms_present    = names(fixef(m))
  interaction_term = paste0(e, ":targetCD-19")
  main_term        = e

  if (!(interaction_term %in% terms_present)) {
    warning("Interaction term not found in model: ", interaction_term)
    return(tidytable())
  }

  # Bootstrap both male (main effect) and interaction term
  boot = bootMer(
    m,
    FUN = function(x) c(
      cd19   = fixef(x)[[main_term]],
      bcma   = fixef(x)[[main_term]] + fixef(x)[[interaction_term]],
      diff   = fixef(x)[[interaction_term]]
    ),
    nsim = 1000,
    use.u = FALSE,
    type = "parametric",
    seed = 2025
  )

  boot_df = as_tibble(boot$t, .name_repair = "minimal") |> setNames(c("cd19", "bcma", "diff"))

  est_cd19   = fixef(m)[[main_term]]
  est_bcma   = est_cd19 + fixef(m)[[interaction_term]]
  est_diff   = fixef(m)[[interaction_term]]

  ci_cd19   = quantile(boot_df$cd19,   probs = c(0.025, 0.975), na.rm = TRUE)
  ci_bcma = quantile(boot_df$bcma, probs = c(0.025, 0.975), na.rm = TRUE)
  ci_diff   = quantile(boot_df$diff,   probs = c(0.025, 0.975), na.rm = TRUE)

  tidytable(
    exposure        = e,
    outcome         = time_v,
    target          = c("CD-19", "BCMA", "CD-19-BCMA"),
    rmst_per_h      = c(est_cd19, est_bcma, est_diff),
    ci_lo           = c(ci_cd19[1], est_bcma[1], ci_diff[1]),
    ci_hi           = c(ci_cd19[2], est_bcma[2], ci_diff[2]),
    type            = c("stratum", "stratum", "interaction")
  )
}

## grid of target × exposure × outcome --------------------------------------------

interaction_grid = expand.grid(
  exposure         = exposures,
  outcome          = names(rmst_sets),
  stringsAsFactors = FALSE
)

rmst_combined = pmap_dfr(
  interaction_grid,
  function(exposure, outcome) {
    vars = rmst_sets[[outcome]]
    rmst_interaction(exposure, vars["time"], vars["status"], covariates, df)
  }
)

plot_rmst_by_target = function(data, outcome_filter) {
  data |>
    filter(outcome == outcome_filter, type == "stratum") |>
    ggplot(aes(x = target, y = rmst_per_h, color = target)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15) +
    facet_wrap(~exposure, scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = paste("RMST Effect Per Hour by target —", outcome_filter),
      y = "Δ RMST per Hour (days)",
      x = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
}

plot_rmst_by_target(rmst_combined, "overall_time")

