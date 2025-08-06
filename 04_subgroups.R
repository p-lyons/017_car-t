
# Unified Subgroup Analysis Script ---------------------------------------------------

## Setup ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(broom.mixed); library(gtsummary); library(patchwork); library(tidytable)
  library(lubridate); library(collapse); library(stringr); library(logistf)
  library(ggplot2); library(EValue); library(arrow); library(here); library(lme4)
  library(margins); library(sandwich); library(survminer); library(survival)
  library(cmprsk); library(coxme); library(scales); library(gt)
  library(dplyr); library(tidyr); library(purrr); library(forcats); library(RColorBrewer)
})

options(dplyr.summarise.inform = FALSE)
today <- format(Sys.Date(), "%Y%m%d")

source(here("code/02_analysis/02_model_functions.R"))

## Load Data ---------------------------------------------------------------------
files  <- list.files(here("clean"), pattern = "^analysis_\\d{6}\\.csv$", full.names = TRUE)
dates  <- str_extract(files, "\\d{6}")
latest <- files[which.max(dates)]

safely_fread <- purrr::possibly(data.table::fread, otherwise = data.frame())
df_main_90 <- safely_fread(latest) |>
  mutate(
    cancer_cat_4 = case_when(
      cancer_cat %in% c("MCL", "Follicular Lymphoma") ~ "Other Lymphomas",
      TRUE ~ cancer_cat
    ),
    cancer_cat_3 = case_when(
      cancer_cat %in% c("DLBCL", "Other Lymphomas") ~ "Lymphoma",
      TRUE ~ cancer_cat
    ),
    cancer_cat_2 = if_else(cancer_cat == "Multiple Myeloma", "Multiple Myeloma", "Leukemia/Lymphoma"),
    target = if_else(cancer_cat_2 == "Multiple Myeloma", "BCMA", "CD-19"),
    Sex = factor(if_else(female_01 == 1, "Female", "Male"), levels = c("Female", "Male")),
    surv_90 = if_else(o_d90_01 == 1, 0L, 1L),
    efs_90 = if_else(o_e90_01 == 1, 0L, 1L)
  )

stopifnot("admin_dttm" %in% names(df_main_90))

df_main_365 <- df_main_90 |>
  filter(admin_dttm <= as.POSIXct("2024-07-01")) |>
  mutate(
    surv_365 = if_else(o_d365_01 == 1, 0L, 1L),
    efs_365  = if_else(o_e365_01 == 1, 0L, 1L)
  )

message("Loaded file: ", latest)
rm(files, dates, latest)
gc()

## Variable Definitions ----------------------------------------------------------

exposures <- c("e_hours_precise", "e_hours_sunrise")

outcomes_90 <- c(
  "surv_90", "efs_90", "o_imv_01", "o_vasoactive_01", "o_infx_any_01",
  "o_crs_01", "o_crs34_01", "o_icans_01", "o_icans34_01",
  "o_icu_01", "tocilizumab_01", "anakinra_01"
)

core_outcomes_90 <- c(
  "o_os_90", "o_efs_90", "o_icu_01", "o_imv_01", "o_vasoactive_01",
  "o_infx_any_01", "o_crs34_01", "o_icans34_01", "tocilizumab_01", "anakinra_01"
)

core_outcomes_365 <- c("o_os_365", "o_efs_365")

covariates <- c(
  "age", "nhw_01", "female_01", "ecog", "vw", "saps2",
  "cancer_cat", "product_cat", "cell_dose_std", 
  "conditioning_flucy_01", "los_pre_cart", "year_cat", 
  "season", "hospital", "ldh"
)

outcome_labels <- c(
  surv_90 = "Overall Survival at 90 Days", surv_365 = "Overall Survival at 365 Days", 
  efs_90 = "Event-Free Survival at 90 Days", efs_365 = "Event-Free Survival at 365 Days",
  o_os_90 = "90-Day OS", o_os_365 = "365-Day OS", 
  o_efs_90 = "90-Day EFS", o_efs_365 = "365-Day EFS",
  o_icu_01 = "ICU Admission", o_imv_01 = "Mechanical Ventilation",
  o_vasoactive_01 = "Vasopressor Use", o_infx_any_01 = "Any Infection",
  o_crs01_01 = "Any CRS", o_crs34_01 = "Grade 3/4 CRS", 
  o_icans01_01 = "Any ICANS", o_icans34_01 = "Grade 3/4 ICANS",
  tocilizumab_01 = "Tocilizumab", anakinra_01 = "Anakinra"
)

get_adjusted_covariates <- function(group_var, covs) {
  setdiff(covs, group_var)
}

## Subgroup Configuration --------------------------------------------------------
subgroup_configs <- list(
  outpatient = list(
    var = "outpt_01",
    labels = c("0" = "Inpatient, n=605", "1" = "Outpatient, n=65"),
    colors = c("Inpatient, n=605" = "purple4", "Outpatient, n=65" = "salmon"),
    data_90 = df_main_90, data_365 = df_main_365
  ),
  hospital = list(
    var = "hospital",
    labels = c("BJH" = "WU, n=363", "OHSU" = "OHSU, n=307"),
    colors = c("WU, n=363" = "#BA0C2F", "OHSU, n=307" = "#0E4D8F"),
    data_90 = df_main_90, data_365 = df_main_365
  ),
  sex = list(
    var = "Sex",
    labels = c("Female" = "Female", "Male" = "Male"),
    colors = c("Female" = "#cd0055", "Male" = "#2975f1"),
    data_90 = df_main_90, data_365 = df_main_365
  ),
  target = list(
    var = "target",
    labels = c("BCMA" = "BCMA", "CD-19" = "CD-19"),
    colors = c("BCMA" = "#2BAE66", "CD-19" = "#EEA47F"),
    data_90 = df_main_90, data_365 = df_main_365
  )
)

## Unified Subgroup Analysis Runner ---------------------------------------------
run_all_subgroup_analyses <- function(subgroup_configs, outcomes, data_key) {
  map_dfr(names(subgroup_configs), function(sg_name) {
    config <- subgroup_configs[[sg_name]]
    data <- config[[data_key]]
    outcomes_to_use <- if (sg_name == "outpatient" && data_key == "data_90") outcomes_90 else outcomes
    
    safely_run_models <- purrr::possibly(run_unified_subgroup_models, otherwise = tibble())
    
    unadj <- safely_run_models(
      data = data,
      group_var = config$var,
      exposures = exposures,
      outcomes = outcomes_to_use,
      model_type = "unadjusted"
    ) |> mutate(model = "Unadjusted", analysis = "Subgroup")
    
    adj_covs <- get_adjusted_covariates(config$var, covariates)
    
    adj <- safely_run_models(
      data = data,
      group_var = config$var,
      exposures = exposures,
      outcomes = outcomes_to_use,
      covariates = adj_covs,
      model_type = "adjusted",
      use_random = TRUE
    )
    
    if (exists("compute_evalue_row", mode = "function")) {
      adj <- adj |> mutate(evals = pmap(list(or_adj, ci_lo, ci_hi), compute_evalue_row)) |> unnest_wider(evals)
    }
    
    bind_rows(unadj, adj) |> mutate(subgroup_type = sg_name)
  })
}

## Run Subgroup Analyses --------------------------------------------------------
all_results_90  <- run_all_subgroup_analyses(subgroup_configs, core_outcomes_90,  "data_90")
all_results_365 <- run_all_subgroup_analyses(subgroup_configs, core_outcomes_365, "data_365")

## Combine All Results ----------------------------------------------------------
all_results <-
  bind_rows(all_results_90, all_results_365) |> 
  mutate(or = coalesce(or_adj, or_unadj)) |>
  select(analysis, subgroup_type, subgroup_var, subgroup_val, model,
         exposure, outcome, or, ci_lo, ci_hi, e_value, e_value_lo)

fwrite(all_results, here("output", paste0("unified_subgroup_results_", today, ".csv")))

## Plotting ----------------------------------------------------------------------

# (1) KM Curves by subgroup using 1500 cutoff
plot_km_by_subgroup <- function(data, time_var, event_var, subgroup_var, cutoff = 1500) {
  data <- data |> mutate(time_group = if_else(e_hours_precise >= cutoff, "Late", "Early"))
  
  surv_obj <- Surv(time = !!sym(time_var), event = !!sym(event_var))
  
  ggsurvplot(
    survfit(surv_obj ~ time_group, data = data),
    data = data,
    pval = TRUE,
    legend.title = "Infusion Time",
    legend.labs = c("Early", "Late"),
    palette = c("#2975f1", "#cd0055"),
    ggtheme = theme_minimal()
  )
}

# (2) AOR dot-whisker plots
plot_dotwhisker <- function(results, subgroup_name, group_labels, colors) {
  results |> filter(subgroup_type == subgroup_name, model == "Adjusted") |>
    mutate(outcome_label = outcome_labels[outcome]) |>
    ggplot(aes(x = or, y = fct_rev(outcome_label), xmin = ci_lo, xmax = ci_hi, color = as.factor(subgroup_val))) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(position = position_dodge(width = 0.5), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(x = "Adjusted Odds Ratio", y = NULL, color = subgroup_name)
}

# (3) Marginal effects plots for sunrise model on OS, with p-value
plot_marginal_effects <- function(data, subgroup_var, subgroup_val, covs) {
  df_sub <- data |> filter(!!sym(subgroup_var) == subgroup_val)
  
  model <- glm(o_os_90 ~ e_hours_sunrise + age + saps2 + vw + ecog + female_01 + nhw_01 + year_cat + hospital,
               data = df_sub, family = binomial())
  
  pval <- summary(model)$coefficients["e_hours_sunrise", "Pr(>|z|)"]
  
  new_data <- tibble(e_hours_sunrise = seq(min(df_sub$e_hours_sunrise, na.rm=TRUE),
                                           max(df_sub$e_hours_sunrise, na.rm=TRUE), length.out = 100))
  for (v in covs) new_data[[v]] <- median(df_sub[[v]], na.rm = TRUE)
  
  new_data$pred <- predict(model, newdata = new_data, type = "response")
  
  ggplot(new_data, aes(x = e_hours_sunrise, y = pred)) +
    geom_line(size = 1.2) +
    labs(
      title = paste("Marginal Effect of Sunrise Time in", subgroup_val),
      subtitle = paste("p-value:", signif(pval, 3)),
      x = "Hours Since Sunrise",
      y = "Predicted Probability of 90-Day OS"
    ) +
    theme_minimal()
}

## Loop through and generate plots ----------------------------------------------
for (sg_name in names(subgroup_configs)) {
  cfg <- subgroup_configs[[sg_name]]
  
  # KM
  km_plot <- plot_km_by_subgroup(cfg$data_90, time_var = "o_os_90", event_var = "surv_90", subgroup_var = cfg$var)
  ggsave(here("figures", paste0("kmcurve_", sg_name, "_", today, ".pdf")), km_plot$plot, width = 8, height = 6)
  
  # Dot-whisker
  dw_plot <- plot_dotwhisker(all_results_90, sg_name, cfg$labels, cfg$colors)
  ggsave(here("figures", paste0("aor_dotwhisker_", sg_name, "_", today, ".pdf")), dw_plot, width = 8, height = 6)
  
  # Marginal effects
  for (sg_val in unique(cfg$data_90[[cfg$var]])) {
    me_plot <- plot_marginal_effects(cfg$data_90, cfg$var, sg_val, covariates)
    ggsave(here("figures", paste0("meplot_", sg_name, "_", sg_val, "_", today, ".pdf")), me_plot, width = 8, height = 6)
  }
}

message("All subgroup analysis outputs completed successfully.")
