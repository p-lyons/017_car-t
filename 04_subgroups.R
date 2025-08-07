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

# Check if model functions file exists
if (file.exists(here("code/02_analysis/02_model_functions.R"))) {
  source(here("code/02_analysis/02_model_functions.R"))
}

## Load Data ---------------------------------------------------------------------
files  <- list.files(here("clean"), pattern = "^analysis_\\d{6}\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No analysis CSV files found in clean/ directory")
}

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
    surv_365 = o_os_365,
    efs_365  = o_efs_365
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
    labels  = c("0" = "Inpatient, n=645", "1" = "Outpatient, n=70"),
    colors  = c("Inpatient, n=645" = "purple4", "Outpatient, n=70" = "salmon"),
    data_90 = df_main_90, data_365 = df_main_365
  ),
  hospital = list(
    var = "hospital",
    labels  = c("BJH" = "WU, n=384", "OHSU" = "OHSU, n=331"),
    colors  = c("WU, n=384" = "#BA0C2F", "OHSU, n=331" = "#0E4D8F"),
    data_90 = df_main_90, data_365 = df_main_365
  ),
  sex = list(
    var = "Sex",
    labels  = c("Female" = "Female, n=233", "Male" = "Male, n=482"),
    colors  = c("Female, n=233" = "#cd0055", "Male, n=482" = "#2975f1"),
    data_90 = df_main_90, data_365 = df_main_365
  ),
  target = list(
    var = "target",
    labels = c("BCMA" = "BCMA, n=164", "CD-19" = "CD-19, n=551"),
    colors = c("BCMA, n=164" = "#2BAE66", "CD-19, n=551" = "#EEA47F"),
    data_90 = df_main_90, data_365 = df_main_365
  )
)

# Unadjusted logistic regression
fn_lr_unadj <- function(data, exposure_var, outcome_var) {
  tryCatch({
    # Check if exposure variable exists and has variation
    if (!exposure_var %in% names(data)) {
      stop("Exposure variable not found: ", exposure_var)
    }
    if (!outcome_var %in% names(data)) {
      stop("Outcome variable not found: ", outcome_var)
    }

    # Check for sufficient variation
    if (length(unique(na.omit(data[[exposure_var]]))) < 2) {
      warning("Insufficient variation in exposure: ", exposure_var)
      return(tibble(outcome = outcome_var, exposure = exposure_var,
                    or_unadj = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_))
    }

    formula_obj <- as.formula(paste(outcome_var, "~", exposure_var))
    model <- glm(formula_obj, data = data, family = binomial())

    results <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE, conf.method = "Wald")
    results |>
      filter(term == exposure_var) |>
      select(
        exposure  = term,
        or_unadj = estimate,
        ci_lo    = conf.low,
        ci_hi    = conf.high,
        p_value  = p.value
      ) |>
      mutate(outcome = outcome_var, .before = exposure)
  }, error = function(err) {
    warning("Error in unadjusted model for ", outcome_var, " ~ ", exposure_var, ": ", err$message)
    tibble(
      outcome  = outcome_var,
      exposure = exposure_var,
      or_unadj = NA_real_,
      ci_lo    = NA_real_,
      ci_hi    = NA_real_,
      p_value  = NA_real_
    )
  })
}

# Adjusted logistic regression - IMPROVED
fn_lr_adj <- function(data, exposure_var, covariate_vars, outcome_var, use_random = TRUE) {
  tryCatch({
    # Validate inputs
    missing_vars <- setdiff(c(exposure_var, outcome_var, covariate_vars), names(data))
    if (length(missing_vars) > 0) {
      stop("Missing variables: ", paste(missing_vars, collapse = ", "))
    }

    # Filter out covariates with no variation in this subset
    valid_covs <- covariate_vars[vapply(covariate_vars, function(v) {
      length(unique(na.omit(data[[v]]))) > 1
    }, logical(1))]

    if (length(valid_covs) == 0) {
      warning("No valid covariates with variation found")
      return(tibble(outcome = outcome_var, exposure = exposure_var,
                    or_adj = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_))
    }

    formula_string <- paste0(
      outcome_var, " ~ ", exposure_var, " + ",
      paste(valid_covs, collapse = " + ")
    )

    if (use_random && "hospital" %in% names(data) && length(unique(data$hospital)) > 1) {
      formula_string <- paste0(formula_string, " + (1 | hospital)")
      model <- glmer(as.formula(formula_string), data = data, family = binomial(),
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      results <- broom.mixed::tidy(model, exponentiate = TRUE, conf.int = TRUE, conf.method = "Wald")
    } else {
      model <- glm(as.formula(formula_string), data = data, family = binomial())
      results <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE, conf.method = "Wald")
    }

    results_sub <- results[results$term == exposure_var, ]
    if (nrow(results_sub) == 0) {
      warning("No coefficient found for ", exposure_var)
      return(tibble(outcome = outcome_var, exposure = exposure_var,
                    or_adj = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_))
    }

    results_sub |>
      mutate(outcome = outcome_var) |>
      select(outcome, exposure = term, or_adj = estimate, ci_lo = conf.low, ci_hi = conf.high, p_value = p.value)

  }, error = function(err) {
    warning("Error in adjusted model for ", outcome_var, " ~ ", exposure_var, ": ", err$message)
    tibble(outcome = outcome_var, exposure = exposure_var,
           or_adj = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_)
  })
}

# E-value computation - IMPROVED
compute_evalue_row <- function(or_est, ci_lo, ci_hi) {
  or_num <- as.numeric(or_est)
  lo_num <- as.numeric(ci_lo)
  hi_num <- as.numeric(ci_hi)

  if (any(is.na(c(or_num, lo_num, hi_num))) || lo_num >= hi_num || or_num <= 0) {
    return(c(e_value = NA_real_, e_value_lo = NA_real_, e_value_hi = NA_real_))
  }

  tryCatch({
    ev_mat <- EValue::evalues.OR(est = or_num, lo = lo_num, hi = hi_num, rare = FALSE)
    c(
      e_value    = ev_mat["E-values", "point"],
      e_value_lo = ev_mat["E-values", "lower"],
      e_value_hi = ev_mat["E-values", "upper"]
    )
  }, error = function(e) {
    c(e_value = NA_real_, e_value_lo = NA_real_, e_value_hi = NA_real_)
  })
}

# Core runner for one subgroup - IMPROVED
run_models_by_group <- function(data, group_var, exposures, outcomes, covariate_vars = NULL, model_fn, use_random_default = TRUE) {
  group_vals <- sort(unique(data[[group_var]]))
  func_name <- deparse(substitute(model_fn))

  map_dfr(group_vals, function(g) {
    df_sub <- data[data[[group_var]] == g, ]

    if (nrow(df_sub) < 10) {
      warning("Very small subgroup size (n=", nrow(df_sub), ") for ", group_var, "=", g)
    }

    combo <- expand.grid(e = exposures, o = outcomes, stringsAsFactors = FALSE)
    res <- pmap_dfr(combo, function(e, o) {
      if (func_name == "fn_lr_adj") {
        current_use_random <- ifelse(group_var == "hospital", FALSE, use_random_default)
        return(model_fn(
          data           = df_sub,
          exposure_var   = e,
          covariate_vars = covariate_vars,
          outcome_var    = o,
          use_random     = current_use_random
        ))
      } else {
        return(model_fn(data = df_sub, exposure_var = e, outcome_var = o))
      }
    })
    mutate(res, subgroup_var = group_var, subgroup_val = as.character(g))
  })
}

## Unified Subgroup Analysis Runner ---------------------------------------------
run_all_subgroup_analyses <- function(subgroup_configs, outcomes, data_key) {
  map_dfr(names(subgroup_configs), function(sg_name) {
    config <- subgroup_configs[[sg_name]]
    data <- config[[data_key]]
    outcomes_to_use <- if (sg_name == "outpatient" && data_key == "data_90") outcomes_90 else outcomes

    safely_run_models <- purrr::possibly(run_models_by_group, otherwise = tibble())

    unadj <- safely_run_models(
      data     = data,
      group_var = config$var,
      exposures = exposures,
      outcomes  = outcomes_to_use,
      model_fn  = fn_lr_unadj
    ) |> mutate(model = "Unadjusted", analysis = "Subgroup")

    adj_covs <- get_adjusted_covariates(config$var, covariates)
    adj <- safely_run_models(
      data           = data,
      group_var      = config$var,
      exposures      = exposures,
      outcomes       = outcomes_to_use,
      covariate_vars = adj_covs,
      model_fn       = fn_lr_adj,
      use_random_default = TRUE
    ) |> mutate(model = "Adjusted", analysis = "Subgroup")

    # Add E-values if function exists
    if (exists("compute_evalue_row", mode = "function")) {
      adj <- adj |>
        rowwise() |>
        mutate(evals = list(compute_evalue_row(or_adj, ci_lo, ci_hi))) |>
        unnest_wider(evals) |>
        ungroup()
    }

    bind_rows(unadj, adj) |> mutate(subgroup_type = sg_name)
  })
}

## Run Subgroup Analyses --------------------------------------------------------
message("Running 90-day analyses...")
all_results_90  <- run_all_subgroup_analyses(subgroup_configs, core_outcomes_90,  "data_90")

message("Running 365-day analyses...")
all_results_365 <- run_all_subgroup_analyses(subgroup_configs, core_outcomes_365, "data_365")

## Combine All Results ----------------------------------------------------------
all_results <-
  bind_rows(all_results_90, all_results_365) |>
  mutate(or = coalesce(or_adj, or_unadj)) |>
  select(analysis, subgroup_type, subgroup_var, subgroup_val, model,
         exposure, outcome, or, ci_lo, ci_hi, e_value, e_value_lo)

# Create output directory if it doesn't exist
if (!dir.exists(here("output"))) {
  dir.create(here("output"), recursive = TRUE)
}

data.table::fwrite(all_results, here("output", paste0("unified_subgroup_results_", today, ".csv")))
message("Results saved to: ", here("output", paste0("unified_subgroup_results_", today, ".csv")))

## Plotting Functions ------------------------------------------------------------

# (1) KM Curves by subgroup using 1500 cutoff - IMPROVED
plot_km_by_subgroup <- function(data, time_var, event_var, subgroup_var, subgroup_val, cutoff = 1500) {
  tryCatch({
    # Filter to specific subgroup first
    data_sub <- data |>
      filter(!!sym(subgroup_var) == subgroup_val) |>
      filter(!is.na(!!sym(time_var)), !is.na(!!sym(event_var)), !is.na(e_hours_precise)) |>
      mutate(time_group = if_else(e_hours_precise >= cutoff, "Late", "Early"))

    # Skip if only one time group present or insufficient data
    time_group_counts <- table(data_sub$time_group)
    if (length(time_group_counts) < 2 || min(time_group_counts) < 5 || nrow(data_sub) < 10) {
      warning("Insufficient data for KM plot in subgroup: ", subgroup_var, " = ", subgroup_val,
              " (n=", nrow(data_sub), ", time groups: ", paste(names(time_group_counts), "=", time_group_counts, collapse = ", "), ")")
      return(NULL)
    }

    surv_obj <- Surv(time = data_sub[[time_var]], event = data_sub[[event_var]])
    fit <- survfit(surv_obj ~ time_group, data = data_sub)

    ggsurvplot(
      fit,
      data = data_sub,
      pval = TRUE,
      conf.int = TRUE,
      legend.title = "Infusion Time",
      legend.labs = c("Early", "Late"),
      palette = c("#2975f1", "#cd0055"),
      ggtheme = theme_minimal(),
      title = paste("Kaplan-Meier Curves:", subgroup_var, "=", subgroup_val)
    )
  }, error = function(e) {
    warning("Error creating KM plot for ", subgroup_var, " = ", subgroup_val, ": ", e$message)
    return(NULL)
  })
}

# (2) AOR dot-whisker plots - IMPROVED
plot_dotwhisker <- function(results, subgroup_name, group_labels, colors) {
  # Better data validation and filtering
  plot_data <- results |>
    filter(subgroup_type == subgroup_name, model == "Adjusted") |>
    # Create the or column if it doesn't exist
    mutate(or = coalesce(or_adj, or_unadj)) |>
    # Filter out rows with missing odds ratios
    filter(!is.na(or), !is.na(ci_lo), !is.na(ci_hi)) |>
    # Add outcome labels
    mutate(
      outcome_label = case_when(
        outcome %in% names(outcome_labels) ~ outcome_labels[outcome],
        TRUE ~ outcome  # fallback to original outcome name
      )
    ) |>
    # Add subgroup labels
    mutate(
      subgroup_label = case_when(
        subgroup_val %in% names(group_labels) ~ group_labels[subgroup_val],
        TRUE ~ subgroup_val
      )
    ) |>
    # Filter out rows with missing labels
    filter(!is.na(outcome_label), !is.na(subgroup_label))

  if (nrow(plot_data) == 0) {
    warning("No valid data to plot for subgroup: ", subgroup_name)
    return(ggplot() +
             theme_void() +
             ggtitle(paste("No data available for", str_to_title(subgroup_name))) +
             theme(plot.title = element_text(hjust = 0.5)))
  }

  ggplot(plot_data, aes(x = or, y = fct_rev(outcome_label),
                        xmin = ci_lo, xmax = ci_hi, color = subgroup_label, group = subgroup_label)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(position = position_dodge(width = 0.5), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7) +
    scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
                  labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10")) +
    scale_color_manual(values = colors, name = str_to_title(subgroup_name)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Adjusted Odds Ratio (log scale)",
      y = NULL,
      title = paste("Adjusted Odds Ratios by", str_to_title(subgroup_name))
    )
}

# (3) FIXED Marginal effects plots with better error handling
plot_marginal_effects <- function(data, subgroup_var, subgroup_val, covs, outcome_var = "o_os_90") {
  tryCatch({
    df_sub <- data |> filter(!!sym(subgroup_var) == subgroup_val)

    if (nrow(df_sub) < 10) {  # Lower threshold
      warning("Insufficient data (n=", nrow(df_sub), ") for marginal effects in ", subgroup_val)
      return(ggplot() +
               theme_void() +
               ggtitle(paste("Insufficient data for", subgroup_val, "(n =", nrow(df_sub), ")")) +
               theme(plot.title = element_text(hjust = 0.5)))
    }

    # Convert to data frame
    df_sub <- as.data.frame(df_sub)

    # Check required variables exist
    required_vars <- c("e_hours_sunrise", outcome_var)
    missing_vars <- setdiff(required_vars, names(df_sub))
    if (length(missing_vars) > 0) {
      warning("Missing required variables: ", paste(missing_vars, collapse = ", "))
      return(ggplot() + theme_void())
    }

    # For very small subgroups or when covariates cause issues, use unadjusted model
    if (nrow(df_sub) < 30 || subgroup_val == "1") {  # Outpatient or small groups
      formula_str <- paste(outcome_var, "~ e_hours_sunrise")
      model_type <- "Unadjusted"
    } else {
      # Try to find covariates with variation
      valid_covs <- c()
      for (v in covs) {
        if (v %in% names(df_sub)) {
          x <- df_sub[[v]]
          if (!is.null(x) && sum(!is.na(x)) >= 5) {
            if (is.factor(x)) {
              # For factors, ensure there are at least 2 levels with observations
              level_counts <- table(x, useNA = "no")
              if (length(level_counts) >= 2 && min(level_counts) >= 1) {
                valid_covs <- c(valid_covs, v)
              }
            } else {
              # For numeric, ensure variation
              if (length(unique(x[!is.na(x)])) > 1) {
                valid_covs <- c(valid_covs, v)
              }
            }
          }
        }
      }

      if (length(valid_covs) == 0) {
        formula_str <- paste(outcome_var, "~ e_hours_sunrise")
        model_type <- "Unadjusted (no valid covariates)"
      } else {
        formula_str <- paste(outcome_var, "~ e_hours_sunrise +", paste(valid_covs, collapse = " + "))
        model_type <- "Adjusted"
      }
    }

    # Fit model
    model <- glm(as.formula(formula_str), data = df_sub, family = binomial())

    # Extract p-value
    coef_summary <- summary(model)$coefficients
    if (!"e_hours_sunrise" %in% rownames(coef_summary)) {
      warning("e_hours_sunrise coefficient not found in model")
      return(ggplot() + theme_void())
    }

    pval <- coef_summary["e_hours_sunrise", "Pr(>|z|)"]

    # Create prediction data
    sunrise_range <- range(df_sub$e_hours_sunrise, na.rm = TRUE)
    if (sunrise_range[1] == sunrise_range[2]) {
      warning("No variation in e_hours_sunrise")
      return(ggplot() + theme_void())
    }

    new_data <- data.frame(
      e_hours_sunrise = seq(sunrise_range[1], sunrise_range[2], length.out = 100)
    )

    # Set covariates to reference values if they exist in the model
    if (length(valid_covs) > 0) {
      for (v in valid_covs) {
        x <- df_sub[[v]]
        if (is.factor(x)) {
          # Use most common level
          most_common <- names(sort(table(x, useNA = "no"), decreasing = TRUE))[1]
          new_data[[v]] <- factor(rep(most_common, 100), levels = levels(x))
        } else if (is.character(x)) {
          most_common <- names(sort(table(x, useNA = "no"), decreasing = TRUE))[1]
          new_data[[v]] <- rep(most_common, 100)
        } else {
          new_data[[v]] <- rep(median(x, na.rm = TRUE), 100)
        }
      }
    }

    # Predict with confidence intervals
    pred_link <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
    new_data$pred <- plogis(pred_link$fit)
    new_data$pred_lwr <- plogis(pred_link$fit - 1.96 * pred_link$se.fit)
    new_data$pred_upr <- plogis(pred_link$fit + 1.96 * pred_link$se.fit)

    # Create plot
    ggplot(new_data, aes(x = e_hours_sunrise, y = pred)) +
      geom_ribbon(aes(ymin = pred_lwr, ymax = pred_upr), alpha = 0.3, fill = "blue") +
      geom_line(size = 1.2, color = "blue") +
      labs(
        title = paste("Marginal Effect:", subgroup_val, paste0("(", model_type, ")")),
        subtitle = paste("p-value:", format.pval(pval, digits = 3), "| n =", nrow(df_sub)),
        x = "Hours Since Sunrise",
        y = paste("Predicted Probability of", outcome_labels[outcome_var] %||% outcome_var)
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 11))

  }, error = function(e) {
    warning("Error in marginal effects plot for ", subgroup_val, ": ", e$message)
    ggplot() +
      theme_void() +
      ggtitle(paste("Error plotting", subgroup_val)) +
      labs(subtitle = e$message) +
      theme(plot.title = element_text(size = 10, hjust = 0.5),
            plot.subtitle = element_text(size = 8, hjust = 0.5))
  })
}

# (4) Unadjusted marginal effects for outpatient comparison
plot_outpatient_marginal_effects <- function(data, outcome_var = "o_os_90") {
  tryCatch({
    # Check if we have both inpatient and outpatient data
    outpt_counts <- table(data$outpt_01)
    if (length(outpt_counts) < 2) {
      warning("Only one outpatient status present in data")
      return(ggplot() + theme_void())
    }

    # Check for sufficient data in each group
    if (any(outpt_counts < 10)) {
      warning("Small sample size in outpatient groups: ", paste(names(outpt_counts), "=", outpt_counts, collapse = ", "))
    }

    # Create separate datasets
    inpt_data <- data |> filter(outpt_01 == 0) |> as.data.frame()
    outpt_data <- data |> filter(outpt_01 == 1) |> as.data.frame()

    plots_list <- list()

    # Inpatient model (unadjusted)
    if (nrow(inpt_data) >= 10) {
      inpt_model <- glm(as.formula(paste(outcome_var, "~ e_hours_sunrise")),
                        data = inpt_data, family = binomial())
      inpt_pval <- summary(inpt_model)$coefficients["e_hours_sunrise", "Pr(>|z|)"]

      sunrise_range_inpt <- range(inpt_data$e_hours_sunrise, na.rm = TRUE)
      new_data_inpt <- data.frame(
        e_hours_sunrise = seq(sunrise_range_inpt[1], sunrise_range_inpt[2], length.out = 100)
      )

      pred_link_inpt <- predict(inpt_model, newdata = new_data_inpt, type = "link", se.fit = TRUE)
      new_data_inpt$pred <- plogis(pred_link_inpt$fit)
      new_data_inpt$pred_lwr <- plogis(pred_link_inpt$fit - 1.96 * pred_link_inpt$se.fit)
      new_data_inpt$pred_upr <- plogis(pred_link_inpt$fit + 1.96 * pred_link_inpt$se.fit)
      new_data_inpt$group <- "Inpatient"

      plots_list[["inpatient"]] <- new_data_inpt
    }

    # Outpatient model (unadjusted)
    if (nrow(outpt_data) >= 5) {  # Lower threshold for outpatient due to small n
      outpt_model <- glm(as.formula(paste(outcome_var, "~ e_hours_sunrise")),
                         data = outpt_data, family = binomial())
      outpt_pval <- summary(outpt_model)$coefficients["e_hours_sunrise", "Pr(>|z|)"]

      sunrise_range_outpt <- range(outpt_data$e_hours_sunrise, na.rm = TRUE)
      new_data_outpt <- data.frame(
        e_hours_sunrise = seq(sunrise_range_outpt[1], sunrise_range_outpt[2], length.out = 100)
      )

      pred_link_outpt <- predict(outpt_model, newdata = new_data_outpt, type = "link", se.fit = TRUE)
      new_data_outpt$pred <- plogis(pred_link_outpt$fit)
      new_data_outpt$pred_lwr <- plogis(pred_link_outpt$fit - 1.96 * pred_link_outpt$se.fit)
      new_data_outpt$pred_upr <- plogis(pred_link_outpt$fit + 1.96 * pred_link_outpt$se.fit)
      new_data_outpt$group <- "Outpatient"

      plots_list[["outpatient"]] <- new_data_outpt
    }

    if (length(plots_list) == 0) {
      return(ggplot() + theme_void() + ggtitle("Insufficient data for outpatient analysis"))
    }

    # Combine data for plotting
    combined_data <- bind_rows(plots_list)

    # Create combined plot
    p <- ggplot(combined_data, aes(x = e_hours_sunrise, y = pred, color = group, fill = group)) +
      geom_ribbon(aes(ymin = pred_lwr, ymax = pred_upr), alpha = 0.2, color = NA) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("Inpatient" = "purple4", "Outpatient" = "salmon")) +
      scale_fill_manual(values = c("Inpatient" = "purple4", "Outpatient" = "salmon")) +
      labs(
        title = "Unadjusted Marginal Effects: Inpatient vs Outpatient",
        subtitle = paste("Outcome:", outcome_labels[outcome_var]),
        x = "Hours Since Sunrise",
        y = "Predicted Probability",
        color = "Setting",
        fill = "Setting"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10)
      )

    # Add p-values as text if both models exist
    if (exists("inpt_pval") && exists("outpt_pval")) {
      p <- p + annotate("text", x = Inf, y = Inf,
                        label = paste("Inpatient p =", format.pval(inpt_pval, digits = 3),
                                      "\nOutpatient p =", format.pval(outpt_pval, digits = 3)),
                        hjust = 1.1, vjust = 1.1, size = 3)
    } else if (exists("inpt_pval")) {
      p <- p + annotate("text", x = Inf, y = Inf,
                        label = paste("Inpatient p =", format.pval(inpt_pval, digits = 3)),
                        hjust = 1.1, vjust = 1.1, size = 3)
    } else if (exists("outpt_pval")) {
      p <- p + annotate("text", x = Inf, y = Inf,
                        label = paste("Outpatient p =", format.pval(outpt_pval, digits = 3)),
                        hjust = 1.1, vjust = 1.1, size = 3)
    }

    return(p)

  }, error = function(e) {
    warning("Error in outpatient marginal effects plot: ", e$message)
    ggplot() +
      theme_void() +
      ggtitle(paste("Error in outpatient analysis:", e$message)) +
      theme(plot.title = element_text(size = 10))
  })
}

## Generate Plots ---------------------------------------------------------------
# Create figures directory if it doesn't exist
if (!dir.exists(here("figures"))) {
  dir.create(here("figures"), recursive = TRUE)
}

message("Generating plots...")

for (sg_name in names(subgroup_configs)) {
  message("Processing subgroup: ", sg_name)
  cfg <- subgroup_configs[[sg_name]]

  # KM plots
  tryCatch({
    km_plot <- plot_km_by_subgroup(cfg$data_90, time_var = "o_os_90", event_var = "surv_90", subgroup_var = cfg$var)
    if (!is.null(km_plot)) {
      ggsave(here("figures", paste0("kmcurve_", sg_name, "_", today, ".pdf")),
             km_plot$plot, width = 10, height = 6, dpi = 300)
    }
  }, error = function(e) {
    warning("Failed to create KM plot for ", sg_name, ": ", e$message)
  })

  # Dot-whisker plots
  tryCatch({
    dw_plot <- plot_dotwhisker(all_results_90, sg_name, cfg$labels, cfg$colors)
    ggsave(here("figures", paste0("aor_dotwhisker_", sg_name, "_", today, ".pdf")),
           dw_plot, width = 10, height = 8, dpi = 300)
  }, error = function(e) {
    warning("Failed to create dot-whisker plot for ", sg_name, ": ", e$message)
  })

  # Marginal effects plots
  if (sg_name == "outpatient") {
    # Special handling for outpatient - create unadjusted comparison plot
    tryCatch({
      outpt_me_plot <- plot_outpatient_marginal_effects(cfg$data_90, outcome_var = "o_os_90")
      ggsave(here("figures", paste0("meplot_outpatient_comparison_", today, ".pdf")),
             outpt_me_plot, width = 10, height = 6, dpi = 300)
    }, error = function(e) {
      warning("Failed to create outpatient comparison plot: ", e$message)
    })
  } else {
    # Regular marginal effects plots for other subgroups
    tryCatch({
      for (sg_val in unique(cfg$data_90[[cfg$var]])) {
        me_plot <- plot_marginal_effects(cfg$data_90, cfg$var, sg_val, covariates)
        ggsave(here("figures", paste0("meplot_", sg_name, "_", sg_val, "_", today, ".pdf")),
               me_plot, width = 8, height = 6, dpi = 300)
      }
    }, error = function(e) {
      warning("Failed to create marginal effects plots for ", sg_name, ": ", e$message)
    })
  }
}

message("All subgroup analysis outputs completed successfully.")

# Print summary
message("\nSummary:")
message("- Results saved to: ", here("output", paste0("unified_subgroup_results_", today, ".csv")))
message("- Figures saved to: ", here("figures"))
message("- Total models run: ", nrow(all_results))
message("- Successful models: ", sum(!is.na(all_results$or)))
