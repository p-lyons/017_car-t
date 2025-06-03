

## libraries
library(survminer)
library(tidytable)
library(collapse)
library(survival)
library(ggplot2)
library(stringr)
library(cmprsk)
library(EValue)
library(coxme)
library(arrow)
library(here)

## helpers
today = format(Sys.Date(), "%y%m%d")

# data -------------------------------------------------------------------------

## load most recent analysis file ----------------------------------------------

files  = list.files(here("clean"), pattern = "^analysis_\\d{6}\\.csv$", full.names = TRUE)
dates  = str_extract(files, "\\d{6}")
latest = files[which.max(dates)]
df     = fread(latest)
message("Loaded file: ", latest)

censor_date = as.Date("2025-04-01", format = "%Y-%m-%d")

## clean up event times --------------------------------------------------------

df_time =
  select(df, mrn, admin_date, event_cat, event_dttm, death_dttm, o_d90_01, o_d365_01) |>
  fmutate(
    admin_date = as.Date(admin_date, format = "%Y-%m-%d"),
    death_date = as.Date(death_dttm, format = "%Y-%m-%d"),
    event_date = as.Date(event_dttm, format = "%Y-%m-%d"),
    event_date = if_else(event_date > censor_date, censor_date, event_date),
    death_date = if_else(event_cat == "death" & is.na(death_date), event_date, death_date)
  )

missing_death_date =
  fsubset(df_time, o_d365_01 == 1 & is.na(death_date)) |>
  select(mrn, admin_date, o_d90_01, o_d365_01, event_dttm, death_date)

mrns = funique(missing_death_date$mrn)

dem =
  arrow::read_parquet(here("../008_clif/rclif/CLIF-2.0.0-dttm-update-R/clif_patient_2016_2024.parquet")) |>
  fsubset(patient_id %in% mrns) |>
  # Some people have deaths recorded in demo that weren't in CAR-T and vice-versa
  # mutate(death_date_2 = as.Date(death_dttm)) |>
  select(mrn  = patient_id, death_dttm) |>
  funique()

dem2 =
  fread(here("pre/patients_250402.csv")) |>
  fsubset(PrimaryMrn %in% mrns) |>
  select(
    mrn          = PrimaryMrn,
    death_date_3 = DeathDate
  ) |>
  funique()

dem =
  join(dem, dem2, how = "full", multiple = F) |>
  fmutate(death_date_2 = as.Date(death_dttm,    format = "%Y-%m-%d")) |>
  fmutate(death_date_3 = as.Date(death_date_3, format = "%Y-%m-%d")) |>
  fmutate(death_date_2 = pmin(death_date_2, death_date_3, na.rm = T)) |>
  select(mrn, death_date_2)

df_time =
  join(df_time, dem, how = "left", multiple = F) |>
  fmutate(death_date = if_else(is.na(death_date), death_date_2, death_date)) |>
  select(mrn, admin_date, event_date, death_date, event_cat, o_d90_01, o_d365_01)

rm(dem, dem2, missing_death_date)
gc()

df_time = df_time |>
  fmutate(
    overall_time   = as.numeric(if_else(is.na(death_date), censor_date - admin_date, death_date - admin_date)),
    overall_status = if_else(event_cat == "death" | !is.na(death_date), 1L, 0L)  # 1 = event (death), 0 = censored
  ) |>
  fmutate(
    efs_event_date = if_else(is.na(event_date), censor_date, event_date),
    efs_time       = as.numeric(efs_event_date - admin_date),
    efs_status     = if_else(event_cat == "none", 0L, 1L)  # 1 = event, 0 = censored
  )

## exposure for KM figure ------------------------------------------------------

df =
  select(df, mrn, e_late_f15_01, hospital) |>
  join(df_time, how = "left", multiple = F)

write_parquet(df, here("clean", paste0("df_survival_", today, ".parquet")))

# KM figures -------------------------------------------------------------------

### Create survival objects for KM analysis
os_surv  = Surv(time = df$overall_time, event = df$overall_status)
efs_surv = Surv(time = df$efs_time,     event = df$efs_status)
os_fit   = survfit(os_surv  ~ e_late_f15_01, data = df)
efs_fit  = survfit(efs_surv ~ e_late_f15_01, data = df)

km_os =
  ggsurvplot(
    os_fit,
    data              = df,
    size              = 1.8,
    risk.table        = F,
    pval              = T,
    conf.int          = F,
    palette           = c("#F1A226", "#800072"),
    #surv.median.line  = "hv",
    xlim              = c(0,730),
    xlab              = "Days Since CAR T-cell Administration",
    ylab              = "Overall Survival Probability",
    break.time.by     = 90,
    legend.labs       = c("CAR T-cells Before 15:00", "CAR T-cells After 15:00"),
   # risk.table.y.text = F,
    ggtheme           = theme_bw(base_size = 14),
    title             = "A. Overall Survival"
  )

km_os
ggsave(here("figures", paste0("km_os_15_", today, ".pdf")))

km_efs =
  ggsurvplot(
    efs_fit,
    data              = df,
    size              = 1.8,
    risk.table        = F,
    pval              = T,
    conf.int          = F,
    palette           = c("#F1A226", "#800072"),
    #surv.median.line  = "hv",
    xlim              = c(0,730),
    xlab              = "Days Since CAR T-cell Administration",
    ylab              = "Event-Free Survival Probability",
    break.time.by     = 90,
    legend.labs       = c("CAR T-cells Before 15:00", "CAR T-cells After 15:00"),
    # risk.table.y.text = F,
    ggtheme           = theme_bw(base_size = 14),
    title             = "B. Event-Free Survival"
  )

ggsave(here("figures", paste0("km_efs_15_", today, ".pdf")))

library(patchwork)
p1 = km_os$plot
p2 = km_efs$plot

p1 =
  p1 +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank()
  )

p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(
  here("figures", paste0("figure_03_survival_stacked_", today, ".pdf")),
  width = 11,
  height = 8,
  units = "in",
  dpi = 600
)

# 3rd color = "#298C8C"

# RMST -------------------------------------------------------------------------

library(survRM2)

# Subset the data to the groups of interest: Earliest and Latest Tertiles.
# Create a binary variable ("group") that is 0 for Earliest and 1 for Latest.
# (The rmst2 function requires a binary arm indicator.)
df_os = df |> fmutate(group = if_else(e_late_f15_01 == 0, 0L, 1L))

# Set tau for OS analysis (e.g., 730 days as per your xlim)
tau_os = 730

# Calculate RMST for Overall Survival using overall_time and overall_status.
rmst_os = rmst2(time = df_os$overall_time,
                 status = df_os$overall_status,
                 arm = df_os$group,
                 tau = tau_os)

# Print RMST results for overall survival.
print(rmst_os)

# cox models -------------------------------------------------------------------

files  = list.files(here("clean"), pattern = "^analysis_\\d{6}\\.csv$", full.names = TRUE)
dates  = str_extract(files, "\\d{6}")
latest = files[which.max(dates)]
df2    = fread(latest)
message("Loaded file: ", latest)

censor_date = as.Date("2025-04-01", format = "%Y-%m-%d")

exposures = c(
  "e_hours_precise",
  "e_hours_sunrise"
)

covs = c(
  "age",
  "female_01",
  "nhw_01",
  "ecog",
  "ecog_23",
  "vw",
  "saps2",
  "saps2_noage",
  "cancer_cat",
  "product_cat",
  "cell_dose_std",
  "conditioning_flucy_01",
  "tki_01",
  "outpt_01",
  "los_pre_cart",
  "year_cat",
  "season"
)

df =
  select(df2, mrn, all_of(exposures), all_of(covs)) |>
  join(df, how = "left", multiple = F) |>
  fmutate(cancer_cat = if_else(cancer_cat  %in% c("MCL", "Follicular Lymphoma"), "Other Lymphoma", cancer_cat)) |>
  fmutate(cancer_cat = qF(cancer_cat))

p95             = fquantile(df$los_pre_cart, 0.95)
df$los_pre_cart = if_else(df$los_pre_cart > p95, p95, df$los_pre_cart)
df$log_pre_cart = log(df$los_pre_cart + 1)


# adjusted RMST ----------------------------------------------------------------

library(broom.mixed)
library(pseudo)
library(lme4)

# 1. compute pseudo‐values at tau = 365
df$phi <- pseudomean(df$overall_time, df$overall_status, tmax = 730)

# 2. linear regression on the pseudo‐values
formula_precise = as.formula(paste("phi ~ e_hours_precise + ", paste(covs, collapse = " + "), "+ (1|hospital)"))
formula_sunrise = as.formula(paste("phi ~ e_hours_sunrise + ", paste(covs, collapse = " + "), "+ (1|hospital)"))
fit_precise     = lme4::lmer(formula_precise, data = df)
fit_sunrise     = lme4::lmer(formula_sunrise, data = df)

# 3. bootstrap RMST coefficient (can't get robust SE with only 2 clusters)
boot_precise =
  bootMer(
    fit_precise,
    FUN   = function(x) fixef(x)["e_hours_precise"],
    nsim  = 1000,
    use.u = FALSE,
    type  = "parametric",
    seed  = 2025
  )

boot_sunrise =
  bootMer(
    fit_sunrise,
    FUN   = function(x) fixef(x)["e_hours_sunrise"],
    nsim  = 1000,
    use.u = FALSE,
    type  = "parametric",
    seed  = 2025
  )

ci_precise  = quantile(boot_precise$t, probs = c(0.025, 0.975), na.rm = TRUE)
ci_sunrise  = quantile(boot_sunrise$t, probs = c(0.025, 0.975), na.rm = TRUE)
est_precise = fixef(fit_precise)["e_hours_precise"]
est_sunrise = fixef(fit_sunrise)["e_hours_sunrise"]

coefs_rmst_os =
  tidytable(
    term       = c("e_hours_precise", "e_hours_sunrise"),
    rmst_per_h = c(est_precise, est_sunrise),
    ci_lo      = c(ci_precise[1], ci_sunrise[1]),
    ci_hi      = c(ci_precise[2], ci_sunrise[2])
  )

coefs_rmst_os

### check binary exposure too
# formula_binary = as.formula(paste("phi ~ e_late_f15_01 + ", paste(covs, collapse = " + ")))
# fit_binary     = lm(formula_binary, data = df)
# coef_binary    = broom::tidy(fit_binary, conf.int = T, vcov = vcovHC(fit_binary, type = "HC1"))
# fsubset(coef_binary, term == "e_late_f15_01")

## repeat for EFS

df$phi <- pseudomean(df$efs_time, df$efs_status, tmax = 730)

formula_precise = as.formula(paste("phi ~ e_hours_precise + ", paste(covs, collapse = " + "), "+ (1|hospital)"))
formula_sunrise = as.formula(paste("phi ~ e_hours_sunrise + ", paste(covs, collapse = " + "), "+ (1|hospital)"))
fit_precise     = lme4::lmer(formula_precise, data = df)
fit_sunrise     = lme4::lmer(formula_sunrise, data = df)

# 3. bootstrap RMST coefficient (can't get robust SE with only 2 clusters)
boot_precise =
  bootMer(
    fit_precise,
    FUN   = function(x) fixef(x)["e_hours_precise"],
    nsim  = 1000,
    use.u = FALSE,
    type  = "parametric",
    seed  = 2025
  )

boot_sunrise =
  bootMer(
    fit_sunrise,
    FUN   = function(x) fixef(x)["e_hours_sunrise"],
    nsim  = 1000,
    use.u = FALSE,
    type  = "parametric",
    seed  = 2025
  )

ci_precise  = quantile(boot_precise$t, probs = c(0.025, 0.975), na.rm = TRUE)
ci_sunrise  = quantile(boot_sunrise$t, probs = c(0.025, 0.975), na.rm = TRUE)
est_precise = fixef(fit_precise)["e_hours_precise"]
est_sunrise = fixef(fit_sunrise)["e_hours_sunrise"]

coefs_rmst_efs =
  tidytable(
    term       = c("e_hours_precise", "e_hours_sunrise"),
    rmst_per_h = c(est_precise, est_sunrise),
    ci_lo      = c(ci_precise[1], ci_sunrise[1]),
    ci_hi      = c(ci_precise[2], ci_sunrise[2])
  )

# eTable 3 ---------------------------------------------------------------------
coefs_rmst_os
coefs_rmst_efs

