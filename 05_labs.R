


# labs of interest


## libraries
library(tidytable)
library(lubridate)
library(collapse)
library(ggplot2)
library(stringr)
library(arrow)
library(here)

## helpers
today = format(Sys.Date(), "%y%m%d")

# data -------------------------------------------------------------------------

## load most recent analysis cohort file ---------------------------------------

files  = list.files(here("clean"), pattern = "^analysis_\\d{6}\\.csv$", full.names = TRUE)
dates  = str_extract(files, "\\d{6}")
latest = files[which.max(dates)]
df     = fread(latest)
df     = select(df, mrn, hospital, product_cat, starts_with("e_"), contains("crs"), contains("icans"), o_d90_01, o_e90_01)
rm(files, dates, latest)
gc()

## bjh max labs ----------------------------------------------------------------

bjh_max =
  fread(here("pre/analysis_250316.csv")) |>
  janitor::clean_names() |>
  fsubset(hospital == "BJH") |>
  fmutate(mrn = str_pad(mrn, width = 8, side = "left", pad = "0")) |>
  select(mrn, starts_with("max")) |>
  janitor::remove_empty()

bjh_max =
  join(bjh_max, df, how = "inner", multiple = F) |>
  fmutate(surv_90 = if_else(o_d90_01 == 1, 0L, 1L)) |>
  fmutate(efs_90  = if_else(o_e90_01 == 1, 0L, 1L)) |>
  mutate(across(
    .cols  = starts_with("max"),
    .fns   = ~if_else(.x == 0, NA, .x)
  )) |>
  mutate(across(
    .cols  = starts_with("max"),
    .fns   = ~log(.x + 1),
    .names = "log_{.col}"
  )) |>
  select(
    mrn,
    e_hours_precise,
    e_late_f15_01,
    contains("max_"),
    ends_with("90"),
    o_crs_01,
    o_crs34_01,
    o_icans_01,
    o_icans34_01
  ) |>
  pivot_longer(
    cols      = contains("max"),
    names_to  = "lab_category",
    values_to = "lab_value"
  )

## ohsu max labs ---------------------------------------------------------------

files  = list.files(here("clean"), pattern = "^labs_inflam_\\d{6}\\.csv$", full.names = TRUE)
dates  = str_extract(files, "\\d{6}")
latest = files[which.max(dates)]

o_max =
  fread(latest) |>
  fmutate(day = as.numeric(difftime(lab_collect_dttm, admin_dttm), "hours")/24) |>
  fmutate(lab_category = if_else(lab_category == "il-2", "il_02", lab_category)) |>
  fmutate(lab_category = if_else(lab_category == "il-6", "il_06", lab_category)) |>
  fsubset(day >= 0 & day <= 14) |>
  fgroup_by(mrn, lab_category) |>
  fsummarize(lab_value = fmax(lab_value)) |>
  fmutate(log_value = log(lab_value + 1))

rm(files, dates, latest)
gc()

logs =
  select(o_max, mrn, lab_category, lab_value = log_value) |>
  fmutate(lab_category = paste0("log_max_", lab_category))

o_max =
  select(o_max, -log_value) |>
  fmutate(lab_category = paste0("max_", lab_category)) |>
  rowbind(logs)

o_max =
  join(o_max, df, how = "left", multiple = T) |>
  fmutate(surv_90 = if_else(o_d90_01 == 1, 0L, 1L)) |>
  fmutate(efs_90  = if_else(o_e90_01 == 1, 0L, 1L)) |>
  select(
    mrn,
    e_hours_precise,
    e_late_f15_01,
    starts_with("lab_"),
    ends_with("90"),
    o_crs_01,
    o_crs34_01,
    o_icans_01,
    o_icans34_01
  )

labs = rowbind(o_max, bjh_max)

ggplot(labs, aes(x = lab_value, fill = lab_category)) +
  geom_density(alpha = 0.6) +
  theme_bw() +
  facet_wrap(~lab_category, scales = "free")

## log scale = ferritin, ifn_gamma, all il-, tnf, ldh
## reg scale = crp

keepers = c(
  "crp",
  "log_max_ferritin",
  "log_max_ldh",
  "log_max_tnf",
  "log_max_ifn_gamma",
  #"log_max_il_01",
  "log_max_il_02",
  #"log_max_il_04",
  "log_max_il_05",
  "log_max_il_06",
  "log_max_il_08",
  "log_max_il_10",
  #"log_max_il_12",
  "log_max_il_13"
  #"log_max_il_17"
)

outcome_labels = c("CRS", "CRS 3/4", "ICANS", "ICANS 3/4")

lab_labels = c(
  "Ferritin",
  "IFN-g",
  "IL-02",
  "IL-05",
  "IL-06",
  "IL-08",
  "IL-10",
  "IL-13",
  "LDH",
  "TNF-a"
)

plot_tod =
  fsubset(labs, lab_category %in% keepers) |>
  pivot_longer(
    cols      = c(o_crs_01, o_crs34_01, o_icans_01, o_icans34_01),
    names_to  = "outcome",
    values_to = "event_01"
  ) |>
  fmutate(lab_cat  = factor(lab_category, labels = lab_labels)) |>
  fmutate(outcome  = factor(outcome,      labels = outcome_labels)) |>
  fmutate(event_01 = factor(event_01,     labels = c("Did Not Experience Outcome", "Experienced Outcome"))) |>
  ggplot(aes(x = e_hours_precise, y = lab_value, color = event_01, fill = event_01)) +
  geom_point(size = 0.4, alpha = 0.3) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(8, 18), breaks = seq(9, 18, 3)) +
  scale_color_manual(values = c("#5ec962", "#3b528b")) +
  scale_fill_manual(values = c("#5ec962", "#3b528b")) +
  labs(
    x     = "Hour of CAR T-cell Infusion",
    y     = "Maximal Value (Log Scale) Through Day 14",
    color = "",
    fill  = ""
  ) +
  facet_grid(outcome~lab_cat, scales = "free_y")

plot_tod

ggsave(
  here("figures", paste0("cytokines_outcomes_", today, ".pdf")),
  width = 11,
  height = 8,
  units = "in",
  dpi = 600
)

# inspect data -----------------------------------------------------------------

skimr::skim(df)

eda_labs_plot =
  select(df, starts_with("log_max_il"), e_hours_precise) |>
  GGally::ggpairs(
    upper   = list(continuous = GGally::wrap("cor", size = 3)),
    lower   = list(continuous = GGally::wrap("points", alpha = 0.6, size = 2)),
    diag    = list(continuous = GGally::wrap("densityDiag", alpha = 0.4))
  ) +
  theme_minimal()

eda_labs_plot

# plot -------------------------------------------------------------------------

## interleukins by time of day and product
select(df, e_hours_precise, product_cat, starts_with("log_max_il"), log_max_ifn_gamma) |>
  pivot_longer(cols = starts_with("log_max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = product_cat)) +
  geom_point() +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)

## interleukins by time of day and outcomes
select(df, e_hours_precise, o_d90_01, starts_with("log_max_il"), log_max_ifn_gamma) |>
  fmutate(o_d90_01 = qF(o_d90_01)) |>
  pivot_longer(cols = starts_with("log_max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = o_d90_01)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)

select(df, e_hours_precise, o_e90_01, starts_with("log_max_il"), log_max_ifn_gamma) |>
  fmutate(o_e90_01 = qF(o_e90_01)) |>
  pivot_longer(cols = starts_with("log_max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = o_e90_01)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)

select(df, e_hours_precise, o_icans34_01, starts_with("log_max_il"), log_max_ifn_gamma) |>
  fmutate(o_icans34_01 = qF(o_icans34_01)) |>
  pivot_longer(cols = starts_with("log_max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = o_icans34_01)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)

select(df, e_hours_precise, o_crs34_01, starts_with("log_max_il"), log_max_ifn_gamma) |>
  fmutate(o_crs34_01 = qF(o_crs34_01)) |>
  pivot_longer(cols = starts_with("log_max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = o_crs34_01)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)


# non interleukins

select(df, e_hours_precise, o_d90_01, max_crp, max_ldh, max_ferritin, max_tnf) |>
  fmutate(o_d90_01 = qF(o_d90_01)) |>
  pivot_longer(cols = starts_with("max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = o_d90_01)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)

select(df, e_hours_precise, o_e90_01, max_crp, max_ldh, max_ferritin, max_tnf) |>
  fmutate(o_e90_01 = qF(o_e90_01)) |>
  pivot_longer(cols = starts_with("max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = o_e90_01)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)

select(df, e_hours_precise, o_icans34_01, max_crp, max_ldh, max_ferritin, max_tnf) |>
  fmutate(o_icans34_01 = qF(o_icans34_01)) |>
  pivot_longer(cols = starts_with("max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = o_icans34_01)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)

select(df, e_hours_precise, o_crs34_01, max_crp, max_ldh, max_ferritin, max_tnf) |>
  fmutate(o_crs34_01 = qF(o_crs34_01)) |>
  pivot_longer(cols = starts_with("max"), names_to = "lab", values_to = "value") |>
  ggplot(aes(x = e_hours_precise, y = value, color = o_crs34_01)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~lab, scales = "free_y") +
  theme_minimal(base_size = 16)


## time from administration to outcomes

labs |>
  fmutate(e_hours_precise = round(e_hours_precise)) |>
  fmutate(e_hours_precise = qF(e_hours_precise)) |>
  mutate(across(
    .cols  = starts_with("i"),
    .fns   = ~log(.x + 1),
    .names = "log_{.col}"
  )) |>
  ggplot(aes(x = lab_time, y = lab_value, color = e_hours_precise)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~lab_category, scales = "free_y") +
  theme_minimal()

# plot alc serially ------------------------------------------------------------

### which tests?
test_names = c("abs_lymphs")

### empty table to fill
labs = tidytable()
mrns = funique(df$mrn)

### list lab files
files = list.files(
  here("../000_data/ohsu/intermediate"),
  pattern    = "labs_intermediate_",
  full.names = TRUE
)

### loop over all files to open
for (file in files) {
  table =
    open_dataset(file) |>
    dplyr::filter(mrn %in% mrns, test_name %in% test_names) |>
    dplyr::select(mrn, CollectionInstant, test_name, value = ValueLCRF) |>
    dplyr::collect()

  labs = rowbind(labs, table, fill = T)
}

rm(table)
gc()

## clean up labs ---------------------------------------------------------------

files  = list.files(here("clean"), pattern = "^analysis_\\d{6}\\.csv$", full.names = TRUE)
dates  = str_extract(files, "\\d{6}")
latest = files[which.max(dates)]
df     = fread(latest)
df     = select(df, mrn, hospital, admin_dttm, starts_with("e_"), contains("crs"), contains("icans"), o_d90_01, o_e90_01)
rm(files, dates, latest)
gc()

labs =
  join(labs, df, how = "inner", multiple = T) |>
  ftransform(test_name = if_else(test_name == "abs_neut", "abs_neutrophil", test_name)) |>
  ftransform(test_name = if_else(test_name == "crp_hisense", "crp",         test_name)) |>
  ftransform(value = str_remove(value, "<")) |>
  ftransform(value = str_remove(value, ">")) |>
  ftransform(value = str_remove(value, " H")) |>
  ftransform(lab_value = as.numeric(value))

labs =
  fsubset(labs, !is.na(lab_value)) |>
  ftransform(lab_collect_dttm = lubridate::ymd_hms(CollectionInstant, tz = "America/Los_Angeles")) |>
  select(mrn, admin_dttm, lab_collect_dttm, lab_category = test_name, lab_value) |>
  fsubset(lab_collect_dttm > admin_dttm - lubridate::ddays(3)) |>
  fgroup_by(mrn, admin_dttm, lab_collect_dttm, lab_category) |>
  fmax()

### 00622597 (2021-10-06) has no labs for some reason: chart review
### zero lab points
### remaining labs:
milabs =
  fread(here("pre", "missing_labs_patient.csv")) |>
  fmutate(mrn              = "00622597") |>
  fmutate(admin_dttm       = "2021-10-06 15:16:00") |>
  fmutate(admin_dttm       = lubridate::ymd_hms(admin_dttm)) |>
  fmutate(lab_collect_dttm = lubridate::mdy_hm(lab_collect_dttm)) |>
  fmutate(lab_value        = as.numeric(lab_value))

labs = rowbind(labs, milabs, fill = T)

### cell counts

# fill in missing to make sure there are values for all patients

cellcount =
  select(labs, -lab_category) |>
  fmutate(day = as.numeric(difftime(lab_collect_dttm, admin_dttm), "hours")/24) |>
  fmutate(day = floor(day)) |>
  fsubset(day %in% c(0, 7, 0, 14)) |>
  fgroup_by(mrn, day) |>
  fsummarize(alc = fmax(lab_value))

bjh =
  fread(here("pre/CircadianClocks_BJH_ALC Day 0,7,10,14.csv")) |>
  janitor::clean_names() |>
  fmutate(mrn = as.character(record_id)) |>
  fmutate(mrn = str_pad(mrn, width = 8, side = "left", pad = "0")) |>
  pivot_longer(
    cols      = starts_with("day"),
    names_to  = "day",
    values_to = "alc"
  ) |>
  fmutate(day = readr::parse_number(day)) |>
  select(mrn, day, alc)

cellcount = rowbind(cellcount, bjh)

df =
  select(df, mrn, e_hours_precise, o_d90_01, o_icans34_01) |>
  join(cellcount, how = "left", multiple = T) |>
  fmutate(alc = as.numeric(alc)) |>
  fmutate(
    tgroup = case_when(
      e_hours_precise < 11 ~ 1,
      e_hours_precise < 14 ~ 2,
      T                    ~ 3
    )) |>
  fmutate(tgroup = qF(tgroup))

# df = fsubset(df, alc < 10) # chart review the high ones mrn 00622597

ggplot(df, aes(x = factor(day), y = alc, fill = tgroup)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3)) +
  xlab("Day") +
  ylab("ALC") +
  labs(fill = "Time Group")

# try a continuous version... --------------------------------------------------

cellcount =
  select(labs, -lab_category) |>
  fmutate(day = as.numeric(difftime(lab_collect_dttm, admin_dttm), "hours")/24) |>
  fmutate(day = floor(day)) |>
  #fsubset(day %in% c(0, 7, 0, 14)) |>
  fgroup_by(mrn, day) |>
  fsummarize(alc = fmax(lab_value)) |>
  fsubset(day >= 0 & day <= 14)

bjh =
  fread(here("pre/CircadianClocks_BJH_ALC Day 0,7,10,14.csv")) |>
  janitor::clean_names() |>
  fmutate(mrn = as.character(record_id)) |>
  fmutate(mrn = str_pad(mrn, width = 8, side = "left", pad = "0")) |>
  pivot_longer(
    cols      = starts_with("day"),
    names_to  = "day",
    values_to = "alc"
  ) |>
  fmutate(day = readr::parse_number(day)) |>
  select(mrn, day, alc)

cellcount = rowbind(cellcount, bjh)

df =
  select(df, mrn, e_hours_precise, o_d90_01, o_icans34_01) |>
  join(cellcount, how = "left", multiple = T) |>
  fmutate(alc = as.numeric(alc)) |>
  fmutate(
    tgroup = case_when(
      e_hours_precise < 11 ~ 1,
      e_hours_precise < 14 ~ 2,
      T                    ~ 3
    )) |>
  fmutate(tgroup = qF(tgroup))

ggplot(df, aes(x = day, y = alc, color = tgroup, fill = tgroup)) +
  geom_smooth(method = "gam", se = T) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3)) +
  xlab("Day") +
  ylab("ALC")