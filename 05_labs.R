

# Examine laboratory values for CAR T-cell time of day project

# setup ------------------------------------------------------------------------

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
df     = select(df, mrn, hospital, product_cat, starts_with("e_"), ldh, contains("crs"), contains("icans"), o_d90_01, o_e90_01)
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

bjh_new =
  fread(here("pre/CircadianClocksInCAR_DATA_LABELS_2025-07-23_1155.csv")) |>
  janitor::clean_names() |>
  fsubset(record_id > 360)

crp_cols = c(paste0("value_day_", seq(0, 10, by = 1), "_", seq(119, 139, by = 2)))
ldh_cols = c(paste0("value_day_", seq(0, 10, by = 1), "_", seq(165, 185, by = 2)))
fer_cols = c(paste0("value_day_", seq(0, 10, by = 1), "_", seq(142, 162, by = 2)))

crp =
  select(bjh_new, mrn = record_id, all_of(crp_cols)) |>
  pivot_longer(-mrn, names_to = "var", values_to = "crp") |>
  fmutate(crp = as.numeric(crp)) |>
  fgroup_by(mrn) |>
  fsummarize(max_crp = fmax(crp))

ldh =
  select(bjh_new, mrn = record_id, all_of(ldh_cols)) |>
  pivot_longer(-mrn, names_to = "var", values_to = "ldh") |>
  fmutate(ldh = as.numeric(ldh)) |>
  fgroup_by(mrn) |>
  fsummarize(max_ldh = fmax(ldh))

fer =
  select(bjh_new, mrn = record_id, all_of(fer_cols)) |>
  pivot_longer(-mrn, names_to = "var", values_to = "ferritin") |>
  fmutate(ferritin = as.numeric(ferritin)) |>
  fgroup_by(mrn) |>
  fsummarize(max_ferritin = fmax(ferritin))

bjh_new =
  join(crp, ldh, how = "full", multiple = F) |>
  join(fer,      how = "full", multiple = F) |>
  mutate(across(
    .cols  = starts_with("max"),
    .fns   = ~if_else(.x == 0, NA, .x)
  )) |>
  mutate(across(
    .cols  = starts_with("max"),
    .fns   = ~log(.x + 1),
    .names = "log_{.col}"
  )) |>
  pivot_longer(-mrn, names_to = "lab_category", values_to = "lab_value") |>
  fmutate(mrn = as.character(mrn)) |>
  fmutate(mrn = str_pad(mrn, 8, "left", "0"))

bjh_new =
  join(bjh_new, df, how = "left", multiple = T) |>
  fmutate(surv_90 = if_else(o_d90_01 == 1, 0L, 1L)) |>
  fmutate(efs_90  = if_else(o_e90_01 == 1, 0L, 1L)) |>
  select(
    mrn,
    e_hours_precise,
    e_late_f15_01,
    surv_90,
    efs_90,
    o_crs_01,
    o_crs34_01,
    o_icans_01,
    o_icans34_01,
    lab_category,
    lab_value
  )

bjh_max = rowbind(bjh_max, bjh_new) |> funique()

rm(crp, fer, ldh, crp_cols, fer_cols, ldh_cols, bjh_new); gc()

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

# table of max values by admin time --------------------------------------------

library(gtsummary)
library(gt)

labtable =
  select(labs, mrn, e_late_f15_01, lab_category, lab_value) |>
  funique() |>
  pivot_wider(names_from = lab_category, values_from = lab_value) |>
  select(-mrn) |>
  tbl_summary(
    by   = e_late_f15_01,
    type = everything() ~ "continuous"
  ) |>
  add_p() |>
  as_gt()

gtsave(labtable, filename = here(paste0("tables/labtable_", today, ".docx")))

# plot max values with admin time and outcomes ---------------------------------

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
  scale_fill_manual(values  = c("#5ec962", "#3b528b")) +
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

# plot alc serially ------------------------------------------------------------

## bjh lab data ----------------------------------------------------------------

blab =
  fread(here("pre/CircadianClocksInCAR_DATA_LABELS_2025-07-23_1155.csv"), colClasses = "character") |>
  janitor::clean_names() |>
  fmutate(admin_dttm = lubridate::mdy_hm(date_and_time_of_infusion)) |>
  select(mrn = record_id, admin_dttm, contains("lymphocyte")) |>
  fmutate(e_late_f15_01 = if_else(hour(admin_dttm) > 15, 1L, 0L)) |>
  pivot_longer(
    cols      = contains("lympho"),
    names_to  = "day",
    values_to = "alc"
  ) |>
  fmutate(day = readr::parse_integer(str_extract(day, "(?<=day_)[0-9]+"))) |>
  fmutate(alc = as.numeric(alc)) |>
  fsubset(!is.na(alc))

## ohsu lab data ---------------------------------------------------------------

### setup ----------------------------------------------------------------------

test_names = c("abs_lymphs")
labs       = tidytable()      # need empty table to fill
mrns       = funique(df$mrn) # only patients of interest

### list lab files -------------------------------------------------------------

files = list.files(
  here("../000_data/ohsu/intermediate"),
  pattern    = "labs_intermediate_",
  full.names = TRUE
)

### loop over all files to open ------------------------------------------------

for (file in files) {
  table =
    open_dataset(file) |>
    dplyr::filter(mrn %in% mrns, test_name %in% test_names) |>
    dplyr::select(mrn, CollectionInstant, alc = ValueLCRF) |>
    dplyr::collect()

  labs = rowbind(labs, table, fill = T)
}

labs$CollectionInstant = ymd_hms(labs$CollectionInstant)

### labs from 2025 -------------------------------------------------------------

labs2 = fread(here("../_clif_data/z_source/labs/labs_250101_250630.csv"))

labs2 =
  fsubset(labs2, PrimaryMrn %in% mrns) |>
  fsubset(SpecimenType == "Blood") |>
  fsubset(LabComponentName %in% c("LYMPHOCYTE #", "LYMPHOCYTE#")) |>
  select(mrn = PrimaryMrn, CollectionInstant, alc = ValueLCRF)

### ohsu labs on relevant days -------------------------------------------------

labs =
  rowbind(labs, labs2) |>
  fmutate(alc = as.numeric(alc)) |>
  fsubset(!is.na(alc)) |>
  funique()

labs =
  join(labs, df, how = "inner", multiple = T) |>
  fmutate(day = as.numeric(difftime(CollectionInstant, admin_dttm), "hours")/24) |>
  fmutate(day = round(day)) |>
  fsubset(day >= -1 & day <= 14) |>
  select(mrn, admin_dttm, e_late_f15_01, day, alc)

### tidy up --------------------------------------------------------------------

rm(table, labs2); gc()

## plotting df -----------------------------------------------------------------

labs =
  rowbind(labs, blab) |>
  fsubset(day %in% c(0, 7, 10, 14)) |>
  fgroup_by(mrn, day) |>
  fmutate(alc = fmax(alc)) |>
  fungroup() |>
  fmutate(tgroup = if_else(e_late_f15_01 == 1, "After 15:00", "Before 15:00")) |>
  fmutate(tgroup = factor(tgroup, levels = c("Before 15:00", "After 15:00")))

## make plot -------------------------------------------------------------------

ggplot(labs, aes(x = factor(day), y = alc, fill = tgroup, color = tgroup)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    alpha    = 0.4,
    color    = "black",
    outliers = F
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width  = 0.5,
      jitter.height = 0.01,
      dodge.width   = 0.75
    ),
    size     = 1.2,
    alpha    = 0.66
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1.6)) +
  xlab("Day") +
  ylab("ALC") +
  labs(fill = "Time Group", color = "Time Group") +
  scale_fill_manual(values  = c("#f89540", "#7e03a8")) +
  scale_color_manual(values = c("#f89540", "#7e03a8"))

## save plot -------------------------------------------------------------------

ggsave(here("figures", paste0("boxplot_alc_days_", today, ".pdf")))