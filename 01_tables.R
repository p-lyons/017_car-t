
# This script makes tables for the CAR T-cell circadian pattern manuscript.

# setup ------------------------------------------------------------------------

## libraries -------------------------------------------------------------------

library(broom.mixed)
library(gtsummary)
library(tidytable)
library(collapse)
library(stringr)
library(arrow)
library(here)
library(lme4)
library(gt)

## helpers ---------------------------------------------------------------------

today  = format(Sys.Date(), "%y%m%d")
cutoff = as.POSIXct("2024-06-30") # which encounters have 1 year obs time?

# load data --------------------------------------------------------------------

## most recent analysis file ---------------------------------------------------

files   = list.files(here("clean"), pattern = "^analysis_\\d{6}\\.csv$", full.names = T)
dates   = str_extract(files, "\\d{6}")
latest  = files[which.max(dates)]

## categorize cancer_cat -------------------------------------------------------

df_main_90 =
  fread(latest) |>
  fmutate(
    cancer_cat_4 = case_when(
      cancer_cat %in% c("MCL", "Follicular Lymphoma") ~ "Other Lymphomas",
      TRUE                                            ~ cancer_cat
    ),
    cancer_cat_3 = case_when(
      cancer_cat %in% c("DLBCL", "Other Lymphomas")   ~ "Lymphoma",
      TRUE                                            ~ cancer_cat
    ),
    cancer_cat_2 = case_when(
      cancer_cat == "Multiple Myeloma"                ~ "Multiple Myeloma",
      TRUE                                            ~ "Leukemia/Lymphoma"
    )
  )

## df for 365-day outcomes -----------------------------------------------------

df_main_365 = fsubset(df_main_90, admin_dttm <= cutoff)

message("Loaded file: ", latest)
rm(files, dates, latest); gc()

## basic overview --------------------------------------------------------------

fnunique(df_main_90$mrn)
fnunique(df_main_365$mrn)
janitor::tabyl(df_main_90$hospital)
janitor::tabyl(df_main_365$hospital)

# table 1 (characteristics and outcomes by hospital) ---------------------------

## define table 1 variables ----------------------------------------------------

t1_vars = c(
  "hospital",
  "year_cat",
  "season",
  "e_hours_precise",
  "e_hours_sunrise",
  "age",
  "female_01",
  "nhw_01",
  "cancer_cat",
  "product_cat",
  "ldh",
  "conditioning_flucy_01",
  "outpt_01",
  "los_pre_cart",
  "ecog",
  "vw",
  "saps2",
  "tocilizumab_01",
  "anakinra_01"
)

## create 90-day table ---------------------------------------------------------

table_01 =
  select(df_main_90, all_of(t1_vars), starts_with("o_")) |>
  select(-o_d90_01, -o_e90_01, -contains("365"), -ends_with("days")) |>
  tbl_summary(
    by        = hospital,
    type      = list(ecog ~ "continuous"),
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})")
  ) |>
  add_p() |>
  as_gt()

## create 365-day table --------------------------------------------------------

table_01_365 =
  select(df_main_365, hospital, o_os_365, o_efs_365) |>
  tbl_summary(by = hospital) |>
  add_p() |>
  as_gt()

## save tables -----------------------------------------------------------------

table_01
table_01_365

gt::gtsave(table_01,     here("tables", paste0("table_01_hospital_",  today, ".docx")))
gt::gtsave(table_01_365, here("tables", paste0("table_01a_hospital_", today, ".docx")))

# table 2 (characteristics by early/late @ 1500) -------------------------------

## create table 2 --------------------------------------------------------------

table_02 =
  select(df_main_90, all_of(t1_vars), e_late_f15_01) |>
  tbl_summary(
    by        = e_late_f15_01,
    type      = list(ecog ~ "continuous"),
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})")
  ) |>
  add_p() |>
  as_gt()

## save table 2 ----------------------------------------------------------------

table_02

gt::gtsave(table_02, here("tables", paste0("table_02_late15f_",  today, ".docx")))

# table 3 (outcomes by early/late @ 1500) --------------------------------------

table_03 =
  select(df_main_90, e_late_f15_01, starts_with("o_"), anakinra_01, tocilizumab_01) |>
  select(-o_d90_01, -o_e90_01, -contains("365"), -ends_with("days")) |>
  tbl_summary(
    by        = e_late_f15_01,
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})")
  ) |>
  add_p() |>
  as_gt()

## create 365-day table --------------------------------------------------------

table_03_365 =
  select(df_main_365, e_late_f15_01, o_os_365, o_efs_365) |>
  tbl_summary(by = e_late_f15_01) |>
  add_p() |>
  as_gt()

## save tables -----------------------------------------------------------------

table_03
table_03_365

gt::gtsave(table_03,     here("tables", paste0("table_03_hospital_",  today, ".docx")))
gt::gtsave(table_03_365, here("tables", paste0("table_03a_hospital_", today, ".docx")))


# table 4 -----------------





