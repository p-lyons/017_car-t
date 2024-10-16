
### This script identifies the cohort of patients who received CAR-T therapy
### at OHSU during the study period

# setup ------------------------------------------------------------------------

## libraries
library(tidytable)
library(collapse)
library(arrow)
library(here)

## helpers
today = format(Sys.Date(), "%y%m%d")

# data -------------------------------------------------------------------------

## load datasets
df  = open_dataset(here("../000_data/ohsu/intermediate/medadmin_ohsu"))
enc = read_parquet(here("../000_data/ohsu/clean/encounters-all_ohsu.parquet"))
dem = read_parquet(here("../000_data/ohsu/clean/demog_all_240409.parquet"))

## find CAR-T administrations
df <-
  df |>
  dplyr::filter(mar_result %in% c("DUE", "GIVEN", "MISSED", "NEW BAG", "RESTARTED")) |>
  dplyr::select(mrn, enc, admin_dttm, med_nm) |>
  dplyr::collect() |>
  funique() |>
  fmutate(
    keep = case_when(
      stringr::str_detect(med_nm, "LEUCEL") ~ 1L,
      stringr::str_detect(med_nm, "CAR-T")  ~ 1L,
      stringr::str_detect(med_nm, "GENE")   ~ 1L,
      TRUE                                  ~ 0L
    )
  ) |>
  fsubset(keep == 1) |>
  select(-keep) |>
  fsubset(med_nm != "TOCILIZUMAB IV INFUSION (CAR-T USE)") |>
  fsubset(med_nm != "DEXRAZOXANE (BRECKENRIDGE GENERIC ONLY) IV INFUSION") |>
  fsubset(med_nm != "ONASEMNOGENE ABEPARVOVEC-XIOI 2 X 10EXP13 VG/ML IV SUSPENSION,KIT") |>
  slice_min(admin_dttm, .by = enc)

## add demographics and encounters
enc <-
  enc |>
  select(
    mrn,
    enc,
    age = age_yrs,
    admit_dttm = admit_dttm_enc,
    disch_dttm = disch_dttm_enc,
    los_hosp_days,
    dischg_dispo
  )

df <-
  df |>
  left_join(enc) |>
  left_join(dem |> select(mrn, female_01, race_cat, ethn_cat)) |>
  fsubset(!is.na(age))

# cleanup ----------------------------------------------------------------------

## clean product names ---------------------------------------------------------

df <-
  df |>
  fmutate(
    product_cat = case_when(
      stringr::str_detect(med_nm, "AXICABTAGENE")     ~ "axi-cel",
      stringr::str_detect(med_nm, "BREXUCABTAGENE")   ~ "brexu-cel",
      stringr::str_detect(med_nm, "CILTACABTAGENE")   ~ "cilta-cel",
      stringr::str_detect(med_nm, "TISAGENLECLEUCEL") ~ "tisa-cel",
      stringr::str_detect(med_nm, "LISOCABTAGENE")    ~ "liso-cel",
      stringr::str_detect(med_nm, "POSOLEUCEL")       ~ "poso-cel",
      stringr::str_detect(med_nm, "LIFILEUCEL")       ~ "lifi-cel",
      stringr::str_detect(med_nm, "TABELECLEUCEL")    ~ "ltab-cel",
      stringr::str_detect(med_nm, "SIPULEUCEL")       ~ "sipu-cel",
      stringr::str_detect(med_nm, "IDECABTAGENE")     ~ "ide-cel",
      TRUE                                          ~ "help"
    ),
    product_cat = qF(product_cat),
    inv_01      = if_else(stringr::str_starts(med_nm, "INV|ZZZ"), 1L, 0L)
  ) |>
  select(-med_nm)

## time of day -----------------------------------------------------------------

df <-
  df |>
  fmutate(
    tod_hrs = lubridate::hour(admin_dttm),
    tod_cat = case_when(
      tod_hrs >= 08 & tod_hrs < 10 ~ "08:00-10:59",
      tod_hrs >= 10 & tod_hrs < 14 ~ "11:00-13:59",
      TRUE                         ~ "14:00-20:00"
    ),
    tod_late_01  = if_else(tod_hrs == "14:00-20:00", 1L, 0L)
  )

## confounders -----------------------------------------------------------------

df <-
  df |>
  fmutate(
    los_pre_days = as.numeric(difftime(admin_dttm, admit_dttm), "days"),
    year         = lubridate::year(admin_dttm)
  )

## outcomes --------------------------------------------------------------------

### for chart review
select(df, mrn, admit_dttm) |> funique() |> fwrite(here("misc", "for_charts.csv"))

ddt = lubridate::ymd(death_date),
d90 = if_else(ddt > time + lubridate::ddays(90) | is.na(ddt), 0L, 1L)

