

# setup ------------------------------------------------------------------------

## libraries
library(tidytable)
library(collapse)
library(stringr)
library(arrow)
library(here)


## helpers
today = format(Sys.Date(), "today")

# data -------------------------------------------------------------------------

df   = read_parquet()
mrns = funique(df$mrn)
meds = fread("Desktop/car-t_other-meds.csv") |> fsubset(PrimaryMrn %in% mrns)

# remove unecessary meds -------------------------------------------------------

routes = c("*Unspecified", "injection", "intravenous", "oral", "subcutaneous")

meds <-
  meds |>
  select(
    mrn        = PrimaryMrn,
    csn        = EncounterEpicCsn,
    admin_dttm = AdministrationInstant,
    Name,
    Dose,
    DoseUnit,
    Route
  ) |>
  janitor::clean_names() |>
  fsubset(dose_unit != "drop") |>
  fsubset(route %in% routes) |>
  fsubset(
    !name %in% c(
      "HEPARIN-LIDOCAINE-HYDROCORTISONE INTRAVESICAL INSTILLATION",
      "HEPATIC ARTERIAL INFUSION (HAIP) MEDTRONIC PUMP - HEPARIN +/- DEXAMETHASONE",
      "INV ROPIVACAINE 0.4%-DEXAMETHASONE 1 MG PPF BLOCK SOLUTION (IRB 26224)",
      "INV ROPIVACAINE 20 MG-DEXAMETHASONE 4 MG OR PLACEBO 5 ML PPF BLOCK SOLUTION (IRB 26224)",
      "ONE-STEP - LIDOCAINE 1%-EPINEPHRINE 1:100K - DEXAMETHASONE 4 MG/ML",
      "HEPATIC ARTERIAL INFUSION (HAIP) MEDTRONIC PUMP - FLOXURIDINE + DEXAMETHASONE",
      "ZZZ INV DEXAMETHASONE/PLACEBO 1 MG/ML ORAL SOLUTION (IRB 8071)",
      "ZZZ INV HYDROCORTISONE SODIUM SUCCINATE (SOLU-CORTEF)/PLACEBO 100 MG INJECTION (IRB 18320)",
      "ZZZ INV METHYLPREDNISOLONE/PLACEBO IV"
    )
  )

# final cleanup ----------------------------------------------------------------

meds <-
  meds |>
  fmutate(
    mrn      = as.integer(mrn),
    csn      = str_remove_all(csn, ","),
    csn      = as.integer(csn),
    dose     = as.numeric(dose),
    med_name = case_when(
      str_detect(name, "TOCILIZ") ~ "tocilizumab",
      str_detect(name, "ANAKINR") ~ "anakinra",
      str_detect(name, "DEXAMET") ~ "dexamethasone",
      str_detect(name, "HYDROCO") ~ "hydrocortisone",
      str_detect(name, "METHYLP") ~ "methylprednisolone",
      str_detect(name, "ISOLONE") ~ "prednisolone",
      TRUE                        ~ "prednisone"
    )
  ) |>
  select(-route, -name) |>
  fsubset(!is.na(mrn)) |>
  funique()



toci     = fsubset(meds, med_name == "tocilizumab") |> select(-dose_unit)
anakinra = fsubset(meds, med_name == "anakinra") |> select(-dose_unit)
meds     = fsubset(meds, !med_name %in% c("anakinra", "tocilizumab"))

