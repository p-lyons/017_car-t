


# Draw KM Curve ----------------------------------------------------------------




df <-
  df |>
  mutate(
    age = as.numeric(difftime(time, birth_date)/365.25),
    ddt = lubridate::ymd(death_date),
    d90 = if_else(ddt > time + lubridate::ddays(90) | is.na(ddt), 0L, 1L),
    h   = lubridate::hour(time),
    tod  = case_when(
      h >= 08 & h < 10 ~ "08:00-10:59",
      h >= 10 & h < 14 ~ "11:00-13:59",
      TRUE             ~ "14:00-20:00"
    ),
    evening_01 = if_else(tod == "14:00-20:00", 1L, 0L)
  )


df <-
  df |>
  mutate(
    today = lubridate::ymd(Sys.Date()),
    days_elapsed = as.numeric(difftime(today, time)),
    status = if_else(is.na(ddt), 1L, 2L)
  )

library("survival")
fit <- survfit(Surv(days_elapsed, status) ~ evening_01, data = df)


library("survminer")
library("ggplot2")

p <- ggsurvplot(fit, data = df)

p +
  labs(x = "Days Elapsed")



mrnlist = funique(df$mrn)

lab =
  read_parquet(here("../000_data/ohsu/intermediate/labs_ohsu_temp_big.parquet"), as_data_frame = F) |>
  dplyr::filter(mrn %in% mrnlist) |>
  dplyr::collect() |>
  fsubset(test_name == "il_06") |>
  select(mrn, result_dttm, il_6 = result) |>
  fmutate(
    il_6 = stringr::str_remove_all(il_6, "<"),
    il_6 = as.numeric(il_6)
  )

df |>
  left_join(lab) |>
  fsubset(result_dttm >= admin_dttm) |>
  fmutate(
    hour = lubridate::hour(admin_dttm),
    tod  = case_when(
      hour >= 08 & hour < 10 ~ "08:00-09:59",
      hour >= 10 & hour < 14 ~ "10:00-13:59",
      TRUE                   ~ "14:00-20:00"
    )
  ) |>
  fgroup_by(tod, mrn) |>
  fsummarize(il_6 = fmax(il_6)) |>
  fgroup_by(tod) |>
  fsummarize(il_6 = fmedian(il_6))
