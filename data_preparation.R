
library(haven)
library(tidyverse)

# data import -----
dt = read_dta("./data/inf_0.csv") |>
  janitor::clean_names()

# function ----
cut_1_23_4 = function(x) {
  
  # cut and produce Q1 ~ 0, Q2-3 ~ 1, Q4 ~ 2
  x_cat = as.factor(cut(x, 
              breaks = quantile(x, c(0, 0.25, 0.75, 1), na.rm = TRUE), 
              include.lowest = TRUE, 
              labels = c("1", "2", "3")))
  
  return (x_cat)

}

# data cleaning ----
dt_slim =
  dt |>
  mutate(
    sex = as.factor(sex), #sex as factor, 0 as female, 1 as male
    ethnicity_cat = as.factor(ethnicity_cat), # ethnicity as factor, 0 as white, 1 as others
    age_cat = as.factor(case_when(
      age < 50 ~ 0,
      age >= 50 & age < 60 ~ 1,
      age >= 60 & age < 70 ~ 2,
      age >= 70 ~3
    )),
    income = as.factor(income),
    smoke_cat = as.factor(smoke_cat),
    bmi = n_21001_0_0,
    nafld_date = s_131670_0_0,
    cirrho_date = s_131666_0_0
  ) |>
  filter(
    !((!is.na(nafld_date)) & nafld_date <= baseline_date),
    !((!is.na(cirrho_date)) & cirrho_date <= baseline_date),
  ) |> # excluding participants having disease at baseline
  mutate(
    nafld_outcome = as.factor(case_when(
      is.na(nafld_date) ~ 0,
      (!is.na(nafld_date)) & nafld_date > baseline_date ~ 1
    )),
    cirrho_outcome = as.factor(case_when(
      is.na(cirrho_date) ~ 0,
      (!is.na(cirrho_date)) & cirrho_date > baseline_date ~ 1
    )) 
  ) |> # assign values to outcome status, 0 as healthy, 1 as diseased
  mutate(
    nafld_date = ifelse(is.na(nafld_date) & !is.na(death_date), death_date, nafld_date),
    nafld_date = ifelse(is.na(nafld_date) & !is.na(lost_date), lost_date, nafld_date),
    nafld_date = ifelse(is.na(nafld_date), as.Date("2022-10-31"), nafld_date),
    cirrho_date = ifelse(is.na(cirrho_date) & !is.na(death_date), death_date, cirrho_date),
    cirrho_date = ifelse(is.na(cirrho_date) & !is.na(lost_date), lost_date, cirrho_date),
    cirrho_date = ifelse(is.na(cirrho_date), as.Date("2022-10-31"), cirrho_date)
  ) |> # replace those dead, lost to follow up or no disease at the end
  mutate(
    nafld_surv_duration = nafld_date - baseline_date,
    cirrho_surv_duration = cirrho_date - baseline_date
  ) |> # calculate survival duration
  mutate(
    b_lympho = n_30120_0_0,
    b_mono = n_30130_0_0,
    b_plate = n_30080_0_0,
    b_neutro = n_30140_0_0,
    b_neutro_percent = n_30200_0_0,
    b_album = n_30600_0_0,
    b_rbc = n_30010_0_0,
    b_wbc = n_30000_0_0,
  ) |> #renaming blood test as exposure
  drop_na(b_lympho:b_album) |> #drop NAs in blood test variables as exposure
  mutate(
    expo_lmr = b_lympho / b_mono,
    expo_sii = b_neutro * b_plate / b_lympho,
    expo_npar = b_neutro_percent / b_album,
      cri_album = ifelse(b_album < 4, 1, 0),
      cri_nlr = ifelse(b_neutro / b_lympho > 3, 1, 0),
      cri_lmr = ifelse(b_lympho / b_mono < 4.44, 1, 0),
      cri_plr = ifelse(b_plate / b_lympho > 150, 1, 0),
    expo_nps = rowSums(across(c(cri_album, cri_nlr, cri_lmr, cri_plr))),
    
    expo_lmr_cat = cut_1_23_4(expo_lmr),
    expo_sii_cat = cut_1_23_4(expo_sii),
    expo_npar_cat = cut_1_23_4(expo_npar),
    expo_nps_cat = case_when(
      expo_nps == 0 ~ 1,
      expo_nps == 1 | expo_nps == 2 ~ 2,
      expo_nps == 3 | expo_nps == 4 ~ 3
    )
  ) |> # calculate exposure values and cut them into categorical values
  select(
    -starts_with("s_"),
    -starts_with("n_"),
    ends_with("_cat"),
    ends_with("_date"),
    ends_with("_outcome"),
    ends_with("_origin"),
    starts_with("b_"),
    starts_with("expo_"),
    starts_with("cri_"),
    n_eid, sex, age, ethnicity_cat, townsend, income, marriage,
    total_met, diet_quality, sleep_hour, smoke_cat, total_alcohol,
    tc, hdl, non_hdl, tg, hba1c, dbp, sbp, bp_cat
  )

write_dta(dt_slim, "./external_data/data_prepared.dta")


# data for msm ----
library(mstate)

dt_msm =
  dt |>
  mutate(
    sex = as.factor(sex), #sex as factor, 0 as female, 1 as male
    ethnicity_cat = as.factor(ethnicity_cat), # ethnicity as factor, 0 as white, 1 as others
    bmi = n_21001_0_0,
    nafld_date = s_131670_0_0,
    cirrho_date = s_131666_0_0
  ) |>
  filter(
    !((!is.na(nafld_date)) & nafld_date <= baseline_date),
    !((!is.na(cirrho_date)) & cirrho_date <= baseline_date),
  ) |> # excluding participants having disease at baseline
  mutate(
    nafld_outcome = as.factor(case_when(
      is.na(nafld_date) ~ 0,
      (!is.na(nafld_date)) & nafld_date > baseline_date ~ 1
    )),
    cirrho_outcome = as.factor(case_when(
      is.na(cirrho_date) ~ 0,
      (!is.na(cirrho_date)) & cirrho_date > baseline_date ~ 1
    )),
    death_outcome = as.factor(case_when(
      is.na(death_date) ~ 0,
      (!is.na(death_date)) & death_date > baseline_date ~ 1
    ))
  ) |> # assign values to outcome status, 0 as healthy, 1 as having
  filter(
    is.na(death_date) | 
      (nafld_date <= death_date & cirrho_date <= death_date)
  ) |> # drop those reprot disease later than death
  mutate(
    nafld_date = ifelse(
      !is.na(nafld_date) & !is.na(cirrho_date) & nafld_date == cirrho_date,
      cirrho_date - 1,
      nafld_date
    )
  ) |> # change nafld date a day forward if participants report having disease on the same day
  mutate(
    nafld_date = ifelse(is.na(nafld_date) & !is.na(lost_date), lost_date, nafld_date),
    nafld_date = ifelse(is.na(nafld_date), as.Date("2022-10-31"), nafld_date),
    cirrho_date = ifelse(is.na(cirrho_date) & !is.na(lost_date), lost_date, cirrho_date),
    cirrho_date = ifelse(is.na(cirrho_date), as.Date("2022-10-31"), cirrho_date),
    death_date = ifelse(is.na(death_date) & !is.na(lost_date), lost_date, death_date),
    death_date = ifelse(is.na(death_date), as.Date("2022-10-31"), death_date),
  ) |> # replace those lost to follow up or no disease at the end
  mutate(
    nafld_surv_duration = nafld_date - baseline_date,
    cirrho_surv_duration = cirrho_date - baseline_date
  ) |> # calculate survival duration
  mutate(
    b_lympho = n_30120_0_0,
    b_mono = n_30130_0_0,
    b_plate = n_30080_0_0,
    b_neutro = n_30140_0_0,
    b_neutro_percent = n_30200_0_0,
    b_album = n_30600_0_0,
    b_rbc = n_30010_0_0,
    b_wbc = n_30000_0_0,
  ) |> #renaming blood test as exposure
  drop_na(b_lympho:b_album) |> #drop NAs in blood test variables as exposure
  mutate(
    expo_lmr = b_lympho / b_mono,
    expo_sii = b_neutro * b_plate / b_lympho,
    expo_npar = b_neutro_percent / b_album,
    cri_album = ifelse(b_album < 4, 1, 0),
    cri_nlr = ifelse(b_neutro / b_lympho > 3, 1, 0),
    cri_lmr = ifelse(b_lympho / b_mono < 4.44, 1, 0),
    cri_plr = ifelse(b_plate / b_lympho > 150, 1, 0),
    expo_nps = rowSums(across(c(cri_album, cri_nlr, cri_lmr, cri_plr))),
    
    expo_lmr_cat = cut_1_23_4(expo_lmr),
    expo_sii_cat = cut_1_23_4(expo_sii),
    expo_npar_cat = cut_1_23_4(expo_npar),
    expo_nps_cat = case_when(
      expo_nps == 0 ~ 1,
      expo_nps == 1 | expo_nps == 2 ~ 2,
      expo_nps == 3 | expo_nps == 4 ~ 3
    )
  ) |> # calculate exposure values and cut them into categorical values
  select(
    -starts_with("s_"),
    -starts_with("n_"),
    ends_with("_cat"),
    ends_with("_date"),
    ends_with("_outcome"),
    ends_with("_origin"),
    starts_with("b_"),
    starts_with("expo_"),
    starts_with("cri_"),
    n_eid, sex, age, ethnicity_cat, townsend, income, marriage,
    total_met, sleep_hour, total_alcohol, tc, hdl, hba1c, dbp, sbp
  ) |>
  
  pivot_longer(
    cols = c(nafld_date, cirrho_date, death_date),
    names_to = "state",
    values_to = "state_date"
  ) |> 
  mutate(
    state_num = case_when(
      state == "nafld_date" ~ 2,
      state == "cirrho_date" ~ 3,
      state == "death_date" ~ 4,
      TRUE ~ 1
    ),
    time_start = baseline_date,
    time_end = state_date    
  )

# transform matrix for msm ----
tmat =
  transMat(list(
    c(2, 3, 4),  # healthy to nafld, cirrhosis, death
    c(3, 4),     # nafld to Cirrhosis, Death
    c(4),        # cirrhosis to death
    c()          # absorbing state
  ), names = c("Healthy", "NAFLD", "Cirrhosis", "Death"))

print(tmat)

dt_msm =
  dt_msm |>
  msprep(
    trans = tmat,
    time = c("time_start", "time_end"),
    status = "state_num",
    keep = c("age", "sex", "expo_lmr_cat")
  )

write_dta("./data/data_msm.dta")


  