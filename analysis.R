
library(tidyverse)
library(haven)
library(survival)
library(tableone)
library(survminer)


dt = read_dta("./external_data/data_prepared.dta") |>
  janitor::clean_names() |>
  mutate(
    expo_lmr_cat = as.factor(expo_lmr_cat),
    expo_sii_cat = as.factor(expo_sii_cat),
    expo_npar_cat = as.factor(expo_npar_cat),
    expo_nps_cat = as.factor(expo_nps_cat),
    sex = as.factor(sex),
    income = as.factor(income),
    marriage = as.factor(marriage),
    smoke_cat = as.factor(smoke_cat)
  )

# baseline ----
vars <- c("sex", "age", "income", "townsend",
          "total_met", "diet_quality", "sleep_hour", "smoke_cat", "total_alcohol",
          "non_hdl", "tg", "bp_cat")

tableone_create = function(dt, expo) {
  
  tableone = 
    dt |>
    CreateTableOne(vars = vars,
                   strata = expo,
                   data = _,
                   factorVars = c("sex", "income", "marriage", "smoke_cat"),
                   test = TRUE,
                   testApprox = chisq.test,
                   testExact = fisher.test,
                   testNonNormal = kruskal.test,
                   testNormal = oneway.test,
                   smd = TRUE,
                   addOverall = TRUE) |>
    print(quote = FALSE, noSpaces = TRUE, printToggle = FALSE) |>
    write.csv(file  = paste0("./csv/", "Baseline_", expo, ".csv"))
    
}

expo_list = c("expo_lmr_cat", "expo_sii_cat", "expo_npar_cat", "expo_nps_cat")
for (expo in expo_list) {
  
  tableone_create(dt, expo) 
  
}


# main survival ----
expo_list = c("expo_lmr_cat", "expo_sii_cat", "expo_npar_cat", "expo_nps_cat")
outcome_list = c("nafld", "cirrho")

#survival cox function
survival_cox = function(dt, model, expo, outcome) {
  
  dt =
    dt |>
    mutate(
      expo = as.factor(.data[[expo]]))

  surv_obj = paste0("Surv(time =", outcome, "_surv_duration, event =", outcome, "_outcome)")
  
  formula = if (model == 1) {
    as.formula(paste0(surv_obj, "~", expo, "+ age + sex"))
  } else {
    as.formula(paste0(surv_obj, "~", expo, "+ age + sex + income + townsend + total_met +
                      diet_quality + sleep_hour + smoke_cat + total_alcohol + non_hdl + tg + bp_cat"))
  }
  
    cox_model = coxph(formula, data = dt)
      
    result = 
      broom::tidy(cox_model) |>
      select(term, estimate, std.error, p.value) |>
      slice(1:2) |>
      mutate(
        HR = round(exp(estimate), 3),
        lower_conf = round(exp(estimate - 1.96 * std.error), 3),
        upper_conf = round(exp(estimate + 1.96 * std.error), 3),
        CI = paste(lower_conf, upper_conf, sep = ", "),
        P = p.value
      ) |>
      select(-estimate, -std.error, -p.value)
    
  return(result)
  
}

# calculte p for trend function
p_for_trend = function(dt, model, outcome, expo) {
  
  dt =
    dt |>
    mutate(
      expo = as.numeric(.data[[expo]]))
  
  surv_obj = paste0("Surv(time =", outcome, "_surv_duration, event =", outcome, "_outcome)")
  
  formula = if (model == 1) {
    as.formula(paste0(surv_obj, "~", expo, "+ age + sex"))
  } else {
    as.formula(paste0(surv_obj, "~", expo, "+ age + sex + income + townsend + total_met +
                      diet_quality + sleep_hour + smoke_cat + total_alcohol + non_hdl + tg + bp_cat"))
  }
  
  cox_model = coxph(formula, data = dt)
  
  result = 
    broom::tidy(cox_model) |>
    select(p.value) |>
    slice(1)
  
  return(result)
}

# table generating function
generate_tables = function(dt, outcome_list, expo_list) {
  results_list = list()
  
  for (outcome in outcome_list) {
    results = data.frame()
    
    for (expo in expo_list) {
      
      res_model_1 = survival_cox(dt, model = 1, expo = expo, outcome = outcome)
      p_trend_1 = p_for_trend(dt, model = 1, outcome = outcome, expo = expo)
      res_model_1 = res_model_1 |>
        mutate(Adjustment = "Model 1", p_for_trend = ifelse(row_number() == 1, p_trend_1, NA))
      
      res_model_2 = survival_cox(dt, model = 2, expo = expo, outcome = outcome)
      p_trend_2 = p_for_trend(dt, model = 2, outcome = outcome, expo = expo)
      res_model_2 = res_model_2 |>
        mutate(Adjustment = "Model 2", p_for_trend = ifelse(row_number() == 1, p_trend_2, NA))
      
      combined_res = rbind(res_model_1, res_model_2)
      results = rbind(results, combined_res)
    
    }
    
    results_list[[outcome]] = results
  }
  
  return(results_list)
}

# using lists to create results lists

results = generate_tables(dt, outcome_list, expo_list)

# pull out the result and tidy the content, exporting as csv
nafld_results = 
  results[["nafld"]] |>
  rename(
    Model = Adjustment,
    Category = term,
    `P for trend` = p_for_trend) |>
  select(Model, everything()) |>
  mutate(across(where(is.list), as.character))

write.csv(nafld_results, file  = "./csv/Main_cox_nafld.csv")
  
cirrho_results =
  results[["cirrho"]] |>
  rename(
    Model = Adjustment,
    Category = term,
    `P for trend` = p_for_trend) |>
  select(Model, everything()) |>
  mutate(across(where(is.list), as.character))

write.csv(cirrho_results, file  = "./csv/Main_cox_cirrhosis.csv")


# cumulative incidence curve ----

for (expo in expo_list) {
  
  for (outcome in outcome_list) {
    
    surv_formula = as.formula(
      paste0("Surv(", outcome, "_surv_duration, ", outcome, "_outcome)~", expo)
    )

    fit = surv_fit(surv_formula, data = dt)
    
    plot = ggsurvplot(fit, data = dt, 
               conf.int = TRUE,
               fun = "cumhaz",
               censor = FALSE,
               pval = TRUE,
               risk.table = "absolute",
               risk.table.col = "strata",
               xlab = "Follow up time(d)", 
               legend = c(0.2,0.8),
               legend.title = "Kaplan-Meier",
               break.x.by = 1000,
               ylim = c(0, 2),
               palette = c("#E7B800", "#2E9FDF", "green4"),
               ggtheme = theme_minimal())
    
    combined_plot =
      arrange_ggsurvplots(
        list(plot),
        print = FALSE
      )
    
    name = paste0("./png/", "cumhaz_", expo , "_", outcome, ".png")
    ggsave(name, plot = combined_plot, width = 16, height = 6, dpi = 300)
    
  }
}
