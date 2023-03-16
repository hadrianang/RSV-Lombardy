#The code here reads in some input data then passes it on
#to a Stan model for the fitting process. 

library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#All paths relative to base directory (assumed to be source file location)
base_dir = getwd()
model_input_dir = file.path(base_dir, "ModelInputFiles")

#Start by reading in the data 
input_data = readRDS(file.path(model_input_dir, "InputData.rds"))

hosp_vals = input_data$hosp_mat
surv_vals = input_data$surv_mat
test_vals = input_data$test_mat
pos_vals = input_data$pos_mat
pop_vals = input_data$pop_mat


n_seasons = 4
n_pop_age = 102 #0-100, including indices for 0-6m and 7-12m 
n_hosp_age = 83 #0, 0.5, 1 to 80, 81+
n_surv_age = 4 #0-4, 5-14, 15-64, 65+

hosp_range_lower = c(1:83) #0 to 81 years-old
hosp_range_upper = c(1:82, n_pop_age) #0 to 80 then 100 years-old
surv_range_lower = c(1, 7, 17, 67) #0, 5, 15, 65 years-old
surv_range_upper = c(6, 16, 66, n_pop_age) #4, 14, 64, 100 years-old

pop_ages = c(0, 0, 1:(n_pop_age - 2))
real_pop_ages = c(0, 0.5, 1:(n_pop_age - 2))

#Parameterisation here based on work by van Boven et al. (see paper for full citation)
n_rho_hosp_group = 5 #five hospitalisation groups: 0-6m, 7-12m, 1, 2, 3+
rho_hosp_lower = c(1, 2, 3, 4, 5)
rho_hosp_upper = c(1, 2, 3, 4, 83)
rho_hosp_priors = c(0.014, 0.014, 0.0015, 0.0015, 0.00070)

n_rho_surv_group = 4 #four surveillance groups: 0-4, 5-14, 15-64, 65+
rho_surv_lower = c(1, 2, 3, 4)
rho_surv_upper = c(1, 2, 3, 4)
rho_surv_priors = c(0.19, 0.17, 0.17, 0.15)

n_rep_groups = 5 #five reporting groups for attack rate: 0-6m, 7-12m, 1-4, 5-14, 15+
rep_groups_lower = c(1, 2, 3, 7, 17)
rep_groups_upper = c(1, 2, 6, 16, n_pop_age)

weight_7_12 = 0.3 #unused
s = 0.77
model_code = 4
F_code = 2
r_code = 1

data = list(
  n_seasons = n_seasons,
  n_pop_age = n_pop_age, 
  n_hosp_age = n_hosp_age, 
  n_surv_age = n_surv_age, 
  hosp_range_lower = hosp_range_lower,
  hosp_range_upper = hosp_range_upper,
  surv_range_lower = surv_range_lower,
  surv_range_upper = surv_range_upper, 
  pop_ages = pop_ages,
  real_pop_ages = real_pop_ages,
  pop_vals = pop_vals,
  hosp_vals = hosp_vals,
  surv_vals = surv_vals,
  pos_vals = pos_vals, 
  test_vals = test_vals,
  n_rho_hosp_group = n_rho_hosp_group,
  n_rho_surv_group = n_rho_surv_group, 
  rho_hosp_lower = rho_hosp_lower,
  rho_hosp_upper = rho_hosp_upper, 
  rho_surv_lower = rho_surv_lower, 
  rho_surv_upper = rho_surv_upper,
  rho_hosp_priors = rho_hosp_priors,
  rho_surv_priors = rho_surv_priors, 
  n_rep_groups = n_rep_groups,
  rep_groups_lower = rep_groups_lower,
  rep_groups_upper = rep_groups_upper,
  s = s,
  weight_7_12 = weight_7_12,
  model_code = model_code,
  F_code = F_code,
  r_code = r_code
)
pars = c("lp__", "beta", "beta_child", "beta_ya", "rho_hosp_a", "rho_surv_a", "rho_test_a", 
                      "Prob_test_summ", "Prob_Inf_Total", "Prob_pInf", "Prob_sInf",
                      "Theta_hosp", "hosp_totals_model", "Sigma_surv", "surv_totals_model", "prop_non_naive_before", 
                      "prop_non_naive_after", "log_lik", "pop_non_naive_before", 
                      "total_non_naive_after", "total_non_naive_before", "trans_rate_inf_c", "trans_rate_inf_ya", "trans_rate_inf_a", 
                      "trans_rate_ili_c", "trans_rate_ili_ya", "trans_rate_ili_a", "r", "attack_rates_rep_groups", "rsv_naive_rep_groups_before",
                      "rsv_naive_rep_groups_after")


model = stan_model("CatalyticModel.stan")
fit_model = sampling(model,
                     data = data,
                     iter = 20000,
                     chains = 4, refresh = 100, verbose = TRUE, pars = pars,  
                     seed = 0)

saveRDS(fit_model, file.path(base_dir, "ModelOutputFiles", "Model_D_F_code=2_R_code=1_20k_v2.rds"))
