/*
Title: Catalytic Model
Date: 10/01/2023

The code here is an implementation of the catalytic model for the
RSV paper. It allows us to estimate key parameters such as seasonal
FOI, and probabilities of hospitalisation and reporting to surveillance
by fitting to hospitalisation and surveillance data. Three age-groups
are taken into account, 0-4, 5-14, and 15+ years. 

We implement 4 variants of the model here
Model A: Uniform FOI across all age-groups
Model B: FOI among 5-14 = FOI among 0-4
Model C: FOI among 5-14 = FOI among 15+
Model D: FOI allowed to vary across all age-groups. 

The model here has been modified to account for maternally derived immunity
in the 0-6 months age-group lasting for 4 months (based on the paper 
by Ochola et al.). We propose three different hypotheses regarding maternally 
derived immunity, reflected in the value of the variable F (which depends on the F_code). Note that
F has been renamed to mu in the paper and supplementary material. 

1. F is a time-dependent value F(y) which is a function of the average 
probability of RSV infection among women of childbearing age (from age a1 to a2). -- not implemented
2. Maternally derived immunity depends only on prior exposure to RSV, and since almost all
individuals of childbearing age would have experienced primary RSV infection, then F = 1.
3. No maternally derived immunity so F = 0. 

In this updated model, we also test two hypotheses regarding r, a factor denoting a reduced
probability of hospitalisation or reporting to surveillance relative to primary infections.

1. r = 1 and post-primary and primary infections have the same probability 
2. r is a calibrated parameter with r ~ U(0,1)

Note that the beta values described here have been renamed lambda in the paper.
*/

functions{
    /*
    Function: beta_sum
    params:
        beta_t: beta_a values per season (usually for 15+ years)
        beta_child_t: beta_child values per season (usually for 0-4 years)
        beta_ya_t: beta_ya values per season (usually for 5-14 years)
        a: age of individual we are computing the probability of infection for
        t: current time 
        T: breakpoint for 1st age-group (usually 4)
        T2: breakpoint for 2nd age-group (usually 14)

    This function gets the sum of appropriate beta values used in the computations for the
    probability of RSV infection. In other words, it collects the beta values (choosing from
    the beta for 0-4, 5-14, and 15+ years) for each year of life of an individual age a in season t from 
    1 to (a-1) years-old. The value of the FOI for the first year of life has to be computed separately due to 
    maternally derived immunity being a factor.
    */
    real beta_sum(vector beta_t, vector beta_child_t, vector beta_ya_t, int a, int t, int T, int T2){
        //vector to collect the beta values we will be summing, from 1 to a-1 years
        vector[a-1] beta_collect; 

        //From 1 years to a-1 years (the final beta at a years is selected separately)
        for(i in 1:(a-1)){
            int curr_age = i;

            //Curent season given individual is a-years-old at time t
            int y = t - a + i;
            real curr_beta;
            
            if(curr_age <= T){
                //If age is <= the first breakpoint then beta_child
                curr_beta = beta_child_t[y]; 
            }else if((curr_age > T) && (curr_age <= T2)){
                //If T < age <= T2 then beta_ya
                curr_beta = beta_ya_t[y]; 
            }else{
                //Otherwise, we use beta_a
                curr_beta = beta_t[y]; 
            }
            beta_collect[i] = curr_beta; 
        }
        //Return the sum of collected beta values 
        return sum(beta_collect); 
    }
}

/*
Data values (note that some indices related to age are offset due to the presence of 0-6m and 7-12m age-groups)
n_seasons: Number of seasons of data to fit to
n_pop_age: Number of age-groups in the population
n_hosp_age: Number of age-groups in hospitalisation data
n_surv_age: Number of age-groups in surveillance data
hosp_range_lower: Array with lower bound of each hospitalisation age-group
hosp_range_upper: Array with upper bound of each hospitalisation age-group
surv_range_lower: Array with lower bound of each surveillance age-group
surv_range_upper: Array with upper bound of each surveillance age-group
pop_ages: Age values in integer type (both 0-6 and 7-12 months are 0)
real_pop_ages: Age values in real type (7-12 months are 0.5)
pop_vals: Matrix containing the population per age-group, per season 
hosp_vals: 2D array containing the number of hospitalisations per age-group, per season
surv_vals: 2D array containing the number of RSV-attributable-ILI cases per age-group, per season
test_vals: 2D array containing the number of RSV tests carried out per age-group, per season
pos_vals: 2D array containing the number of positive RSV tests per age-group, per season

n_rho_hosp_group: Number of age-groups for probability of hospitalisation
n_rho_surv_group: Number of age-groups for probability of reporting to surveillance
rho_hosp_lower: Array with lower bound of each probability of hospitalisation age-group
rho_hosp_upper: Array with upper bound of each probability of hospitalisation age-group
rho_surv_lower: Array with lower bound of each probability of reporting to surveillance age-group
rho_surv_upper: Array with upper bound of each probability of reporting to surveillance age-group

rho_hosp_priors: Mean for the prior of the probability of hospitalisation
rho_surv_priors: Mean for the prior of the probability of reporting to surveillance

n_rep_groups: Number of reporting groups for attack rate and proportion RSV naive
rep_groups_lower: Array with lower bound of each reporting group for attack rate and proportion RSV naive 
rep_groups_upper: Array with upper bound of each reporting group for attack rate and proportion RSV naive 

weight_7_12: Weight of the previous season for probability of infection of 7-12 months -- no longer used
s: Factor denoting a reduced susceptibility to post-primary infection
model_code: Integer denoting the model variant to run (1 = A, 2 = B, 3 = C, 4 = D)
F_code: Integer denoting the hypothesis to run regarding maternally derived immunity (F) 
        (1: hypothesis 1 (F depends on infections in mothers during previous season) --not implemented, 
        2: hypothesis 2 (F = 1), 3: hypothesis 3 (F = 0))
r_code: Integer denoting the hypothesis to run regarding the probability of hospitalisation
        due to post-primary infection relative to primary infection (r) [1: r = 1, 2: r is calibrated]
*/

data{
    int n_seasons; 
    int n_pop_age;
    int n_hosp_age; 
    int n_surv_age; 

    array[n_hosp_age] int hosp_range_lower;
    array[n_hosp_age] int hosp_range_upper;
    array[n_surv_age] int surv_range_lower;
    array[n_surv_age] int surv_range_upper;
    
    array[n_pop_age] int pop_ages; 
    vector[n_pop_age] real_pop_ages;

    matrix<lower = 0>[n_pop_age, n_seasons] pop_vals; 
    array[n_hosp_age, n_seasons] int<lower = 0> hosp_vals;
    array[n_surv_age, n_seasons] int<lower = 0> surv_vals; 
    
    array[n_surv_age, n_seasons] int<lower = 0> test_vals; 
    array[n_surv_age, n_seasons] int<lower = 0> pos_vals; 

    int n_rho_hosp_group;
    int n_rho_surv_group;
    array[n_rho_hosp_group] int rho_hosp_lower;
    array[n_rho_hosp_group] int rho_hosp_upper;
    array[n_rho_surv_group] int rho_surv_lower;
    array[n_rho_surv_group] int rho_surv_upper;

    vector[n_rho_hosp_group] rho_hosp_priors;
    vector[n_rho_surv_group] rho_surv_priors;

    int n_rep_groups; 
    array[n_rep_groups] int rep_groups_lower;
    array[n_rep_groups] int rep_groups_upper;

    //Weight is not used anymore due to maternally derived immunity
    real weight_7_12; 
    
    real s;
    int model_code;
    int F_code; 
    int r_code;
}

transformed data {
    //Breakpoints for age-groups
    int T = 4;
    int T2 = 14; 

    //Compute total number of hospitalisations and RSV-attributable-ILI cases per season
    array[n_seasons] int hosp_totals_data;
    array[n_seasons] int surv_totals_data;

    for(i in 1:n_seasons){
        hosp_totals_data[i] = sum(hosp_vals[ ,i]);
        surv_totals_data[i] = sum(surv_vals[ ,i]);
    }

    //Index for the start of model time depending on the number of age-groups modelled
    int t_start = n_pop_age - 1; 
    int n_beta = n_seasons;
    
    //These matrices contain the size of the population for each 
    //hospitalisation and surveillance age-group
    matrix[n_hosp_age, n_seasons] pop_hosp_age; 
    matrix[n_surv_age, n_seasons] pop_surv_age;
    matrix[n_rep_groups, n_seasons] pop_rep_groups; 

    //Proportion of each population age in each reporting group per season
    //Ex: proportion of those 7 years in the 5-14 age-group

    //This is used to compute Pi, the (population-weighted) average probability of 
    //infection for an individual in a given age-group
    matrix<lower = 0, upper = 1>[n_pop_age, n_seasons] hosp_group_props; 
    matrix<lower = 0, upper = 1>[n_pop_age, n_seasons] surv_group_props;
    matrix<lower = 0, upper = 1>[n_pop_age, n_seasons] rep_group_props;

    //We compute the values for the matrices above
    for(j in 1:n_seasons){
        for(i in 1:n_hosp_age){
            int curr_lower = hosp_range_lower[i];  
            int curr_upper = hosp_range_upper[i]; 
            //First, we get the sum of the population for each hospitalisation group
            pop_hosp_age[i, j] = sum(pop_vals[curr_lower:curr_upper, j]);

            //We compute the proportion of the hospitalisation group in each age
            //by dividing the population per age by the total size of the group 
            hosp_group_props[curr_lower:curr_upper, j] = pop_vals[curr_lower:curr_upper, j] ./ pop_hosp_age[i,j];
        }
    }

    //We mirror the operations for hospitalisations above and apply them for RSV-attributable ILI
    for(j in 1:n_seasons){
        for(i in 1:n_surv_age){
            int curr_lower = surv_range_lower[i];
            int curr_upper = surv_range_upper[i];
            pop_surv_age[i,j] = sum(pop_vals[curr_lower:curr_upper, j]);
            surv_group_props[curr_lower:curr_upper, j] = pop_vals[curr_lower:curr_upper, j] ./ pop_surv_age[i,j];
        }
    }

    //We mirror the operations again for the reporting groups (used for attack rates and proportion RSV naive)
    for(j in 1:n_seasons){
        for(i in 1:n_rep_groups){
            int curr_lower = rep_groups_lower[i];
            int curr_upper = rep_groups_upper[i];
            pop_rep_groups[i,j] = sum(pop_vals[curr_lower:curr_upper, j]);
            rep_group_props[curr_lower:curr_upper, j] = pop_vals[curr_lower:curr_upper, j] ./ pop_rep_groups[i,j];
        }
    }
}

parameters{
    vector<lower = 0> [n_beta] beta; 
    //We make use of the ? operator and zero-sized vectors to allow for optional parameters
    //We do not calibrate beta_child if we are using Model A (uniform FOI across all groups)
    //We only calibrate beta_ya if we are using Model D (FOI varies across all groups)
    vector<lower = 0> [model_code == 1 ? 0:n_beta] beta_child; 
    vector<lower = 0> [model_code == 4 ? n_beta:0] beta_ya;

    //We also calibrate probabilities of hospitalisation, reporting to surveillance, and the test
    //positive probability
    vector<lower = 0, upper = 1>[n_rho_hosp_group] rho_hosp_a;
    vector<lower = 0, upper = 1>[n_rho_surv_group] rho_surv_a; 
    vector<lower = 0, upper = 1>[n_rho_surv_group] rho_test_a; 

    //If r_code is 2 then we calibrate r, otherwise it is just set to 1. Note that this is not season dependent. 
    vector<lower = 0, upper = 1>[r_code == 2 ? 1:0] r; 
}

transformed parameters {
    //We start by declaring all the variables that we need for computations 
    //The number of time-dependent beta values = number of seasons + oldest age (number of ages - 2 to account for 0-6 and 7-12 months)
    //These values are the operational beta values per year. They are filled later on depending on 
    //the model code (ex: beta_child_t = beta_t for model A and we do not calibrate a separate beta for those 0-4 years)
    vector[n_seasons + (n_pop_age - 2)] beta_t; 
    vector[n_seasons + (n_pop_age - 2)] beta_child_t; 
    vector[n_seasons + (n_pop_age - 2)] beta_ya_t;

    //Probability of primary and post-primary infection and total
    matrix<lower = 0, upper = 1>[n_pop_age, n_seasons] Prob_pInf;
    matrix<lower = 0, upper = 1>[n_pop_age, n_seasons] Prob_sInf; 
    matrix<lower = 0, upper = 1>[n_pop_age, n_seasons] Prob_Inf_Total;

    //Probability of infection per age-group Pi (divided into primary, post-primary, and total)
    //Average probability of infection for hospitalisation age-groups (0-6m, 7-12m, 1, 2, ..., 81+)
    matrix<lower = 0, upper = 1>[n_hosp_age, n_seasons] Pi_y_g_prime; 
    matrix<lower = 0, upper = 1>[n_hosp_age, n_seasons] Pi_y_g_post; 
    matrix<lower = 0, upper = 1>[n_hosp_age, n_seasons] Pi_y_g;

    //Average probability of infection for surveillance age-groups (0-4, 5-14, 15-64, 65+)
    matrix<lower = 0, upper = 1>[n_surv_age, n_seasons] Pi_y_j_prime; 
    matrix<lower = 0, upper = 1>[n_surv_age, n_seasons] Pi_y_j_post; 
    matrix<lower = 0, upper = 1>[n_surv_age, n_seasons] Pi_y_j; 

    //Expected number of hospitalisations & RSV-attributable ILI (split into primary, post-primary, and total)
    matrix[n_hosp_age, n_seasons] Lambda_hosp_prime; 
    matrix[n_hosp_age, n_seasons] Lambda_hosp_post; 
    matrix[n_hosp_age, n_seasons] Lambda_hosp; 

    matrix[n_surv_age, n_seasons] Lambda_surv_prime; 
    matrix[n_surv_age, n_seasons] Lambda_surv_post; 
    matrix[n_surv_age, n_seasons] Lambda_surv; 

    //Age-distribution of hospitalisations and RSV-attributable ILI cases 
    matrix[n_hosp_age, n_seasons] Theta_hosp; 
    matrix[n_surv_age, n_seasons] Sigma_surv; 

    //Value used as the probabiity in the binomial likelihood
    //for test positivity rate: test-positive probability * probability of infection
    matrix<lower = 0, upper = 1>[n_surv_age, n_seasons] Prob_test_summ; 

    //Vectors for reporting rates
    vector[n_hosp_age] x_a_hosp; 
    vector[n_surv_age] x_a_surv;
    vector[n_surv_age] x_a_test;

    //Set the F_val based on F_code
    real F; 
    if(F_code == 1){
        //If F_code = 1 then F becomes a time-dependent value
        //Not yet implemented
    }else if(F_code == 2){
        F = 1; 
    }else if(F_code == 3){
        F = 0; 
    }

    //Set the r_val based on r_code 
    real r_val; 
    if(r_code == 1){
        //If r_code is 1 then r = 1
        r_val = 1; 
    }else if(r_code == 2){
        //If r_code is 2 then we calibrate
        r_val = r[1]; 
    }


    //------------------Start of Computations--------------------------

    //Set the time-dependent beta values, first the main values then all the ones prior
    //to the study period, which are assumed to be equal to the pre-COVID season.

    //The value of beta_a does not depend on the model type
    beta_t[t_start:rows(beta_t)] = beta;
    beta_t[1:t_start] = rep_vector(beta[1], t_start);

    //If we are running Model A (uniform FOI) then beta_child = beta_a
    //Otherwise, we calibrate beta_child separately
    if(model_code == 1){
        beta_child_t[t_start:rows(beta_child_t)] = beta; 
        beta_child_t[1:t_start] = rep_vector(beta[1], t_start);
    }else{
        beta_child_t[t_start:rows(beta_child_t)] = beta_child; 
        beta_child_t[1:t_start] = rep_vector(beta_child[1], t_start); 
    }

    //If we are running Model A or C then beta_ya = beta_a
    if(model_code == 1 || model_code == 3){
        beta_ya_t[t_start:rows(beta_ya_t)] = beta; 
        beta_ya_t[1:t_start] = rep_vector(beta[1], t_start);
    }else if(model_code == 2){
        //If we are running Model B then beta_ya = beta_child
        beta_ya_t[t_start:rows(beta_ya_t)] = beta_child; 
        beta_ya_t[1:t_start] = rep_vector(beta_child[1], t_start);
    }else{
        //If we are running Model D then we fit beta_ya separately
        beta_ya_t[t_start:rows(beta_ya_t)] = beta_ya; 
        beta_ya_t[1:t_start] = rep_vector(beta_ya[1], t_start);
    }

    //Fill up vectors for the reporting rates
    //We do this because of the different age-groups between data and reporting
    //hospitalisation data is from 0-81+, but rates are for 0 (0-6m), 0.5 (7-12m), 1, 2, 3+. 
    for(i in 1:n_rho_hosp_group){
        int curr_lower = rho_hosp_lower[i]; 
        int curr_upper = rho_hosp_upper[i]; 
        x_a_hosp[curr_lower:curr_upper] = rep_vector(rho_hosp_a[i], (curr_upper - curr_lower) + 1); 
    }

    //Fill up the surveillance reporting rate vectors (special case here where the age-groups)
    //of the data and reporting rates match (0-4, 5-14, 15-64, 65+)
    x_a_surv = rho_surv_a; 
    x_a_test = rho_test_a; 

    //Computing probabilities of primary and post-primary infection

    //Special case for 0-6 months (assuming 4 months maternally derived immunity)
    Prob_pInf[1, ] = to_row_vector((1 - exp(-beta_child_t[t_start:rows(beta_child_t)])) * (1 - (2*F)/3)); 
    Prob_sInf[1, ] = rep_row_vector(0, n_seasons);

    //Special case for 7-12 months (similar to those 0-6 months, but no maternal immunity)
    Prob_pInf[2, ] = to_row_vector(1 - exp(-beta_child_t[t_start:rows(beta_child_t)]));
    Prob_sInf[2, ] = rep_row_vector(0, n_seasons);

    //Main probabilities
    for(j in 1:n_seasons){
        int t = t_start + (j-1);
        
        //Starting at index i=3 since those 0-6m and 7-12m have already been computed
        for(i in 3:n_pop_age){
            //Retrieve the integer age value
            int a = pop_ages[i]; 
            
            //We use the beta_sum function to get the sum of the FOI values an individual experiences.
            //This goes from age 1 to (a-1) 
            real val_sum_beta = beta_sum(beta_t, beta_child_t, beta_ya_t, a, t, T, T2);

            //fin_beta then contains the beta value for the current season
            real fin_beta; 
            if(a <= T){
                fin_beta = beta_child_t[t]; 
            }else if((a > T) && (a <= T2)){
                fin_beta = beta_ya_t[t]; 
            }else{
                fin_beta = beta_t[t];
            }

            //Compute probability of escaping primary infection in first year of life for individual age a in time t (age 0 at t-a)
            //This always makes use of beta_child_t (since it's for 0 years-old)
            real Prob_pInf_bar = 1 - ((1 - exp(-beta_child_t[t - a])) * (1 - F/3));

            //We then compute the probability of primary and post-primary infection
            Prob_pInf[i,j] = Prob_pInf_bar * exp(-val_sum_beta) * (1 - exp(-fin_beta)); 
            Prob_sInf[i,j] = (1 - (Prob_pInf_bar * exp(-val_sum_beta))) * (1 - exp(-s * fin_beta)); 
        }
    }
    //Overall probability of infection is just the sum of the probabilities for 
    //primary and post-primary
    Prob_Inf_Total = Prob_pInf + Prob_sInf; 

    //After computing the probability of infection, we then
    //compute (population-weighted) average probabilities of infection per age-group from the data
    // (broken down into primary and post-primary infection probabilities)
    for(j in 1:n_seasons){
        for(i in 1:n_hosp_age){
            int curr_lower = hosp_range_lower[i]; 
            int curr_upper = hosp_range_upper[i];

            Pi_y_g_prime[i,j] = (sum(Prob_pInf[curr_lower:curr_upper, j] .* hosp_group_props[curr_lower:curr_upper, j]));
            Pi_y_g_post[i,j] = (sum(Prob_sInf[curr_lower:curr_upper, j] .* hosp_group_props[curr_lower:curr_upper, j]));
        }
    }
    
    for(j in 1:n_seasons){
        for(i in 1:n_surv_age){
            int curr_lower = surv_range_lower[i]; 
            int curr_upper = surv_range_upper[i];

            Pi_y_j_prime[i,j] = (sum(Prob_pInf[curr_lower:curr_upper, j] .* surv_group_props[curr_lower:curr_upper, j]));
            Pi_y_j_post[i,j] = (sum(Prob_sInf[curr_lower:curr_upper, j] .* surv_group_props[curr_lower:curr_upper, j]));
        }
    }

    Pi_y_g = Pi_y_g_prime + Pi_y_g_post; 
    Pi_y_j = Pi_y_j_prime + Pi_y_j_post;

    //From the average probabilities of infection, we apply hospitalisation / reporting rates
    //to compute the expected number of hospitalisations / RSV-attributable ILI cases
    //accounting for the value of r in the case of post-primary infections

    //Expected number of hospitalisations (we use the matrix form of the calculations)
    //Probability of hospitalisation * probability of infection * population
    Lambda_hosp_prime = diag_pre_multiply(x_a_hosp, (Pi_y_g_prime .* pop_hosp_age));
    Lambda_hosp_post = diag_pre_multiply((x_a_hosp * r_val), (Pi_y_g_post .* pop_hosp_age));
    Lambda_hosp = Lambda_hosp_prime + Lambda_hosp_post; 

    //Expected number of RSV-attributable-ILI cases (we use the matrix form of the calculations)
    //Probability of reporting to surveillance * probability of infection * population
    Lambda_surv_prime = diag_pre_multiply(x_a_surv, (Pi_y_j_prime .* pop_surv_age)); 
    Lambda_surv_post = diag_pre_multiply((x_a_surv * r_val), (Pi_y_j_post .* pop_surv_age)); 
    Lambda_surv = Lambda_surv_prime + Lambda_surv_post; 

    //Modelling proportion testing positive for RSV
    //Test positive probability * probability of infection 
    Prob_test_summ = diag_pre_multiply(x_a_test, Pi_y_j);

    //Compute season totals for model-derived hospitalisations and RSV-attributable-ILI cases
    vector[n_seasons] hosp_totals_model;
    vector[n_seasons] surv_totals_model;
    for(i in 1:n_seasons){
        hosp_totals_model[i] = sum(Lambda_hosp[,i]); 
        surv_totals_model[i] = sum(Lambda_surv[,i]); 
    }

    //Compute age-distributions of hospitalisations and RSV-attributable-ILI
    //We then use these values for the multinomials in the likelihood
    for(j in 1:n_seasons){
        for(i in 1:n_hosp_age){
            Theta_hosp[i,j] = Lambda_hosp[i,j] / hosp_totals_model[j]; 
        }
        for(i in 1:n_surv_age){
            Sigma_surv[i,j] = Lambda_surv[i,j] / surv_totals_model[j]; 
        }
    }
}

model{
    //Setting priors for beta values depends on the model variant
    //being run 
    for(i in 1:n_beta){
        //We set a special prior for the 2020-2021 season due to expected COVID impact (as shown in other countries)
        if(i == n_beta-1){
            beta[i] ~ exponential(1000);
            if(model_code == 4){
                //If Model D then we calibrate beta_ya
                beta_ya[i] ~ exponential(1000); 
            }
            if(model_code != 1){
                //We calibrate beta_child in all cases except Model A
                beta_child[i] ~ exponential(1000);
            }
        }else{
            beta[i] ~ normal(0.673977,0.2); 
            if(model_code == 4){
                //If Model D then we calibrate beta_ya
                beta_ya[i] ~ normal(0.673977,0.2); 
            }
            if(model_code != 1){
                //We calibrate beta_child in all cases except Model A
                beta_child[i] ~ normal(0.673977,0.2); 
            }

        }
    }

    for(i in 1:n_rho_hosp_group){
        //For priors here, we apply calibrated probabilities
        //of hospitalisation from van Boven et al. 
        rho_hosp_a[i] ~ normal(rho_hosp_priors[i], 0.2); 
    }

    for(i in 1:n_rho_surv_group){
        //For priors here, we apply calibrated probabilities
        //of GP appointment from van Boven et al. (same priors for
        //reporting to surveillance and test positive probability)
        rho_surv_a[i] ~ normal(rho_surv_priors[i], 0.2);
        rho_test_a[i] ~ normal(rho_surv_priors[i], 0.2);
    }

    for(i in 1:n_seasons){
        //can modify here to exclude the multinomial for 2020-2021
        //hosp_vals[, i] ~ multinomial(Theta_hosp[, i]);

        //Version where we do not fit the multinomial to 2020-2021
         if(i != (n_seasons-1)){
            hosp_vals[, i] ~ multinomial(Theta_hosp[, i]); 
        }
        hosp_totals_data[i] ~ poisson(hosp_totals_model[i]);  
    }

    for(i in 1:n_seasons){
        //We do not fit the multinomial to 2020-2021
        if(i != (n_seasons-1)){
            surv_vals[, i] ~ multinomial(Sigma_surv[, i]); 
        }
        surv_totals_data[i] ~ poisson(surv_totals_model[i]); 
    }

    //Fit to RSV test positivity rate 
    for(j in 1:n_seasons){
        for(i in 1:n_surv_age){
            pos_vals[i,j] ~ binomial(test_vals[i,j], Prob_test_summ[i,j]); 
        }
    }

    //Set a prior for r depending on the r_code
    if(r_code == 2){
        r ~ uniform(0,1);
    }
}

generated quantities {
    //We need a vector log lik with dimensions equal to the number of data points
    //we are fitting to. This pointwise log-likelihood will be used later on to compute the DIC. 

    //We compute pointwise log-likelihoods separately for hospitalisations and surveillance
    //then combine them together in the end.

    //Likelihood for hospitalisation data
    vector[n_seasons - 1] log_lik_hosp_mult; 
    vector[n_seasons] log_lik_hosp_pois; 
    for(i in 1:n_seasons){
        //We skip the 2020-2021 season for the multinomial
        if(i != n_seasons - 1){
            if(i == n_seasons){
                //Put value for last season in correct index
                log_lik_hosp_mult[i-1] = multinomial_lpmf(hosp_vals[, i] | Theta_hosp[, i]);
            }else{
                log_lik_hosp_mult[i] = multinomial_lpmf(hosp_vals[, i] | Theta_hosp[, i]); 
            }  
        }
        //log_lik_hosp_mult[i] = multinomial_lpmf(hosp_vals[,i] | Theta_hosp[,i]);
        log_lik_hosp_pois[i] = poisson_lpmf(hosp_totals_data[i] | hosp_totals_model[i]); 
    }

    //Likelihood for surveillance data
    vector[n_seasons - 1] log_lik_surv_mult; 
    vector[n_seasons] log_lik_surv_pois; 
    for(i in 1:n_seasons){
        //We skip the 2020-2021 season for the multinomial
        if(i != n_seasons - 1){
            if(i == n_seasons){
                //Put value for last season in correct index
                log_lik_surv_mult[i-1] = multinomial_lpmf(surv_vals[, i] | Sigma_surv[, i]);
            }else{
                log_lik_surv_mult[i] = multinomial_lpmf(surv_vals[, i] | Sigma_surv[, i]); 
            }  
        }
        log_lik_surv_pois[i] = poisson_lpmf(surv_totals_data[i] | surv_totals_model[i]);
    }

    //Likelihood for testing
    matrix[n_surv_age, n_seasons] log_lik_test; 
    for(j in 1:n_seasons){
        for(i in 1:n_surv_age){
            log_lik_test[i,j] = binomial_lpmf(pos_vals[i,j] | test_vals[i,j], Prob_test_summ[i,j]); 
        }
    }

    //We then combine all the log-likelihoods into a single vector, which will then be used in 
    //the DIC computations. 
    vector[(2*n_seasons - 1) + (2*n_seasons - 1) + (n_seasons * n_surv_age)] log_lik = append_row(append_row(append_row(log_lik_hosp_mult, log_lik_hosp_pois), append_row(log_lik_surv_mult, log_lik_surv_pois)), to_vector(log_lik_test'));


    //We also compute the proportion non-naive to RSV before and after each season
    matrix[n_pop_age, n_seasons] prop_non_naive_after = rep_matrix(0, n_pop_age, n_seasons);
    matrix[n_pop_age, n_seasons] prop_non_naive_before = rep_matrix(0, n_pop_age, n_seasons);

    //We compute the proportion non-naive to RSV but we have to account for those who had maternally derived immunity.
    //These individuals were immune to RSV but they are counted as RSV naive. 

    //Special case for 0-year-olds (we compute the proportion who are infected after each season)
    prop_non_naive_before[1, ] = rep_row_vector(0, n_seasons);
    prop_non_naive_after[1, ] = to_row_vector((1 - exp(-beta_child_t[t_start:rows(beta_child_t)])) * (1 - (2*F)/3));

    //Special case for 7-12 months old
    prop_non_naive_before[2, ] = rep_row_vector(0, n_seasons);
    prop_non_naive_after[2, ] = to_row_vector((1 - exp(-beta_child_t[t_start:rows(beta_child_t)])));

    for(j in 1:n_seasons){
        int t = t_start + (j - 1); 
        for(i in 3:n_pop_age){
            int a = pop_ages[i]; 
            
            real val_sum_beta = beta_sum(beta_t, beta_child_t, beta_ya_t, a, t, T, T2); 
            
            //We retrieve the value of beta for the current season, which is then used for
            //computing the "after" state
            real fin_beta; 
            if(a <= T){
                fin_beta = beta_child_t[t]; 
            }else if((a > T) && (a <= T2)){
                fin_beta = beta_ya_t[t]; 
            }else{
                fin_beta = beta_t[t]; 
            }

            real Prob_pInf_bar = 1 - ((1 - exp(-beta_child_t[t - a])) * (1 - F/3));
            prop_non_naive_before[i,j] = (1 - Prob_pInf_bar * exp(-val_sum_beta)); 
            prop_non_naive_after[i,j] = (1 - Prob_pInf_bar * exp(-(val_sum_beta + fin_beta))); 
        }
    }
    //To get the population non-naive, we multiply proportions by the size of the population
    matrix[n_pop_age, n_seasons] pop_non_naive_after = prop_non_naive_after .* pop_vals; 
    matrix[n_pop_age, n_seasons] pop_non_naive_before = prop_non_naive_before .* pop_vals;

    //Total proportion of the population non-naive to RSV before and after each season
    vector[n_seasons] total_non_naive_after = ((pop_non_naive_after' * rep_vector(1, n_pop_age)) ./ (pop_vals' * rep_vector(1, n_pop_age)));
    vector[n_seasons] total_non_naive_before = ((pop_non_naive_before' * rep_vector(1, n_pop_age)) ./ (pop_vals' * rep_vector(1, n_pop_age)));

    //We break down the FOI values into transmission rates. This is done in two ways: 
    //(1) using the total number of RSV-attributable ILI cases 
    //and (2) using the total model estimated number of infections. 

    //First, we compute the total number of RSV-attributable ILI cases from the data. This is just a column sum of the surv_vals matrix
    //Second, we do the same but this time compute the number of model estimated RSV infections by multiplying
    //the probability of infection by the population 
    matrix[n_pop_age, n_seasons] num_infections = Prob_Inf_Total .* pop_vals;
    vector[n_seasons] total_rsv_infections;
    vector[n_seasons] total_rsv_att_ili; 
    for(i in 1:n_seasons){
        total_rsv_infections[i] = sum(num_infections[,i]); 
        total_rsv_att_ili[i] = sum(surv_vals[,i]); 
    }

    //Next we divide our beta estimates with these values. Note that there are 0
    //RSV-attributable ILI cases in the 2020-2021 season. To prevent a division by 0,
    //we skip this season.

    //Note that the output here depends on the beta values thus they are linked to the model code
    vector[model_code == 1 ? 0:n_seasons] trans_rate_inf_c; 
    vector[model_code == 4 ? n_seasons:0] trans_rate_inf_ya;
    vector[n_seasons] trans_rate_inf_a; 


    vector[model_code == 1 ? 0:(n_seasons-1)] trans_rate_ili_c; 
    vector[model_code == 4 ? (n_seasons-1):0] trans_rate_ili_ya;
    vector[n_seasons - 1] trans_rate_ili_a; 

    for(i in 1:n_seasons){
        //Compute the transmission rate for 0-4, 5-14, and 15+ years using model estimated infections 

        if(model_code != 1){
            trans_rate_inf_c[i] = beta_child[i] / total_rsv_infections[i]; 
        }
        if(model_code == 4){
            trans_rate_inf_ya[i] = beta_ya[i] / total_rsv_infections[i]; 
        }
        trans_rate_inf_a[i] = beta[i] / total_rsv_infections[i]; 

        //Compute the transmission rate for 0-4, 5-14, and 15+ years using RSV-attributable ILI
        if(i == (n_seasons-1)){
            continue; 
        }else{
            int ind = i; 

            if(i == n_seasons){
                ind = i - 1; 
            }

            if(model_code != 1){
                trans_rate_ili_c[ind] = beta_child[i] / total_rsv_att_ili[i]; 
            }
            if(model_code == 4){
                trans_rate_ili_ya[ind] = beta_ya[i] / total_rsv_att_ili[i]; 
            }
            trans_rate_ili_a[ind] = beta[i] / total_rsv_att_ili[i]; 
        }
    }

    //We generate population weighted average attack rate and proportion RSV naive values for the reporting groups
    //Attack rates are based on the probability of infection (Prob_Inf_Total)
    //While proportion RSV naive are based on prop_non_naive_before and prop_non_naive_after
    matrix[n_rep_groups, n_seasons] attack_rates_rep_groups;
    matrix[n_rep_groups, n_seasons] rsv_naive_rep_groups_before; 
    matrix[n_rep_groups, n_seasons] rsv_naive_rep_groups_after;

    for(j in 1:n_seasons){
        for(i in 1:n_rep_groups){
            int curr_lower = rep_groups_lower[i];
            int curr_upper = rep_groups_upper[i];

            //[n_pop_age, n_seasons] Prob_Inf_Total;
            attack_rates_rep_groups[i, j] = sum(Prob_Inf_Total[curr_lower:curr_upper, j] .* rep_group_props[curr_lower:curr_upper, j]); 
            rsv_naive_rep_groups_before[i, j] = sum(prop_non_naive_before[curr_lower:curr_upper, j] .* rep_group_props[curr_lower:curr_upper, j]);
            rsv_naive_rep_groups_after[i, j] = sum(prop_non_naive_after[curr_lower:curr_upper, j] .* rep_group_props[curr_lower:curr_upper, j]);
        }
    }

}
