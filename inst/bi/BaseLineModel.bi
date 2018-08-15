/**
 * Baseline TB and BCG vaccination model
 */
model Baseline {
  
  // Model dimensions
  const e_bcg = 2 // 0 = unvaccinated, 1 = vaccinated
  const e_age = 15 // 0,..14= 5 year age groups (i.e 0-4), and 15 = 70-89

  dim bcg(e_bcg)
  dim age(e_age)
  dim age2(e_age)
  
  //Age at vaccination
  const vac_scheme = 0 //0 = vaccination at school age, 1 = vaccination at birth, 2 = no vaccination.
  const e_d_of_p = 6 // Duration of protection (must be at least 1)
  
  dim d_of_p(e_d_of_p)

  //Control model
  const const_pop = 0 //Set to 1 for constant population (i.e births == deaths)
  const no_age = 0 //Set to 1 to turn off ageing
  const no_disease = 0 //Set to 1 to prevent disease from being initialised / importation
  // Time dimensions
  const ScaleTime = 1 / 12 // Scale model over a year 
  //const ScaleTime = 1 // Scale model over a month
  
  //Initialise model
  const init_pop = 37359045 //Estimated intial population - http://www.visionofbritain.org.uk/census/table/EW1931COU1_M3
  const init_P_cases = 49798 // TB cases in England (and Wales)
  const init_E_cases = 16084 // https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/554455/TB_case_notifications_1913_to_2015.pdf

  //Disease model parameters
  
  // Effective contact rate
  param c_eff
  
  // Historic effective contact rate
  param c_hist
  
  // Protection from infection due to prior latent infection
  param delta
  
  // Transition from high risk latent disease to active disease
  param epsilon_h_0_4 //Age specific parameters
  param epsilon_h_5_14
  param epsilon_h_15_89
  param epsilon_h[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Transition to low risk latent disease from high risk latent disease
  param kappa_0_4 //Age specific parameters
  param kappa_5_14
  param kappa_15_89
  param kappa[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Transition from low risk latent disease to active disease
  param epsilon_l_0_4 //Age specific parameters
  param epsilon_l_5_14
  param epsilon_l_15_89
  param epsilon_l[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Rate of succesful treatment completion
  param phi_0_14 //Age specific parameters
  param phi_15_59
  param phi_60_89
  param phi[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Proportion of cases that have pulmonary TB
  param Upsilon_0_14 //Age specific parameters
  param Upsilon_15_59
  param Upsilon_60_89
  param Upsilon[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Proportion of cases that have pulmonary smear postive TB
  param rho_0_14 //Age specific parameters
  param rho_15_59
  param rho_60_89
  param rho[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Rate of starting treatment - pulmonary/extra-pulmonary
  // Pulmonary
  param nu_p_0_14 //Age specific parameters
  param nu_p_15_89
  param nu_p[age](has_output = 0, has_input = 0) // Dummy model parameter
  // Extra-pulmonary
  param nu_e_0_14 //Age specific parameters
  param nu_e_15_89
  param nu_e[age](has_output = 0, has_input = 0) // Dummy model parameter
  // Average rate of starting treatment   
  param avg_nu_p[age](has_output = 0, has_input = 0)
  
  // Rate loss to follow up - pulmonary/extra-pulmonary
  // Pulmonary
  param zeta_p_0_14 //Age specific parameters
  param zeta_p_15_59
  param zeta_p_60_89
  param zeta_p[age](has_output = 0, has_input = 0) // Dummy model parameter
  // Extra-pulmonary
  param zeta_e_0_14 //Age specific parameters
  param zeta_e_15_59
  param zeta_e_60_89
  param zeta_e[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Rate of TB mortality
  param mu_p_0_14 //Age specific parameters
  param mu_p_15_59
  param mu_p_60_89
  param mu_p[age](has_output = 0, has_input = 0) // Dummy model parameter
  // Extra-pulmonary
  param mu_e_0_14 //Age specific parameters
  param mu_e_15_59
  param mu_e_60_89
  param mu_e[age](has_output = 0, has_input = 0) // Dummy model parameter
        
  
  //BCG vaccination parameters
  
  // Age specific protection from infection conferred by BCG vaccination
  param chi_init
  // Protection from active disease due to BCG vaccination
  param alpha_t[d_of_p]
  
  //Demographic model parameters
  //Ageing
  param theta[age](has_output = 0, has_input = 0)
  
  //Noise parameters
  noise CNoise[age, age2](has_output = 0, has_input = 0) // Sampled contact rate
  noise coverage(has_output = 0, has_input = 0) //Estimated coverage of BCG vaccination in the UK born
    
  // Time varying parameter states
  state CSample[age, age2](has_output = 0, has_input = 0) // Sampled contact rate (symmetric)
  state SelfContacts[age](has_output = 0, has_input = 0) //Contacts within age group
  state TotalContacts[age](has_output = 0, has_input = 0) // Average number of contacts (across age groups)
  state beta[age](has_output = 0, has_input = 0) //Probability of transmission
  state foi[age] // force of infection
  
  //Calculation parameters
  param I_age[age](has_output = 0, has_input = 0)
  param I_bcg[bcg](has_output = 0, has_input = 0)
  
  // Time varing parameter states
  //Demographic
  state births //All births
  state mu_all[age] // All cause natural mortality
  state mu[age](has_output = 0, has_input = 0) //All cause mortality excluding TB
  
  //Vaccination
  state age_at_vac(has_output = 0, has_input = 0) //Age at vaccination
  state chi[age](has_output = 0, has_input = 0) //Protection from initial infection due to BCG vaccination
  state alpha[age](has_output = 0, has_input = 0) //Protection from active disease due to BCG vaccination
  state gamma[age](has_output = 0, has_input = 0) //Coverage of the vaccination program by age
  
  //Population states
  state S[bcg, age] // susceptible
  state H[bcg, age] // high risk latent
  state L[bcg, age] // low risk latent
  state P[bcg, age] // pulmonary TB
  state E[bcg, age] // extra-pulmonary TB only
  state T_P[bcg, age] // pulmonary TB on treatment
  state T_E[bcg, age] // extra-pulmonary TB on treatment
  state N[bcg, age] // Overall population
  state NAge[age](has_output = 0, has_input = 0) //Age summed population
  state NSum[age](has_input = 0, has_output = 0) // Sum of population (vector but repeating values)
  state death_sum[age](has_output = 0, has_input = 0) //Used to estimate deaths
      
  //Accumalator states
  //state PulCases[bcg, age] // monthly pulmonary cases starting treatment
  //state EPulCases[bcg, age] // monthly extra-pulmonary cases starting treatment
  state YearlyPulCases[bcg, age] // yearly pulmonary cases starting treatment
  state YearlyEPulCases[bcg, age] // yearly extra-pulmonary cases starting treatment
  //state PulDeaths[bcg, age] // Pulmonary TB deaths
  //state EPulDeaths[bcg, age] // Extra-pulmonary TB deaths
  state YearlyPulDeaths[bcg, age] // Pulmonary TB deaths (yearly)
  state YearlyEPulDeaths[bcg, age] // Extra-pulmonary TB deaths (yearly)
  
  // Reporting states
  state YearlyAgeCases[age]
  state YearlyCasesCumSum[age](has_input = 0, has_output = 0)
  state YearlyCases
  
  //Noise variables
  noise births_sample(has_output = 0, has_input = 0) //Sampled noisy births
  noise mu_all_sample[age](has_output = 0, has_input = 0) //Sampled noisy deaths
     
  //Input
  input births_input //Births (time varying)
  input pop_dist[age] //Population distribution (average from 2000 to 2015).
  input exp_life_span[age] //Expected life span (time varying)
  input polymod[age, age2] //Polymod contact matrix
  input polymod_sd[age, age2] //Polymod SD contact matrix
  
  //Observations
  obs YearlyInc // Yearly overall incidence
  obs YearlyAgeInc[age] // Yearly incidence by age group
  
      sub parameter {
        
        // Parameter scales
        inline dscale = 12 / 365.25 * ScaleTime 
        inline mscale = 1 * ScaleTime
        inline yscale = 12 * ScaleTime
        
        // Priors for BCG vaccination
        //Protection from infection at vaccination
        chi_init ~ truncated_gaussian(mean = 0.185, std = 0.0536, lower = 0, upper = 1)
        //Protection from active TB (decaying from time since vaccination)
        alpha_t[0] ~ gaussian(mean = -1.86, std = 0.22)
        alpha_t[1] ~ gaussian(mean = -1.19, std = 0.24)
        alpha_t[2] ~ gaussian(mean = -0.84, std = 0.22)
        alpha_t[3] ~ gaussian(mean = -0.84, std = 0.2)
        alpha_t[4] ~ gaussian(mean = -0.28, std = 0.19)
        alpha_t[5] ~ gaussian(mean = -0.23, std = 0.29)
        alpha_t[d_of_p] <- 1 - exp(alpha_t[d_of_p]) // Previous log transformed

        //Disease priors
        c_eff ~ uniform(0, 5)
        c_hist ~ uniform(10, 15)
        delta ~ truncated_gaussian(mean = 0.78, std = 0.0408, lower = 0, upper = 1)
        
        // Transition from high risk latent to active TB
        epsilon_h_0_4 ~ truncated_gaussian(mean = 0.00695, std =  0.0013, lower = 0)
        epsilon_h_5_14 ~ truncated_gaussian(mean = 0.0028, std =  0.000561, lower = 0)
        epsilon_h_15_89 ~ truncated_gaussian(mean = 0.000335, std =  0.0000893, lower = 0)
        epsilon_h_0_4 <- dscale / (epsilon_h_0_4)
        epsilon_h_5_14 <- dscale / (epsilon_h_5_14)
        epsilon_h_15_89 <- dscale / (epsilon_h_15_89)
        
        epsilon_h[0] <- epsilon_h_0_4
        epsilon_h[age = 1:2] <-  epsilon_h_5_14
        epsilon_h[age = 3:(e_age - 1)] <- epsilon_h_15_89
        epsilon_h <- 1 / epsilon_h
        
        // Rate of transition from high risk to low risk latents
        kappa_0_4 ~ truncated_gaussian(mean = 0.0133, std = 0.00242, lower = 0)
        kappa_5_14 ~ truncated_gaussian(mean = 0.012, std = 0.00207, lower = 0)
        kappa_15_89 ~ truncated_gaussian(mean = 0.00725, std = 0.00191, lower = 0)
        kappa_0_4 <- dscale / (kappa_0_4)
        kappa_5_14 <- dscale / (kappa_5_14)
        kappa_15_89 <- dscale / (kappa_15_89)
        kappa[0] <- kappa_0_4
        kappa[age = 1:2] <- kappa_5_14
        kappa[age = 3:(e_age - 1)] <- kappa_15_89
        kappa <- 1 / kappa
        
        // Rate of transition for low risk latent to active TB
        epsilon_l_0_4 ~ truncated_gaussian(mean = 0.000008, std =  0.00000408, lower = 0)
        epsilon_l_5_14 ~ truncated_gaussian(mean = 0.00000984, std =  0.00000467, lower = 0)
        epsilon_l_15_89 ~ truncated_gaussian(mean = 0.00000595, std =  0.00000207, lower = 0)
        epsilon_l_0_4 <- dscale / (epsilon_l_0_4)
        epsilon_l_5_14 <- dscale / (epsilon_l_5_14)
        epsilon_l_15_89 <- dscale / (epsilon_l_15_89)
        epsilon_l[0] <- epsilon_l_0_4
        epsilon_l[age = 1:2] <-  epsilon_l_5_14
        epsilon_l[age = 3:(e_age - 1)] <- epsilon_l_15_89
        epsilon_l <- 1 / epsilon_l
        
        // Rate of successful treatment
        phi_0_14 ~ gamma(shape = 9.86, scale = 0.061)
        phi_15_59 ~ gamma(shape = 7.73, scale = 0.0837)
        phi_60_89 ~ gamma(shape = 8.46, scale = 0.0734)
        phi_0_14 <- phi_0_14 * yscale
        phi_15_59 <-  phi_15_59 * yscale
        phi_60_89 <- phi_60_89 * yscale
        phi[age = 0:2] <- phi_0_14
        phi[age = 3:11] <-  phi_15_59
        phi[age = 12:(e_age - 1)] <- phi_60_89
        phi <- 1 / phi

        
        // Proportion of TB cases with pulmonary TB
        Upsilon_0_14  ~ truncated_gaussian(mean = 0.629, std = 0.0101, lower = 0, upper = 1)
        Upsilon_15_59 ~ truncated_gaussian(mean = 0.706, std = 0.00411, lower = 0, upper = 1)
        Upsilon_60_89 ~ truncated_gaussian(mean = 0.75, std = 0.00569, lower = 0, upper = 1)
        Upsilon[age = 0:2] <- Upsilon_0_14
        Upsilon[age = 3:11] <-  Upsilon_15_59
        Upsilon[age = 12:(e_age - 1)] <- Upsilon_60_89
          
        // Propotion of pulmonary TB cases that are smear positive
        rho_0_14  ~ truncated_gaussian(mean = 0.302, std = 0.0189, lower = 0, upper = 1)
        rho_15_59 ~ truncated_gaussian(mean = 0.652, std = 0.00518, lower = 0, upper = 1)
        rho_60_89 ~ truncated_gaussian(mean = 0.536, std = 0.00845, lower = 0, upper = 1)
        rho[age = 0:2] <- rho_0_14
        rho[age = 3:11] <-  rho_15_59
        rho[age = 12:(e_age - 1)] <- rho_60_89
        
        // Rate of starting treatment - pulmonary/extra-pulmonary
        // Pulmonary
        nu_p_0_14  ~ gamma(shape = 0.878, scale = 0.206)
        nu_p_15_89 ~ gamma(shape = 1.1, scale = 0.3)
        nu_p_0_14  <- nu_p_0_14 * yscale
        nu_p_15_89 <- nu_p_15_89 * yscale
        nu_p[age = 0:2] <- nu_p_0_14
        nu_p[age = 3:(e_age - 1)] <- nu_p_15_89
        nu_p <- 1 / nu_p
        // Extra-pulmonary
        nu_e_0_14 ~ gamma(shape = 0.686, scale = 0.446)
        nu_e_15_89 ~ gamma(shape = 0.897, scale = 0.536)
        nu_e_0_14  <- nu_e_0_14 * yscale
        nu_e_15_89 <- nu_e_15_89 * yscale
        nu_e[age = 0:2] <- nu_e_0_14
        nu_e[age = 3:(e_age - 1)] <- nu_e_15_89
        nu_e <- 1 / nu_e
  
        // Rate loss to follow up - pulmonary/extra-pulmonary
        // Pulmonary
        zeta_p_0_14  ~ truncated_gaussian(mean = 0.0107, std = 0.0225, lower = 0)
        zeta_p_15_59 ~ truncated_gaussian(mean = 0.0356, std = 0.00977, lower = 0)
        zeta_p_60_89 ~ truncated_gaussian(mean = 0.00847, std = 0.0147, lower = 0)
        zeta_p_0_14  <- yscale / (zeta_p_0_14)
        zeta_p_15_59 <- yscale / (zeta_p_15_59)
        zeta_p_60_89 <- yscale / (zeta_p_60_89)
        zeta_p[age = 0:2] <- zeta_p_0_14
        zeta_p[age = 3:11] <-  zeta_p_15_59
        zeta_p[age = 12:(e_age - 1)] <- zeta_p_60_89
        zeta_p <- 1 / zeta_p
        // Extra-pulmonary
        zeta_e_0_14  ~ truncated_gaussian(mean = 0.00807, std = 0.0298, lower = 0)
        zeta_e_15_59 ~ truncated_gaussian(mean = 0.0281, std = 0.0151, lower = 0)
        zeta_e_60_89 ~ truncated_gaussian(mean = 0.00754, std = 0.025, lower = 0)
        zeta_e_0_14  <- yscale / (zeta_e_0_14)
        zeta_e_15_59 <- yscale / (zeta_e_15_59)
        zeta_e_60_89 <- yscale / (zeta_e_60_89)
        zeta_e[age = 0:2] <- zeta_e_0_14
        zeta_e[age = 3:11] <-  zeta_e_15_59
        zeta_e[age = 12:(e_age - 1)] <- zeta_e_60_89
        zeta_e <- 1 / zeta_e
            
        // Rate of TB mortality
        mu_p_0_14  ~ truncated_gaussian(mean = 0.00413, std = 0.0227, lower = 0)
        mu_p_15_59 ~ truncated_gaussian(mean = 0.0225, std = 0.0101, lower = 0)
        mu_p_60_89 ~ truncated_gaussian(mean = 0.112, std = 0.0151, lower = 0)
        mu_p_0_14  <- yscale / (mu_p_0_14)
        mu_p_15_59 <- yscale / (mu_p_15_59)
        mu_p_60_89 <- yscale / (mu_p_60_89)
        mu_p[age = 0:2] <- mu_p_0_14
        mu_p[age = 3:11] <-  mu_p_15_59
        mu_p[age = 12:(e_age - 1)] <- mu_p_60_89
        mu_p <- 1 / mu_p
        // Extra-pulmonary
        mu_e_0_14  ~ truncated_gaussian(mean = 0.00363, std = 0.0301, lower = 0)
        mu_e_15_59 ~ truncated_gaussian(mean = 0.00516, std = 0.0156, lower = 0)
        mu_e_60_89 ~ truncated_gaussian(mean = 0.0438, std = 0.0258, lower = 0)
        mu_e_0_14  <- yscale / (mu_e_0_14)
        mu_e_15_59 <- yscale / (mu_e_15_59)
        mu_e_60_89 <- yscale / (mu_e_60_89)
        mu_e[age = 0:2] <- mu_e_0_14
        mu_e[age = 3:11] <-  mu_e_15_59
        mu_e[age = 12:(e_age - 1)] <- mu_e_60_89
        mu_e <- 1 / mu_e
        
        //Calculation parameters
        I_age <- 1
        I_bcg <- 1
        
        // Avg period infectious for pulmonary cases
        avg_nu_p <- inclusive_scan(nu_p)
        avg_nu_p[age] <- avg_nu_p[e_age - 1] / e_age
        
        //Demographic model parameters
        theta[age=0:14] <- (no_age == 0 ? 1 / (5 * yscale) : 0)
        theta[15] <- (no_age == 0 ? 1 / (20 * yscale) : 0)
      }
    
    sub initial {
      S[0, age] ~  truncated_gaussian(mean = init_pop * pop_dist[age], std = 0.05 * init_pop * pop_dist[age], lower = 0) // susceptible
      S[1, age] <- 0 // BCG vaccinated susceptibles
      H[bcg, age] <- 0 // high risk latent
      L[bcg, age] <- 0 // low risk latent
      P[0, age] <- (no_disease == 0 ? init_P_cases * pop_dist[age]: 0)  // pulmonary TB
      P[1, age] <- 0 //BCG vaccinated pulmonary TB
      E[0, age] <- (no_disease == 0 ? init_E_cases * pop_dist[age] : 0) // extra-pulmonary TB only
      E[1, age] <- 0 // BCG extra-pulmonary TB only
      T_P[bcg, age] <- 0// pulmonary TB on treatment
      T_E[bcg, age] <- 0 // extra-pulmonary TB on treatment
    }
    
    sub transition {
      
      // Reset accumalator variables
      //PulCases[bcg, age] <- 0
      //EPulCases[bcg, age] <- 0
      //PulDeaths[bcg, age] <- 0
      //EPulDeaths[bcg, age] <- 0
      inline yr_reset = 12 * ScaleTime
      YearlyPulCases[bcg, age] <- (t_now % 1 < yr_reset ? 0 : YearlyPulCases[bcg, age])
      YearlyEPulCases[bcg, age] <- (t_now % 1 < yr_reset ? 0 : YearlyEPulCases[bcg, age])
      YearlyPulDeaths[bcg, age] <- (t_now % 1 < yr_reset ? 0 : YearlyPulDeaths[bcg, age])
      YearlyEPulDeaths[bcg, age] <- (t_now % 1 < yr_reset ? 0 : YearlyEPulDeaths[bcg, age])
      //Apply BCG vaccination to correct populations
      inline policy_change = 74 * 12 * ScaleTime // Assume policy switch occurred in 2005
      //Set up age at vaccination
      age_at_vac <- (t_now >= policy_change ? (vac_scheme == 0 ? 3 : (vac_scheme == 1 ? 0 : -1)) : 3)
      
      // Apply linear transform for protection against initial infection
      chi[age] <- (age_at_vac < 0 ? 0 : (age_at_vac > age ? 0 : (age >= (age_at_vac + e_d_of_p) ? 0 : alpha_t[age - age_at_vac] * chi_init / alpha_t[0])))
      //Back calculate protection from latent disease based on initial protection and overall protection
      alpha[age] <- (age_at_vac < 0 ? 0 : (age_at_vac > age ? 0 : (age >= (age_at_vac + e_d_of_p) ? 0 : (alpha_t[age - age_at_vac] - chi[age])/ (1 - chi[age]))))
      //Apply coverage of vac program to correct population
      coverage ~ truncated_gaussian(mean = 0.8, std = 0.05, lower = 0, upper = 1)
      // Set vaccination to begin in 1953
      inline vac_start = 22 * 12 * ScaleTime
      gamma[age] <- (age_at_vac == age ? (t_now > vac_start ? coverage : 0) : 0)
      
      //Contact rate
      CNoise[age, age2] ~  truncated_gaussian(mean = polymod[age, age2], std = polymod_sd[age, age2], lower = 0)
      CSample <- CNoise
      CSample[age, age2] <- CSample[age2, age]
      SelfContacts[age] <- CSample[age, age]
      TotalContacts <- CSample * I_age + SelfContacts
      TotalContacts <- inclusive_scan(TotalContacts)
      TotalContacts[age] <- TotalContacts[e_age - 1] / 2
      
      // Population
      N <- S + H + L + P + E + T_P + T_E
      NAge[age] <- N[0, age] + N[1, age]
      NSum <- inclusive_scan(NAge)
      NSum[age] <- NSum[e_age - 1]
      
      // Estimate force of infection - start with probability of transmission
      inline curr_start = 59 * 12 * ScaleTime // Modern day is 1990 with a baseline date of 1931
      beta <- avg_nu_p * (c_eff + c_hist * ((t_now > curr_start  ? 0 : (curr_start - t_now) / curr_start))) 
      beta <- beta ./ TotalContacts
      
      //Now build force of infection
      foi <- transpose(P) * I_bcg
      foi[age] <- rho[age] * foi[age]
      foi <- CSample * foi
      foi[age] <- beta[age] * foi[age] / NSum[age]
      
      //All-cause mortality excluding TB
      mu_all_sample[age] ~ truncated_gaussian(mean = exp_life_span[age], std = exp_life_span[age]*0.05, lower = 0)
      mu_all[age] <- mu_all_sample[age]
      mu[age] <- 1 / mu_all[age] - (mu_p[age] * (P[0, age] + P[1, age] + T_P[0, age] + T_P[1, age]) + mu_e[age] * (E[0, age] + E[1, age] + T_E[0, age] + T_E[1, age])) / NAge[age]
      // All used to fix births to deaths for testing 
      death_sum[age] <-  NAge[age] / mu_all[age]
      death_sum <- inclusive_scan(death_sum)
      births_sample ~ gaussian(mean = births_input, std = 0.05 * births_input)
      births <- (const_pop == 0 ? births_sample : death_sum[e_age - 1] + theta[e_age - 1] * (N[0, e_age - 1] +  N[1, e_age - 1])) //Use to fix births to deaths
      
      ode {
        
        //Year of treatment becoming available
        inline treat_start = 12 * ScaleTime * 21 //Treatment first becomes available in 1952
        inline today_treat = 12 * ScaleTime * 59 //Treatment is as available in 1990 as it is today
        inline curr_treat = (t_now > today_treat ? 1 : (t_now < treat_start ? 0 : (t_now - treat_start) / (today_treat - treat_start))) //Scale treatment rate linearly based on time from 1952 to 1990
        // Model equations
        
        dS[bcg, age]/dt = 
          // Disease model updates
          - (1 - (bcg == 1 ? chi[age] : 0)) * foi[age] * S[bcg, age]
          // Demographic model updates
          + (age == 0 ? (bcg == 1 ? gamma[age] * births : (1 - gamma[age]) * births) : 0) //Births
          + (age == 0 ? 0 : (bcg == 1 ? 
          (gamma[age] * theta[age - 1] *  S[0, age - 1] + theta[age - 1] *  S[1, age - 1]) :
          (1 - gamma[age]) * theta[age - 1] *  S[0, age - 1])) //Ageing into bucket
          - theta[age] * S[bcg, age] //Ageing out of bucket
          - mu[age] * S[bcg, age] //All cause (excluding TB) mortality
          dH[bcg, age]/dt = 
          // Disease model updates
          + (1 - (bcg == 1 ? chi[age] : 0)) * foi[age] * S[bcg, age] 
          + (1 - delta) * foi[age] * L[bcg, age] 
          - (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
          - kappa[age] * H[bcg, age]
          // Demographic model updates
          + (age == 0 ? 0 : theta[age - 1] *  H[bcg, age - 1]) //Ageing into bucket
          - theta[age] * H[bcg, age] //Ageing out of bucket
          - mu[age] * H[bcg, age] //All cause (excluding TB) mortality
          dL[bcg, age]/dt = 
          // Disease model updates
          + kappa[age] * H[bcg, age]
          - (1 - delta) * foi[age] * L[bcg, age] 
          - (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age] 
          + phi[age] * T_P[bcg, age] + phi[age] * T_E[bcg, age]
          // Demographic model updates
          + (age == 0 ? 0 : theta[age - 1] *  L[bcg, age - 1]) //Ageing into bucket
          - theta[age] * L[bcg, age] //Ageing out of bucket
          - mu[age] * L[bcg, age] //All cause (excluding TB) mortality
          dP[bcg, age]/dt = 
          // Disease model updates
          + Upsilon[age] * (
              (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
        + (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
          ) 
        + zeta_p[age] * T_P[bcg, age] 
        - curr_treat * nu_p[age] * P[bcg, age]  
        - mu_p[age] * P[bcg, age]
        // Demographic model updates
        + (age == 0 ? 0 : theta[age - 1] *  P[bcg, age - 1]) //Ageing into bucket
        - theta[age] * P[bcg, age] //Ageing out of bucket
        - mu[age] * P[bcg, age] //All cause (excluding TB) mortality
        dE[bcg, age]/dt = 
        // Disease model updates
        + (1 - Upsilon[age]) * (
            (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
        +   (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
        ) 
        + zeta_e[age] * T_E[bcg, age] 
        - curr_treat * nu_e[age] * E[bcg, age]
        - mu_e[age] * E[bcg, age]
        // Demographic model updates
        + (age == 0 ? 0 : theta[age - 1] *  E[bcg, age - 1]) //Ageing into bucket
        - theta[age] * E[bcg, age] //Ageing out of bucket
        - mu[age] * E[bcg, age] //All cause (excluding TB) mortality
        dT_P[bcg, age]/dt = 
        // Disease model updates
        + curr_treat * nu_p[age] * P[bcg, age]
        - zeta_p[age] * T_P[bcg, age] 
        - phi[age] * T_P[bcg, age] 
        - mu_p[age] * T_P[bcg, age]
        // Demographic model updates
        + (age == 0 ? 0 : theta[age - 1] *  T_P[bcg, age - 1]) //Ageing into bucket
        - theta[age] * T_P[bcg, age] //Ageing out of bucket
        - mu[age] * T_P[bcg, age] //All cause (excluding TB) mortality
        dT_E[bcg, age]/dt = 
        // Disease model updates
        + curr_treat * nu_e[age] * E[bcg, age]
        - zeta_e[age] * T_E[bcg, age] 
        - phi[age] * T_E[bcg, age]
        - mu_e[age] * T_E[bcg, age]
        // Demographic model updates
        + (age == 0 ? 0 : theta[age - 1] *  T_E[bcg, age - 1]) //Ageing into bucket
        - theta[age] * T_E[bcg, age] //Ageing out of bucket
        - mu[age] * T_E[bcg, age] //All cause (excluding TB) mortality
        
        
        
        //Accumalator states
        //dPulCases[bcg, age]/dt = nu_p[age] * P[bcg, age]
        dYearlyPulCases[bcg, age]/dt = curr_treat * nu_p[age] * P[bcg, age]
        //dEPulCases[bcg, age]/dt = nu_e[age] * E[bcg, age]
        dYearlyEPulCases[bcg, age]/dt = curr_treat * nu_e[age] * E[bcg, age]  
        //dPulDeaths[bcg, age]/dt = mu_p[age] * T_P[bcg, age] + mu_p[age] * P[bcg, age]
        //dEPulDeaths[bcg, age]/dt =  mu_e[age] * T_E[bcg, age] + mu_e[age] * E[bcg, age]
        dYearlyPulDeaths[bcg, age]/dt = mu_p[age] * T_P[bcg, age] + mu_p[age] * P[bcg, age]
        dYearlyEPulDeaths[bcg, age]/dt = mu_e[age] * T_E[bcg, age] + mu_e[age] * E[bcg, age]
        
      }
      
      // Reporting states
      //By year all summarised
      YearlyAgeCases[age] <- YearlyPulCases[0, age] + YearlyEPulCases[0, age] + YearlyPulCases[1, age] + YearlyEPulCases[1, age]
      YearlyCasesCumSum <- inclusive_scan(YearlyAgeCases)
      YearlyCases <- YearlyCasesCumSum[e_age - 1]
      
      
    }
    
    sub observation {
      
      YearlyInc ~ poisson(rate = YearlyCases)
      YearlyAgeInc[age] ~ poisson(rate = YearlyAgeCases[age])
      
    }
    
}