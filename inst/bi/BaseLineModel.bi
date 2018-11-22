/**
 * Baseline TB and BCG vaccination model
 */
model Baseline {
  
  // Model dimensions
  const e_bcg = 2 // 0 = unvaccinated, 1 = vaccinated
  const e_age = 12 // 0,..9 = 5 year age groups (i.e 0-4), 10 = 50-69 and 11 = 70-89
  const e_AgeGroup = 3 // Age groups for parameters

  dim bcg(e_bcg)
  dim age(e_age)
  dim age2(e_age)
  dim AgeGroup(e_AgeGroup)
  
  //Age at vaccination
  const vac_scheme = 0 //0 = vaccination at school age, 1 = vaccination at birth, 2 = no vaccination.
  const e_d_of_p = 6 // Duration of protection (must be at least 1)
  
  dim d_of_p(e_d_of_p)

  //Control model
  const const_pop = 0 //Set to 1 for constant population (i.e births == deaths)
  const no_age = 0 //Set to 1 to turn off ageing
  const no_disease = 0 //Set to 1 to prevent disease from being initialised / importation
  const noise_switch = 1 // Set noise to 1 to include process noise, and 0 to exclude.
  const scale_rate_treat = 1 //Scale up rate of starting treatment over time (0 to turn off)
  const non_uk_born_scaling = 1 // Scale up of non-UK born cases (from 1960 to 2000). 1 = linear, 2 = log, 3+ = linear
  const beta_df = 1 //Degrees of freedom for transmission prob. 1 = constant across age groups, 2 = modified for children, 3 = modified for children and older adults.
  const mix_type = 1 //Amount of mixing between non UK born and UK born. 1=Homogeneous mixing (0.5 - 1), 2 = hetergeneous mixing (0 - 0.5)ı
  // Time dimensions
  const ScaleTime = 1 / 12 // Scale model over a year 
  //const ScaleTime = 1 // Scale model over a month
  
  // Parameter scales
  const dscale = 12 / 365.25 * ScaleTime 
  const mscale = 1 * ScaleTime
  const yscale = 12 * ScaleTime
  
  //Initialise model
  const init_pop = 37359045 //Estimated intial population - http://www.visionofbritain.org.uk/census/table/EW1931COU1_M3
  const init_P_cases = 49798 // TB cases in England (and Wales)
  const init_E_cases = 16084 // https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/554455/TB_case_notifications_1913_to_2015.pdf

  //Disease model parameters
  
  // Non-UK born mixing
  param M
  // Effective contact rate
  param c_eff
  
  // Historic effective contact rate
  param c_hist
  
  // Modifier for transmission probability
  param beta_child_mod
  param beta_older_adult_mod
  
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
  param phi_15_69
  param phi_70_89
  param phi[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Proportion of cases that have pulmonary TB
  param Upsilon_0_14 //Age specific parameters
  param Upsilon_15_69
  param Upsilon_70_89
  param Upsilon[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Proportion of cases that have pulmonary smear postive TB
  param rho_0_14 //Age specific parameters
  param rho_15_69
  param rho_70_89
  param rho[age](has_output = 0, has_input = 0) // Dummy model parameter
  
  // Rate of starting treatment - pulmonary/extra-pulmonary
  // Pulmonary
  param nu_p_0_14 //Age specific parameters
  param nu_p_15_89
  // Extra-pulmonary
  param nu_e_0_14 //Age specific parameters
  param nu_e_15_89
  
  // Rate loss to follow up - pulmonary/extra-pulmonary
  param zeta_0_14 //Age specific parameters
  param zeta_15_69
  param zeta_70_89
  param zeta[age](has_output = 0, has_input = 0) // Dummy model parameter
    
  // Rate of TB mortality
  param mu_t_0_14 //Age specific parameters
  param mu_t_15_69
  param mu_t_70_89
  param mu_t[age](has_output = 0, has_input = 0) // Dummy model parameter

  //BCG vaccination parameters
  
  // Age specific protection from infection conferred by BCG vaccination
  param chi_init
  // Protection from active disease due to BCG vaccination
  param alpha_t_init // Initial effectiveness of BCG
  param alpha_t_decay //Linear decay observed in suscptibility to disease for those vaccination
  param alpha_t[d_of_p](has_output = 0, has_input = 0) //Estimation effectiveness of BCG vaccine by age group
  
  //Demographic model parameters
  //Ageing
  param theta[age](has_output = 0, has_input = 0)
  
  //Noise parameters
  noise CNoise[age, age2](has_output = 0, has_input = 0) // Sampled contact rate
  noise coverage_sample(has_output = 0, has_input = 0) //Estimated coverage of BCG vaccination in the UK born
    
  // Time varying parameter states
  state CSample[age, age2](has_output = 0, has_input = 0) // Sampled contact rate (symmetric)
  state TotalContacts[age](has_output = 0, has_input = 0) // Average number of contacts (across age groups)
  state beta[age](has_output = 0, has_input = 0) //Probability of transmission
  state foi[age] // force of infection
  state nu_p[age](has_output = 0, has_input = 0) // Dummy model parameter
  state nu_e[age](has_output = 0, has_input = 0) // Dummy model parameter
  // Average rate of starting treatment   
  state avg_nu_p[age](has_output = 0, has_input = 0)
  
  //Calculation parameters
  param I_age[age](has_output = 0, has_input = 0)
  param I_bcg[bcg](has_output = 0, has_input = 0)
    
  //Observational parameters
  param HistMeasError
    
  // Time varing parameter states
  //Demographic
  state births(has_output = 0, has_input = 0) //All births
  state mu_all[age](has_output = 0, has_input = 0) // All cause natural mortality
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
  state T[bcg, age] // TB on treatment
  state N[bcg, age](has_output = 0, has_input = 0) // Overall population
  state NAge[age](has_output = 0, has_input = 0) //Age summed population
  state NSum[age](has_input = 0, has_output = 0) // Sum of population (vector but repeating values)
  state death_sum[age](has_output = 0, has_input = 0) //Used to estimate deaths
  //Temporary population states for discrete demographic model updates
  state d_S[bcg, age](has_output = 0, has_input = 0) // susceptible
  state d_H[bcg, age](has_output = 0, has_input = 0) // high risk latent
  state d_L[bcg, age](has_output = 0, has_input = 0) // low risk latent
  state d_P[bcg, age](has_output = 0, has_input = 0) // pulmonary TB
  state d_E[bcg, age](has_output = 0, has_input = 0) // extra-pulmonary TB only
  state d_T[bcg, age](has_output = 0, has_input = 0) // TB on treatment

  //Accumalator states
  state YearlyPulCases[bcg, age] // yearly pulmonary cases starting treatment
  state YearlyEPulCases[bcg, age] // yearly extra-pulmonary cases starting treatment
  state YearlyDeaths[bcg, age] // TB deaths (yearly)
  
  // Reporting states
  state YearlyPAgeCases[age] 
  state YearlyPCasesCumSum[age](has_input = 0, has_output = 0)
  state YearlyPCases 
  state YearlyEAgeCases[age]
  state YearlyECasesCumSum[age](has_input = 0, has_output = 0)
  state YearlyECases 
    
  state YearlyAgeCases[age]
  state YearlyCases
    
  // States for tracking non UK born in model
  state EstNUKCases[age](has_output = 0, has_input = 0)
  state NoiseNUKCases[age](has_output = 0, has_input = 0)
  state NonUKBornCum[age](has_input = 0, has_output = 0)
  state YearlyNonUKborn
 
    
  //Noise variables
  noise NoiseNUKCasesSample[age](has_output = 0, has_input = 0)
  noise births_sample(has_output = 0, has_input = 0) //Sampled noisy births
  noise mu_all_sample[age](has_output = 0, has_input = 0) //Sampled noisy deaths
    
  //Input
  input births_input //Births (time varying)
  input pop_dist[age] //Population distribution (average from 2000 to 2015).
  input exp_life_span[age] //Expected life span (time varying)
  input polymod[age, age2] //Polymod contact matrix
  input polymod_sd[age, age2] //Polymod SD contact matrix
  input NonUKBornPCases[age] //Non UK born Pulmonary cases (time varying)
  input NUKCases2000[age] //Non UK born cases in 2000.
    
  //Observations
  obs YearlyHistPInc //Historic yearly incidence (pulmonary)
  obs YearlyInc // Yearly overall incidence
  obs YearlyAgeInc[age] // Yearly incidence by age group
  obs YearlyChildInc //Yearly incidence in children
  obs YearlyAdultInc //Yearly incidence in adults
  obs YearlyOlderAdultInc //Yearly incidence in older adults
  
      sub parameter {
        
        // Priors for BCG vaccination
        //Protection from infection at vaccination
        chi_init ~ truncated_gaussian(mean = 0.185, std = 0.0536, lower = 0, upper = 1)
        //Protection from active TB (decaying from time since vaccination)
        alpha_t_init ~ log_gaussian(mean = -1.86, std = 0.22)
        alpha_t_init <- 1 - alpha_t_init
        alpha_t_decay ~ gaussian(mean = -0.134, std = 0.0513)
        alpha_t[0] <- alpha_t_init
        alpha_t[1] <- alpha_t_init + alpha_t_decay 
        alpha_t[2] <- alpha_t_init + 2 * alpha_t_decay 
        alpha_t[3] <- alpha_t_init + 3 * alpha_t_decay 
        alpha_t[4] <- alpha_t_init + 4 * alpha_t_decay 
        alpha_t[5] <- alpha_t_init + 5 * alpha_t_decay 
        
        //Disease priors
        M ~ uniform(0.5, 1)
        M <- (mix_type == 1 ? M : M - 0.5)
        c_eff ~ uniform(0, 5)
        c_hist ~ uniform(10, 15)
        delta ~ truncated_gaussian(mean = 0.78, std = 0.0408, lower = 0, upper = 1)
        
        //Modification of transmission probability by age.
        beta_child_mod ~ truncated_gaussian(mean = 1, std = 0.5, lower = 0)
        beta_older_adult_mod ~ truncated_gaussian(mean = 1, std = 0.5, lower = 0)
        
        // Transition from high risk latent to active TB
        epsilon_h_0_4 ~ truncated_gaussian(mean = (dscale) / 0.00695, std =  (dscale) / 0.0013, lower = 0)
        epsilon_h_5_14  ~ truncated_gaussian(mean = (dscale) / 0.0028, std =  (dscale) / 0.000561, lower = 0)
        epsilon_h_15_89 ~ truncated_gaussian(mean = (dscale) / 0.000335, std =  (dscale) / 0.0000893, lower = 0)
        
        epsilon_h[0] <- epsilon_h_0_4
        epsilon_h[age = 1:2] <-   epsilon_h_5_14
        epsilon_h[age = 3:(e_age - 1)] <- epsilon_h_15_89
        epsilon_h <- 1 / epsilon_h
        
        // Rate of transition from high risk to low risk latents
        kappa_0_4  ~ truncated_gaussian(mean = (dscale) / 0.0133, std = (dscale) / 0.00242, lower = 0)
        kappa_5_14 ~ truncated_gaussian(mean = (dscale) / 0.012, std = (dscale) / 0.00207, lower = 0)
        kappa_15_89 ~ truncated_gaussian(mean = (dscale) / 0.00725, std = (dscale) / 0.00191, lower = 0)
        kappa[0] <- kappa_0_4
        kappa[age = 1:2] <- kappa_5_14
        kappa[age = 3:(e_age - 1)] <- kappa_15_89
        kappa <- 1 / kappa
        
        // Rate of transition for low risk latent to active TB
        epsilon_l_0_4 ~ truncated_gaussian(mean = (dscale) / 0.000008, std =  (dscale) / 0.00000408, lower = 0)
        epsilon_l_5_14 ~ truncated_gaussian(mean = (dscale) / 0.00000984, std =  (dscale) / 0.00000467, lower = 0)
        epsilon_l_15_89  ~ truncated_gaussian(mean = (dscale) / 0.00000595, std =  (dscale) / 0.00000207, lower = 0)
        epsilon_l[0] <- epsilon_l_0_4
        epsilon_l[age = 1:2] <-  epsilon_l_5_14
        epsilon_l[age = 3:(e_age - 1)] <- epsilon_l_15_89
        epsilon_l <- 1 / epsilon_l
        
        // Rate of successful treatment
        phi_0_14 ~ gamma(shape = 9.86, scale = 0.061)
        phi_15_69 ~ gamma(shape = 7.80, scale = 0.0827)
        phi_70_89 ~ gamma(shape = 8.52, scale = 0.0724)
        phi_0_14 <- phi_0_14 * (yscale)
        phi_15_69 <-  phi_15_69 * (yscale)
        phi_70_89 <- phi_70_89 * (yscale)
        phi[age = 0:2] <- phi_0_14
        phi[age = 3:(e_age - 2)] <-  phi_15_69
        phi[(e_age - 1)] <- phi_70_89
        phi <- 1 / phi

        
        // Proportion of TB cases with pulmonary TB
        Upsilon_0_14  ~ truncated_gaussian(mean = 0.629, std = 0.0101, lower = 0, upper = 1)
        Upsilon_15_69 ~ truncated_gaussian(mean = 0.713, std = 0.00377, lower = 0, upper = 1)
        Upsilon_70_89 ~ truncated_gaussian(mean = 0.748, std = 0.00718, lower = 0, upper = 1)
        Upsilon[age = 0:2] <- Upsilon_0_14
        Upsilon[age = 3:(e_age - 2)] <-  Upsilon_15_69
        Upsilon[(e_age - 1)] <- Upsilon_70_89
          
        // Propotion of pulmonary TB cases that are smear positive
        rho_0_14  ~ truncated_gaussian(mean = 0.302, std = 0.0189, lower = 0, upper = 1)
        rho_15_69 ~ truncated_gaussian(mean = 0.637, std = 0.00487, lower = 0, upper = 1)
        rho_70_89 ~ truncated_gaussian(mean = 0.531, std = 0.0107, lower = 0, upper = 1)
        rho[age = 0:2] <- rho_0_14
        rho[age = 3:(e_age - 2)] <-  rho_15_69
        rho[(e_age - 1)] <- rho_70_89
        
        // Rate of starting treatment - pulmonary/extra-pulmonary
        // Pulmonary
        nu_p_0_14  ~ gamma(shape = 0.878, scale = 0.206)
        nu_p_15_89 ~ gamma(shape = 1.1, scale = 0.3)
        nu_p_0_14  <- nu_p_0_14 * (yscale)
        nu_p_15_89 <- nu_p_15_89 * (yscale)
    
        // Extra-pulmonary
        nu_e_0_14 ~ gamma(shape = 0.686, scale = 0.446)
        nu_e_15_89 ~ gamma(shape = 0.897, scale = 0.536)
        nu_e_0_14  <- nu_e_0_14 * (yscale)
        nu_e_15_89 <- nu_e_15_89 * (yscale)
  
        // Rate loss to follow up - pulmonary/extra-pulmonary
        // Extra-pulmonary
        zeta_0_14  ~ truncated_gaussian(mean = (yscale) / 0.00976, std = (yscale) / 0.0179, lower = 0)
        zeta_15_69 ~ truncated_gaussian(mean = (yscale) / 0.0304, std = (yscale) / 0.00764, lower = 0)
        zeta_70_89 ~ truncated_gaussian(mean = (yscale) / 0.00614, std = (yscale) / 0.0159, lower = 0)
        zeta[age = 0:2] <- zeta_0_14
        zeta[age = 3:(e_age - 2)] <-  zeta_15_69
        zeta[(e_age - 1)] <- zeta_70_89
        zeta <- 1 / zeta
            
        // Rate of TB mortality
        mu_t_0_14  ~ truncated_gaussian(mean = (yscale) / 0.00390, std = (yscale) / 0.0180, lower = 0)
        mu_t_15_69 ~ truncated_gaussian(mean = (yscale) / 0.0226, std = (yscale) / 0.00787, lower = 0)
        mu_t_70_89 ~ truncated_gaussian(mean = (yscale) / 0.119, std = (yscale) / 0.0165, lower = 0)
        mu_t[age = 0:2] <- mu_t_0_14
        mu_t[age = 3:(e_age - 2)] <-  mu_t_15_69
        mu_t[(e_age - 1)] <- mu_t_70_89
        mu_t <- 1 / mu_t
        
        //Calculation parameters
        I_age <- 1
        I_bcg <- 1
        
        //Demographic model parameters
        theta[age=0:(e_age - 3)] <- (no_age == 0 ? 1 / (5 * yscale) : 0)
        theta[age=(e_age - 2):(e_age - 1)] <- (no_age == 0 ? 1 / (20 * yscale) : 0)
        
        //Historic measurement error
        HistMeasError ~ truncated_gaussian(mean = 1, std = 0.2, lower = 0)
        
      }
    
    sub initial {
      
      //Initial states
      S[0, age] ~  truncated_gaussian(mean = init_pop * pop_dist[age], std = 0.05 * init_pop * pop_dist[age], lower = 0) // susceptible
      S[0, age] <-  (noise_switch == 0 ? init_pop * pop_dist[age] :  S[0, age])
      S[1, age] <- 0 // BCG vaccinated susceptibles
      H[0, age] ~ truncated_gaussian(mean = (init_E_cases + init_P_cases) * pop_dist[age], std = 0.05 * (init_E_cases + init_P_cases) * pop_dist[age], lower = 0) // high risk latents 
      H[0, age] <- (no_disease == 0 ? (noise_switch == 0 ? (init_E_cases + init_P_cases) * pop_dist[age]: H[0, age]) : 0) 
      H[1, age] <- 0 // BCG high risk latent
      L[0, age] ~ truncated_gaussian(mean = 9 * (init_E_cases + init_P_cases) * pop_dist[age], std = 0.05 * 9 * (init_E_cases + init_P_cases) * pop_dist[age], lower = 0) // high risk latents 
      L[0, age] <- (no_disease == 0 ? (noise_switch == 0 ? 9 * (init_E_cases + init_P_cases) * pop_dist[age] : L[0, age]) : 0) 
      L[1, age] <- 0 // BCG low risk latent
      P[0, age] ~  truncated_gaussian(mean = init_P_cases * pop_dist[age], std = 0.05 * init_P_cases * pop_dist[age], lower = 0) // inital pulmonary cases
      P[0, age] <- (no_disease == 0 ? (noise_switch == 0 ? init_P_cases * pop_dist[age] : P[0, age]) : 0) 
      P[1, age] <- 0 //BCG vaccinated pulmonary TB
      E[0, age] ~  truncated_gaussian(mean = init_E_cases * pop_dist[age], std = 0.05 * init_E_cases * pop_dist[age], lower = 0) // inital pulmonary cases
      E[0, age] <- (no_disease == 0 ? (noise_switch == 0 ? init_E_cases * pop_dist[age] : E[0, age]) : 0) 
      E[1, age] <- 0 // BCG extra-pulmonary TB only
      T[bcg, age] <- 0// TB on treatment
      
      //Priors samples without updating against data
    }
    
    sub transition {
      
      // Reset accumalator variables
      inline yr_reset = yscale
      YearlyPulCases[bcg, age] <- (t_now % 1 < yr_reset ? 0 : YearlyPulCases[bcg, age])
      YearlyEPulCases[bcg, age] <- (t_now % 1 < yr_reset ? 0 : YearlyEPulCases[bcg, age])
      YearlyDeaths[bcg, age] <- (t_now % 1 < yr_reset ? 0 : YearlyDeaths[bcg, age])

      //Apply BCG vaccination to correct populations
      inline policy_change = 74 * yscale // Assume policy switch occurred in 2005
      //Set up age at vaccination
      age_at_vac <- (t_now >= policy_change ? (vac_scheme == 0 ? 3 : (vac_scheme == 1 ? 0 : -1)) : 3)
      
      // Apply linear transform for protection against initial infection
      chi[age] <- (age_at_vac < 0 ? 0 : (age_at_vac > age ? 0 : (age >= (age_at_vac + e_d_of_p) ? 0 : alpha_t[age - age_at_vac] * chi_init / alpha_t[0])))
      //Back calculate protection from latent disease based on initial protection and overall protection
      alpha[age] <- (age_at_vac < 0 ? 0 : (age_at_vac > age ? 0 : (age >= (age_at_vac + e_d_of_p) ? 0 : (alpha_t[age - age_at_vac] - chi[age])/ (1 - chi[age]))))
      //Apply coverage of vac program to correct population
      inline coverage_est = 0.8
      coverage_sample ~ truncated_gaussian(mean = coverage_est, std = 0.05, lower = 0, upper = 1)
      
      // Set vaccination to begin in 1953
      inline vac_start = 22 * yscale
      gamma[age] <- (age_at_vac == age ? (t_now > vac_start ? (noise_switch == 1 ? coverage_sample : coverage_est) : 0) : 0)
      
      //Set time from active symptoms to treatment - adjust based on modern standards and log distribution
      inline treat_start = 21 * yscale //Treatment first becomes available in 1952
      inline modern_treat = 59 * yscale //Treatment reaches modern levels in 1990
      inline scale_infectious_time = (t_now < treat_start ? 0 : (t_now >= modern_treat ? 1 : (scale_rate_treat == 0 ? 1 : log(2 + t_now - treat_start) / log(2 + modern_treat - treat_start))))
      
      nu_p[age = 0:2] <- nu_p_0_14
      nu_p[age = 3:(e_age - 1)] <-  nu_p_15_89
      nu_p <- scale_infectious_time / nu_p
      nu_e[age = 0:2] <-  nu_e_0_14
      nu_e[age = 3:(e_age - 1)] <- nu_e_15_89
      nu_e <- scale_infectious_time / nu_e
      
      // Avg period infectious for pulmonary cases
      avg_nu_p <- inclusive_scan(nu_p)
      avg_nu_p[age] <- avg_nu_p[e_age - 1] / e_age
      
      //Contact rate
      CNoise[age, age2] ~  truncated_gaussian(mean = polymod[age, age2], std = polymod_sd[age, age2], lower = 0)
      CSample <- (noise_switch == 1 ? CNoise : polymod)
      CSample[age, age2] <- CSample[age2, age]
      TotalContacts <- CSample * I_age
      TotalContacts <- inclusive_scan(TotalContacts)
      TotalContacts[age] <- TotalContacts[e_age - 1] / e_age
      
      // Population
      N <- S + H + L + P + E + T
      NAge[age] <- N[0, age] + N[1, age]
      NSum <- inclusive_scan(NAge)
      NSum[age] <- NSum[e_age - 1]
      
      // Estimate force of infection - start with probability of transmission
      inline modern_contacts = 59 * yscale // Modern day is 1990 with a baseline date of 1931
      inline scale_historic_contacts = (t_now > modern_contacts ? 0 : 1 - log(t_now + 1) / log(modern_treat + 1))
      beta <- avg_nu_p * (c_eff + c_hist * scale_historic_contacts) 
      beta <- beta ./ TotalContacts
      
      beta[age = 0:2] <- (beta_df == 1 ? beta[age] : beta[age] * beta_child_mod)
      beta[11] <-  (beta_df < 3 ? beta[11] : beta[11] * beta_older_adult_mod)
      // Estimate the number of nonuk born cases
      inline nuk_start = (29 * yscale) //Start introducing non-UK born cases from 1960
      inline nuk_data = (69 * yscale)  //Start using data from 2000
      EstNUKCases[age] <- (t_now <  nuk_data ?  
                                (t_now < nuk_start ? 0 : 
                                   (non_uk_born_scaling == 1 ?   (t_now - nuk_start) / nuk_data * NUKCases2000[age] / (yscale) : //Estimate cases using a linear relationship based on cases in 2000:
                                      (non_uk_born_scaling == 2 ?  log(2 + t_now - nuk_start) / log(2 + nuk_data - nuk_start) * NUKCases2000[age] : //Assume that cases increased using a log link
                                         NUKCases2000[age]))) : //Assume nonUK born cases were constant from 1960 to 2000)))
                                  NonUKBornPCases[age])
      NoiseNUKCasesSample[age] ~ truncated_gaussian(mean = EstNUKCases[age], std =  0.05 * EstNUKCases[age], lower = 0)
      NoiseNUKCases <- (noise_switch == 1 ? NoiseNUKCasesSample : EstNUKCases)
      
      //Now build force of infection
      foi <- transpose(P) * I_bcg
      foi <- rho .* foi + M * NoiseNUKCases ./ nu_p
      foi <- CSample * foi
      foi <- beta .* foi ./ NSum
      
      //All-cause mortality excluding TB
      mu_all_sample[age] ~ truncated_gaussian(mean = exp_life_span[age], std = exp_life_span[age]*0.05, lower = 0)
      
      mu_all <- (noise_switch == 1 ? mu_all_sample : exp_life_span)
      
      mu[age] <- 1 / mu_all[age] - (mu_t[age] * (P[0, age] + P[1, age] + E[0, age] + E[1, age] + T[0, age] + T[1, age])) / NAge[age]
      // All used to fix births to deaths for testing 
      death_sum <-  NAge ./ mu_all
      death_sum <- inclusive_scan(death_sum)
      births_sample ~ gaussian(mean = births_input, std = 0.05 * births_input)
      births <- (const_pop == 0 ? (noise_switch == 1 ? births_sample : births_input) : death_sum[e_age - 1] + theta[e_age - 1] * (N[0, e_age - 1] +  N[1, e_age - 1])) //Use to fix births to deaths
      
      // Demographic model updates (discrete at start of year)
      d_S[bcg, age] <-  S[bcg, age]         
        // Demographic model updates
        + (age == 0 ? (bcg == 1 ? gamma[age] * births : (1 - gamma[age]) * births) : 0) //Births
        + (age == 0 ? 0 : (bcg == 1 ? 
        (gamma[age] * theta[age - 1] *  S[0, age - 1] + theta[age - 1] *  S[1, age - 1]) :
                             (1 - gamma[age]) * theta[age - 1] *  S[0, age - 1])) //Ageing into bucket
      - theta[age] * S[bcg, age] //Ageing out of bucket
      - mu[age] * S[bcg, age] //All cause (excluding TB) mortality
      d_H[bcg, age] <- H[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  H[bcg, age - 1]) //Ageing into bucket
      - theta[age] * H[bcg, age] //Ageing out of bucket
      - mu[age] * H[bcg, age] //All cause (excluding TB) mortality
      d_L[bcg, age] <- L[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  L[bcg, age - 1]) //Ageing into bucket
      - theta[age] * L[bcg, age] //Ageing out of bucket
      - mu[age] * L[bcg, age] //All cause (excluding TB) mortality
      d_P[bcg, age] <- P[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  P[bcg, age - 1]) //Ageing into bucket
      - theta[age] * P[bcg, age] //Ageing out of bucket
      - mu[age] * P[bcg, age] //All cause (excluding TB) mortality
      d_E[bcg, age] <- E[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  E[bcg, age - 1]) //Ageing into bucket
      - theta[age] * E[bcg, age] //Ageing out of bucket
      - mu[age] * E[bcg, age] //All cause (excluding TB) mortality
      d_T[bcg, age] <- T[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  T[bcg, age - 1]) //Ageing into bucket
      - theta[age] * T[bcg, age] //Ageing out of bucket
      - mu[age] * T[bcg, age] //All cause (excluding TB) mortality

      //Update all states and force to be above 0
      S[bcg, age] <- (d_S[bcg, age] < 0 ? 0 : d_S[bcg, age])
      H[bcg, age] <- (d_H[bcg, age] < 0 ? 0 : d_H[bcg, age])
      L[bcg, age] <- (d_L[bcg, age] < 0 ? 0 : d_L[bcg, age])
      P[bcg, age] <- (d_P[bcg, age] < 0 ? 0 : d_P[bcg, age])
      E[bcg, age] <- (d_E[bcg, age] < 0 ? 0 : d_E[bcg, age])
      T[bcg, age] <- (d_T[bcg, age] < 0 ? 0 : d_T[bcg, age])
      
      
      ode(alg='RK4(3)', h=1e-1, atoler=1e-2, rtoler=1e-5) {
        
        // Model equations
        dS[bcg, age]/dt = 
          // Disease model updates
          - (1 - (bcg == 1 ? chi[age] : 0)) * foi[age] * S[bcg, age]
          dH[bcg, age]/dt = 
          // Disease model updates
          + (1 - (bcg == 1 ? chi[age] : 0)) * foi[age] * S[bcg, age] 
          + (1 - delta) * foi[age] * L[bcg, age] 
          - (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
          - kappa[age] * H[bcg, age]
          dL[bcg, age]/dt = 
          // Disease model updates
          + kappa[age] * H[bcg, age]
          - (1 - delta) * foi[age] * L[bcg, age] 
          - (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age] 
          + phi[age] * T[bcg, age]
          dP[bcg, age]/dt = 
          // Disease model updates
          + Upsilon[age] * (
              (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
        + (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
        + zeta[age] * T[bcg, age] 
          ) 
        
        - nu_p[age] * P[bcg, age]  
        - mu_t[age] * P[bcg, age]
        dE[bcg, age]/dt = 
        // Disease model updates
        + (1 - Upsilon[age]) * (
            (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
        +   (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
        + zeta[age] * T[bcg, age]
        ) 
        - nu_e[age] * E[bcg, age]
        - mu_t[age] * E[bcg, age]
        dT[bcg, age]/dt = 
        // Disease model updates
        + nu_p[age] * P[bcg, age]
        + nu_e[age] * E[bcg, age]
        - zeta[age] * T[bcg, age] 
        - phi[age] * T[bcg, age] 
        - mu_t[age] * T[bcg, age]
        
        //Accumalator states
        dYearlyPulCases[bcg, age]/dt = nu_p[age] * P[bcg, age]
        dYearlyEPulCases[bcg, age]/dt = nu_e[age] * E[bcg, age]  
        dYearlyDeaths[bcg, age]/dt = mu_t[age] * T[bcg, age] + mu_t[age] * (P[bcg, age] + E[bcg, age])
        
      }
      //Enforce reporting states to be above 0
      YearlyPulCases[bcg, age] <- (YearlyPulCases[bcg, age] < 0 ? 0 : YearlyPulCases[bcg, age])
      YearlyEPulCases[bcg, age]<- (YearlyEPulCases[bcg, age] < 0 ? 0 : YearlyEPulCases[bcg, age])

      // Reporting states
      //By year all summarised reporting states
      YearlyPAgeCases[age] <-  YearlyPulCases[0, age] + YearlyPulCases[1, age]
      YearlyPCasesCumSum <- inclusive_scan(YearlyPAgeCases)
      YearlyPCases <- YearlyPCasesCumSum[e_age - 1]
      
      YearlyEAgeCases[age] <-  YearlyEPulCases[0, age] + YearlyEPulCases[1, age]
      YearlyECasesCumSum <- inclusive_scan(YearlyEAgeCases)
      YearlyECases <- YearlyECasesCumSum[e_age - 1]
      
      YearlyAgeCases <- YearlyPAgeCases + YearlyEAgeCases
      YearlyCases <- YearlyPCases + YearlyECases
      
      //Non UK born cases
      NonUKBornCum <- inclusive_scan(EstNUKCases)
      YearlyNonUKborn <- NonUKBornCum[e_age - 1]
    }
    
    sub observation {
      
      YearlyHistPInc ~ poisson(rate = HistMeasError * (YearlyPCases + YearlyNonUKborn))
      YearlyInc ~ poisson(rate = YearlyCases)
      YearlyAgeInc[age] ~ poisson(rate = YearlyAgeCases[age])
      YearlyChildInc ~ poisson(rate = YearlyAgeCases[0] + YearlyAgeCases[1] + YearlyAgeCases[2])
      YearlyAdultInc ~ poisson(rate = YearlyAgeCases[3] + YearlyAgeCases[4] + YearlyAgeCases[5]
                                    + YearlyAgeCases[6] + YearlyAgeCases[7] + YearlyAgeCases[8]
                                    + YearlyAgeCases[9] + YearlyAgeCases[10])
      YearlyOlderAdultInc ~ poisson(rate = YearlyAgeCases[11])
      
    }
}