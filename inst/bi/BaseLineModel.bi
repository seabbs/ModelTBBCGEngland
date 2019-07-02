/**
 * Baseline TB and BCG vaccination model
 */
model Baseline {
  
  //Timestep
  const timestep = 1 //timestep (if using)
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
  const initial_uncertainty_switch = 1 //Set to 1 to include initial state and parameter uncertainty. 0 o exclude.
  const noise_switch = 1 // Set noise to 1 to include process noise, and 0 to exclude.
  const scale_rate_treat = 1 //Scale up rate of starting treatment over time (0 to turn off)
  const non_uk_born_scaling = 1 // Scale up of non-UK born cases (from 1960 to 2000). 1 = linear, 2 = log, 3+ = linear
  const beta_df = 1 //Degrees of freedom for transmission prob. 1 = constant across age groups, 2 = modified for children, 3 = modified for children and older adults.
  const M_df = 1 //Degrees of freedom for non-UK born mixing. 1 = constant across age groups, 2 = modified for children, 3 = modified for children and older adults.
  const measurement_model = 1 // Should the measurement model be given non-constant priors
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
  const initial_infectious_period = 2 // Assume that cases are initially infectious for 2 years prior to detection and isolation/treastment by health services
  //Disease model parameters
  
  // Non-UK born mixing
  param M
  
  //Non-UK born scaling from 1960 to 1999
  param NonUKScaling 
  // Effective contact rate
  param c_eff
  
  // Historic effective contact rate
  param c_hist
  
  // Historic effective contact rate half life
  param c_hist_half

  // Modifier for transmission probability
  param beta_child
  param beta_older_adult
  
  // Modifier for non-UK born mixing
  param M_child
  param M_older_adult
  
  
  // Transition from low risk latent disease to active disease
  param epsilon_l_0_4 //Age specific parameters
  param epsilon_l_5_14
  param epsilon_l_15_69
  param epsilon_l_70_89
  
  
  // Protection from infection due to prior latent infection
  state delta(has_output = 0, has_input = 0) 
  
  // Transition from high risk latent disease to active disease
  state epsilon_h_0_4(has_output = 0, has_input = 0)  //Age specific parameters
  state epsilon_h_5_14(has_output = 0, has_input = 0) 
  state epsilon_h_15_89(has_output = 0, has_input = 0) 

  // Transition to low risk latent disease from high risk latent disease
  state kappa_0_4(has_output = 0, has_input = 0)  //Age specific parameters
  state kappa_5_14(has_output = 0, has_input = 0) 
  state kappa_15_89(has_output = 0, has_input = 0) 

  // Rate of succesful treatment completion
  state phi_0_14(has_output = 0, has_input = 0)  //Age specific parameters
  state phi_15_69(has_output = 0, has_input = 0) 
  state phi_70_89(has_output = 0, has_input = 0) 
    
  // Proportion of cases that have pulmonary TB
  state Upsilon_0_14(has_output = 0, has_input = 0)  //Age specific parameters
  state Upsilon_15_69(has_output = 0, has_input = 0) 
  state Upsilon_70_89(has_output = 0, has_input = 0) 

  // Proportion of cases that have pulmonary smear postive TB
  state rho_0_14(has_output = 0, has_input = 0)  //Age specific parameters
  state rho_15_69(has_output = 0, has_input = 0) 
  state rho_70_89(has_output = 0, has_input = 0) 

  // Rate of starting treatment - pulmonary/extra-pulmonary
  // Pulmonary
  state nu_p_0_14(has_output = 0, has_input = 0) //Age specific parameters
  state nu_p_15_89(has_output = 0, has_input = 0)
  state scaled_nu_p_0_14(has_output = 0, has_input = 0) //Age specific parameters
  state scaled_nu_p_15_89(has_output = 0, has_input = 0)
    
  // Extra-pulmonary
  state nu_e_0_14(has_output = 0, has_input = 0) //Age specific parameters
  state nu_e_15_89(has_output = 0, has_input = 0)
  state scaled_nu_e_0_14(has_output = 0, has_input = 0) //Age specific parameters
  state scaled_nu_e_15_89(has_output = 0, has_input = 0)
    
  // Rate loss to follow up - pulmonary/extra-pulmonary
  state zeta_0_14(has_output = 0, has_input = 0)  //Age specific parameters
  state zeta_15_69(has_output = 0, has_input = 0) 
  state zeta_70_89(has_output = 0, has_input = 0) 

  // Rate of TB mortality
  state mu_t_0_14(has_output = 0, has_input = 0)  //Age specific parameters
  state mu_t_15_69(has_output = 0, has_input = 0) 
  state mu_t_70_89(has_output = 0, has_input = 0) 
    
  //BCG vaccination parameters
  
  // Age specific protection from infection conferred by BCG vaccination
  state chi_init(has_output = 0, has_input = 0) 
  // Protection from active disease due to BCG vaccination
  state alpha_t[d_of_p](has_output = 0, has_input = 0) //Estimation effectiveness of BCG vaccine by age group

    
  //Demographic model parameters
  //Ageing
  param theta[age](has_output = 0, has_input = 0)
  
  //Noise parameters
  noise CSample[age, age2](has_output = 0, has_input = 0) // Sampled contact rate
  noise gamma[age](has_output = 0, has_input = 0)   //Coverage of the vaccination program by age
    
  // Time varying parameter states
  state foi[age](has_output = 0, has_input = 0)// force of infection

  //Calculation parameters
  param I_age[age](has_output = 0, has_input = 0)
  param I_bcg[bcg](has_output = 0, has_input = 0)
  
  //Latent case distribution
  param latent_dist[age](has_output = 0, has_input = 0)
    
  //Observational parameters
  param MeasError
  param MeasStd
    
  // Time varing parameter states
  //Demographic
  state mu[age](has_output = 0, has_input = 0) //All cause mortality excluding TB
  
  //Average rate out of pulmonary TB population
  state avg_rate_from_pulmonary(has_output = 0, has_input = 0)
    
  //Vaccination
  state age_at_vac(has_output = 0, has_input = 0) //Age at vaccination
  state chi[age](has_output = 0, has_input = 0) //Protection from initial infection due to BCG vaccination
  state alpha[age](has_output = 0, has_input = 0) //Protection from active disease due to BCG vaccination

  //Population states
  state S[bcg, age] // susceptible
  state H[bcg, age] // high risk latent
  state L[bcg, age] // low risk latent
  state P[bcg, age] // pulmonary TB
  state E[bcg, age] // extra-pulmonary TB only
  state T_E[bcg, age] // TB on treatment (extra-pulmonary)
  state T_P[bcg, age] // TB on treatment (pulmonary)
  state N[bcg, age](has_output = 0, has_input = 0) // Overall population
  state death_sum[age](has_output = 0, has_input = 0) //Used to estimate deaths

  //Accumalator states
  state YearlyPulCases[bcg, age] // yearly pulmonary cases starting treatment
  state YearlyEPulCases[bcg, age] // yearly extra-pulmonary cases starting treatment
  state YearlyDeaths[bcg, age] // TB deaths (yearly)
  
  // Reporting states
  state YearlyPCases 
  state YearlyECases 
    
  state YearlyAgeCases[age]
    
  // States for tracking non UK born in model
  state EstNUKCases[age](has_output = 0, has_input = 0)
    
  //Noise variables
  noise NoiseNUKCases[age]
  noise births(has_output = 0, has_input = 0) //Sampled noisy births
  noise mu_all[age](has_output = 0, has_input = 0) // All cause natural mortality
    
  //Input
  input births_input //Births (time varying)
  input pop_dist[age] //Population distribution (average from 2000 to 2015).
  input exp_life_span[age] //Expected life span (time varying)
  input polymod[age, age2] //Polymod contact matrix
  input polymod_sd[age, age2] //Polymod SD contact matrix
  input avg_contacts  // Average number of contacts (across age groups)
  input NonUKBornPCases[age] //Non UK born Pulmonary cases (time varying)
  input NUKCases2000[age] //Non UK born cases in 2000.
  input DistUKCases2000[age] //Distribution of UK born cases in 2000
  
  //Observations
  obs YearlyHistPInc //Historic yearly incidence (pulmonary)
  obs YearlyInc // Yearly overall incidence
  obs YearlyAgeInc[age] // Yearly incidence by age group
  obs YearlyChildInc //Yearly incidence in children
  obs YearlyAdultInc //Yearly incidence in adults
  obs YearlyOlderAdultInc //Yearly incidence in older adults

  sub parameter {
        
        //Disease priors
        M  ~ uniform(0, 20)
        c_eff ~ truncated_gaussian(mean = 1, std = 1, lower = 0)
        c_hist ~ uniform(10, 20)
    
    
        // Half life of historic effect contact rate.
        c_hist_half ~ truncated_gaussian(mean = 5, std = 2, lower = 0)
    
        //Modification of transmission probability by age.
        beta_child ~ truncated_gaussian(mean = 0.5, 
                                            std = (beta_df == 1 ? 0 : 0.25), 
                                            lower = 0,
                                            upper = 1)
        beta_older_adult ~ truncated_gaussian(mean = 0.5, std = (beta_df < 3 ? 0 : 0.25), 
                                              lower = 0,
                                              upper = 1)
        
        
        // Modification of non-UK born mixing by age.
        M_child ~ truncated_gaussian(mean = 0.5, 
                                     std = (M_df == 1 ? 0 : 0.25), 
                                     lower = 0,
                                     upper = 1)
          
        M_older_adult ~ truncated_gaussian(mean = 0.5, 
                                           std = (M_df == 1 ? 0 : 0.25), 
                                           lower = 0,
                                           upper = 1)
    
        
        //Non-UK born scaling
        NonUKScaling ~ truncated_gaussian(mean = 0, std = 20)
    
        //Measurement error
        MeasError ~ truncated_gaussian(mean = 0.9, 
                                       std = (measurement_model == 0 ? 0 :  0.05), 
                                       lower = 0.8,
                                       upper = 1.0)
    
        // Rate of transition for low risk latent to active TB
        epsilon_l_0_4 ~ truncated_gaussian(mean = 0.000008 / dscale, 
                                           std = 0.00000408 / dscale, 
                                           lower = 0)
        epsilon_l_5_14 ~ truncated_gaussian(mean =  0.00000984 / dscale, 
                                            std = 0.00000467 / dscale,
                                            lower = 0)
        epsilon_l_15_69  ~ truncated_gaussian(mean = 0.00000595 / (2 * dscale), 
                                              std = 0.00000207 / dscale,
                                              lower = 0)
        
        epsilon_l_70_89  ~ truncated_gaussian(mean = 0.00000595 * 2 / dscale, 
                                              std = 0.00000207 / dscale,
                                              lower = 0)
        
        // Prior on measurement Std
        MeasStd ~ uniform(0, 0.05)
    
        //Calculation parameters
        I_age <- 1
        I_bcg <- 1
        
        //Demographic model parameters
        theta[age=0:(e_age - 3)] <- (no_age == 0 ? 1 / (5 * yscale) : 0)
        theta[age=(e_age - 2):(e_age - 1)] <- (no_age == 0 ? 1 / (20 * yscale) : 0)
    
        //Latent distribution
        latent_dist <- exclusive_scan(DistUKCases2000)
          
      }
  
  sub proposal_parameter {
    //Proposal at 5% of prior SD or range
    inline proposal_scaling = 2
    //Disease priors
    M ~ truncated_gaussian(mean = M,
                           std = 2 / proposal_scaling,
                           lower = 0,
                           upper = 20)
    c_eff ~ truncated_gaussian(mean = c_eff,
                               std = 0.1 / proposal_scaling,
                               lower = 0)
    c_hist ~ truncated_gaussian(mean = c_hist,
                                std = 1 / proposal_scaling,
                                lower = 10,
                                upper = 20)
    
    c_hist_half ~ truncated_gaussian(mean = c_hist_half,
                                std = 0.5 / proposal_scaling,
                                lower = 0)
    
    //Modification of transmission probability by age.
    beta_child ~ truncated_gaussian(mean = beta_child, 
                                    std = (beta_df == 1 ? 0 : 0.025 / proposal_scaling), 
                                    lower = 0,
                                    upper = 1)
    beta_older_adult ~  truncated_gaussian(mean = beta_older_adult, 
                                          std = (beta_df < 3 ? 0 : 0.025 / proposal_scaling), 
                                          lower = 0,
                                          upper = 1)
    
    
    //Modification of non-UK born mixing by age.
    M_child ~ truncated_gaussian(mean = M_child, 
                                 std = (M_df == 1 ? 0 : 0.025 / proposal_scaling), 
                                 lower = 0,
                                 upper = 1)
    M_older_adult ~  truncated_gaussian(mean = M_older_adult, 
                                        std = (M_df < 3 ? 0 : 0.025 / proposal_scaling), 
                                        lower = 0,
                                        upper = 1)
    
    
    //Non-UK born scaling
    NonUKScaling ~ truncated_gaussian(mean = NonUKScaling,
                                      std = 2 / proposal_scaling)
    
    //Measurement error
    MeasError ~ truncated_gaussian(mean = MeasError,
                                   std =  (measurement_model == 0 ? 0 :  0.005 / proposal_scaling),
                                   lower = 0)
    
    // Prior on measurement Std
    MeasStd ~ truncated_gaussian(mean = MeasStd,
                                 std =  (measurement_model == 0 ? 0 : 0.005 / proposal_scaling), 
                                 lower = 0,
                                 upper = 0.05)
    
  }
  
    sub initial {
    
      //Priors samples without updating against data
      // Priors for BCG vaccination
      //Protection from infection at vaccination
      chi_init ~ truncated_gaussian(mean = 0.185, std = 0.0536, lower = 0, upper = 1)
      
      chi_init <- (initial_uncertainty_switch == 0 ? 0.185 : chi_init)
      
      //Protection from active TB
      alpha_t[0] ~ log_gaussian(mean = -1.86, 
                                std = (initial_uncertainty_switch == 0 ? 0 : 0.22))
      alpha_t[1] ~ log_gaussian(mean = -1.19,
                                std = (initial_uncertainty_switch == 0 ? 0 : 0.24))
      alpha_t[2] ~ log_gaussian(mean = -0.84, 
                                std = (initial_uncertainty_switch == 0 ? 0 : 0.22))
      alpha_t[3] ~ log_gaussian(mean = -0.84, 
                                std = (initial_uncertainty_switch == 0 ? 0 : 0.2))
      alpha_t[4] ~ log_gaussian(mean = -0.28, 
                                std = (initial_uncertainty_switch == 0 ? 0 : 0.19))
      alpha_t[5] ~ log_gaussian(mean = -0.23, 
                                std = (initial_uncertainty_switch == 0 ? 0 : 0.29))
      
      alpha_t <- 1 - alpha_t
      
      //Disease priors
      delta ~ truncated_gaussian(mean = 0.78, std = (initial_uncertainty_switch == 0 ? 0 : 0.0408), 
                                 lower = 0, upper = 1)
      
      // Transition from high risk latent to active TB
      epsilon_h_0_4 ~ truncated_gaussian(mean = 0.00695,
                                         std =  (initial_uncertainty_switch == 0 ? 0 : 0.0013), lower = 0)
      epsilon_h_5_14  ~ truncated_gaussian(mean = 0.0028,
                                           std =  (initial_uncertainty_switch == 0 ? 0 : 0.000561), lower = 0)
      epsilon_h_15_89 ~ truncated_gaussian(mean = 0.000335, 
                                           std =  (initial_uncertainty_switch == 0 ? 0 : 0.0000893), lower = 0)
      
      epsilon_h_0_4 <-  epsilon_h_0_4 / dscale
      epsilon_h_5_14 <- epsilon_h_5_14 / dscale
      epsilon_h_15_89 <- epsilon_h_15_89 / dscale
      
      // Rate of transition from high risk to low risk latents
      kappa_0_4  ~ truncated_gaussian(mean = 0.0133, 
                                      std =  (initial_uncertainty_switch == 0 ? 0 : 0.00242), lower = 0)
      kappa_5_14 ~ truncated_gaussian(mean =  0.012,
                                      std =  (initial_uncertainty_switch == 0 ? 0 : 0.00207), lower = 0)
      kappa_15_89 ~ truncated_gaussian(mean = 0.00725, 
                                       std = (initial_uncertainty_switch == 0 ? 0 : 0.00191), lower = 0)
      
      kappa_0_4 <- kappa_0_4 / dscale
      kappa_5_14 <- kappa_5_14 / dscale
      kappa_15_89 <- kappa_15_89 / dscale
      
      // Rate of successful treatment
      phi_0_14 ~ truncated_gaussian(mean = yscale * 0.606, 
                                    std =  (initial_uncertainty_switch == 0 ? 0 : yscale * 0.237),
                                    lower = 4 / 12)
      phi_15_69 ~ truncated_gaussian(mean = yscale * 0.645, 
                                     std =  (initial_uncertainty_switch == 0 ? 0 : yscale * 0.290), 
                                     lower = 4 / 12)
      phi_70_89 ~ truncated_gaussian(mean = yscale * 0.616, 
                                     std =  (initial_uncertainty_switch == 0 ? 0 :  yscale * 0.265),
                                     lower = 4 / 12)
      
      phi_0_14 <- 1 /  phi_0_14
      phi_15_69 <-  1 /  phi_15_69
      phi_70_89 <- 1 /  phi_70_89
      
      // Rate of starting treatment - pulmonary/extra-pulmonary
      // Pulmonary
      nu_p_0_14  ~ truncated_gaussian(mean = yscale * 0.181, 
                                      std =  (initial_uncertainty_switch == 0 ? 0 :  yscale * 0.310),
                                      lower = 0)
      nu_p_15_89 ~ truncated_gaussian(mean = yscale * 0.328,
                                      std =   (initial_uncertainty_switch == 0 ? 0 :  yscale * 0.447),
                                      lower = 0)
      nu_p_0_14  <- 1 / nu_p_0_14 
      nu_p_15_89 <- 1 / nu_p_15_89
      
      // Extra-pulmonary
      nu_e_0_14 ~ truncated_gaussian(mean = yscale * 0.306,
                                     std =   (initial_uncertainty_switch == 0 ? 0 :  yscale * 0.602),
                                     lower = 0)
      nu_e_15_89 ~ truncated_gaussian(mean = yscale * 0.480,
                                      std =  (initial_uncertainty_switch == 0 ? 0 : yscale * 0.866),
                                      lower = 0)
      nu_e_0_14  <- 1 / nu_e_0_14 
      nu_e_15_89 <- 1 / nu_e_15_89 
      
      
      // Rate loss to follow up - pulmonary/extra-pulmonary
      // Extra-pulmonary
      zeta_0_14 ~ truncated_gaussian(mean = yscale * 0.00976,
                                     std = (initial_uncertainty_switch == 0 ? 0 :  yscale * 0.0179),
                                     lower = 0)
      zeta_15_69 ~ truncated_gaussian(mean = yscale * 0.0304,
                                      std = (initial_uncertainty_switch == 0 ? 0 : yscale * 0.00764),
                                      lower = 0)
      zeta_70_89 ~ truncated_gaussian(mean = yscale * 0.00614,
                                      std = (initial_uncertainty_switch == 0 ? 0 : yscale * 0.0159),
                                      lower = 0)
      
      // Rate of TB mortality
      mu_t_0_14 ~ truncated_gaussian(mean = yscale * 0.00390, 
                                     std = (initial_uncertainty_switch == 0 ? 0 : yscale * 0.0180),
                                     lower = 0)
      mu_t_15_69 ~ truncated_gaussian(mean = yscale * 0.0226,
                                      std = (initial_uncertainty_switch == 0 ? 0 : yscale * 0.00787),
                                      lower = 0)
      mu_t_70_89 ~ truncated_gaussian(mean = yscale * 0.117,
                                      std = (initial_uncertainty_switch == 0 ? 0 : yscale * 0.0165),
                                      lower = 0)
      
      // Proportion of TB cases with pulmonary TB
      Upsilon_0_14 ~ truncated_gaussian(mean = 0.629, 
                                        std = (initial_uncertainty_switch == 0 ? 0 : 0.0101), 
                                        lower = 0, upper = 1)
      Upsilon_15_69 ~ truncated_gaussian(mean = 0.713, 
                                         std = (initial_uncertainty_switch == 0 ? 0 : 0.00377), 
                                         lower = 0, upper = 1)
      Upsilon_70_89 ~ truncated_gaussian(mean = 0.748, 
                                         std = (initial_uncertainty_switch == 0 ? 0 : 0.00718), 
                                         lower = 0, upper = 1)
      
      // Propotion of pulmonary TB cases that are smear positive
      rho_0_14 ~ truncated_gaussian(mean = 0.302, 
                                    std = (initial_uncertainty_switch == 0 ? 0 : 0.0189), 
                                    lower = 0, upper = 1)
      rho_15_69 ~ truncated_gaussian(mean = 0.637, 
                                     std = (initial_uncertainty_switch == 0 ? 0 : 0.00487),
                                     lower = 0, upper = 1)
      rho_70_89 ~ truncated_gaussian(mean = 0.531, 
                                     std = (initial_uncertainty_switch == 0 ? 0 : 0.0107), 
                                     lower = 0, upper = 1)
      
      
      // Estimate the proprotion that go from high risk latency to active disease
      inline avg_rate_high_disease = (epsilon_h_0_4 + 2 * epsilon_h_5_14 + 9 * epsilon_h_15_89) / 12
      inline avg_rate_high_low =  (kappa_0_4 + 2 * kappa_5_14 + 9 * kappa_15_89) / 12
      inline avg_prop_high_disease = avg_rate_high_disease / (avg_rate_high_disease + avg_rate_high_low)
      
      //Scale the initial disease states + summaries
      inline scaled_pul = init_P_cases 
      inline scaled_epul = init_E_cases
      inline overall_cases = scaled_pul + scaled_epul
      inline initial_high_risk = 1 / avg_rate_high_disease * 1 / avg_prop_high_disease * overall_cases / initial_infectious_period
      inline initial_latent = initial_high_risk / 2
      
      //Initial states
      S[0, age] ~  truncated_gaussian(mean = init_pop * pop_dist[age], std = (initial_uncertainty_switch == 0 ? 0 : 0.05 * init_pop * pop_dist[age]), lower = 0) // susceptible
      S[1, age] <- 0 // BCG vaccinated susceptibles
      H[0, age] ~ truncated_gaussian(mean = initial_high_risk * DistUKCases2000[age], std = (initial_uncertainty_switch == 0 ? 0 : 0.05 * initial_high_risk * DistUKCases2000[age]), lower = 0) // high risk latents 
      H[1, age] <- 0 // BCG high risk latent
      L[0, age] ~ truncated_gaussian(mean = initial_latent * latent_dist[age], std = (initial_uncertainty_switch == 0 ? 0 : 0.05 * initial_latent * latent_dist[age]), lower = 0) // low risk latents 
      S[0, age] <-  S[0, age]  - L[0, age]
      L[1, age] <- 0 // BCG low risk latent
      P[0, age] ~  truncated_gaussian(mean = scaled_pul * DistUKCases2000[age], std = (initial_uncertainty_switch == 0 ? 0 : 0.05 * scaled_pul * DistUKCases2000[age]), lower = 0) // inital pulmonary cases
      P[1, age] <- 0 //BCG vaccinated pulmonary TB
      E[0, age] ~  truncated_gaussian(mean = scaled_epul * DistUKCases2000[age], std = (initial_uncertainty_switch == 0 ? 0 : 0.05 * scaled_epul * DistUKCases2000[age]), lower = 0) // inital pulmonary cases
      E[1, age] <- 0 // BCG extra-pulmonary TB only
      T_E[bcg, age] <- 0// TB on treatment (extra-pulmonary)
      T_P[bcg, age] <- 0// TB on treatment (pulmonary)
      
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
      
      // Set vaccination to begin in 1953
      inline vac_start = 22 * yscale
      
      //Apply coverage of vac program to correct population
      gamma[age] ~ truncated_gaussian(mean = (age_at_vac == age ? 
                                                (t_now > vac_start ? 
                                                0.8 : 0) : 0), 
                                      std = (age_at_vac == age ? 
                                               (t_now > vac_start ? 
                                               (noise_switch == 0 ? 0 : 0.05) : 0) : 0), 
                                      lower = 0, upper = 1)
      
      
      //Set time from active symptoms to treatment - adjust based on modern standards
      //Linear scaled between 0 and 1
      inline treat_start = 21 * yscale //Treatment first becomes available in 1952
      inline modern_treat = 59 * yscale //Treatment reaches modern levels in 1990
      inline scale_infectious_time = (t_now < treat_start ? 0 : 
                                        (t_now >= modern_treat ? 1 : 
                                           (scale_rate_treat == 0 ? 1 :
                                              (t_now - treat_start) / (modern_treat - treat_start))))
      
      // Restrict treatment time to a maximum of 2 years (unless they die in the meantime)
      scaled_nu_p_0_14 <- 1 / initial_infectious_period + scale_infectious_time * (nu_p_0_14 - 1 / initial_infectious_period)
      scaled_nu_p_15_89 <- 1 / initial_infectious_period + scale_infectious_time * (nu_p_15_89 - 1 /  initial_infectious_period)
      scaled_nu_e_0_14 <- 1 / initial_infectious_period + scale_infectious_time * (nu_e_0_14 - 1 / initial_infectious_period)
      scaled_nu_e_15_89 <- 1 / initial_infectious_period + scale_infectious_time * (nu_e_15_89 - 1 / initial_infectious_period)
  
      //Contact rate - sample
      CSample[age, age2] ~  truncated_gaussian(mean = polymod[age, age2],
                                               std = (noise_switch == 0 ? 0 : polymod_sd[age, age2]), 
                                               lower = 0)

      
      // Population
      N <- S + H + L + P + E + T_E + T_P
      
      //All-cause mortality excluding TB
      mu_all[age] ~ truncated_gaussian(mean = exp_life_span[age],
                                       std = (noise_switch == 0 ? 0 : exp_life_span[age]*0.05), 
                                       lower = 0)

      
      mu[age] <- 1 / mu_all[age] - ((age < 3 ? mu_t_0_14 : (age < 11 ?  mu_t_15_69 : mu_t_70_89)) *
        (P[0, age] + P[1, age] + 
        E[0, age] + E[1, age] + T_E[0, age] + T_E[1, age]+ T_P[0, age] + 
        T_P[1, age])) / 
        (N[0, age] + N[1, age])
      mu[age] <-(mu[age] < 0 ? 0 : mu[age])
      
      //Average rate out of pulmonary
      avg_rate_from_pulmonary <- ((mu[0] + mu_t_0_14 + scaled_nu_p_0_14) * (P[0, 0] + P[1, 0]) + 
        (mu[1] + mu_t_0_14 + scaled_nu_p_0_14) * (P[0, 1] + P[1, 1]) +
        (mu[2] + mu_t_0_14 + scaled_nu_p_0_14) * (P[0, 2] + P[1, 2]) +
        (mu[3] + mu_t_15_69 + scaled_nu_p_15_89) * (P[0, 3] + P[1, 3]) +
        (mu[4] + mu_t_15_69 + scaled_nu_p_15_89) * (P[0, 4] + P[1, 4]) +
        (mu[5] + mu_t_15_69 + scaled_nu_p_15_89) * (P[0, 5] + P[1, 5]) +
        (mu[6] + mu_t_15_69 + scaled_nu_p_15_89) * (P[0, 6] + P[1, 6]) +
        (mu[7] + mu_t_15_69 + scaled_nu_p_15_89) * (P[0, 7] + P[1, 7]) +
        (mu[8] + mu_t_15_69 + scaled_nu_p_15_89) * (P[0, 8] + P[1, 8]) +
        (mu[9] + mu_t_15_69 + scaled_nu_p_15_89) * (P[0, 9] + P[1, 9]) +
        (mu[10] + mu_t_15_69 + scaled_nu_p_15_89) * (P[0, 10] + P[1, 10]) +
        (mu[11] + mu_t_70_89 + scaled_nu_p_15_89) * (P[0, 11] + P[1, 11])) / (
            (P[0, 1] + P[1, 1]) +
              (P[0, 2] + P[1, 2]) +
              (P[0, 3] + P[1, 3]) +
              (P[0, 4] + P[1, 4]) +
              (P[0, 5] + P[1, 5]) +
              (P[0, 6] + P[1, 6]) +
              (P[0, 7] + P[1, 7]) +
              (P[0, 8] + P[1, 8]) +
              (P[0, 9] + P[1, 9]) +
              (P[0, 10] + P[1, 10]) +
              (P[0, 11] + P[1, 11])
        )

      // Estimate the number of nonuk born cases
      inline nuk_start = (29 * yscale) //Start introducing non-UK born cases from 1960
      inline nuk_data = (69 * yscale)  //Start using data from 2000
      EstNUKCases[age] <- (t_now <  nuk_data ?  
                                (t_now < nuk_start ? 0 : 
                                   (NonUKScaling == 0 ? 0 : 
                                   (exp((t_now - nuk_start) / (log(2) * NonUKScaling)) - 1) / 
                                     (exp((nuk_data - nuk_start) / (log(2) * NonUKScaling)) - 1) *
                                     NUKCases2000[age]
                                         //Fit cases increase between exponential increase and log constrained increase
                                   )) : 
                             NonUKBornPCases[age])
      
      NoiseNUKCases[age] ~ truncated_gaussian(mean =  EstNUKCases[age] / MeasError, 
                                              std = (noise_switch == 0 ? 0 : MeasStd * EstNUKCases[age]  / MeasError), 
                                              lower = 0)
      
      // Estimate force of infection - start with probability of transmission
      inline decay_start = 4 * yscale //Assuming historic contacts start decaying from 1935.
      inline decay_end = 49 * yscale //Assuming that modern infection starts from 1980.
      inline actual_contact = (t_now > decay_end ?  c_eff : 
                                          (t_now < decay_start ? c_hist :
                                          c_eff + (c_hist - c_eff) * pow(1/2, (t_now - decay_start) / c_hist_half)
                                                                                   )
                                          )
      

      //Now build force of infection
      foi <- transpose(P) * I_bcg
      foi[age] <- (age < 3 ? rho_0_14 : (age < 11 ?  rho_15_69 : rho_70_89))
      * (foi[age] + 
        (age < 3 ? M * M_child : (age < 11 ? M : M * M_older_adult)) * NoiseNUKCases[age] 
           / (age < 3 ? scaled_nu_p_0_14 : scaled_nu_p_15_89))
      foi <- CSample * foi
      // i.e beta * foi / N
      foi[age] <- (avg_rate_from_pulmonary *  actual_contact / avg_contacts) * 
        foi[age] / (N[0, age] + N[1, age])
        
      // Account for age based adjustment if present.
      foi[age = 0:2] <- foi[age] * beta_child
      foi[11] <-  foi[11] * beta_older_adult
      
      //Births
      // All used to fix births to deaths for testing (uncomment for this functionality)
      //death_sum[age] <-  N[1, age] + N[0, age]) / mu_all[age]
     // death_sum <- inclusive_scan(death_sum)
      births ~ gaussian(mean = births_input, std = (noise_switch == 1 ? 0 : 0.05 * births_input))
      //Use to fix births to deaths
      //births <- (const_pop == 0 ? (noise_switch == 1 ? births_sample : births_input) : death_sum[e_age - 1] + theta[e_age - 1] * (N[0, e_age - 1] +  N[1, e_age - 1])) 

          
        ode(alg="RK4(3)", h = 1.0, atoler = 1e-3, rtoler = 1e-3) { 
            // Model equations
            dS[bcg, age]/dt = (
            // Disease model updates
            - (1 - (bcg == 1 ? chi[age] : 0)) * foi[age] * S[bcg, age]
            // Demographic model updates
            + (age == 0 ? (bcg == 1 ? gamma[age] * births : (1 - gamma[age]) * births) : 0) //Births
            + (age == 0 ? 0 : (bcg == 1 ? 
            (gamma[age] * theta[age - 1] *  S[0, age - 1] + theta[age - 1] *  S[1, age - 1]) :
                                 (1 - gamma[age]) * theta[age - 1] *  S[0, age - 1])) //Ageing into bucket
            - theta[age] * S[bcg, age] //Ageing out of bucket
            - mu[age] * S[bcg, age] //All cause (excluding TB) mortality
            ) 
   
            dH[bcg, age]/dt = (
            // Disease model updates
            + (1 - (bcg == 1 ? chi[age] : 0)) * foi[age] * S[bcg, age] 
            + (1 - delta) * foi[age] * L[bcg, age] 
            - (1 - (bcg == 1 ? alpha[age] : 0)) * 
            (age < 1 ? epsilon_h_0_4 : (age < 3 ?  epsilon_h_5_14 : epsilon_h_15_89)) * H[bcg, age] 
            - (age < 1 ?  kappa_0_4 : (age < 3 ?  kappa_5_14 :  kappa_15_89)) * H[bcg, age]
            // Demographic model updates
            + (age == 0 ? 0 : theta[age - 1] *  H[bcg, age - 1]) //Ageing into bucket
            - theta[age] * H[bcg, age] //Ageing out of bucket
            - mu[age] * H[bcg, age] //All cause (excluding TB) mortality
            ) 
          
            dL[bcg, age]/dt = (
            // Disease model updates
            + (age < 1 ?  kappa_0_4 : (age < 3 ?  kappa_5_14 :  kappa_15_89)) * H[bcg, age]
            - (1 - delta) * foi[age] * L[bcg, age] 
            - (1 - (bcg == 1 ? alpha[age] : 0)) * 
            (age < 1 ? epsilon_l_0_4 : (age < 3 ?  epsilon_l_5_14 : 
                                          (age < 11 ? epsilon_l_15_69 : epsilon_l_70_89)
                                          )) * L[bcg, age] 
            +  (age < 3 ? phi_0_14 : (age < 11 ?  phi_15_69 : phi_70_89)) * (T_E[bcg, age] + T_P[bcg, age])
            // Demographic model updates
            + (age == 0 ? 0 : theta[age - 1] *  L[bcg, age - 1]) //Ageing into bucket
            - theta[age] * L[bcg, age] //Ageing out of bucket
            - mu[age] * L[bcg, age] //All cause (excluding TB) mortality
            ) 
            
            dP[bcg, age]/dt = (
            // Disease model updates
            +  (age < 3 ? Upsilon_0_14 : (age < 11 ?  Upsilon_15_69 : Upsilon_70_89)) *
              (
                  (1 - (bcg == 1 ? alpha[age] : 0)) *  
                  (age < 1 ? epsilon_h_0_4 : (age < 3 ?  epsilon_h_5_14 :   epsilon_h_15_89)) * H[bcg, age] 
            +     (1 - (bcg == 1 ? alpha[age] : 0)) *
              (age < 1 ? epsilon_l_0_4 : (age < 3 ?  epsilon_l_5_14 : 
                                            (age < 11 ? epsilon_l_15_69 : epsilon_l_70_89)
              )) * L[bcg, age]
          ) 
            + (age < 3 ? zeta_0_14 : (age < 11 ?  zeta_15_69 : zeta_70_89)) * T_P[bcg, age] 
            - (age < 3 ? scaled_nu_p_0_14 : scaled_nu_p_15_89) * P[bcg, age]  
            - (age < 3 ? mu_t_0_14 : (age < 11 ?  mu_t_15_69 : mu_t_70_89)) * P[bcg, age]
            // Demographic model updates
            + (age == 0 ? 0 : theta[age - 1] *  P[bcg, age - 1]) //Ageing into bucket
            - theta[age] * P[bcg, age] //Ageing out of bucket
            - mu[age] * P[bcg, age] //All cause (excluding TB) mortality
            )
            
            dE[bcg, age]/dt = (
              // Disease model updates
              + (1 - (age < 3 ? Upsilon_0_14 : (age < 11 ?  Upsilon_15_69 : Upsilon_70_89))) *
                (
                  (1 - (bcg == 1 ? alpha[age] : 0)) * 
                  (age < 1 ? epsilon_h_0_4 : (age < 3 ?  epsilon_h_5_14 :   epsilon_h_15_89)) * H[bcg, age] 
            +     (1 - (bcg == 1 ? alpha[age] : 0)) * 
              (age < 1 ? epsilon_l_0_4 : (age < 3 ?  epsilon_l_5_14 : 
                                            (age < 11 ? epsilon_l_15_69 : epsilon_l_70_89)
              )) * L[bcg, age]
                )
            + (age < 3 ? zeta_0_14 : (age < 11 ?  zeta_15_69 : zeta_70_89)) * T_E[bcg, age]
            - (age < 3 ? scaled_nu_e_0_14 : scaled_nu_e_15_89) * E[bcg, age]
            - (age < 3 ? mu_t_0_14 : (age < 11 ?  mu_t_15_69 : mu_t_70_89)) * E[bcg, age]
            // Demographic model updates
            + (age == 0 ? 0 : theta[age - 1] *  E[bcg, age - 1]) //Ageing into bucket
            - theta[age] * E[bcg, age] //Ageing out of bucket
            - mu[age] * E[bcg, age] //All cause (excluding TB) mortality
            ) 
            
            dT_E[bcg, age]/dt = (
              // Disease model updates
            + (age < 3 ? scaled_nu_e_0_14 : scaled_nu_e_15_89) * E[bcg, age]
            - (age < 3 ? zeta_0_14 : (age < 11 ?  zeta_15_69 : zeta_70_89)) * T_E[bcg, age] 
            - (age < 3 ? phi_0_14 : (age < 11 ?  phi_15_69 : phi_70_89)) * T_E[bcg, age] 
            - (age < 3 ? mu_t_0_14 : (age < 11 ?  mu_t_15_69 : mu_t_70_89)) * T_E[bcg, age]
            // Demographic model updates
            + (age == 0 ? 0 : theta[age - 1] *  T_E[bcg, age - 1]) //Ageing into bucket
            - theta[age] * T_E[bcg, age] //Ageing out of bucket
            - mu[age] * T_E[bcg, age] //All cause (excluding TB) mortality
            ) 
            
            dT_P[bcg, age]/dt = (
            // Disease model updates
            + (age < 3 ? scaled_nu_p_0_14 : scaled_nu_p_15_89) * P[bcg, age]
            - (age < 3 ? zeta_0_14 : (age < 11 ?  zeta_15_69 : zeta_70_89)) * T_P[bcg, age] 
            - (age < 3 ? phi_0_14 : (age < 11 ?  phi_15_69 : phi_70_89)) * T_P[bcg, age] 
            - (age < 3 ? mu_t_0_14 : (age < 11 ?  mu_t_15_69 : mu_t_70_89)) * T_P[bcg, age]
            // Demographic model updates
            + (age == 0 ? 0 : theta[age - 1] *  T_P[bcg, age - 1]) //Ageing into bucket
            - theta[age] * T_P[bcg, age] //Ageing out of bucket
            - mu[age] * T_P[bcg, age] //All cause (excluding TB) mortality
            ) 
            //Accumalator states
            dYearlyPulCases[bcg, age]/dt = ((age < 3 ? scaled_nu_p_0_14 : scaled_nu_p_15_89) * P[bcg, age])
            dYearlyEPulCases[bcg, age]/dt = ((age < 3 ? scaled_nu_e_0_14 : scaled_nu_e_15_89) * E[bcg, age]) 
            dYearlyDeaths[bcg, age]/dt = (age < 3 ? mu_t_0_14 : (age < 11 ?  mu_t_15_69 : mu_t_70_89)) *
              (T_E[bcg, age] + T_P[bcg, age] + P[bcg, age] + E[bcg, age])
            
          }

      // Reporting states
      //By year all summarised reporting states
      YearlyPCases <- YearlyPulCases[0, 0] + YearlyPulCases[0, 1] + YearlyPulCases[0, 2] + 
        YearlyPulCases[0, 3] + YearlyPulCases[0, 4] + YearlyPulCases[0, 5] + 
        YearlyPulCases[0, 6] + YearlyPulCases[0, 7] + YearlyPulCases[0, 8] + 
        YearlyPulCases[0, 9] + YearlyPulCases[0, 10] + YearlyPulCases[0, 11] +
        YearlyPulCases[1, 0] + YearlyPulCases[1, 1] + YearlyPulCases[1, 2] + 
        YearlyPulCases[1, 3] + YearlyPulCases[1, 4] + YearlyPulCases[1, 5] + 
        YearlyPulCases[1, 6] + YearlyPulCases[1, 7] + YearlyPulCases[1, 8] + 
        YearlyPulCases[1, 9] + YearlyPulCases[1, 10] + YearlyPulCases[1, 11]
      
      YearlyECases <- YearlyEPulCases[0, 0] + YearlyEPulCases[0, 1] + YearlyEPulCases[0, 2] + 
        YearlyEPulCases[0, 3] + YearlyEPulCases[0, 4] + YearlyEPulCases[0, 5] + 
        YearlyEPulCases[0, 6] + YearlyEPulCases[0, 7] + YearlyEPulCases[0, 8] + 
        YearlyEPulCases[0, 9] + YearlyEPulCases[0, 10] + YearlyEPulCases[0, 11] +
        YearlyEPulCases[1, 0] + YearlyEPulCases[1, 1] + YearlyEPulCases[1, 2] + 
        YearlyEPulCases[1, 3] + YearlyEPulCases[1, 4] + YearlyEPulCases[1, 5] + 
        YearlyEPulCases[1, 6] + YearlyEPulCases[1, 7] + YearlyEPulCases[1, 8] + 
        YearlyEPulCases[1, 9] + YearlyEPulCases[1, 10] + YearlyEPulCases[1, 11]
      
      YearlyAgeCases[age] <- YearlyPulCases[0, age] + YearlyPulCases[1, age] +
        YearlyEPulCases[0, age] + YearlyEPulCases[1, age]
      
    }
    
    sub observation {
      
      YearlyHistPInc ~ truncated_gaussian(MeasError * (YearlyPCases + 
        NoiseNUKCases[0] + NoiseNUKCases[1] + NoiseNUKCases[2] + 
        NoiseNUKCases[3] + NoiseNUKCases[4] + NoiseNUKCases[5] + 
        NoiseNUKCases[6] + NoiseNUKCases[7] + NoiseNUKCases[8] + 
        NoiseNUKCases[9] + NoiseNUKCases[10] + NoiseNUKCases[11]),
                                          2 * MeasStd * (YearlyPCases +
                                            NoiseNUKCases[0] + NoiseNUKCases[1] + NoiseNUKCases[2] + 
                                            NoiseNUKCases[3] + NoiseNUKCases[4] + NoiseNUKCases[5] + 
                                            NoiseNUKCases[6] + NoiseNUKCases[7] + NoiseNUKCases[8] + 
                                            NoiseNUKCases[9] + NoiseNUKCases[10] + NoiseNUKCases[11]),
                                            0)
      YearlyInc ~ truncated_gaussian(MeasError * (YearlyPCases + YearlyECases),
                                     MeasStd *  (YearlyPCases + YearlyECases),
                                     0)
      YearlyAgeInc[age] ~ truncated_gaussian(MeasError * YearlyAgeCases[age], 
                                             MeasStd * YearlyAgeCases[age],
                                             0)
      YearlyChildInc ~ truncated_gaussian(MeasError * (YearlyAgeCases[0] + YearlyAgeCases[1] + YearlyAgeCases[2]), 
                                          MeasStd * (YearlyAgeCases[0] + YearlyAgeCases[1] + YearlyAgeCases[2]), 
                                          0)
      YearlyAdultInc ~ truncated_gaussian(MeasError * (YearlyAgeCases[3] + YearlyAgeCases[4] + YearlyAgeCases[5]
                                    + YearlyAgeCases[6] + YearlyAgeCases[7] + YearlyAgeCases[8]
                                    + YearlyAgeCases[9] + YearlyAgeCases[10]),
                                           MeasStd * (YearlyAgeCases[3] + YearlyAgeCases[4] + YearlyAgeCases[5]
                                                 + YearlyAgeCases[6] + YearlyAgeCases[7] + YearlyAgeCases[8]
                                                 + YearlyAgeCases[9] + YearlyAgeCases[10]),
                                          0)
      YearlyOlderAdultInc ~ truncated_gaussian(MeasError * YearlyAgeCases[11], 
                                               MeasStd * YearlyAgeCases[11], 
                                               0)
      
    }
}
