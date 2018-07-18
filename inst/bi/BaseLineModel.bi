/**
 * Baseline TB and BCG vaccination model
 */
model Baseline {

  // Model dimensions
  const e_bcg = 1 // 0 = unvaccinated, 1 = vaccinated
  const e_age = 1 // 0,..14= 5 year age groups (i.e 0-4), and 15 = 70-89
  
  dim bcg(e_bcg)
  dim age(e_age)
  
  //Age at vaccination
  const age_at_vac = 1 // Age group that are vaccinated
  const e_d_of_p = 1 // Duration of protection (must be at least 1)
  
  dim d_of_p(e_d_of_p)
  
  // Placeholders
  
  // Contact matrix
  param C[age, age]
  
  //Disease model parameters
  
  // Age specific protection from infection conferred by BCG vaccination
  param chi[age]
  
  // Protection from infection due to prior latent infection
  param delta
  
  // Protection from active disease due to BCG vaccination
  param alpha[age]
  
  // Transition from high risk latent disease to active disease
  param epsilon_h[age]
  
  // Transition to low risk latent disease from high risk latent disease
  param kappa[age]
  
  // Transition from low risk latent disease to active disease
  param epsilon_l[age]
  
  // Rate of succesful treatment completion
  param phi[age]
  
  // Proportion of cases that have pulmonary TB
  param Upsilon[age]
  
  // Proportion of cases that have pulmonary smear postive TB
  param rho[age]
  
  // Rate of starting treatment - pulmonary/extra-pulmonary
  param nu_p[age]
  param nu_e[age]
  
  // Rate loss to follow up - pulmonary/extra-pulmonary
  param zeta_p[age]
  param zeta_e[age]
  
  // Rate of TB mortality
  param mu_p[age]
  param mu_e[age]
  
  //BCG vaccination parameters
  
  //Demographic model parameters
  
  //Noise parameters
  
  //Population states
  state S[bcg, age] // susceptible
  state H[bcg, age] // high risk latent
  state L[bcg, age] // low risk latent
  state P[bcg, age] // pulmonary TB
  state E[bcg, age] // extra-pulmonary TB only
  state T_P[bcg, age] // pulmonary TB on treatment
  state T_E[bcg, age] // extra-pulmonary TB on treatment
  
  //Calculation states

  //Accumalator states
  state PulCases[bcg, age] // monthly pulmonary cases starting treatment
  state EPulCases[bcg, age] // monthly extra-pulmonary cases starting treatment
  state YearlyPulCases[bcg, age] // yearly pulmonary cases starting treatment
  state YearlyEPulCases[bcg, age] // yearly extra-pulmonary cases starting treatment
  state PulDeaths[bcg, age] // Pulmonary TB deaths
  state EPulDeaths[bcg, age] // Extra-pulmonary TB deaths
  state YearlyPulDeaths[bcg, age] // Pulmonary TB deaths (yearly)
  state YearlyEPulDeaths[bcg, age] // Extra-pulmonary TB deaths (yearly)
  
  //Observations
  obs YearlyInc // Yearly overall incidence
  obs AgeInc // Invidence stratified by age
  
  sub parameter {

    // Place hold polymod values
    C[age, age] <- 1
    
    // Parameter scales
    inline dscale = 365.25/12
    inline mscale = 1
    inline yscale = 1/12
    
    // Priors for BCG vaccination + transforms
    alpha[age] ~ 1 - exp(gaussian(mean = -1.86, std = 0.22)) // Previous log transformed
    chi[age] ~ truncated_gaussian(mean = 0.185, std = 0.0536, lower = 0, upper = 1)
    
    //Disease priors
    delta ~ truncated_gaussian(mean = 0.78, std = 0.0408, lower = 0, upper = 1)
    epsilon_h[age] ~ truncated_gaussian(mean = 0.00695 * dscale, std =  0.0013 * dscale, lower = 0)
    kappa[age] ~ truncated_gaussian(mean = 0.0133 * dscale, std = 0.00242 * dscale, lower = 0)
    epsilon_l[age] ~ truncated_gaussian(mean = 0.000008 * dscale, std =  0.00000408 * dscale, lower = 0)
    phi[age] ~ pow(gamma(shape = 9.86, scale = 16.28) * yscale, -1)
    // Proportion of TB cases with pulmonary TB
    Upsilon[age] ~ truncated_gaussian(mean = 0.629, std = 0.0101, lower = 0, upper = 1)
    // Propotion of pulmonary TB cases that are smear positive
    rho[age] ~ truncated_gaussian(mean = 0.302, std = 0.0189, lower = 0, upper = 1)
    // Rate of starting treatment - pulmonary/extra-pulmonary
    nu_p[age] ~ pow(gamma(shape = 0.986, scale = 4.57) * yscale, -1)
    nu_e[age] ~ pow(gamma(shape = 0.697, scale = 2.21) * yscale, -1)
    // Rate loss to follow up - pulmonary/extra-pulmonary
    zeta_p[age] ~ truncated_gaussian(mean = 0.00807 * yscale, std = 0.0298 * yscale, lower = 0)
    zeta_e[age] ~ truncated_gaussian(mean = 0.00807 * yscale, std = 0.0298 * yscale, lower = 0)
    // Rate of TB mortality
    mu_p[age] ~ truncated_gaussian(mean = 0.00413 * yscale, std = 0.0227 * yscale, lower = 0)
    mu_e[age] ~ truncated_gaussian(mean = 0.00363 * yscale, std = 0.0301 * yscale, lower = 0)
      
  }

  sub initial {
    S[bcg, age] <- 1000 // susceptible
    H[bcg, age] <- 0 // high risk latent
    L[bcg, age] <- 0 // low risk latent
    P[bcg, age] <- 10 // pulmonary TB
    E[bcg, age] <- 0 // extra-pulmonary TB only
    T_P[bcg, age] <- 0// pulmonary TB on treatment
    T_E[bcg, age] <- 0 // extra-pulmonary TB on treatment
  }

  sub transition {


    // Reset accumalator variables
    PulCases[bcg, age] <- 0
    EPulCases[bcg, age] <- 0
    PulDeaths[bcg, age] <- 0
    EPulDeaths[bcg, age] <- 0
    YearlyPulCases[bcg, age] <- (t_now % 12 == 0 ? 0 : YearlyPulCases[bcg, age])
    YearlyEPulCases[bcg, age] <- (t_now % 12 == 0 ? 0 : YearlyEPulCases[bcg, age])
    YearlyPulDeaths[bcg, age] <- (t_now % 12 == 0 ? 0 : YearlyPulDeaths[bcg, age])
    YearlyEPulDeaths[bcg, age] <- (t_now % 12 == 0 ? 0 : YearlyEPulDeaths[bcg, age])
    
    ode {
      
       // Force of infection
       inline N[bcg, age] <- S[bcg, age] + H[bcg, age] + L[bcg, age] + P[bcg, age] + E[bcg, age] + T_P[bcg, age] + T_E[bcg, age]
       
       inline lambda[age] = beta[age] / inclusive_scan(inclusive_scan(N[,age])) * inclusive_scan(rho * C[age,] * (inclusive_scan(P[,age] + 1))
       
      //Disease model updates
      
      inline nl_S[bcg, age] = (1 - (bcg == 1 ? chi[age] : 0)) * lambda[age] * S[bcg, age]
      inline nl_L[bcg, age] = (1 - delta) * lambda[age] * L[bcg, age]
      inline hlc[bcg, age] = (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] H[bcg, age]
      inline hlll[bcg, age] = kappa[bcg] * H[bcg, age]
      inline llc[bcg, age] = (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
      inline tp_l[bcg, age] = phi[age] * T_P[bcg, age]
      inline te_l[bcg, age] = phi[age] * T_E[bcg, age]
      inline ncp[bcg, age] = Upsilon[age] * (hlc[bcg, age] +  llc[bcg, age])
      inline nce[bcg, age] = (1 - Upsilon[age]) * (hlc[bcg, age] +  llc[bcg, age])
      inline tpp[bcg, age] = zeta_p[age] * T_P[bcg, age]
      inline ptp[bcg, age] = nu_p[age] * P[bcg, age]
      inline tee[bcg, age] = zeta_e[age] * E_P[bcg, age]
      inline ete[bcg, age] = nu_e[age] * E[bcg, age]
                                 
      
      // Disease model equations
      
      inline S_d[bcg, age] = -nl_S[bcg, age]
      inline H_d[bcg, age] = +nl_S[bcg, age] + nl_L[bcg, age] - hlc[bcg, age] - hlll[bcg, age]
      inline L_d[bcg, age] = +hlll[bcg, age] - nl_L[bcg, age] - llc[bcg, age] + tp_l[bcg, age] + te_l[bcg, age]
      inline P_d[bcg, age] = +ncp[bcg, age] + tpp[bcg, age] - ptp[bcg, age] - mu_p[age] * P[bcg, age]
      inline E_d[bcg, age] = +nce[bcg, age] + tee[bcg, age] - ete[bcg, age] -  mu_e[age] * E[bcg, age]
      inline T_P_d[bcg, age] = ptp[bcg, age] - tpp[bcg, age] - tp_l[bcg, age] - mu_p[age] * T_P[bcg, age]
      inline T_E_d[bcg, age] = ete[bcg, age] - tee[bcg, age] - te_l[bcg, age] - mu_e[age] * T_E[bcg, age]
      
      // Demographic model updates
      
      
      // Demographic model equations
      
      //Update final states
      dS[bcg, age]/dt = S_d[bcg, age]
      dH[bcg, age]/dt = H_d[bcg, age]
      dL[bcg, age]/dt = L_d[bcg, age]
      dP[bcg, age]/dt = P_d[bcg, age]
      dE[bcg, age]/dt = E_d[bcg, age]
      dT_P[bcg, age]/dt = T_P_d[bcg, age]
      dT_E[bcg, age]/dt = T_E_d[bcg, age]
      
      //Accumalator states
      dPulCases[bcg, age]/dt = ptp[bcg, age]
      dYearlyPulCases[bcg, age]/dt = ptp[bcg, age]
      dEPulCases[bcg, age]/dt = ete[bcg, age]
      dYearlyEPulCases[bcg, age]/dt = ete[bcg, age]
      dPulDeaths[bcg, age]/dt = mu_p[age] * T_P[bcg, age] + mu_p[age] * P[bcg, age]
      dEPulDeaths[bcg, age]/dt =  mu_e[age] * T_E[bcg, age] + mu_e[age] * E[bcg, age]
      dYearlyPulDeaths[bcg, age]/dt = mu_p[age] * T_P[bcg, age] + mu_p[age] * P[bcg, age]
      dYearlyEPulDeaths[bcg, age]/dt = mu_e[age] * T_E[bcg, age] + mu_e[age] * E[bcg, age]

    }
  }

  sub observation {
   
   YearlyInc ~ poisson(rate = inclusive_scan(YearlyPulCases) + inclusive_scan(YearlyEPulCases))
   AgeInc[age] ~ poisson(rate = inclusive_scan(YearlyPulCases[,age]) + inclusive_scan(YearlyEPulCases[,age]))
    
  }

}