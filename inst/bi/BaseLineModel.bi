/**
 * Baseline TB and BCG vaccination model
 */
model Baseline {
  
  // Model dimensions
  const e_bcg = 2 // 0 = unvaccinated, 1 = vaccinated
  const e_age = 1 // 0,..14= 5 year age groups (i.e 0-4), and 15 = 70-89
  
  dim bcg(e_bcg)
  dim age(e_age)
  dim age2(e_age)
  
  //Age at vaccination
  const age_at_vac = 1 // Age group that are vaccinated
  const e_d_of_p = 1 // Duration of protection (must be at least 1)
  
  dim d_of_p(e_d_of_p)
  
  // Placeholders
  
  // Contact matrix
  param C[age, age2]
  
  
  //Disease model parameters
  
  // Effective contact rate
  param c_eff
  
  // Historic effective contact rate
  param c_hist
  
  // Protection from infection due to prior latent infection
  param delta
  
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
  param avg_nu_p[age]
  
  // Rate loss to follow up - pulmonary/extra-pulmonary
  param zeta_p[age]
  param zeta_e[age]
  
  // Rate of TB mortality
  param mu_p[age]
  param mu_e[age]
  
  //BCG vaccination parameters
  
  // Age specific protection from infection conferred by BCG vaccination
  param chi[age]
  
  // Protection from active disease due to BCG vaccination
  param alpha[age]
  
  //Demographic model parameters
  
  //Noise parameters
  noise CNoise[age, age2] // Sampled contact rate
  
  // Time varying parameter states
  state CSample[age, age2] // Sampled contact rate (symmetric)
  state TotalContacts[age] // Average number of contacts (across age groups)
  state beta[age] //Probability of transmission
  state foi[age] // force of infection
  
  //Calculation parameters
  param I_age[age]
  param I_bcg[bcg]
  
  //Population states
  state S[bcg, age] // susceptible
  state H[bcg, age] // high risk latent
  state L[bcg, age] // low risk latent
  state P[bcg, age] // pulmonary TB
  state E[bcg, age] // extra-pulmonary TB only
  state T_P[bcg, age] // pulmonary TB on treatment
  state T_E[bcg, age] // extra-pulmonary TB on treatment
  state N[bcg, age] // Overall population
  state NSum[age] // Sum of population (vector but repeating values)
  
  //Accumalator states
  state PulCases[bcg, age] // monthly pulmonary cases starting treatment
  state EPulCases[bcg, age] // monthly extra-pulmonary cases starting treatment
  state YearlyPulCases[bcg, age] // yearly pulmonary cases starting treatment
  state YearlyEPulCases[bcg, age] // yearly extra-pulmonary cases starting treatment
  state PulDeaths[bcg, age] // Pulmonary TB deaths
  state EPulDeaths[bcg, age] // Extra-pulmonary TB deaths
  state YearlyPulDeaths[bcg, age] // Pulmonary TB deaths (yearly)
  state YearlyEPulDeaths[bcg, age] // Extra-pulmonary TB deaths (yearly)
  
  // Reporting states
  state YearlyAgeCases[age](has_input = 0, has_output = 0)  
    state YearlyCasesCumSum[age](has_input = 0, has_output = 0)
      state YearlyCases
      
      
      //Observations
      obs YearlyInc // Yearly overall incidence
      
      sub parameter {
        
        // Place hold polymod values
        C[age, age2] <- 1000
        
        // Parameter scales
        inline dscale = 12/365.25
        inline mscale = 1
        inline yscale = 12
        
        // Priors for BCG vaccination + transforms
        alpha[age] ~ gaussian(mean = -1.86, std = 0.22)
        alpha[age] <- 1 - exp(alpha[age]) // Previous log transformed
        chi[age] ~ truncated_gaussian(mean = 0.185, std = 0.0536, lower = 0, upper = 1)
        
        //Disease priors
        c_eff ~ uniform(0, 5)
        c_hist ~ uniform(10, 15)
        delta ~ truncated_gaussian(mean = 0.78, std = 0.0408, lower = 0, upper = 1)
        epsilon_h[age] ~ truncated_gaussian(mean = 0.00695 * dscale, std =  0.0013 * dscale, lower = 0)
        kappa[age] ~ truncated_gaussian(mean = 0.0133 * dscale, std = 0.00242 * dscale, lower = 0)
        epsilon_l[age] ~ truncated_gaussian(mean = 0.000008 * dscale, std =  0.00000408 * dscale, lower = 0)
        phi[age] ~ inverse_gamma(shape = 9.86, scale = 16.28)
        phi[age] <- phi[age] * yscale
        // Proportion of TB cases with pulmonary TB
        Upsilon[age] ~ truncated_gaussian(mean = 0.629, std = 0.0101, lower = 0, upper = 1)
        // Propotion of pulmonary TB cases that are smear positive
        rho[age] ~ truncated_gaussian(mean = 0.302, std = 0.0189, lower = 0, upper = 1)
        // Rate of starting treatment - pulmonary/extra-pulmonary
        nu_p[age] ~ inverse_gamma(shape = 0.986, scale = 4.57)
        nu_p[age] <- nu_p[age] * yscale
        nu_e[age] ~ inverse_gamma(shape = 0.697, scale = 2.21)
        nu_e[age] <- nu_e[age] * yscale
        // Rate loss to follow up - pulmonary/extra-pulmonary
        zeta_p[age] ~ truncated_gaussian(mean = 0.00807 * yscale, std = 0.0298 * yscale, lower = 0)
        zeta_e[age] ~ truncated_gaussian(mean = 0.00807 * yscale, std = 0.0298 * yscale, lower = 0)
        // Rate of TB mortality
        mu_p[age] ~ truncated_gaussian(mean = 0.00413 * yscale, std = 0.0227 * yscale, lower = 0)
        mu_e[age] ~ truncated_gaussian(mean = 0.00363 * yscale, std = 0.0301 * yscale, lower = 0)
        
        
        //Calculation parameters
        I_age <- 1
        I_bcg <- 1
        
        // Avg period infectious for pulmonary cases
        avg_nu_p <- inclusive_scan(nu_p)
        avg_nu_p[age] <- avg_nu_p[e_age -1] / e_age
        
      }
    
    sub initial {
      S[bcg, age] <- 100000 // susceptible
      H[bcg, age] <- 1000 // high risk latent
      L[bcg, age] <- 10000 // low risk latent
      P[bcg, age] <- 100 // pulmonary TB
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
      
      //Contact rate
      CNoise[age, age2] ~  truncated_gaussian(mean = C[age, age2], std = C[age, age2]/ 10, lower = 0)
      CSample <- CNoise
      CSample[age, age2] <- CSample[age2, age]
      TotalContacts <- CSample * I_age
      TotalContacts <- inclusive_scan(TotalContacts)
      TotalContacts[age] <- TotalContacts[e_age - 1] / e_age
      
      // Population
      N <- S + H + L + P + E + T_P + T_E
      NSum[age] <- N[0, age] + N[1, age]
      NSum <- inclusive_scan(NSum)
      NSum[age] <- NSum[e_age - 1]
      
      // Estimate force of infection - start with probability of transmission
      inline curr_start = 59 * 12 // Modern day is 1990 with a baseline date of 1931
      beta <- avg_nu_p * (c_eff + c_hist * ((t_now > curr_start  ? 0 : (curr_start - t_now) / curr_start))) 
      beta <- beta ./ TotalContacts
      
      //Now build force of infection
      beta <- beta ./ NSum
      foi <- transpose(P) * I_bcg
      foi[age] <- rho[age] * foi[age]
      foi <- CSample * foi
      foi[age] <- beta[age] * foi[age]
      
      ode {
        
        // Model equations
        
        dS[bcg, age]/dt = 
          // Disease model updates
          - (1 - (bcg == 1 ? chi[age] : 0)) * foi[age] * S[bcg, age]
          // Demographic model updates
          dH[bcg, age]/dt = 
          // Disease model updates
          + (1 - (bcg == 1 ? chi[age] : 0)) * foi[age] * S[bcg, age] 
          + (1 - delta) * foi[age] * L[bcg, age] 
          - (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
          - kappa[age] * H[bcg, age]
          // Demographic model updates
          dL[bcg, age]/dt = 
          // Disease model updates
          + kappa[age] * H[bcg, age]
          - (1 - delta) * foi[age] * L[bcg, age] 
          - (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age] 
          + phi[age] * T_P[bcg, age] + phi[age] * T_E[bcg, age]
          // Demographic model updates
          dP[bcg, age]/dt = 
          // Disease model updates
          + Upsilon[age] * (
              (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
        + (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
          ) 
        + zeta_p[age] * T_P[bcg, age] 
        - nu_p[age] * P[bcg, age] 
        - mu_p[age] * P[bcg, age]
        // Demographic model updates
        dE[bcg, age]/dt = 
        // Disease model updates
        + (1 - Upsilon[age]) * (
            (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
        +   (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
        ) 
        + zeta_e[age] * T_E[bcg, age] 
        - nu_e[age] * E[bcg, age] 
        - mu_e[age] * E[bcg, age]
        // Demographic model updates
        dT_P[bcg, age]/dt = 
        // Disease model updates
        + nu_p[age] * P[bcg, age] 
        - zeta_p[age] * T_P[bcg, age] 
        - phi[age] * T_P[bcg, age] 
        - mu_p[age] * T_P[bcg, age]
        // Demographic model updates
        dT_E[bcg, age]/dt = 
        // Disease model updates
        + nu_e[age] * E[bcg, age] 
        - zeta_e[age] * T_E[bcg, age] 
        - phi[age] * T_E[bcg, age]
        - mu_e[age] * T_E[bcg, age]
        // Demographic model updates
        
        
        
        //Accumalator states
        dPulCases[bcg, age]/dt = nu_p[age] * P[bcg, age]
        dYearlyPulCases[bcg, age]/dt = nu_p[age] * P[bcg, age]
        dEPulCases[bcg, age]/dt = nu_e[age] * E[bcg, age]
        dYearlyEPulCases[bcg, age]/dt = nu_e[age] * E[bcg, age]
        dPulDeaths[bcg, age]/dt = mu_p[age] * T_P[bcg, age] + mu_p[age] * P[bcg, age]
        dEPulDeaths[bcg, age]/dt =  mu_e[age] * T_E[bcg, age] + mu_e[age] * E[bcg, age]
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
      //AgeInc[age] ~ poisson(rate = YearlyAgeCases[age])
      
    }
    
}