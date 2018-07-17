/**
 * Baseline TB and BCG vaccination model
 */
model Basline {

  // Model dimensions
  const e_bcg = 1 // 0 = unvaccinated, 1 = vaccinated
  const e_age = 1 // 0,..14= 5 year age groups (i.e 0-4), and 15 = 70-89
  
  dim bcg(e_bcg)
  dim age(e_age)
  
  //Age at vaccination
  const age_at_vac = 1 // Age group that are vaccinated
  const e_d_of_p = 1 // Duration of protection (must be at least 1)
  
  dim d_of_p(e_d_of_p)
  
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
  
  //Accumalator states
  state Inc[bcg, age] // monthly incidence
  state IncYearly[bcg, age] / yearly incidence

  //Observations
  obs Cases[bcg, age] // monthly cases starting treatment
  obs YearlyCases[bcg, age] // yearly cases starting treatment
  obs PulDeaths[bcg, age] // Pulmonary TB deaths
  obs EPulDeaths[bcg, age] // Extra-pulmonary TB deaths
  
  sub parameter {

    // Parameter scales
    inline dscale = 365.25/12
    inline mscale = 1
    inline yscale = 1/12
    
    // Priors for BCG vaccination + transforms
    alpha[age] ~ 1 - exp(gaussian(mean = -1.86, sd = 0.22)) // Previous log transformed
    chi[age] ~ truncated_gaussian(mean = 0.185, std = 0.0536, lower = 0, upper = 1)
    
    //Disease priors
    delta ~ truncated_gaussian(mean = 0.78, std = 0.0408, lower = 0, upper = 1)
    epsilon_h[age] ~ truncated_gaussian(mean = 0.00695 * dscale, sd =  0.0013 * dscale, lower = 0)
    kappa[age] ~ truncated_gaussian(mean = 0.0133 * dscale, sd = 0.00242 * dscale, lower = 0)
    epsilon_l[age] ~ truncated_gaussian(mean = 0.00695 * dscale, sd =  0.0013 * dscale, lower = 0)
    
      
      
  }

  sub initial {

  }

  sub transition {

    \\ Update the force of infection
    inline lambda[age] = 

    ode {
      

      //Disease model updates
      
      inline nl_S[bcg, age] = (1 - (bcg == 1 ? chi[age] : 0)) * lambda[age] * S[bcg, age]
      inline nl_L[bcg, age] = (1 - delta) * lambda[age] * L[bcg, age]
      inline hlc[bcg, age] = (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] H[bcg, age]
      inline hlll[bcg, age] = kappa[bcg] * H[bcg, age]
      inline llc[bcg, age] = (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
      inline tp_l[bcg, age] = phi_p[age] * T_P[bcg, age]
      inline te_l[bcg, age] = phi_e[age] * T_P[bcg, age]
                  
                                 
      
      // Disease model equations
      
      inline S_d[bcg, age] = -nl_S[bcg, age]
      inline H_d[bcg, age] = +nl_S[bcg, age] + nl_L[bcg, age] - hlc[bcg, age] - hlll[bcg, age]
      inline L_d[bcg, age] = +hlll[bcg, age] - nl_L[bcg, age] - llc[bcg, age] + tp_l[bcg, age] + te_l[bcg, age]
      inline P_d[bcg, age] = 
      inline E_d[bcg, age] = 
      inline T_E_d[bcg, age] = 
      inline T_P_d[bcg, age] = 
      
      //Update final states
      dS[bcg, age]/dt = S_d[bcg, age]
      dH[bcg, age]/dt = H_d[bcg, age]
      dL[bcg, age]/dt = L_d[bcg, age]
      dP[bcg, age]/dt = P_d[bcg, age]
      dE[bcg, age]/dt = E_d[bcg, age]
      dT_P[bcg, age]/dt = T_E_d[bcg, age]
      dT_E[bcg, age]/dt = T_P_d[bcg, age]
      
      //Accumalator states
      dInc[bcg, age]/dt =
      dIncYearly[bcg, age]/dt  = 

    }
  }

  sub observation {
   
  }


}