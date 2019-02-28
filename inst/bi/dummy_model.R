//ode(alg='RK4', h=1) {
  ode { 
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
      - (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
      - kappa[age] * H[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  H[bcg, age - 1]) //Ageing into bucket
      - theta[age] * H[bcg, age] //Ageing out of bucket
      - mu[age] * H[bcg, age] //All cause (excluding TB) mortality
    ) 
    
    dL[bcg, age]/dt = (
      // Disease model updates
      + kappa[age] * H[bcg, age]
      - (1 - delta) * foi[age] * L[bcg, age] 
      - (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age] 
      + phi[age] * (T_E[bcg, age] + T_P[bcg, age])
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  L[bcg, age - 1]) //Ageing into bucket
      - theta[age] * L[bcg, age] //Ageing out of bucket
      - mu[age] * L[bcg, age] //All cause (excluding TB) mortality
    ) 
    
    dP[bcg, age]/dt = (
      // Disease model updates
      + Upsilon[age] * (
        (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
      + (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
      )   
      + zeta[age] * T_P[bcg, age] 
      - nu_p[age] * P[bcg, age]  
      - mu_t[age] * P[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  P[bcg, age - 1]) //Ageing into bucket
      - theta[age] * P[bcg, age] //Ageing out of bucket
      - mu[age] * P[bcg, age] //All cause (excluding TB) mortality
    )
    
    dE[bcg, age]/dt = (
      // Disease model updates
      + (1 - Upsilon[age]) * (
        (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_h[age] * H[bcg, age] 
    +   (1 - (bcg == 1 ? alpha[age] : 0)) * epsilon_l[age] * L[bcg, age]
      ) 
      + zeta[age] * T_E[bcg, age]
      - nu_e[age] * E[bcg, age]
      - mu_t[age] * E[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  E[bcg, age - 1]) //Ageing into bucket
      - theta[age] * E[bcg, age] //Ageing out of bucket
      - mu[age] * E[bcg, age] //All cause (excluding TB) mortality
    ) 
    
    dT_E[bcg, age]/dt = (
      // Disease model updates
      + nu_e[age] * E[bcg, age]
      - zeta[age] * T_E[bcg, age] 
      - phi[age] * T_E[bcg, age] 
      - mu_t[age] * T_E[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  T_E[bcg, age - 1]) //Ageing into bucket
      - theta[age] * T_E[bcg, age] //Ageing out of bucket
      - mu[age] * T_E[bcg, age] //All cause (excluding TB) mortality
    ) 
    
    dT_P[bcg, age]/dt = (
      // Disease model updates
      + nu_p[age] * P[bcg, age]
      - zeta[age] * T_P[bcg, age] 
      - phi[age] * T_P[bcg, age] 
      - mu_t[age] * T_P[bcg, age]
      // Demographic model updates
      + (age == 0 ? 0 : theta[age - 1] *  T_P[bcg, age - 1]) //Ageing into bucket
      - theta[age] * T_P[bcg, age] //Ageing out of bucket
      - mu[age] * T_P[bcg, age] //All cause (excluding TB) mortality
    ) 
    //Accumalator states
    dYearlyPulCases[bcg, age]/dt = (nu_p[age] * P[bcg, age])
    dYearlyEPulCases[bcg, age]/dt = (nu_e[age] * E[bcg, age]) 
    dYearlyDeaths[bcg, age]/dt = (mu_t[age] * (T_E[bcg, age] + T_P[bcg, age]) + mu_t[age] * (P[bcg, age] + E[bcg, age]))
    
  }