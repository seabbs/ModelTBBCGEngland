/**
 * Baseline TB and BCG vaccination model
 */
model Basline {

  // Model dimensions
  const e_bcg = 1 // 0 = unvaccinated, 1 = vaccinated
  const e_age = 1 // 0,..14= 5 year age groups (i.e 0-4), and 15 = 70-89
  
  dim bcg(e_bcg)
  dim age(e_age)
  
  //Disease model parameters
  
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
  obs PulDeaths[bcg,age] // Pulmonary TB deaths
  obs EPulDeaths[bcg,age] // Extra-pulmonary TB deaths
  
  sub parameter {

  }

  sub initial {

  }

  sub transition {


    ode {
      
      inline lambda[bcg,age] = 
        
      //Disease model updates
      inline S_d[bcg, age] = -(1 - chi[bcg,age])
      inline H_d[bcg, age] =   
      inline L_d[bcg, age] = 
      inline P_d[bcg, age] = 
      inline E_d[bcg, age] = 
      inline T_E_d[bcg, age] = 
      inline T_P_d[bcg, age] = 
      
      //Update states
      dS[bcg, age]/dt = S_d[bcg, age]
      dH[bcg, age]/dt = 
      dL[bcg, age]/dt = 
      dP[bcg, age]/dt = 
      dE[bcg, age]/dt = 
      dT_P[bcg, age]/dt =
      dT_E[bcg, age]/dt = 
      
      //Accumalator states
      dInc[bcg, age]/dt =
      dIncYearly[bcg, age]/dt  = 

    }
  }

  sub observation {
   
  }


}