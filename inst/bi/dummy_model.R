sub proposal_parameter {
  
  //Disease priors
  //M ~ truncated_gaussian(M, 0.005, 0, 1)
  //c_eff ~ truncated_gaussian(c_eff, 0.05, 0, 5)
  //c_hist ~ truncated_gaussian(c_hist, 0.05, 10, 15)
  
  //Modification of transmission probability by age.
  //beta_child_mod ~ truncated_gaussian(beta_child_mod, 0.05, 0)
  //beta_older_adult_mod ~ truncated_gaussian(beta_child_mod, 0.05, 0)
  
  //Historic measurement error
  //HistMeasError ~ truncated_gaussian(HistMeasError, 0.01, 0)
}