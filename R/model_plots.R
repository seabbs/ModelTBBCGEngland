#' Model Summary Plots
#' 
#' @param model A LibBi model object
#' @param uncertainty Logical, defaults to \code{TRUE}. 
#' @param prior A libBi model object. Containing priors to be compared to. Defaults to \code{NULL}
#' @param verbose Logical, defaults to \code{TRUE}. Should progress messages be printed.
#' @return A list of model summary plots
#' @export
#' @seealso plot_param plot_state
#' @examples
#' 
#' 
#' ## Code 
#' model_plots
model_plots <- function(model = NULL, uncertainty = TRUE, prior = NULL, verbose = TRUE) {
  
  plots <- list()
  
  plots[[1]] <- ModelTBBCGEngland::plot_param(model, prior_params = prior)
  
  if (verbose){
    message("Posterior plot complete")
  }
  
  plots[[2]] <- ModelTBBCGEngland::plot_state(model, summarise = TRUE, 
                                              plot_uncert = uncertainty)
  
  if (verbose){
    message("Overall trajectory plot complete")
  }
  
  plots[[3]] <- ModelTBBCGEngland::plot_state(model, summarise = TRUE, 
                                              start_time = 40, 
                                              plot_uncert = uncertainty)
  
  if (verbose){
    message("Zoomed trajectory plot complete")
  }
  
  plots[[4]] <- ModelTBBCGEngland::plot_state(model, summarise = TRUE,
                                              start_time = 60, 
                                              plot_uncert = uncertainty)
  
  if (verbose){
    message("Further zoomed trajectory plot complete")
  }
  
  plots[[5]] <-ModelTBBCGEngland::plot_state(model, summarise = FALSE,
                                             states = "YearlyAgeInc",
                                             start_time = 60,
                                             plot_uncert = uncertainty)
  
  if (verbose){
    message("Age distributed cases plot complete")
  }
  
  plots[[6]] <- ModelTBBCGEngland::plot_state(model, summarise = FALSE, 
                                              states = "YearlyAgeInc",
                                              start_time = 20,
                                              plot_uncert = uncertainty)
  
  if (verbose){
    message("Zoomed age distributed cases plot complete")
  }
  
  plots[[7]] <- ModelTBBCGEngland::plot_state(model, summarise = FALSE, 
                                              states = "L", start_time = 10,
                                              plot_uncert = uncertainty)
  
  if (verbose){
    message("latent trajectory plot complete")
  }
  
  return(plots)
}