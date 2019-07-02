#' Model Summary Plots
#' 
#' @param model A LibBi model object
#' @param uncertainty Logical, defaults to \code{TRUE}. 
#' @param prior A libBi model object. Containing priors to be compared to. Defaults to \code{NULL}
#'
#' @return A list of model summary plots
#' @export
#' @seealso plot_param plot_state
#' @examples
#' 
model_plots <- function(model, uncertainty = TRUE, prior = NULL) {
  
  plots <- list()
  
  plots[[1]] <- ModelTBBCGEngland::plot_param(model, prior_params = prior)
  
  plots[[2]] <- ModelTBBCGEngland::plot_state(model, summarise = TRUE, 
                                              plot_uncert = uncertainty)
  
  plots[[3]] <- ModelTBBCGEngland::plot_state(model, summarise = TRUE, 
                                              start_time = 40, 
                                              plot_uncert = uncertainty)
  
  plots[[4]] <- ModelTBBCGEngland::plot_state(model, summarise = TRUE,
                                              start_time = 60, 
                                              plot_uncert = uncertainty)
  
  plots[[5]] <-ModelTBBCGEngland::plot_state(model, summarise = FALSE,
                                             states = "YearlyAgeInc",
                                             start_time = 60,
                                             plot_uncert = uncertainty)
  
  plots[[6]] <- ModelTBBCGEngland::plot_state(model, summarise = FALSE, 
                                              states = "YearlyAgeInc",
                                              start_time = 20,
                                              plot_uncert = uncertainty)
  
  plots[[7]] <- ModelTBBCGEngland::plot_state(prior, summarise = FALSE, 
                                              states = "L", start_time = 10,
                                              plot_uncert = uncertainty)
  
  return(plots)
}