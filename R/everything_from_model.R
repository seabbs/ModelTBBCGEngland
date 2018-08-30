#' Adapt LibBi model with all Variables in Output
#'
#' @description Some \code{rbi.helper} functions require all variables to have outputs. This function modifies a LibBi model 
#' so that this is the case.
#' @param model A LibBi model, loaded using \code{bi_model}.
#'
#' @return A LibBi model.
#' @export
#'
#' @importFrom rbi replace_all
#' @examples
#' library(rbi)
#' 
#' model <- system.file(package="rbi", "SIR.bi") # get full file name from package
#' model <- bi_model(model) # load model
#' model <- replace_all(model, "noise n_recovery", "noise n_recovery\\(has_output = 0\\)")
#' 
#' model
#' 
#' everything_from_model(model)
everything_from_model <- function(model) {
  
 model <- replace_all(model, "\\(has_output = 0, has_input = 0\\)", "")
 model <- replace_all(model, "\\(has_input = 0, has_output = 0\\)", "")
 model <- replace_all(model, "\\(has_output = 0\\)", "")
 
 return(model)
}