#' Plot a Subset of States From a LiBBi Model
#'
#' @description Plot either the distribution or trace of a subset of parameters split out by their dimension
#' @param libbi_model  LiBBi model object
#' @param state A character vector containing the complete names of the variables to plot.
#' @param log_scale Logical, defaults to \code{TRUE}. Should the x axis be plotted on a log scale?
#' @param sqrt_scale Logical, defaults to \code{TRUE}. Should the x axis be plotted on a square root scale?
#'
#' @return A plot of the specified states.
#' @export
#' @import ggplot2
#' @import rbi bi_read
#' @importFrom scales comma
#' @importFrom dplyr group_by mutate group_by summarise ungroup
#' @importFrom rbi bi_read 
#' @importFrom purrr
#' @importFrom tibble as_tibble
#' @examples
#' 
#' 
plot_state <- function(libbi_model = NULL, states = NULL, 
                       log_scale = FALSE, sqrt_scale = FALSE,
                       summarise = FALSE,
                       summarise_by = NULL) {
  
  ## Read in data
  data <- bi_read(libbi_model)
  
  ## Filter out for required states/states
  data <- data[states] %>% 
    map(as_tibble)
  
  summarise_state <- function(df, summarise = summarise, summarise_by = summarise_by) {
    if (summarise) {
      df <- df %>% 
        group_by(.dots = c(summarise_by, "np", "time")) %>% 
        summarise(value = sum(value, na.rm = TRUE)) %>% 
        ungroup
    }
    
    df <- df %>% 
      group_by(.dots = setdiff(colnames(df), c("np", "value")))
    return(df)
  }
  
    data <- map_dfr(data, summarise_state, 
                summarise = TRUE, summarise_by = summarise_by, 
                .id = "state")

  
}
  

