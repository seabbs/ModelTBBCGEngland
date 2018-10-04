## dplyr arrange pull mutate
## ggplot2 
#' Plot output sensitivty to parameters
#'
#' @param sensitivity A dataframe of sensitivity data as produced by \code{test_sensitivity}.
#'
#' @return A plot summarising model sensitivity to parameter changes.
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr arrange pull mutate bind_rows 
#' 
#' @examples
#' 
#' ## Code
#' plot_sensitivity
plot_sensitivity <- function(sensitivity = NULL) {
  
  ranked_params <- sensitivity %>% 
    arrange(abs(gamma)) %>% 
    pull(Parameter)
  
  sens <- sensitivity %>% 
    mutate(Correlation = ifelse(gamma < 0, "Negative", "Positive")) %>% 
    mutate(Parameter = Parameter %>% 
             factor(levels = ranked_params))
  
  
  plot <- sens %>% 
    ggplot(aes(x = abs(gamma), y = Parameter, col = Correlation)) +
    geom_point(size = 3) +
    geom_line(data = sens %>% 
                bind_rows(sens %>% 
                            mutate(gamma = 0)),
              aes(group = Parameter), size = 1.2, alpha = 0.6) +
    scale_x_continuous(minor_breaks = NULL) +
    scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = "Correlation") 
  
  return(plot)
}
