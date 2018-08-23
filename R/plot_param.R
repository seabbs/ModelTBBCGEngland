#' Plot a Subset of Parameters From a LiBBi Model
#'
#' @description Plot either the distribution or trace of a subset of parameters split out by their dimension
#' @param libbi_model  LiBBi model object
#' @param parameter A character vector containing the complete names of the variables to plot.
#' @param log_scale Logical, defaults to \code{TRUE}. Should the x axis be plotted on a log scale?
#' @param sqrt_scale Logical, defaults to \code{TRUE}. Should the x axis be plotted on a square root scale?
#' @param plot_type A character string indicating the type of plot to show. \code{"dist"} plots a density and rug plot, whilst \code{"trace"}
#' plots a trace plot.
#' @param plot_data Logical, defaults to \code{TRUE}. Should the summarised data be plotted.
#' @return A plot of the specified parameters.
#' @export
#' @import ggplot2
#' @import rbi.helpers
#' @importFrom scales comma
#' @importFrom dplyr group_by mutate
#' @examples
#' 
#' ##Show function code
#' plot_param
plot_param <- function(libbi_model = NULL, parameter = NULL,
                       log_scale = FALSE, sqrt_scale = FALSE, 
                       plot_type = "dist", plot_data = TRUE) {
  

    p <- suppressWarnings(plot(libbi_model, param = parameter, type = "param", plot = FALSE))
    
    data <- p$data$params %>% 
      group_by(np, distribution, parameter) %>% 
      mutate(length = 1:n())

 if (plot_data) {
   if (plot_type == "dist") {
     plot <- ggplot(data, aes(x = value, col = distribution, fill = distribution)) +
       geom_density(alpha = 0.6) +
       geom_rug() 
   }else if (plot_type =="trace") {
     plot <- ggplot(data, aes(x = np, y = value, col = distribution)) +
       geom_line(alpha = 0.8) 
   }
   
   plot <- plot +
     facet_grid(length~parameter) +
     theme_minimal() +
     theme(legend.position = "top") +
     scale_x_continuous(labels = comma)
   
   if (log_scale) {
     plot <- plot +
       scale_x_log10(labels = comma)
   }
   
   if (sqrt_scale) {
     plot <- plot +
       scale_x_sqrt(labels = comma)
   }
 }else{
   plot <- data
 }
 
 return(plot)
}