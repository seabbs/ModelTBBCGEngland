#' Plot a Subset of States From a LiBBi Model
#'
#' @description Plot a subset of the states from a LiBBi model and summarising multiple runs. 
#' Multi dimension states can be plotted using stratification and facetting.
#' @param libbi_model  LiBBi model object
#' @param states A character vector containing the complete names of the variables to plot.
#' @param burn_in Numeric, indicating the burn in period.
#' @param scales A character string indicating the axis scaling to use for facets. Defaults to 
#' "free_y".
#' @param plot_uncert Logical, defaults to \code{TRUE}. Should simulation uncertainty be plotted.
#' @param plot_data Logical, defaults to \code{TRUE}. Should the summarised data be plotted.
#' @return A plot of the specified states.
#' @export
#' @import ggplot2
#' @importFrom rbi bi_read
#' @importFrom scales comma
#' @importFrom dplyr group_by mutate summarise ungroup filter mutate_at funs
#' @importFrom purrr map map_dfr
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @examples
#' 
#' ## Show function code
#' plot_state
plot_state <- function(libbi_model = NULL, 
                       states = NULL, 
                       summarise = FALSE,
                       summarise_by = NULL,
                       strat_var = "state",
                       burn_in = 0,
                       scales = "free_y",
                       plot_uncert = TRUE,
                       plot_data = TRUE) {
  
  ## Read in data
  data <- bi_read(libbi_model)
  
  ## Filter out for required states/states
  data <- data[states] %>% 
    map(as_tibble) %>% 
    map(~filter(., time > 0))
  
  summarise_state <- function(df, summarise = summarise, summarise_by = summarise_by) {
    if (summarise) {
      df <- df %>% 
        group_by(.dots = c(summarise_by, "np", "time")) %>% 
        summarise(value = sum(value, na.rm = TRUE)) %>% 
        ungroup
    }
  ##  df <- df %>% 
   ##   group_by(.dots = setdiff(colnames(df), c("np", "value")))
    return(df)
  }
  
  ## Summarise by vectorisation if required.
    data <- map_dfr(data, summarise_state, 
                summarise = summarise, summarise_by = summarise_by, 
                .id = "state")
    
    if (!(strat_var %in% colnames(data)) || is.null(strat_var)) {
      strat_var <- "state"
    }
    
   ## Filter for burn in
    data <- data %>% 
      filter(time >= burn_in)
    
   ## Summarise model runs
    sum_data <- data %>% 
      group_by(.dots = setdiff(colnames(data), c("np", "value"))) %>% 
      summarise(
        mean = mean(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        ll = quantile(value, 0.25, na.rm = TRUE),
        lll = quantile(value, 0.025, na.rm = TRUE),
        hh = quantile(value, 0.75, na.rm = TRUE),
        hhh = quantile(value, 0.975, na.rm = TRUE)) %>% 
      ungroup %>% 
      gather("Average", "Count", mean, median) %>%
      mutate_at(.vars = c(strat_var), .funs = funs(as.factor(.)))
   
    if(plot_data) {
      ## Plot model runs 
      plot <- sum_data %>% 
        ggplot(aes_string(x = "time", y = "Count", linetype = "Average", col = strat_var, fill = strat_var)) +
        geom_line(size = 1.2, alpha = 0.8)
      
      if (plot_uncert) {
        plot <- plot + 
          geom_ribbon(aes(ymin = lll, ymax = hhh, col = NULL), alpha = 0.1) +
          geom_ribbon(aes(ymin = ll, ymax = hh, col = NULL), alpha = 0.2)
          
      }
      
      plot <- plot +
        scale_fill_viridis_d(end = 0.8) +
        scale_color_viridis_d(end = 0.8) +
        scale_y_continuous(labels = comma) + 
        theme_minimal() +
        theme(legend.position = "top") +
        labs(x = "Time")
      
      facet_var <- setdiff(colnames(sum_data), c("Average", "Count", "lll", "ll", "hh", "hhh", "time", strat_var))
      
      if (length(facet_var) != 0) {
        plot <- plot +
          facet_wrap(facet_var, scales = scales)
      }
    }else{
      plot <- sum_data
    }

    return(plot)
  
}
  

