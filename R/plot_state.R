#' Plot a Subset of States From a LiBBi Model
#'
#' @description Plot a subset of the states from a LiBBi model and summarising multiple runs. 
#' Multi dimension states can be plotted using stratification and facetting.
#' @param libbi_model  LiBBi model object (or a list of named libbi models).
#' @param model_paths A named character vector of paths to models to plot.
#' @param states A character vector containing the complete names of the variables to plot. Defaults to all states
#' @param plot_obs Logical, defaults to \code{TRUE}. Should the observed data be plotted.
#' @param obs List of dataframes containing observed data. Defaults to generating observed
#'  data using \code{setup_model_obs}.
#' @param burn_in Numeric, indicating the burn in period (samples).
#' @param start_time Numeric, indicating the time to start plotting from.
#' @param start_time_label Numeric, the label to apply to the time variable. Defaults to 1931.
#' @param end_time Numeric, defaults to the maximum value present in the data. The final time to plot.
#' @param scenarios_start Numeric, start time for which to plot alternative scenarios. Defaults to 2005.
#' @param scales A character string indicating the axis scaling to use for facets. Defaults to 
#' "free_y".
#' @param plot_uncert Logical, defaults to \code{TRUE}. Should simulation uncertainty be plotted.
#' @param plot_data Logical, defaults to \code{TRUE}. Should the summarised data be plotted.
#' @param show_mean Logical, defaults to \code{TRUE}. Should the mean value as well as the median be plotted.
#' @param use_comma_formatting Logicial, defaults to \code{TRUE}. Should comma formating be used on the y axis.
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
                       model_paths = NULL,
                       states = NULL, 
                       plot_obs = TRUE,
                       obs = NULL,
                       summarise = FALSE,
                       summarise_by = NULL,
                       strat_var = "state",
                       burn_in = 0,
                       start_time = 0,
                       end_time = NULL,
                       start_time_label = 1931,
                       scenarios_start = 2005,
                       scales = "free_y",
                       plot_uncert = TRUE,
                       plot_data = TRUE,
                       show_mean = FALSE,
                       use_comma_formatting = TRUE) {

  if (!is.null(libbi_model) && !is.null(model_paths)) {
    stop("Both a libbi model and a model path has been passed. Only one of these options is allowed.")
  }else if (!is.null(model_paths) && is.null(libbi_model)) {
    libbi_model <- model_paths
  }
  ## Summarise states - generic function
  summarise_state <- function(df, summarise = summarise, summarise_by = summarise_by) {
    if (summarise) {
      df <- df %>% 
        group_by(.dots = c(summarise_by, "time")) %>% 
        summarise(value = sum(value, na.rm = TRUE)) %>% 
        ungroup
    }
    ##  df <- df %>% 
    ##   group_by(.dots = setdiff(colnames(df), c("np", "value")))
    return(df)
  }
  
  ## Summarise model runs - generic function
  summarise_model_runs <- function(data) {
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
  }
  
  
  load_data <- function(libbi_model) {
    
    if (!is.null(model_paths)) {
      libbi_model <- read_libbi(libbi_model)
    }
    
    ## Read in data
    data <- bi_read(libbi_model, type = c("state", "obs", "noise"))
    
    ## Use all states if none set
    if (is.null(states)) {
      states <- setdiff(names(data), c("clock", "logprior", "loglikelihood"))
    }
    
    ## Filter out for required states/states
    data <- data[states] %>% 
      map(as_tibble) %>% 
      map(~filter(., time > 0,
                  time >= start_time,
                  np >= burn_in))

    if (!is.null(end_time)) {
      data <- data %>% 
        map(~filter(., time <= end_time))
    }
    
    ## Summarise by vectorisation if required.
    data <- map_dfr(data, summarise_state, 
                    summarise = summarise, summarise_by = c(summarise_by, "np"), 
                    .id = "state")
    
    if (!(strat_var %in% colnames(data)) || is.null(strat_var)) {
      strat_var <- "state"
    }
    
    sum_data <- summarise_model_runs(data) %>% 
      mutate(time = time + start_time_label)
    
    return(sum_data)
  }

  if (class(libbi_model) %in% "list") {
    sum_data <- map_dfr(libbi_model, load_data, .id = "Scenario") %>% 
      filter(any(time > scenarios_start - start_time_label, Scenario == names(libbi_model)[1])) %>% 
      mutate(Scenario = factor(Scenario, levels = names(libbi_model)))
    
    strat_var <- "Scenario"
  }else{
    sum_data <- load_data(libbi_model)
  }
  
  ## Use all states if none set
  if (is.null(states)) {
    states <- setdiff(names(sum_data), c("clock", "logprior", "loglikelihood"))
  }
  
  if (!(strat_var %in% colnames(sum_data)) || is.null(strat_var)) {
    strat_var <- "state"
  }
  
   
   ## Use observational data default if not supplied
   if (is.null(obs)) {
     obs <- ModelTBBCGEngland::setup_model_obs(years_of_age = 2000:2015, age_groups = 0:11,
                                               con_age_groups = c("children", "adults", "older adults"))
   }
   
   obs <- map_dfr(obs, summarise_state, 
                  summarise = summarise, summarise_by = summarise_by, 
                  .id = "state")
   
   obs <- obs %>% 
     filter(state %in% unique(as.character(sum_data$state))) %>% 
     filter(time >= start_time) %>% 
     mutate(bcg = NA) %>% 
     mutate(time = time + start_time_label)
   
    if (!show_mean) {
      sum_data <- sum_data %>% 
        filter(!(Average %in% "mean"))
    }
   
    if(plot_data) {
      ## Plot model runs 
      plot <- sum_data %>%
        ggplot(aes_string(x = "time", y = "Count", linetype = "Average", col = strat_var, fill = strat_var)) +
        geom_line(size = 1.2, alpha = 0.6)
      
      if (plot_obs) {
        plot <- plot + 
          geom_point(data = obs, aes(x = time, y = value,
                                     linetype = NULL,
                                     col = NULL, fill = NULL), alpha = 0.6)
      }

      
      if (plot_uncert) {
        plot <- plot + 
          geom_ribbon(aes(ymin = lll, ymax = hhh, col = NULL), alpha = 0.05) +
          geom_ribbon(aes(ymin = ll, ymax = hh, col = NULL), alpha = 0.1)
          
      }
      
      plot <- plot +
        scale_fill_viridis_d(end = 0.8) +
        scale_color_viridis_d(end = 0.8) +
        theme_minimal() +
        theme(legend.position = "top") +
        labs(x = "Time")
      
      if (use_comma_formatting) {
        plot <- plot + 
          scale_y_continuous(labels = scales::comma)
      }
      
      if (strat_var %in% "state") {
        plot <- plot + 
          guides(col = FALSE, fill = FALSE, linetype = show_mean)
      }
      
      facet_var <- setdiff(colnames(sum_data), c("Average", "Count", "lll", "ll", "hh", "hhh", "time", "Scenario"))
      
      if (length(facet_var) != 0) {
        plot <- plot +
          facet_wrap(facet_var, scales = scales)
      }
    }else{
      plot <- sum_data
    }

    return(plot)
  
}
  

