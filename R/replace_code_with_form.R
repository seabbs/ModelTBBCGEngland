#' Replace Coded Parameter names with latex.
#'
#' @param param A vector of coded parameters.
#'
#' @return A vector of parameters formated with latex
#' @export
#' @importFrom stringr str_replace
#' @examples
#' 
#' 
replace_code_with_form <- function(param) {
  param %>%
    str_replace("beta_young_adult", "$\\\\beta_{\\\\text{young adult}}$") %>%
    str_replace("M_young_adult", "$M_{\\\\text{young adult}}$") %>% 
    str_replace("M", "$M$") %>%
    str_replace("\\$M\\$_", "M_") %>% 
    str_replace("meas_error", "$E_{\\\\text{syst}}$") %>% 
    str_replace("meas_std", "$E_{\\\\text{noise}}$") %>% 
    str_replace("c_eff", "$c_{\\\\text{eff}}$") %>% 
    str_replace("c_hist_half", "$c^{\\\\text{hist}}_{\\\\text{half}}$") %>% 
    str_replace("c_hist", "$c^{\\\\text{hist}}_{\\\\text{eff}}$") %>% 
    str_replace("non_uk_scaling", "$\\\\iota_{\\\\text{scale}}$") %>% 
    str_replace("older_adult_activation_scaling", "$\\\\epsilon^{\\\\text{older-adult}}_L$")
}