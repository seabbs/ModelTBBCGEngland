#' @title Read a LiBbi Model
#'  
#' @description
#' Read results of a \code{LibBi} run from a folder. This completely reconstructs the saved \code{LibBi} object
#' This reads all options, files and outputs of a \code{LibBi} run from a specified folder
#'
#' @param folder Name of the folder containing the Libbi output as formated by \code{save_libbi}.
#' @param ... any extra options to pass to \code{\link{read_libbi}} when creating the new object
#' @return a \code{\link{libbi}} object#'
#' @importFrom rbi attach_data bi_write get_name
#' @importFrom purrr map
#' @importFrom stringr str_replace
#' @export
#' @examples
#' 
#'## Function code
#' ModelTBBCGEngland::read_libbi
read_libbi <- function(folder, ...) {
  if (missing(folder)) {
    stop("Need to specify a folder to read")
  }
  
  files <- list.files(folder)
  
  read_obj <- map(files, function(x) {
    if (x == "output") {
      files <- list.files(file.path(folder, x))
      file <- map(files, ~readRDS(file.path(folder, x, .)))
      names(file) <-  str_replace(files,".rds", "")
                          
    }else{
      file <- readRDS(file.path(folder, x))
      file <- file[[1]]
    }
    
    return(file)
  })
  
  names(read_obj) <- str_replace(files,".rds", "")
  
  libbi_options <- list(...)
  pass_options <- c("model", "dims", "time_dim", "coord_dims", 
                    "thin", "output-every", "init", "input", "obs")
  for (option in pass_options) {
    if (!(option %in% names(libbi_options)) && option %in% 
        names(read_obj)) {
      libbi_options[[option]] <- read_obj[[option]]
    }
  }
  if ("options" %in% names(read_obj)) {
    for (option in names(read_obj[["options"]])) {
      if (!(option %in% names(libbi_options))) {
        libbi_options[[option]] <- read_obj[["options"]][[option]]
      }
    }
  }
  new_obj <- do.call(libbi, libbi_options)
  new_obj <- attach_data(new_obj, file = "output", data = read_obj$output, 
                         time_dim = libbi_options$time_dim)
  if ("log" %in% names(read_obj)) {
    writeLines(read_obj[["log"]], new_obj$log_file_name)
  }
  new_obj$supplement <- read_obj$supplement
  
  return(new_obj)
}