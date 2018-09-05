#' @title Save a LiBbi Model
#' @description
#' Write results of a \code{LibBi} run to a folder as a series of RDS files.
#' This saves all options, files and outputs of a \code{LibBi} run to a specified folder.
#'
#' @param x a \code{\link{libbi}} object
#' @param folder A character string indicating the folder name under which to save the model.
#' @param supplement any supplementary data to save
#' @importFrom rbi bi_read
#' @importFrom purrr walk2
#' @export
#' @examples 
#'
#' ## Function code 
#' ModelTBBCGEngland::save_libbi
save_libbi <- function(x, folder, supplement) {
  if (missing(folder)) {
    stop("Need to specify a folder name")
  }
  
  if (!dir.exists(folder)) {
    dir.create(folder)
    dir.create(file.path(folder, "output"))
  }
  
  save_obj <- list(model=x$model,
                   dims=x$dims,
                   time_dim=x$time_dim,
                   coord_dims=x$coord_dims,
                   thin=1,
                   supplement=x$supplement
  )
  
  
  options <- x$options
  
  for (file_type in c("init", "input", "obs")) {
    file_option <- paste(file_type, "file", sep="-")
    if (file_option %in% names(x$options)) {
      save_obj[[file_type]] <- bi_read(x, file=file_type)
      options[[file_option]] <- NULL
    }
  }
  
  save_obj[["options"]] <- options
  
  if (!missing(supplement)) save_obj[["supplement"]] <- supplement
  
  walk2(save_obj, names(save_obj), ~ saveRDS(list(.x), file.path(folder, paste0(.y, ".rds"))))
  
  output <- bi_read(x)
  
  walk2(output, names(output), ~ saveRDS(.x, file.path(folder, "output", paste0(.y, ".rds"))))
  
}