library(ravedash)
# global variables for the module

# Stores global variables. These are required
module_id <- "fragility_explorer"
pipeline <- raveio::pipeline(
  pipeline_name = "fragility_explorer",
  settings_file = "settings.yaml",
  paths = "./modules")
debug <- TRUE

#' Function to check whether data is loaded.
#' @param first_time whether this function is run for the first time
#' @details The function will be called whenever \code{data_changed} event is
#' triggered. This function should only return either \code{TRUE} or
#' \code{FALSE} indicating the check results. If \code{TRUE} is returned,
#' \code{module_html} will be called, and module 'UI' should be displayed.
#' If \code{FALSE} is returned, \code{open_loader} event will be dispatched,
#' resulting in calling function \code{loader_html}.
#' @return Logical variable of length one.
check_data_loaded <- function(first_time = FALSE){

  # Make sure uncomment this in production
  # if( first_time ) { return(FALSE) }

  repo <- pipeline$read("repository")

  if( inherits(repo, "rave_prepare_subject_voltage_with_epoch") ) {
    ravedash::fire_rave_event('loader_message', repo$subject$subject_id)
    return( TRUE )
  } else {
    ravedash::fire_rave_event('loader_message', NULL)
    return(FALSE)
  }

}



# ----------- Initial configurations -----------

# Change the logger level when `debug` is enabled
if(exists('debug', inherits = FALSE) && isTRUE(get('debug'))){
  ravedash::logger_threshold("trace", module_id = module_id)
} else {
  ravedash::logger_threshold("info", module_id = module_id)
}



