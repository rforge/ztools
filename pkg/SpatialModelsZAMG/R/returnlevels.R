#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: returnlevels -------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

returnlevels = function(GEVparam, q) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: returnlevels(GEVparam = sd_GEVparam,  q = 100)
  # output: a vector with the return levels for each location
  #
  # --- this function calculates the return levels according to the given GEV parameters 'GEVparam' and
  #           a return period of 'q' years
  #
  # --- input:
  #           1. 'GEVparam': a named matrix or vector with the GEV parameters
  #                 each row corresponds to one location
  #                 columns are 'loc', 'scale' and 'shape'
  #           2. 'q': the return period for the calculation of return levels -- must be a number greater than 1
  # --- output:
  #           1. a vector with the return levels for each location: 'rl'
  #---------------------------------------------------------------------------------------------------------------#
  
  # Check required packages and input parameters ------------------------------------------------------------------
  
  # # load required package 'SpatialExtremes'
  # if (inherits(try(library(SpatialExtremes, warn.conflicts = FALSE, quietly = TRUE), 
  #                  silent = TRUE), "try-error")) {
  #   message("required package 'SpatialExtremes' is not installed yet -- trying to install package")
  #   install.packages("SpatialExtremes", quiet = TRUE)
  #   if (inherits(try(library(SpatialExtremes, warn.conflicts = FALSE, quietly = TRUE), 
  #                    silent = TRUE), "try-error")) {
  #     stop("package 'SpatialExtremes' couldn't be installed")
  #   } else {
  #     message("package successfully installed and loaded")
  #   }
  # }
  
  # check whether 'GEVparam' is a vector -- if it is, transform it to a matrix
  if (is.vector(GEVparam)) {
    GEVparam = as.matrix(t(GEVparam))
  }
  
  # columns of 'GEVparam' must be 'loc', 'scale' and 'shape'
  if (!(all(c("loc","scale","shape") %in% colnames(GEVparam)))) {
    stop("columns/names of 'GEVparam' must be 'loc', 'scale' and 'shape'")
  }
  
  # 'q' has to be a number greater than 1
  if (q <= 1) {
    stop(sprintf("'q' has to be a number greater than 1 -- '%s' is not allowed",q))
  }
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Perform calculations ------------------------------------------------------------------------------------------
  
  # predefine return level vector 'rl'
  rl = rep(NA,times = nrow(GEVparam))
  
  # calculate return levels
  p  = 1/q
  for (k in 1:nrow(GEVparam)) {
    rl[k] = qgev(1-p, loc = GEVparam[k,"loc"], scale = GEVparam[k,"scale"], 
                 shape = GEVparam[k,"shape"])
  }
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Define function output ----------------------------------------------------------------------------------------
  
  return(rl)
  
  #---------------------------------------------------------------------------------------------------------------#
}