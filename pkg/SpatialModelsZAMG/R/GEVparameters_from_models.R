#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: GEVparameters_from_models ------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

GEVparameters_from_models = function(covariables, coefficients) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: GEVparameters_from_models(covariables = covariables, coefficients = sd_coeff)
  # output: a matrix or a vector with the GEV parameters (see output details below)
  #
  # --- this function calculates the GEV parameters according to the linear models specified with 'coefficients'
  #
  # --- input:
  #           1. 'covariables': a named matrix or vector with the covariables
  #                 each row corresponds to one location, 
  #                 columns are the corresponding covariables (names of the coefficients) from the linear models
  #                     the intercept is added seperately and thus doesn't have to be included
  #           2. 'coefficients": a list with the named coefficients of the linear models including
  #                 'loccoeff' (coefficients of location parameter) -- must be named
  #                 'scalecoeff' (coefficients of scale parameter) -- must be named
  #                 'shapecoeff' (coefficients of shape parameter) -- must be named
  # --- output:
  #           1. a matrix or vector with the GEV parameters: 'GEVparam'
  #                 each row corresponds to one location
  #                 columns are 'loc', 'scale' and 'shape'
  #---------------------------------------------------------------------------------------------------------------#
  
  # Check required packages and input parameters ------------------------------------------------------------------
  
  # 'coefficients' must be a list with names: 'loccoeff', 'scalecoeff' and 'shapecoeff'
  if (!(is.list(coefficients) && all(c("loccoeff","scalecoeff","shapecoeff") %in% names(coefficients)))) {
    stop("'coefficients' must be a list with names: 'loccoeff', 'scalecoeff' and 'shapecoeff'")
  }
  
  # 'loccoeff', 'scalecoeff' and 'shapecoeff' must be named
  if (length(names(coefficients$loccoeff)) != length(coefficients$loccoeff) || 
      length(names(coefficients$scalecoeff)) != length(coefficients$scalecoeff) ||
      length(names(coefficients$shapecoeff)) != length(coefficients$shapecoeff)) {
    stop("'loccoeff', 'scalecoeff' and 'shapecoeff' must be named")
  }
  
  # check whether 'covariables' is a vector -- if it is, transform it to a matrix
  if (is.vector(covariables)) {
    covariables = as.matrix(t(covariables))
  }
  
  # columns of 'covariables' have to be named and include coefficient-names of models
  if (length(colnames(covariables)) != ncol(covariables)) {
    stop("columns of 'covariables' have to be named")
  }
  if (!(all(names(c(coefficients$loccoeff[-1],coefficients$scalecoeff[-1],coefficients$shapecoeff[-1])) %in% 
            c(colnames(covariables),"loc","scale")))) {
    stop("colnames of 'covariables' have to include coefficient-names of models")
  }
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Perform calculations ------------------------------------------------------------------------------------------
  
  COV   = cbind("(Intercept)" = rep(1,times = nrow(covariables)),covariables)
  loc   = COV[,names(coefficients$loccoeff)]%*%coefficients$loccoeff
  COV   = cbind(COV,"loc" = as.numeric(loc))
  scale = abs(COV[,names(coefficients$scalecoeff)]%*%coefficients$scalecoeff)
  COV   = cbind(COV,"scale" = as.numeric(scale))
  shape = COV[,names(coefficients$shapecoeff)]%*%coefficients$shapecoeff
  COV   = cbind(COV,"shape" = as.numeric(shape))
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Define function output ----------------------------------------------------------------------------------------
  
  GEVparam = COV[,c("loc","scale","shape")]
  
  return(GEVparam)
  
  #---------------------------------------------------------------------------------------------------------------#
}