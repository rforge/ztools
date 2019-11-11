model_selection_MEV = function(fitted_par, covariables, plot_station_distr = FALSE) {
  # Information ---------------------------------------------------------------------------------------------------
  # output: a list with 'fitted_par', 'covariables' and 'models'
  #           (see output details below)
  #
  # --- this function finds the best linear models according to Akaike Information Criterion (AIC)
  #
  # --- input:
  #           1. 'fitted_par': a dataframe with the MEV parameters w,C,n as columns
  #                 each row corresponds to one station,
  #           2. 'covariables': a dataframe with the covariables for each station
  #                 each row corresponds to one station,
  #                 columns should include at least 'lon' (longitude), 'lat' (latitude) and 'alt' (altitude)
  # --- optional input:
  #           3.  'plot_station_distribution': logical value
  #                 if TRUE, a plot with the distribution of the used stations for the model selection is generated
  #                 default (if this input is missing) is FALSE
  # --- output: a list with
  #           1. the given 'fitted_par' matrix
  #           2. the given 'covariables' matrix
  #           3. a list of lm-class objects with the best fitted linear models for the MEV parameters: 'models'
  #                 scale_model: scale ~ ...,
  #                 shape_model: shape ~ ...
  #---------------------------------------------------------------------------------------------------------------#

  # Check required packages and input parameters ------------------------------------------------------------------



  # number of stations in 'fitted_par' matrix has to be the same as in 'covariables' matrix
  if (nrow(fitted_par) != nrow(covariables)) {
    stop(sprintf(c("number of stations (%i) in 'fitted_par' matrix and number of stations (%i) \n  in",
                   " 'covariables' matrix don't match up"), nrow(fitted_par), nrow(covariables)))
  }

  # columns of 'covariables' matrix have to be named and include at least 'lon', 'lat' and 'alt'
  if (length(colnames(covariables)) == 0) {
    stop("columns of 'covariables' matrix have to be named")
  }
  if (length(colnames(covariables)) != length(unique(colnames(covariables)))) {
    stop("colnames of 'covariables' matrix have to be unique")
  }
  if (!all(c("lon","lat","alt") %in% colnames(covariables))) {
    stop("'covariables' matrix has to include at least 'lon', 'lat' and 'alt'")
  }

  # 'plot_station_distr' has to be a logical value
  if (!(is.logical(plot_station_distr))) {
    stop(sprintf("'plot_station_distr' has to be TRUE or FALSE -- '%s' is not allowed",plot_station_distr))
  }

  #---------------------------------------------------------------------------------------------------------------#

  # Perform calculations ------------------------------------------------------------------------------------------

  # plot the station distribution
  if (plot_station_distr) {
    layout(matrix(c(1,2), 1, 2, byrow = TRUE))
    plot(covariables$lon, covariables$lat, xlab = "lon", ylab = "lat")
    plot(covariables$alt, xlab = "stations", ylab = "alt")
  }



  # stepwise model selection for the GEV parameters
  n=as.numeric(fitted_par$n)
  # scale parameter
  scale        = as.numeric(fitted_par$C)
  full_scale   = glm(scale ~ . + n, data = as.data.frame(covariables))
  null_scale   = glm(scale ~ 1, data = as.data.frame(covariables))
  scale_models = list(step(full_scale, trace = FALSE),
                      step(glm(scale ~ lon + lat + alt, data = as.data.frame(covariables)),
                           scope = list(upper = full_scale, lower = null_scale), trace = FALSE),
                      step(null_scale, scope = list(upper = full_scale), direction = "forward", trace = FALSE))
  AICs         = c(scale_models[[1]]$aic, scale_models[[2]]$aic, scale_models[[3]]$aic)
  scale_model  = scale_models[[which.min(AICs)]]

  # shape parameter
  shape        = as.numeric(fitted_par$w)
  full_shape   = glm(shape ~ . + n + scale, data = as.data.frame(covariables))
  null_shape   = glm(shape ~ 1, data = as.data.frame(covariables))
  shape_models = list(step(full_shape, trace = FALSE),
                      step(glm(shape ~ lon + lat + alt, data = as.data.frame(covariables)),
                           scope = list(upper = full_shape, lower = null_shape), trace = FALSE),
                      step(null_shape, scope = list(upper = full_shape), direction = "forward", trace = FALSE))
  AICs         = c(shape_models[[1]]$aic, shape_models[[2]]$aic, shape_models[[3]]$aic)
  shape_model  = shape_models[[which.min(AICs)]]

  # update linear models to all stations

  A = as.data.frame(cbind(covariables, "scale" = scale))

  scale_model = glm(scale_model$formula, data = A)
  shape_model = glm(shape_model$formula, data = A)

  # save the models in a list
  models = list(scale_model = scale_model, shape_model = shape_model)

  #---------------------------------------------------------------------------------------------------------------#

  # Define function output ----------------------------------------------------------------------------------------

  ans = list(fitted_par=fitted_par, covariables = covariables,
             models = models)

  return(ans)

  #---------------------------------------------------------------------------------------------------------------#
}
