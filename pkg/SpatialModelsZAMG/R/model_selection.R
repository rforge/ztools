model_selection = function(max_data, covariables, at_least_data = 30, plot_station_distr = FALSE) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: model_selection(max_data = get_data$sd_max_data, covariables = sd_covariables, at_least_data = 30)
  # output: a list with 'max_data', 'covariables', 'point_est', 'models' and 'used_for_model_selection'
  #           (see output details below)
  #
  # --- this function calculates the pointwise GEV parameter estimates for each station and
  #           finds the best linear models according to Akaike Information Criterion (AIC)
  #
  # --- input:
  #           1. 'max_data': a matrix with the yearly maxima of snow depth (sd) or snow water equivalent (swe)
  #                 each row corresponds to one station,
  #                 columns are the corresponding years
  #                 matrix might contain NA's
  #           2. 'covariables': a matrix with the covariables for each station
  #                 each row corresponds to one station,
  #                 columns should include at least 'lon' (longitude), 'lat' (latitude) and 'alt' (altitude)
  # --- optional input:
  #           3. 'at_least_data': how many measurements does each station has to have at least for model selection
  #                 default (if this input is missing) is 30
  #           4. 'plot_station_distribution': logical value
  #                 if TRUE, a plot with the distribution of the used stations for the model selection is generated
  #                 default (if this input is missing) is FALSE
  # --- output: a list with
  #           1. the given 'max_data' matrix
  #           2. the given 'covariables' matrix
  #           3. a matrix with the pointwise GEV parameter estimates: 'point_est'
  #                 each row corresponds to one station,
  #                 columns are 'loc' (location parameter), 'scale' (scale parameter) and 'shape' (shape parameter)
  #           4. a list of lm-class objects with the best fitted linear models for the GEV parameters: 'models'
  #                 loc_model:   loc   ~ ...,
  #                 scale_model: scale ~ ...,
  #                 shape_model: shape ~ ...
  #           5. the indices of the stations which were used for the model selection: 'used_for_model_selection'
  #                 you can use it for example as 'max_data[used_for_model_selection,]'
  #                     to get the maxima data of all stations used for the model selection
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

  # number of stations in 'max_data' matrix has to be the same as in 'covariables' matrix
  if (nrow(max_data) != nrow(covariables)) {
    stop(sprintf(c("number of stations (%i) in 'max_data' matrix and number of stations (%i) \n  in",
                   " 'covariables' matrix don't match up"), nrow(max_data), nrow(covariables)))
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

  # 'at_least_data' has to be a positive integer
  if (at_least_data < 0 || (at_least_data != as.integer(at_least_data))) {
    stop(sprintf("'at_least_data' has to be a positive integer -- '%s' is not allowed",at_least_data))
  }

  # 'plot_station_distr' has to be a logical value
  if (!(is.logical(plot_station_distr))) {
    stop(sprintf("'plot_station_distr' has to be TRUE or FALSE -- '%s' is not allowed",plot_station_distr))
  }

  #---------------------------------------------------------------------------------------------------------------#

  # Perform calculations ------------------------------------------------------------------------------------------

  # get all stations with at least s measurements to find the best covariables
  nrow_data = apply(!is.na(max_data),1,sum)
  s   = at_least_data
  ind = which(nrow_data >= s)
  if (length(ind) == 0){
    stop(sprintf("no station with at least %i measurements -- 'at_least_data' being %i is too high",s,s))
  }
  if (length(ind) < 100){
    if (plot_station_distr) {
      warning(sprintf(c("there are less than 100 (only %i) stations with at least %i measurements --",
                        " results might be bad \n  check station distribution in the plot"),length(ind),s))
    } else {
      warning(sprintf(c("there are less than 100 (only %i) stations with at least %i measurements --",
                        " results might be bad \n  check station distribution by setting the optional input",
                        " argument 'plot_station_distr' to TRUE"),length(ind),s))
    }
  }

  # plot the station distribution
  if (plot_station_distr) {
    layout(matrix(c(1,2), 1, 2, byrow = TRUE))
    plot(covariables[ind,"lon"], covariables[ind,"lat"], xlab = "lon", ylab = "lat")
    plot(covariables[ind,"alt"], xlab = "stations", ylab = "alt")
  }

  # calculate pointwise estimators for GEV parameters at each station
  GEVparam = matrix(NA, nrow = nrow(max_data), ncol = 3,
                    dimnames = list(rownames(max_data), c("loc","scale","shape")))
  for (k in 1:nrow(max_data)) {
    z = max_data[k,!is.na(max_data[k,])]
    GEVparam[k,] = gevmle(z)
  }

  # stepwise model selection for the GEV parameters
  # loc parameter
  loc        = as.numeric(GEVparam[ind,"loc"])
  full_loc   = glm(loc ~ ., data = as.data.frame(covariables[ind,]))
  null_loc   = glm(loc ~ 1, data = as.data.frame(covariables[ind,]))
  loc_models = list(step(full_loc, trace = FALSE),
                    step(glm(loc ~ lon + lat + alt, data = as.data.frame(covariables[ind,])),
                         scope = list(upper = full_loc, lower = null_loc), trace = FALSE),
                    step(null_loc, scope = list(upper = full_loc), direction = "forward", trace = FALSE))
  AICs       = c(loc_models[[1]]$aic, loc_models[[2]]$aic, loc_models[[3]]$aic)
  loc_model  = loc_models[[which.min(AICs)]]

  # scale parameter
  scale        = as.numeric(GEVparam[ind,"scale"])
  full_scale   = glm(scale ~ . + loc, data = as.data.frame(covariables[ind,]))
  null_scale   = glm(scale ~ 1, data = as.data.frame(covariables[ind,]))
  scale_models = list(step(full_scale, trace = FALSE),
                      step(glm(scale ~ lon + lat + alt, data = as.data.frame(covariables[ind,])),
                           scope = list(upper = full_scale, lower = null_scale), trace = FALSE),
                      step(null_scale, scope = list(upper = full_scale), direction = "forward", trace = FALSE))
  AICs         = c(scale_models[[1]]$aic, scale_models[[2]]$aic, scale_models[[3]]$aic)
  scale_model  = scale_models[[which.min(AICs)]]

  # shape parameter
  shape        = as.numeric(GEVparam[ind,"shape"])
  full_shape   = glm(shape ~ . + loc + scale, data = as.data.frame(covariables[ind,]))
  null_shape   = glm(shape ~ 1, data = as.data.frame(covariables[ind,]))
  shape_models = list(step(full_shape, trace = FALSE),
                      step(glm(shape ~ lon + lat + alt, data = as.data.frame(covariables[ind,])),
                           scope = list(upper = full_shape, lower = null_shape), trace = FALSE),
                      step(null_shape, scope = list(upper = full_shape), direction = "forward", trace = FALSE))
  AICs         = c(shape_models[[1]]$aic, shape_models[[2]]$aic, shape_models[[3]]$aic)
  shape_model  = shape_models[[which.min(AICs)]]

  # update linear models to all stations
  loc   = as.numeric(GEVparam[,"loc"])
  scale = as.numeric(GEVparam[,"scale"])
  shape = as.numeric(GEVparam[,"shape"])

  A = as.data.frame(cbind(covariables,"loc" = loc, "scale" = scale))

  loc_model   = glm(loc_model$formula,   data = A)
  scale_model = glm(scale_model$formula, data = A)
  shape_model = glm(shape_model$formula, data = A)

  # save the models in a list
  models = list(loc_model = loc_model, scale_model = scale_model, shape_model = shape_model)

  #---------------------------------------------------------------------------------------------------------------#

  # Define function output ----------------------------------------------------------------------------------------

  ans = list(max_data = max_data, covariables = covariables, point_est = GEVparam,
             models = models, used_for_model_selection = ind)

  return(ans)

  #---------------------------------------------------------------------------------------------------------------#
}
