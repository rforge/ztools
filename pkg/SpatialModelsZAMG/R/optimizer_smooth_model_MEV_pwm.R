optimizer_smooth_model_MEV_pwm = function(m_select, data, method = c("nlminb","BFGS","ucminf","Nelder-Mead"), follow.on = FALSE, 
                                      itnmax = NULL, printParam = FALSE){
  # Information ---------------------------------------------------------------------------------------------------
  # example: optimizer_smooth_model(m_select = selected_models_MEV,data=data_aut, method = "nlminb", save_name = "optim_smooth_model")
  # output: a list with 'summary' and 'coefficients' (see output details below)
  #
  # --- this function optimizes the coefficients of the best fitted linear models 
  #           (from the function 'model_selection_MEV') via probability weighted moments optimization
  #
  # --- input:
  #           1. 'm_select': this input should be a list including 'fitted_param', 'covariables' and 'models' 
  #                 as in the output of the function 'model_selection_MEV'
  #           2. 'data': list whose elements are vectors including all observed daily values at one station,
  #                   the stations have to be the same and used in the same order as the stations used for model selection
  # --- optional input:
  #           3. 'method': optimization method(s) for external function 'optimx', this can also be a vector
  #                 possible methods are: 'Nelder-Mead','BFGS','CG','L-BFGS-B','nlm','nlminb','spg',
  #                     'ucminf','newuoa','bobyqa','nmkb','hjkb','Rcgmin','Rvmmin'
  #                 default (if this input is missing) is 'c("nlminb","BFGS","ucminf","Nelder-Mead")'
  #           4. 'follow.on': logical value 
  #                 if TRUE, and there are multiple methods, then the last set of coefficients from one method 
  #                     is used as the starting set for the next
  #                 default (if this input is missing) is FALSE
  #           5. 'itnmax': if provided as a vector of the same length as the length of 'method', 
  #                 this gives the maximum number of iterations or function values for the corresponding method 
  #                 if a single number is provided, this will be used for all methods
  #           6. 'printParam': logical value
  #                 if TRUE, the GEV parameters during the optimization are printed
  #                 this might be useful to check the proper functioning of the optimization 
  #                     (shape parameter should be approximately between -0.5 and 0.5)
  #                 default (if this input is missing) is FALSE
  #           7. 'save_name': a character string defining the name of the function output to be saved
  #                 default (if this input is missing) is 'optim_smooth_model'
  # --- output: a list with
  #           1. a summary of the optimization results, including an information message whether the 
  #                 optimization was successful or not and which method delivered the best coefficients: 'summary'
  #           2. a list with the optimized coefficients: 'coefficients'
  #                 containing: 'scalecoeff', 'shapecoeff' and 'all_coeff'
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
  # 
  # # load required package 'optimx'
  # if (inherits(try(library(optimx, warn.conflicts = FALSE, quietly = TRUE), silent = TRUE), "try-error")) {
  #   message("required package 'optimx' is not installed yet -- trying to install package")
  #   install.packages("optimx", quiet = TRUE)
  #   if (inherits(try(library(optimx, warn.conflicts = FALSE, quietly = TRUE), silent = TRUE), "try-error")) {
  #     stop("package 'optimx' couldn't be installed")
  #   } else {
  #     message("package successfully installed and loaded")
  #   }
  # }
  
  # 'm_select' should be a list including 'max_data', 'covariables' and
  #     'models' as in the output of the function 'model_selection'
  if (!(is.list(m_select) && all(c("covariables","models") %in% names(m_select)))) {
    stop("please use the output of the function 'model_selection' as input 'm_select'")
  }
  
  # 'method' has to be a vector of possible methods from the external function 'optimx'
  if (!(all(method %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B","nlm","nlminb","spg",
                          "ucminf","newuoa","bobyqa","nmkb","hjkb","Rcgmin","Rvmmin")))) {
    if (length(method) == 1) {
      stop("the used optimization method is unnown -- please check possible methods of function 'optimx'")
    } else {
      x = sum(!(method %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B","nlm","nlminb","spg","ucminf",
                              "newuoa","bobyqa","nmkb","hjkb","Rcgmin","Rvmmin")))
      if (x == 1) {
        stop(c("one of the used optimization methods is unnown -- please check possible methods", 
               " of function 'optimx'"))
      } else {
        stop(sprintf(
          c("%i of the used optimization methods are unnown -- please check possible methods",
            " of function 'optimx'"),x))
      }
    }
  }
  
  # 'follow.on' has to be a logical value
  if (!(is.logical(follow.on))) {
    stop(sprintf("'follow.on' has to be TRUE or FALSE -- '%s' is not allowed",follow.on))
  }
  
  # 'itnmax' has to be a positive integer and of length 1 or of same length as 'method'
  if (!missing(itnmax)) {
    if (length(itnmax) == 1){
      if ((itnmax < 0) || (itnmax != as.integer(itnmax))) {
        stop(sprintf("'itnmax' has to be a positive integer -- '%s' is not allowed",itnmax))
      }
      if (length(method) > 1) {
        message(sprintf(c("since a single number is provided for the 'itnmax' input, this will be used",
                          " for all the %i chosen methods"),length(method)))
      }
    } else if (length(itnmax) != length(method)) {
      stop(sprintf(
        "'itnmax' has to be a vector of either the same length as 'method' (%i) or of length 1",length(method)))
    } else {
      if (any(itnmax < 0) || any(itnmax != as.integer(itnmax))) {
        stop(sprintf("'itnmax' has to be a vector (possibly of length 1) of positive integer(s)"))
      }
    }
  }
  
  # 'printParam' has to be a logical value
  if (!(is.logical(printParam))) {
    stop(sprintf("'printParam' has to be TRUE or FALSE -- '%s' is not allowed",printParam))
  }
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Perform calculations ------------------------------------------------------------------------------------------
  
  # define max_data, cov and models from input argument 'm_select'
  cov      = as.matrix(m_select$covariables)
  models   = m_select$models
  par_n    = m_select$fitted_par$n
  
  scale_model = models$scale_model
  shape_model = models$shape_model
  
  # define number of stations
  n.site = length(data)
  
  # define total number of coefficients to be optimized for each MEV parameter
  n.scale = length(scale_model$coefficients)
  n.shape = length(shape_model$coefficients)
  
  # define starting values for the optimization as the coefficients from the model selection
  start = c(scale_model$coefficients, shape_model$coefficients)
  
  # rename coefficients
  scale_names  = paste0("scale_", c("Int",names(scale_model$coefficients)[-1]), sep="")
  shape_names  = paste0("shape_", c("Int",names(shape_model$coefficients)[-1]), sep="")
  coeff_names  = c(scale_names, shape_names)
  names(start) = coeff_names
  
  # add intercept column to covariable matrix
  C = cbind("(Intercept)" = rep(1,times = n.site), cov)
  
  # calculate sample moments
  M0hat=rep(0,length(data))
  M1hat=rep(0,length(data))
  for(i in 1:length(data)){
    M0hat[i]=PWM_weibull_estimated(k=0,data=data[[i]])
    M1hat[i]=PWM_weibull_estimated(k=1,data=data[[i]])
  }
  
  
  # calculate the MEV parameters for given coefficients in every iteration
  dsgnmat2Param = function(scalecoeff, shapecoeff) {
    
    scales = as.numeric(abs(C[,names(scale_model$coefficients)]%*%scalecoeff))
    C_new  = cbind(C,"scale" = scales)
    shapes = as.numeric(C_new[,names(shape_model$coefficients)]%*%shapecoeff)
    
    ans = cbind("scales" = scales,"shapes" = shapes)
    
    return(ans)
  }
  
  # rescale the starting values for more stable optimization results
  parscale = 10^(floor(log10(abs(start))))
  start    = start/parscale
  
  # define the function to be minimized
  minfun = function(p) {
    
    # scale the coefficients back to original scaling
    p = p*parscale
    
    # define scale and shape coefficients
    scalecoeff = p[1:n.scale]
    shapecoeff = p[(n.scale+1):(n.scale+n.shape)]
    
    # calculate the MEV parameters according to the given coefficients
    mev    = dsgnmat2Param(scalecoeff,shapecoeff)
    scales = mev[,"scales"]
    shapes = mev[,"shapes"]
    if (printParam) {
      print(summary(mev))
    }
    
    # calculate analytical weibull moments
    M0=rep(0,length(data))
    M1=rep(0,length(data))
    for(i in 1:length(data)){
      M0[i]=PWM_weibull(k=0,C=scales[i],w=shapes[i])
      M1[i]=PWM_weibull(k=1,C=scales[i],w=shapes[i])
    }
    
    # calculate the function to be minimized
    sum = 0
    for (k in 1:n.site) {
      summand1 = ((M0hat[k]-M0[k])/M0hat[k])^2
      summand2 = ((M1hat[k]-M1[k])/M1hat[k])^2
      sum     = sum + summand1 + summand2
    }
    # in order to maximize we take the negative value
    return(sum)
  }
  
  
  # perform the actual optimization
  # method 'nlminb' has an additional control parameter 'eval.max'
  if (any(method == "nlminb") && !follow.on && !missing(itnmax)) {
    if (length(method) > 1) {
      pos = which(method == "nlminb")
      optimizer1 = optimx(start, minfun, method = method[-pos], itnmax = itnmax, 
                          control = list(follow.on = follow.on, kkt = FALSE))
      optimizer2 = optimx(start, minfun, method = method[pos], itnmax = itnmax, 
                          control = list(follow.on = follow.on, eval.max = itnmax, kkt = FALSE))
      optimizer = rbind(optimizer1, optimizer2)
      optimizer = optimizer[match(method, rownames(optimizer)),]
    } else {
      optimizer = optimx(start, minfun, method = method, itnmax = itnmax, 
                         control = list(follow.on = follow.on, eval.max = itnmax, kkt = FALSE))
    }
  } else {
    optimizer = optimx(start, minfun, method = method, itnmax = itnmax, 
                       control = list(follow.on = follow.on, kkt = FALSE))
  }
  
  # scale the coefficients back to original scaling
  for (l in 1:nrow(optimizer)) {
    optimizer[l,coeff_names] = optimizer[l,coeff_names]*parscale
    optimizer[l,"value"]     = -optimizer[l,"value"]
  }
  
  # check whether the optimization was successful according to the convergence codes and
  #       find the method with the maximal log likelihood
  if (any(optimizer$convcode %in% c(0,1))) {
    ind = which(optimizer[,"value"] == max(optimizer[which(optimizer$convcode %in% c(0,1)),"value"]))
    # if there are more than one method with the same (best) result choose the one with a convergence code of 0
    #     or the last one of them
    if (any(optimizer$convcode[ind] == 0)) {
      ind = ind[max(which(optimizer$convcode[ind] == 0))]
      message = sprintf(
        "optimization was successful, coefficients were chosen from method '%s'",row.names(optimizer)[ind])
    } else {
      ind = max(ind)
      warning(sprintf(c("optimization method with the best result (method '%s') reached iteration limit;",
                        "\n  coefficients were chosen from this method, but results might be bad;",
                        "\n  compare with other methods or increase the optional input argument 'itnmax'!"),
                      row.names(optimizer)[ind]))
      message = sprintf(c("optimization method with the best result (method '%s') reached iteration limit;",
                          "coefficients were chosen from this method, but results might be bad;",
                          "compare with other methods or increase the optional input argument 'itnmax'!"),
                        row.names(optimizer)[ind])
    }
    coeff = as.numeric(optimizer[ind,coeff_names])
    if (all(abs(coeff - start*parscale) < 10^-3)) {
      warning(
        "optimized coefficients and starting values are (almost) the same -- difference smaller than 10^-3")
    }
  } else {
    print(optimizer)
    stop(c("optimization was NOT successful -- try other method(s)",
           "\n  you can use the optional input argument 'method = c(...)'"))
  }
  
  names(coeff) = coeff_names
  
  # define the coefficients for each MEV parameter separately
  scalecoeff = coeff[1:n.scale]
  shapecoeff = coeff[(n.scale+1):(n.scale+n.shape)]
  
  names(scalecoeff) = names(scale_model$coefficients)
  names(shapecoeff) = names(shape_model$coefficients)
  
  # summarize optimization result and information message
  summary = list(message = message, method = row.names(optimizer)[ind], result = optimizer)
  print(summary)
  coefficients = list("scalecoeff" = scalecoeff,
                      "shapecoeff" = shapecoeff, "all_coeff" = coeff)
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Define function output ----------------------------------------------------------------------------------------
  ans = list(summary = summary, coefficients = coefficients)
  
  return(ans)
  
  #---------------------------------------------------------------------------------------------------------------#
}


#k-th sample probability weighted moment
PWM_weibull <- function(k,C,w){
  Mk=C*gamma(1+1/w)/((1+k)^(1+1/w))
  return(Mk)
}

#estimates the k-th sample probability weighted moment from vector data(k=1,2)
PWM_weibull_estimated <- function(k,data){
  
  data=sort(data)
  
  if(k==0){
    Mk=mean(data)
  } else if(k==1){
    M1hat  = 0
    N      = length(data) # sample size
    for (i in 1:(N-1)){
      M1hat   = M1hat + data[i]*(N - i)
    }
    Mk=M1hat/(N*(N-1))
  } else{
    stop(cat("function only implemented for the first two moments, i.e. k=0,1\n"))
  }
  
  return(Mk)
}