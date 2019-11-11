optimizer_smooth_model_MEV_mle = function(m_select, max_data, method = c("nlminb","BFGS","ucminf","Nelder-Mead"), follow.on = FALSE,
                                      itnmax = NULL, printParam = FALSE){


  # Information ---------------------------------------------------------------------------------------------------
  #
  # --- this function optimizes the coefficients of the best fitted linear models
  #           (from the function 'model_selection_MEV') via smooth modeling and with maximum likelihood estimation
  #
  # --- input:
  #           1. 'm_select': this input should be a list including 'max_data', 'covariables' and 'models'
  #                 as in the output of the function 'model_selection'
  #           2. `max_data`: matrix in format of max_data_aut, i.e. a matrix with the annual maximum of year j
  #               at station i at entry (i,j), meaning stations in rows, years in columns
  # --- optional input:
  #           2. 'method': optimization method(s) for external function 'optimx', this can also be a vector
  #                 possible methods are: 'Nelder-Mead','BFGS','CG','L-BFGS-B','nlm','nlminb','spg',
  #                     'ucminf','newuoa','bobyqa','nmkb','hjkb','Rcgmin','Rvmmin'
  #                 default (if this input is missing) is 'c("nlminb","BFGS","ucminf","Nelder-Mead")'
  #           3. 'follow.on': logical value
  #                 if TRUE, and there are multiple methods, then the last set of coefficients from one method
  #                     is used as the starting set for the next
  #                 default (if this input is missing) is FALSE
  #           4. 'itnmax': if provided as a vector of the same length as the length of 'method',
  #                 this gives the maximum number of iterations or function values for the corresponding method
  #                 if a single number is provided, this will be used for all methods
  #           5. 'printParam': logical value
  #                 if TRUE, the GEV parameters during the optimization are printed
  #                 this might be useful to check the proper functioning of the optimization
  #                     (shape parameter should be approximately between -0.5 and 0.5)
  #                 default (if this input is missing) is FALSE
  #           6. 'save_name': a character string defining the name of the function output to be saved
  #                 default (if this input is missing) is 'optim_smooth_model'
  # --- output: a list with
  #           1. a summary of the optimization results, including an information message whether the
  #                 optimization was successful or not and which method delivered the best coefficients: 'summary'
  #           2. a list with the optimized coefficients: 'coefficients'
  #                 containing: 'scalecoeff', 'shapecoeff' and 'all_coeff'
  #---------------------------------------------------------------------------------------------------------------#

  # Check required packages and input parameters ------------------------------------------------------------------


  # 'm_select' should be a list including 'covariables' and
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

  # define cov and models from input argument 'm_select'
  cov      = as.matrix(m_select$covariables)
  models   = m_select$models
  par_n    = m_select$fitted_par$n

  scale_model = models$scale_model
  shape_model = models$shape_model

  # define number of stations
  n.site = nrow(max_data)

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

  # calculate the GEV parameters for given coefficients
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

  # define the logarithmized likelihood function
  loglikeli = function(p) {

    # scale the coefficients back to original scaling
    p = p*parscale

    # define scale and shape coefficients
    scalecoeff = p[1:n.scale]
    shapecoeff = p[(n.scale+1):(n.scale+n.shape)]

    # calculate the GEV parameters according to the given coefficients
    mev    = dsgnmat2Param(scalecoeff,shapecoeff)
    scales = mev[,"scales"]
    shapes = mev[,"shapes"]
    if (printParam) {
      print(summary(mev))
    }

    # calculate the log likelihood
    sum = 0
    for (k in 1:n.site) {
      z       = as.numeric(max_data[k,!is.na(max_data[k,])])
      #summand = sum(log(dMEV_weibull(z, C = scales[k], w = shapes[k], n =par_n[k])))
      summand = sum(log(dmev(z, w = shapes[k], c = scales[k], n =par_n[k])))
      sum     = sum + summand
    }
    # in order to maximize we take the negative value
    return(-sum)
  }


  optloglik=loglikeli


  # perform the actual optimization
  # method 'nlminb' has an additional control parameter 'eval.max'
  if (any(method == "nlminb") && !follow.on && !missing(itnmax)) {
    if (length(method) > 1) {
      pos = which(method == "nlminb")
      optimizer1 = optimx(start, optloglik, method = method[-pos], itnmax = itnmax,
                          control = list(follow.on = follow.on, kkt = FALSE))
      optimizer2 = optimx(start, optloglik, method = method[pos], itnmax = itnmax,
                          control = list(follow.on = follow.on, eval.max = itnmax, kkt = FALSE))
      optimizer = rbind(optimizer1, optimizer2)
      optimizer = optimizer[match(method, rownames(optimizer)),]
    } else {
      optimizer = optimx(start, optloglik, method = method, itnmax = itnmax,
                         control = list(follow.on = follow.on, eval.max = itnmax, kkt = FALSE))
    }
  } else {
    optimizer = optimx(start, optloglik, method = method, itnmax = itnmax,
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
