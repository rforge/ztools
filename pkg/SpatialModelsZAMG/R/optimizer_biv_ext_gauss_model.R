#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: optimizer_biv_ext_gauss_model --------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

optimizer_biv_ext_gauss_model = function(sd_m_select, swe_m_select, method = "ucminf", follow.on = FALSE,
                                         itnmax = NULL, printParam = FALSE){
  # Information ---------------------------------------------------------------------------------------------------
  # example: optimizer_biv_ext_gauss_model(sd_m_select = sd_m_select, swe_m_select = swe_m_select)
  # output: a list with 'complete_sd_max_data', 'complete_swe_max_data', 'summary' and 'coefficients'
  #           (see output details below)
  #
  # --- this function optimizes the coefficients of the best fitted linear models
  #           (from the function 'model_selection') via bivariate Extremal-Gaussian Model (= ext-t with fixed nu = 1)
  #           and with composite likelihood inference
  #
  # --- input:
  #           1. 'sd_m_select': this input should be a list including 'max_data', 'covariables' and 'models'
  #                 as in the output of the function 'model_selection'
  #           2. 'swe_m_select': this input should be a list including 'max_data', 'covariables' and 'models'
  #                 as in the output of the function 'model_selection'
  # --- optional input:
  #           3. 'method': optimization method(s) for external function 'optimx', this can also be a vector,
  #                     but optimization might take several hours for one method!
  #                 possible methods are: 'Nelder-Mead','BFGS','CG','L-BFGS-B','nlm','nlminb','spg',
  #                     'ucminf','newuoa','bobyqa','nmkb','hjkb','Rcgmin','Rvmmin'
  #                 default (if this input is missing) is 'ucminf'
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
  #                 default (if this input is missing) is 'optim_biv_ext-gauss'
  # --- output: a list with
  #           1. a completed submatrix of the 'max_data' matrix for sd: 'complete_sd_max_data'
  #                 each row corresponds to one station (unchanged)
  #                 columns (years) were chosen such that less than 50% of the original column-entries were NA's
  #           2. a completed submatrix of the 'max_data' matrix for swe: 'complete_swe_max_data'
  #                 each row corresponds to one station (unchanged)
  #                 columns (years) were chosen such that less than 50% of the original column-entries were NA's
  #           3. a summary of the optimization results, including an information message whether the
  #                 optimization was successful or not and which method delivered the best coefficients: 'summary'
  #           4. a list with the optimized coefficients: 'coefficients'
  #                 containing: 'sd_coeff': 'sd_loccoeff', 'sd_scalecoeff', 'sd_shapecoeff',
  #                             'swe_coeff': 'swe_loccoeff', 'swe_scalecoeff', 'swe_shapecoeff',
  #                             'cor_coeff': 'alpha', 'sd_kappa', 'swe_kappa', 'rho' and
  #                             'all_coeff'
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
  #
  # # load required package 'mice'
  # if (inherits(try(library(mice, warn.conflicts = FALSE, quietly = TRUE), silent = TRUE), "try-error")) {
  #   message("required package 'mice' is not installed yet -- trying to install package")
  #   install.packages("mice", quiet = TRUE)
  #   if (inherits(try(library(mice, warn.conflicts = FALSE, quietly = TRUE), silent = TRUE), "try-error")) {
  #     stop("package 'mice' couldn't be installed")
  #   } else {
  #     message("package successfully installed and loaded")
  #   }
  # }

  # 'sd_m_select' should be a list including 'max_data', 'covariables' and
  #     'models' as in the output of the function 'model_selection'
  if (!(is.list(sd_m_select) && all(c("max_data","covariables","models") %in% names(sd_m_select)))) {
    stop("please use the output of the function 'model_selection' as input 'sd_m_select'")
  }

  # 'swe_m_select' should be a list including 'max_data', 'covariables' and
  #     'models' as in the output of the function 'model_selection'
  if (!(is.list(swe_m_select) && all(c("max_data","covariables","models") %in% names(swe_m_select)))) {
    stop("please use the output of the function 'model_selection' as input 'swe_m_select'")
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

  message("optimization might take several hours")

  #---------------------------------------------------------------------------------------------------------------#

  # Perform calculations ------------------------------------------------------------------------------------------

  # define max_data, covariables and models from input arguments 'sd_m_select' and 'swe_m_select'
  z_sd  = sd_m_select$max_data
  z_swe = swe_m_select$max_data
  sd_models   = sd_m_select$models
  swe_models  = swe_m_select$models
  covariables = cbind(sd_m_select$covariables,swe_m_select$covariables)
  covariables = covariables[,unique(colnames(covariables))]

  # complete 'max_data' matrices in case NA's are present
  if (any(is.na(z_sd))) {
    nas_sd = apply(is.na(z_sd),2,sum)
    # chose only columns where there are less than 50% NA's
    col.chose_sd = which(nas_sd < 0.5*nrow(z_sd))
    z_sd = z_sd[,col.chose_sd]
    z_sd = complete(mice(data.frame(z_sd), seed = 50, printFlag = FALSE))
    z_sd = as.matrix(z_sd)
  }
  if (any(is.na(z_swe))) {
    nas_swe = apply(is.na(z_swe),2,sum)
    # chose only columns where there are less than 50% NA's
    col.chose_swe = which(nas_swe < 0.5*nrow(z_swe))
    z_swe = z_swe[,col.chose_swe]
    z_swe = complete(mice(data.frame(z_swe), seed = 50, printFlag = FALSE))
    z_swe = as.matrix(z_swe)
  }

  message(sprintf(
    "optimization is performed with %i stations and %i years -- this submatrix was chosen to be completed",
    nrow(z_sd),ncol(z_sd)))
  # warn if there are less than 30% of the years left
  if (ncol(z_sd) < 0.3*ncol(sd_m_select$max_data)) {
    warning(sprintf(
      "optimization was performed with %i stations and %i years -- this submatrix was chosen to be completed",
      nrow(z_sd),ncol(z_sd)))
  }

  # define number of stations 'n.site' and number of observations 'n.obs'
  n.site = nrow(z_sd)
  n.obs  = ncol(z_sd)

  # define linear models
  sd_loc_model    = sd_models$loc_model
  sd_scale_model  = sd_models$scale_model
  sd_shape_model  = sd_models$shape_model
  swe_loc_model   = swe_models$loc_model
  swe_scale_model = swe_models$scale_model
  swe_shape_model = swe_models$shape_model

  # update pointwise GEV parameter estimators to completed 'max_data' matrices
  gev1 = matrix(NA, nrow = n.site, ncol = 3, dimnames = list(NULL, c("locs","scales","shapes")))
  gev2 = matrix(NA, nrow = n.site, ncol = 3, dimnames = list(NULL, c("locs","scales","shapes")))
  for (k in 1:n.site) {
    gev1[k,] = gevmle(z_sd[k,])
    gev2[k,] = gevmle(z_swe[k,])
  }

  # update linear models to new GEV parameters
  # sd
  loc   = as.numeric(gev1[,"locs"])
  scale = as.numeric(gev1[,"scales"])
  shape = as.numeric(gev1[,"shapes"])

  A = as.data.frame(cbind(covariables,"loc" = loc, "scale" = scale))

  sd_loc_model   = glm(sd_loc_model$formula,   data = A)
  sd_scale_model = glm(sd_scale_model$formula, data = A)
  sd_shape_model = glm(sd_shape_model$formula, data = A)

  # swe
  loc   = as.numeric(gev2[,"locs"])
  scale = as.numeric(gev2[,"scales"])
  shape = as.numeric(gev2[,"shapes"])

  A = as.data.frame(cbind(covariables,"loc" = loc, "scale" = scale))

  swe_loc_model   = glm(swe_loc_model$formula,   data = A)
  swe_scale_model = glm(swe_scale_model$formula, data = A)
  swe_shape_model = glm(swe_shape_model$formula, data = A)

  # define the correlation function for the Extremal-Gaussian model
  matern_h = function(h,alpha,kappa) {
    2^(1-kappa)/gamma(kappa)*(h/alpha)^kappa*besselK(h/alpha, kappa)
  }

  # needed values for function loglikeli1
  z1 = as.vector(z_sd)
  z2 = as.vector(z_swe)
  T1val = 1 + rep(gev1[,"shapes"],times = n.obs)*(z1 - rep(gev1[,"locs"],times = n.obs))/
    rep(gev1[,"scales"],times = n.obs)
  T2val = 1 + rep(gev2[,"shapes"],times = n.obs)*(z2 - rep(gev2[,"locs"],times = n.obs))/
    rep(gev2[,"scales"],times = n.obs)
  pot1  = 1/rep(gev1[,"shapes"],times = n.obs)
  pot2  = 1/rep(gev2[,"shapes"],times = n.obs)
  T1  = T1val^(pot1)
  T2  = T2val^(pot2)
  dT1 = 1/rep(gev1[,"scales"],times = n.obs)*T1val^(pot1 - 1)
  dT2 = 1/rep(gev2[,"scales"],times = n.obs)*T2val^(pot2 - 1)

  # define the logarithmized likelihood function (for the correlation parameters)
  loglikeli1 = function(p) {

    # scale the coefficients back to original scaling
    p = p*parscale1

    # define the correlation parameters
    alpha     = p[1]
    sd_kappa  = p[2]
    swe_kappa = p[3]
    rho12     = p[4]

    # parameters must be in certain intervals
    if (alpha <= 0 || sd_kappa <= 0 || swe_kappa <= 0 || abs(rho12) > 1 ||
        abs(rho12) > sqrt(sd_kappa*swe_kappa)/(1/2*(sd_kappa + swe_kappa))) {
      sum = -Inf
    } else {
      # calculate log likelihood
      sum = 0
      for (j in 1:n.site) {
        h = sqrt(apply(t(t(covariables) - covariables[j,])^2,1,sum))
        if (j < n.site) {
          matern1 = matern_h(h[(j+1):n.site],alpha,sd_kappa)
          matern2 = matern_h(h[(j+1):n.site],alpha,swe_kappa)
          x1 = rep(T1[j+(0:(n.obs-1))*n.site],each = length((j+1):n.site))
          y1 = T1[rep((j+1):n.site,times = n.obs) + rep((0:(n.obs-1))*n.site,each = length((j+1):n.site))]
          x2 = rep(T2[j+(0:(n.obs-1))*n.site],each = length((j+1):n.site))
          y2 = T2[rep((j+1):n.site,times = n.obs) + rep((0:(n.obs-1))*n.site,each = length((j+1):n.site))]
          help1a = (y1/x1)
          help1b = (x1/y1)
          help2a = (y2/x2)
          help2b = (x2/y2)
          helpmat1 = (2/(1 - rep(matern1,times = n.obs)^2))^(1/2)
          helpmat2 = (2/(1 - rep(matern2,times = n.obs)^2))^(1/2)
          g1xy = helpmat1*(help1a - rep(matern1,times = n.obs))
          g1yx = helpmat1*(help1b - rep(matern1,times = n.obs))
          g2xy = helpmat2*(help2a - rep(matern2,times = n.obs))
          g2yx = helpmat2*(help2b - rep(matern2,times = n.obs))

          dx1 = rep(dT1[j+(0:(n.obs-1))*n.site],each = length((j+1):n.site))
          dy1 = dT1[rep((j+1):n.site,times = n.obs) + rep((0:(n.obs-1))*n.site,each = length((j+1):n.site))]
          helpdx1 = 1/x1*helpmat1*help1a*dt(g1xy, 2)
          helpdy1 = 1/y1*helpmat1*help1b*dt(g1yx, 2)

          dx2 = rep(dT2[j+(0:(n.obs-1))*n.site],each = length((j+1):n.site))
          dy2 = dT2[rep((j+1):n.site,times = n.obs) + rep((0:(n.obs-1))*n.site,each = length((j+1):n.site))]
          helpdx2 = 1/x2*helpmat2*help2a*dt(g2xy, 2)
          helpdy2 = 1/y2*helpmat2*help2b*dt(g2yx, 2)

          helpdd1 = gamma(3/2)*3/(sqrt(2*pi)*(1 - rep(matern1,times = n.obs)^2))
          helpdd2 = gamma(3/2)*3/(sqrt(2*pi)*(1 - rep(matern2,times = n.obs)^2))

          ptg1xy = 1/x1*pt(g1xy, 2)
          ptg1yx = 1/y1*pt(g1yx, 2)
          ptg2xy = 1/x2*pt(g2xy, 2)
          ptg2yx = 1/y2*pt(g2yx, 2)

          val1 = ((dy1/y1*(ptg1yx + helpdy1 - helpdx1))*
                    (dx1/x1*(ptg1xy + helpdx1 - helpdy1)) +
                    dx1*dy1/x1/y1*(helpdx1*2 - 1/x1*(helpdd1*(1 + g1xy^2/2)^(-5/2)*g1xy*help1a^2) +
                                     helpdy1*2 - 1/y1*(helpdd1*(1 + g1yx^2/2)^(-5/2)*g1yx*help1b^2)))

          val2 = ((dy2/y2*(ptg2yx + helpdy2 - helpdx2))*
                    (dx2/x2*(ptg2xy + helpdx2 - helpdy2)) +
                    dx2*dy2/x2/y2*(helpdx2*2 - 1/x2*(helpdd2*(1 + g2xy^2/2)^(-5/2)*g2xy*help2a^2) +
                                     helpdy2*2 - 1/y2*(helpdd2*(1 + g2yx^2/2)^(-5/2)*g2yx*help2b^2)))

          if (any(val1 <= 0) || any(val2 <= 0)) {
            sum = -Inf
            break
          } else {
            sum = sum + sum(-ptg1xy - ptg1yx + log(val1)) +
              sum(-ptg2xy - ptg2yx + log(val2))
          }
        }

        matern12    = rho12*matern_h(h,alpha,1/2*(sd_kappa + swe_kappa))
        matern12[j] = rho12
        x1 = rep(T1[j+(0:(n.obs-1))*n.site],each = n.site)
        y2 = T2
        helpa   = (y2/x1)
        helpb   = (x1/y2)
        helpmat = (2/(1 - rep(matern12,times = n.obs)^2))^(1/2)

        gxy = helpmat*(helpa - rep(matern12,times = n.obs))
        gyx = helpmat*(helpb - rep(matern12,times = n.obs))

        dx1 = rep(dT1[j+(0:(n.obs-1))*n.site],each = n.site)
        dy2 = dT2
        helpdx1 = 1/x1*helpmat*helpa*dt(gxy, 2)
        helpdy2 = 1/y2*helpmat*helpb*dt(gyx, 2)

        helpdd = gamma(3/2)*3/(sqrt(2*pi)*(1 - rep(matern12,times = n.obs)^2))

        ptgxy = 1/x1*pt(gxy, 2)
        ptgyx = 1/y2*pt(gyx, 2)

        val = ((dy2/y2*(ptgyx + helpdy2 - helpdx1))*
                 (dx1/x1*(ptgxy + helpdx1 - helpdy2)) +
                 dx1*dy2/x1/y2*(helpdx1*2 - 1/x1*(helpdd*(1 + gxy^2/2)^(-5/2)*gxy*helpa^2) +
                                  helpdy2*2 - 1/y2*(helpdd*(1 + gyx^2/2)^(-5/2)*gyx*helpb^2)))

        if (any(val <= 0)) {
          sum = -Inf
          break
        } else {
          sum = sum + sum(-ptgxy - ptgyx + log(val))
        }
      }
    }
    # in order to maximize we take the negative value
    return(-sum)
  }

  # find good starting values for the correlation parameters by perfoming an optimization of the log likelihood
  start1        = c(100,0.1,0.1,0.9)
  names(start1) = c("alpha","sd_kappa","swe_kappa","rho12")

  # rescale the starting values for more stable optimization results
  parscale1 = 10^(floor(log10(abs(start1))))
  start1    = start1/parscale1

  # perform optimization for the correlation parameters and scale them back
  optimizer1 = optimx(start1, loglikeli1, method = "ucminf", control = list(kkt = FALSE))
  for (l in 1:nrow(optimizer1)) {
    optimizer1[l,names(start1)] = optimizer1[l,names(start1)]*parscale1
    optimizer1[l,"value"]       = -optimizer1[l,"value"]
  }

  # define total number of coefficients to be optimized for each GEV parameter
  n.loc1   = length(sd_loc_model$coefficients)
  n.scale1 = length(sd_scale_model$coefficients)
  n.shape1 = length(sd_shape_model$coefficients)
  n.loc2   = length(swe_loc_model$coefficients)
  n.scale2 = length(swe_scale_model$coefficients)
  n.shape2 = length(swe_shape_model$coefficients)

  # define starting values for the optimization as the optimized correlation parameters and
  #     the coefficients from the model selection
  start = c(as.numeric(optimizer1[1,names(start1)]),
            sd_loc_model$coefficients, sd_scale_model$coefficients, sd_shape_model$coefficients,
            swe_loc_model$coefficients, swe_scale_model$coefficients, swe_shape_model$coefficients)

  # rename the coefficients
  loc_names_1   = paste0("loc1_", c("Int",names(sd_loc_model$coefficients)[-1]), sep="")
  scale_names_1 = paste0("scale1_", c("Int",names(sd_scale_model$coefficients)[-1]), sep="")
  shape_names_1 = paste0("shape1_", c("Int",names(sd_shape_model$coefficients)[-1]), sep="")
  loc_names_2   = paste0("loc2_", c("Int",names(swe_loc_model$coefficients)[-1]), sep="")
  scale_names_2 = paste0("scale2_", c("Int",names(swe_scale_model$coefficients)[-1]), sep="")
  shape_names_2 = paste0("shape2_", c("Int",names(swe_shape_model$coefficients)[-1]), sep="")
  coeff_names   = c(names(start1), loc_names_1, scale_names_1, shape_names_1,
                    loc_names_2, scale_names_2, shape_names_2)
  names(start)  = coeff_names

  # add intercept column to covariable matrix
  C = cbind("(Intercept)" = rep(1,times = n.site), covariables)

  # calculate the GEV parameters for given coefficients
  dsgnmat2Param = function(loccoeff, scalecoeff, shapecoeff, model) {

    if (model == 1) {
      locs   = as.numeric(C[,names(sd_loc_model$coefficients)]%*%loccoeff)
      C_new  = cbind(C,"loc" = locs)
      scales = as.numeric(abs(C_new[,names(sd_scale_model$coefficients)]%*%scalecoeff))
      C_new  = cbind(C_new,"scale" = scales)
      shapes = as.numeric(C_new[,names(sd_shape_model$coefficients)]%*%shapecoeff)
    } else if (model == 2) {
      locs   = as.numeric(C[,names(swe_loc_model$coefficients)]%*%loccoeff)
      C_new  = cbind(C,"loc" = locs)
      scales = as.numeric(abs(C_new[,names(swe_scale_model$coefficients)]%*%scalecoeff))
      C_new  = cbind(C_new,"scale" = scales)
      shapes = as.numeric(C_new[,names(swe_shape_model$coefficients)]%*%shapecoeff)
    }

    ans = cbind("locs" = locs,"scales" = scales,"shapes" = shapes)

    return(ans)
  }

  # rescale the starting values for more stable optimization results
  parscale = 10^(floor(log10(abs(start))))
  start    = start/parscale

  # define the logarithmized likelihood function
  loglikeli = function(p) {

    # scale the coefficients back to original scaling
    p = p*parscale

    # define the correlation parameters and location, scale and shape coefficients
    alpha     = p[1]
    sd_kappa  = p[2]
    swe_kappa = p[3]
    rho12     = p[4]
    loccoeff1   = p[5:(n.loc1+4)]
    scalecoeff1 = p[(n.loc1+4+1):(n.loc1+4+n.scale1)]
    shapecoeff1 = p[(n.loc1+4+n.scale1+1):(n.loc1+4+n.scale1+n.shape1)]
    loccoeff2   = p[(n.loc1+4+n.scale1+n.shape1+1):(n.loc1+4+n.scale1+n.shape1+n.loc2)]
    scalecoeff2 = p[(n.loc1+4+n.scale1+n.shape1+n.loc2+1):(n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2)]
    shapecoeff2 = p[(n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2+1):
                      (n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2+n.shape2)]

    # parameters must be in certain intervals
    if (alpha <= 0 || sd_kappa <= 0 || swe_kappa <= 0 || abs(rho12) > 1 ||
        abs(rho12) > sqrt(sd_kappa*swe_kappa)/(1/2*(sd_kappa + swe_kappa))) {
      sum = -Inf
    } else {
      # calculate log likelihood
      sum  = 0
      gev1 = dsgnmat2Param(loccoeff1,scalecoeff1,shapecoeff1,model = 1)
      gev2 = dsgnmat2Param(loccoeff2,scalecoeff2,shapecoeff2,model = 2)
      if (printParam) {
        print(summary(gev1))
      }

      z1 = as.vector(z_sd)
      z2 = as.vector(z_swe)
      T1val = 1 + rep(gev1[,"shapes"],times = n.obs)*(z1 - rep(gev1[,"locs"],times = n.obs))/
        rep(gev1[,"scales"],times = n.obs)
      T2val = 1 + rep(gev2[,"shapes"],times = n.obs)*(z2 - rep(gev2[,"locs"],times = n.obs))/
        rep(gev2[,"scales"],times = n.obs)
      if(any(T1val <= 0) || any(T2val <= 0)) {
        sum = -Inf
      } else {
        pot1 = 1/rep(gev1[,"shapes"],times = n.obs)
        pot2 = 1/rep(gev2[,"shapes"],times = n.obs)
        T1   = T1val^(pot1)
        T2   = T2val^(pot2)
        dT1  = 1/rep(gev1[,"scales"],times = n.obs)*T1val^(pot1 - 1)
        dT2  = 1/rep(gev2[,"scales"],times = n.obs)*T2val^(pot2 - 1)

        for (j in 1:n.site) {
          h = sqrt(apply(t(t(covariables) - covariables[j,])^2,1,sum))
          if (j < n.site) {
            matern1 = matern_h(h[(j+1):n.site],alpha,sd_kappa)
            matern2 = matern_h(h[(j+1):n.site],alpha,swe_kappa)
            x1 = rep(T1[j+(0:(n.obs-1))*n.site],each = length((j+1):n.site))
            y1 = T1[rep((j+1):n.site,times = n.obs) + rep((0:(n.obs-1))*n.site,each = length((j+1):n.site))]
            x2 = rep(T2[j+(0:(n.obs-1))*n.site],each = length((j+1):n.site))
            y2 = T2[rep((j+1):n.site,times = n.obs) + rep((0:(n.obs-1))*n.site,each = length((j+1):n.site))]
            help1a = (y1/x1)
            help1b = (x1/y1)
            help2a = (y2/x2)
            help2b = (x2/y2)
            helpmat1 = (2/(1 - rep(matern1,times = n.obs)^2))^(1/2)
            helpmat2 = (2/(1 - rep(matern2,times = n.obs)^2))^(1/2)
            g1xy = helpmat1*(help1a - rep(matern1,times = n.obs))
            g1yx = helpmat1*(help1b - rep(matern1,times = n.obs))
            g2xy = helpmat2*(help2a - rep(matern2,times = n.obs))
            g2yx = helpmat2*(help2b - rep(matern2,times = n.obs))

            dx1 = rep(dT1[j+(0:(n.obs-1))*n.site],each = length((j+1):n.site))
            dy1 = dT1[rep((j+1):n.site,times = n.obs) + rep((0:(n.obs-1))*n.site,each = length((j+1):n.site))]
            helpdx1 = 1/x1*helpmat1*help1a*dt(g1xy, 2)
            helpdy1 = 1/y1*helpmat1*help1b*dt(g1yx, 2)

            dx2 = rep(dT2[j+(0:(n.obs-1))*n.site],each = length((j+1):n.site))
            dy2 = dT2[rep((j+1):n.site,times = n.obs) + rep((0:(n.obs-1))*n.site,each = length((j+1):n.site))]
            helpdx2 = 1/x2*helpmat2*help2a*dt(g2xy, 2)
            helpdy2 = 1/y2*helpmat2*help2b*dt(g2yx, 2)

            helpdd1 = gamma(3/2)*3/(sqrt(2*pi)*(1 - rep(matern1,times = n.obs)^2))
            helpdd2 = gamma(3/2)*3/(sqrt(2*pi)*(1 - rep(matern2,times = n.obs)^2))

            ptg1xy = 1/x1*pt(g1xy, 2)
            ptg1yx = 1/y1*pt(g1yx, 2)
            ptg2xy = 1/x2*pt(g2xy, 2)
            ptg2yx = 1/y2*pt(g2yx, 2)

            val1 = ((dy1/y1*(ptg1yx + helpdy1 - helpdx1))*
                      (dx1/x1*(ptg1xy + helpdx1 - helpdy1)) +
                      dx1*dy1/x1/y1*(helpdx1*2 - 1/x1*(helpdd1*(1 + g1xy^2/2)^(-5/2)*g1xy*help1a^2) +
                                       helpdy1*2 - 1/y1*(helpdd1*(1 + g1yx^2/2)^(-5/2)*g1yx*help1b^2)))

            val2 = ((dy2/y2*(ptg2yx + helpdy2 - helpdx2))*
                      (dx2/x2*(ptg2xy + helpdx2 - helpdy2)) +
                      dx2*dy2/x2/y2*(helpdx2*2 - 1/x2*(helpdd2*(1 + g2xy^2/2)^(-5/2)*g2xy*help2a^2) +
                                       helpdy2*2 - 1/y2*(helpdd2*(1 + g2yx^2/2)^(-5/2)*g2yx*help2b^2)))

            if (any(val1 <= 0) || any(val2 <= 0)) {
              sum = -Inf
              break
            } else {
              sum = sum + sum(-ptg1xy - ptg1yx + log(val1)) +
                sum(-ptg2xy - ptg2yx + log(val2))
            }
          }

          matern12    = rho12*matern_h(h,alpha,1/2*(sd_kappa + swe_kappa))
          matern12[j] = rho12
          x1 = rep(T1[j+(0:(n.obs-1))*n.site],each = n.site)
          y2 = T2
          helpa   = (y2/x1)
          helpb   = (x1/y2)
          helpmat = (2/(1 - rep(matern12,times = n.obs)^2))^(1/2)

          gxy = helpmat*(helpa - rep(matern12,times = n.obs))
          gyx = helpmat*(helpb - rep(matern12,times = n.obs))

          dx1 = rep(dT1[j+(0:(n.obs-1))*n.site],each = n.site)
          dy2 = dT2
          helpdx1 = 1/x1*helpmat*helpa*dt(gxy, 2)
          helpdy2 = 1/y2*helpmat*helpb*dt(gyx, 2)

          helpdd = gamma(3/2)*3/(sqrt(2*pi)*(1 - rep(matern12,times = n.obs)^2))

          ptgxy = 1/x1*pt(gxy, 2)
          ptgyx = 1/y2*pt(gyx, 2)

          val = ((dy2/y2*(ptgyx + helpdy2 - helpdx1))*
                   (dx1/x1*(ptgxy + helpdx1 - helpdy2)) +
                   dx1*dy2/x1/y2*(helpdx1*2 - 1/x1*(helpdd*(1 + gxy^2/2)^(-5/2)*gxy*helpa^2) +
                                    helpdy2*2 - 1/y2*(helpdd*(1 + gyx^2/2)^(-5/2)*gyx*helpb^2)))

          if (any(val <= 0)) {
            sum = -Inf
            break
          } else {
            sum = sum + sum(-ptgxy - ptgyx + log(val))
          }
        }
      }
    }
    # in order to maximize we take the negative value
    return(-sum)
  }

  # if loglikeli(start) is infinite, try to eliminate the responsible stations
  if (is.infinite(loglikeli(start))) {
    # define the coefficients
    coeff     = start*parscale
    alpha     = coeff[1]
    sd_kappa  = coeff[2]
    swe_kappa = coeff[3]
    rho12     = coeff[4]
    loccoeff1   = coeff[5:(n.loc1+4)]
    scalecoeff1 = coeff[(n.loc1+4+1):(n.loc1+4+n.scale1)]
    shapecoeff1 = coeff[(n.loc1+4+n.scale1+1):(n.loc1+4+n.scale1+n.shape1)]
    loccoeff2   = coeff[(n.loc1+4+n.scale1+n.shape1+1):(n.loc1+4+n.scale1+n.shape1+n.loc2)]
    scalecoeff2 = coeff[(n.loc1+4+n.scale1+n.shape1+n.loc2+1):(n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2)]
    shapecoeff2 = coeff[(n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2+1):
                          (n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2+n.shape2)]

    # calculate the GEV parameters according to the given coefficients
    gev1 = dsgnmat2Param(loccoeff1,scalecoeff1,shapecoeff1,model = 1)
    gev2 = dsgnmat2Param(loccoeff2,scalecoeff2,shapecoeff2,model = 2)

    # identify the stations that cause the infinity (-> bad stations)
    bad = NULL
    j = 0
    for (k in 1:n.site) {
      z1   = z_sd[k,]
      z2   = z_swe[k,]
      # the terms val1 and val2 should be positive!
      val1 = 1 + gev1[k,"shapes"]*((z1 - gev1[k,"locs"])/gev1[k,"scales"])
      val2 = 1 + gev2[k,"shapes"]*((z2 - gev2[k,"locs"])/gev2[k,"scales"])
      if (any(val1 <= 0) || any(val2 <=0)) {
        j = j + 1
        bad[j] = k
      }
    }

    # warn if there are more than 10% of bad stations
    if (length(bad) >= floor(n.site*0.1)) {
      warning(sprintf(c("more than 10 %% of all stations (%i) caused infinite log likelihood",
                        " with starting values", "\n  all those stations were removed"),length(bad)))
    }

    # drop the bad stations
    z_sd   = z_sd[-bad,]
    z_swe  = z_swe[-bad,]
    n.site = nrow(z_sd)
    covariables = covariables[-bad,]
    C = cbind("(Intercept)" = rep(1,times = n.site), covariables)
  }

  # if loglikeli(start) is still infinite, there is another problem
  if (is.infinite(loglikeli(start))) {
    stop("starting value delivers infinite loglikelihood")
  } else {
    optloglik = loglikeli
  }

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

  # define the correlation parameters and coefficients for each GEV parameter separately
  alpha     = coeff[1]
  sd_kappa  = coeff[2]
  swe_kappa = coeff[3]
  rho12     = coeff[4]
  sd_loccoeff    = coeff[5:(n.loc1+4)]
  sd_scalecoeff  = coeff[(n.loc1+4+1):(n.loc1+4+n.scale1)]
  sd_shapecoeff  = coeff[(n.loc1+4+n.scale1+1):(n.loc1+4+n.scale1+n.shape1)]
  swe_loccoeff   = coeff[(n.loc1+4+n.scale1+n.shape1+1):(n.loc1+4+n.scale1+n.shape1+n.loc2)]
  swe_scalecoeff = coeff[(n.loc1+4+n.scale1+n.shape1+n.loc2+1):(n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2)]
  swe_shapecoeff = coeff[(n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2+1):
                           (n.loc1+4+n.scale1+n.shape1+n.loc2+n.scale2+n.shape2)]

  names(sd_loccoeff)    = names(sd_loc_model$coefficients)
  names(sd_scalecoeff)  = names(sd_scale_model$coefficients)
  names(sd_shapecoeff)  = names(sd_shape_model$coefficients)
  names(swe_loccoeff)   = names(swe_loc_model$coefficients)
  names(swe_scalecoeff) = names(swe_scale_model$coefficients)
  names(swe_shapecoeff) = names(swe_shape_model$coefficients)

  # summarize optimization result and information message
  summary = list(message = message, method = row.names(optimizer)[ind], result = optimizer)
  print(summary)
  sd_coeff  = list("loccoeff" = sd_loccoeff, "scalecoeff" = sd_scalecoeff,
                   "shapecoeff" = sd_shapecoeff)
  swe_coeff = list("loccoeff" = swe_loccoeff, "scalecoeff" = swe_scalecoeff,
                   "shapecoeff" = swe_shapecoeff)
  cor_coeff = c(alpha, sd_kappa, swe_kappa, rho12)
  coefficients = list(sd_coeff = sd_coeff, swe_coeff = swe_coeff, cor_coeff = cor_coeff, all_coeff = coeff)

  #---------------------------------------------------------------------------------------------------------------#

  # Define function output ----------------------------------------------------------------------------------------

  ans = list(complete_sd_max_data = z_sd, complete_swe_max_data = z_swe,
             summary = summary, coefficients = coefficients)
  #filename = paste0(save_name, ".robj", sep = "")
  #message(sprintf("optimization results were saved as '%s' to directory: \n  '%s'",filename,getwd()))
  #assign(save_name, ans)
  #save(list = save_name, file = filename)

  return(ans)

  #---------------------------------------------------------------------------------------------------------------#
}
