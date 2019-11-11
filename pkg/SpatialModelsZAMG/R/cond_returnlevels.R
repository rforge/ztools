#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: cond_returnlevels --------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

cond_returnlevels = function(locations, GEVparam, q, cond_locations, cond_GEVparam, same_var,
                             cor_coeff, cond_B, cond_B_2 = Inf, var = NULL, model = "ext-t", 
                             printObjectives = FALSE) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: cond_returnlevels(locations = covariables, GEVparam = sd_GEVparam, q = 100,
  #                            cond_locations = covariables, cond_GEVparam = swe_GEVparam, same_var = FALSE,
  #                            cond_B = swe_rl, cor_coeff = cor_coeff, model = "hr", printObjectives = TRUE)
  # output: a vector with the conditional return levels
  #
  # --- this function calculates the conditional return levels
  #     the q-year conditional return level of variable Z1 at location x1 given variable Z2 at location x2 is
  #           defined as the threshold B, such that the conditional probability that Z1(x1) exceeds this threshold,
  #           given that Z2(x2) is in the interval ('cond_B','cond_B_2'), is 1/q : 
  #           Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)] = 1/q
  #
  # --- input:
  #           1. 'locations': a matrix or vector with certain location characteristics as columns/entries
  #                 each row corresponds to one location
  #                 columns are for example: longitude, latitude and altitude
  #           2. 'GEVparam': a named matrix or vector with the GEV parameters of the variable for which
  #                   conditional return levels are calculated
  #                 each row corresponds to one location (same ones as in 'locations')
  #                 columns are 'loc', 'scale' and 'shape'
  #           3. 'q': the return period for the calculation of the conditional return levels --
  #                 must be a number greater than 1
  #           4. 'cond_locations': a matrix or a vector with the characteristics of the conditioned locations
  #                 each row corresponds to one location
  #                 columns should be the same as in 'locations'
  #           5. 'cond_GEVparam': a named matrix or vector with the GEV parameters of the conditioned variable
  #                 each row corresponds to one location (same ones as in 'cond_locations')
  #                 columns are 'loc', 'scale' and 'shape'
  #           6. 'same_var': logical value
  #                 if TRUE, the conditioned variable is the same as the unconditioned variable
  #                 this has to be known in order to use the right correlation function
  #           7. 'cor_coeff': a named vector with the correlation parameters
  #                 for the Huesler-Reiss model: 'alpha', 'kappa' and 'lambda12'
  #                 for the Extremal-Gaussian model: 'alpha', 'sd_kappa', 'swe_kappa' and 'rho12'
  #                 for the Extremal-t model: 'alpha', 'sd_kappa', 'swe_kappa', 'rho12' and 'nu'
  #           8. 'cond_B': a vector of real numbers as the lower barriers of the conditioned variable 
  #                 (e.g. the return levels)
  # --- optional input:
  #           9. 'cond_B_2': a vector of real numbers as the upper barriers of the conditioned variable
  #                 default (if this input is missing) is Inf
  #          10. 'var': a character string being 'sd' or 'swe' if 'same_var' is TRUE and model is
  #                 'ext-t' or 'ext-gauss'
  #                   that is, if the conditioned variable is the same as the unconditioned variable,
  #                   for the Extremal-Gaussian and Extremal-t model it has to be known which variable it is
  #          11. 'model': a character string
  #                 chose which bivariate max-stable model should be used to calculate the conditional
  #                   return levels; either 'hr' for the Huesler-Reiss, 'ext-gauss' for the Extremal-Gaussian or
  #                   'ext-t' for the Extremal-t model
  #                 default (if this input is missing) is 'ext-t'
  #          12. 'printObjectives': logical value
  #                 if TRUE, a summary of the values
  #                       f(B) = abs(1/q - Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)])
  #                 is printed, where B is the found conditional return level and 
  #                 Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)] is the conditional probability that Z1(x1) 
  #                 (eg. sd or swe) exceeds B, given that Z2(x2) (eg. swe or sd) is in between cond_B and cond_B_2;
  #                 by the minimaization of the function f we find the 1/q quantile of this conditional 
  #                 distribution, thus we want small values
  #                 default (if this input is missing) is FALSE
  # --- output:
  #           1. a vector with the conditional return levels: 'cond_rl'
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

  # check whether 'locations' is a vector -- if it is, transform it to a matrix
  if (is.vector(locations)) {
    locations = as.matrix(t(locations))
  }

  # check whether 'GEVparam' is a vector -- if it is, transform it to a matrix
  if (is.vector(GEVparam)) {
    GEVparam = as.matrix(t(GEVparam))
  }

  # check whether the number of rows in 'locations' and 'GEVparam' coincide
  if (nrow(locations) != nrow(GEVparam)) {
    stop(sprintf("number of rows (%i) in 'locations' and number of rows (%i) in 'GEVparam' don't match up",
                 nrow(locations), nrow(GEVparam)))
  }

  # columns of 'GEVparam' must be 'loc', 'scale' and 'shape'
  if (!(all(c("loc","scale","shape") %in% colnames(GEVparam)))) {
    stop("columns/names of 'GEVparam' must be 'loc', 'scale' and 'shape'")
  }

  # 'q' has to be a number greater than 1
  if (q <= 1) {
    stop(sprintf("'q' has to be a number greater than 1 -- '%s' is not allowed",q))
  }

  # check whether 'cond_locations' is a vector -- if it is, transform it to a matrix
  if (is.vector(cond_locations)) {
    cond_locations = as.matrix(t(cond_locations))
  }

  # check whether 'cond_GEVparam' is a vector -- if it is, transform it to a matrix
  if (is.vector(cond_GEVparam)) {
    cond_GEVparam = as.matrix(t(cond_GEVparam))
  }

  # 'cond_B' has to be a vector of real numbers
  if (!is.vector(cond_B) || !is.numeric(cond_B)) {
    stop("'cond_B' has to be a vector of real numbers")
  }
  
  # 'cond_B_2' has to be a vector of real numbers
  if (!is.vector(cond_B_2) || !is.numeric(cond_B_2)) {
    stop("'cond_B_2' has to be a vector of real numbers")
  }

  # the number of rows in 'cond_locations' has to be either the same as the number of rows in 'locations' or 1
  if (nrow(cond_locations) != nrow(locations) && nrow(cond_locations) != 1) {
    stop(sprintf(c("number of rows in 'cond_locations' (%i) has to be either the same",
                   "\n  as the number of rows in 'locations' (%i) or 1"),
                 c(nrow(cond_locations),nrow(locations))))
  }

  # check whether the number of columns in 'locations' and 'cond_locations' coincide
  if (ncol(locations) != ncol(cond_locations)) {
    stop(sprintf("number of columns in 'locations' (%i) and 'cond_locations' (%i) don't match up",
                 ncol(locations),ncol(cond_locations)))
  }

  # check whether the number of rows in 'cond_locations' and 'cond_GEVparam' coincide
  if (nrow(cond_locations) != nrow(cond_GEVparam)) {
    stop(sprintf(
      "number of rows (%i) in 'cond_locations' and number of rows (%i) in 'cond_GEVparam' don't match up",
      nrow(cond_locations),nrow(cond_GEVparam)))
  }

  # columns of 'cond_GEVparam' must be 'loc', 'scale' and 'shape'
  if (!(all(c("loc","scale","shape") %in% colnames(cond_GEVparam)))) {
    stop("columns/names of 'cond_GEVparam' must be 'loc', 'scale' and 'shape'")
  }

  # check whether the number of rows in 'cond_locations' and the length of 'cond_B' coincide
  if (nrow(cond_locations) != length(cond_B)) {
    stop(sprintf("number of rows (%i) in 'cond_locations' and length of 'cond_B' (%i) don't match up",
                 nrow(cond_locations),length(cond_B)))
  }
  
  # check whether the length of 'cond_B' and the length of 'cond_B_2' coincide
  if (length(cond_B) != length(cond_B_2)) {
    if (length(cond_B_2) == 1 && is.infinite(cond_B_2)) {
      cond_B_2 = rep(Inf,times = length(cond_B))
    } else {
      stop(sprintf("length of 'cond_B' (%i) and length of 'cond_B_2' (%i) don't match up",
                   length(cond_B),length(cond_B_2)))
    }
  }
  
  # 'cond_B_2' has to be greater than 'cond_B'
  for (k in 1:length(cond_B)) {
    if (cond_B_2[k] <= cond_B[k]) {
      stop(sprintf(c("each entry of 'cond_B_2' has to be greater than its respective entry in 'cond_B'",
             "\n  check e.g. entry %i"),k))
    }
  }

  # 'model' has to be either 'hr', 'ext-gauss' or 'ext-t'
  if (!(model %in% c("hr","ext-gauss","ext-t"))) {
    stop(sprintf("'model' has to be either 'hr', 'ext-gauss' or 'ext-t' -- '%s' is not allowed",model))
  }

  # 'cor_coeff' must be a named vector with the correlation coefficients according to the chosen model
  if (model == "hr" && !all(c("alpha","kappa","lambda12") %in% names(cor_coeff))) {
    stop(c("for the Huesler-Reiss model 'cor_coeff' must be a named vector with the correlation coefficients",
           "\n  'alpha', 'kappa' and 'lambda12' -- you might change input argument 'model'"))
  } else if (model == "ext-gauss" && !all(c("alpha","sd_kappa","swe_kappa","rho12") %in% names(cor_coeff))) {
    stop(c("for the Extremal-Gaussian model 'cor_coeff' must be a named vector with the correlation coefficients",
           "\n  'alpha', 'sd_kappa', 'swe_kappa' and 'rho12' -- you might change input argument 'model'"))
  } else if (model == "ext-t" && !all(c("alpha","sd_kappa","swe_kappa","rho12","nu") %in% names(cor_coeff))) {
    stop(c("for the Extremal-t model 'cor_coeff' must be a named vector with the correlation coefficients",
           "\n  'alpha', 'sd_kappa', 'swe_kappa', 'rho12' and 'nu' -- you might change input argument 'model'"))
  }

  # 'printObjectives' has to be a logical value
  if (!(is.logical(printObjectives))) {
    stop(sprintf("'printObjectives' has to be TRUE or FALSE -- '%s' is not allowed",printObjectives))
  }

  # 'same_var' has to be a logical value
  if (!(is.logical(same_var))) {
    stop(sprintf("'same_var' has to be TRUE or FALSE -- '%s' is not allowed",same_var))
  }

  # if 'same_var' is TRUE, the variable has to be named for the 'ext-gauss' or 'ext-t' model,
  # this can be done with the input argument 'var' being 'sd' or 'swe'
  if (same_var && (model %in% c("ext-gauss","ext-t"))) {
    if (missing(var)) {
      stop(c("if 'same_var' is TRUE, the variable has to be named",
             "\n  this can be done with the input argument 'var' being 'sd' or 'swe'"))
    } else{
      if (!(var %in% c("sd","swe"))) {
        stop(sprintf("'var' has to be either 'sd' or 'swe' -- '%s' is not allowed",var))
      }
    }
  }

  #---------------------------------------------------------------------------------------------------------------#

  # Perform calculations ------------------------------------------------------------------------------------------

  # define number of locations
  n.site = nrow(locations)

  # calculate the spatial lag between locations
  if (nrow(cond_locations) == 1) {
    h = sqrt(apply((locations - matrix(rep(cond_locations,times = n.site), nrow = n.site, byrow = TRUE))^2,1,sum))
  } else {
    h = sqrt(apply((locations - cond_locations)^2,1,sum))
  }

  # if 'same_var' is TRUE, 'locations' and 'cond_locations' can't be the same (h != 0)
  if (same_var) {
    if (any(h == 0)) {
      x = which(h == 0)
      if (length(x) > 1) {
        stop(c("same variable ('same_var = TRUE') can't be conditioned on same locations",
               "\n  either condition on different locations or use different variables"))
      } else
      stop(sprintf(c("entry %i of 'locations' and 'cond_locations' coincide",
                     "\n  you can't condition the same variable ('same_var = TRUE') on the same location",
                     "\n  drop that entry in 'locations'"),x))
    }
  }

  # transform conditional barriers 'cond_B' and 'cond_B_2' to unit Frechet
  if (length(cond_B) == 1) {
    b1 = gev2frech(cond_B, loc = as.numeric(cond_GEVparam[,"loc"]), scale = as.numeric(cond_GEVparam[,"scale"]),
                   shape = as.numeric(cond_GEVparam[,"shape"]))
    b2 = gev2frech(cond_B_2, loc = as.numeric(cond_GEVparam[,"loc"]), scale = as.numeric(cond_GEVparam[,"scale"]),
                   shape = as.numeric(cond_GEVparam[,"shape"]))
    b1 = rep(b1,times = n.site)
    b2 = rep(b2,times = n.site)
  } else {
    b1 = rep(NA,times = n.site)
    b2 = rep(NA,times = n.site)
    for (k in 1:n.site) {
      b1[k] = gev2frech(cond_B[k], loc = as.numeric(cond_GEVparam[k,"loc"]),
                        scale = as.numeric(cond_GEVparam[k,"scale"]), 
                        shape = as.numeric(cond_GEVparam[k,"shape"]))
      b2[k] = gev2frech(cond_B_2[k], loc = as.numeric(cond_GEVparam[k,"loc"]),
                        scale = as.numeric(cond_GEVparam[k,"scale"]), 
                        shape = as.numeric(cond_GEVparam[k,"shape"]))
    }
  }
  
  # define starting values for the optimization to find the conditional return levels (see below)
  start_val = returnlevels(GEVparam = GEVparam, q = q)
  
  # depending on the model calculate the conditional return levels
  if (model == "hr") {
    # define the correlation function for the Huesler-Reiss model
    lambda_h = function(h,alpha,kappa) {
      as.numeric(sqrt((h/alpha)^kappa))
    }
    lambda12_h = function(h,alpha,kappa,lambda12) {
      as.numeric(sqrt(lambda12^2 + lambda_h(h,alpha,kappa)^2))
    }

    # define the correlation parameters
    alpha    = cor_coeff["alpha"]
    kappa    = cor_coeff["kappa"]
    lambda12 = cor_coeff["lambda12"]

    # calculate correlations or cross-correlations, depending on 'same_var'
    if (same_var) {
      lh = lambda_h(h,alpha,kappa)
    } else {
      lh = lambda12_h(h, alpha, kappa, lambda12)
    }

    # predefine conditional return level and objectives vector
    cond_rl    = rep(NA,times = n.site)
    objectives = rep(NA,times = n.site)

    # find the conditional return level for every location
    for (k in 1:n.site) {
      # define function f calculating abs(1/q - the conditional probability),
      #     which will be minimized in order to find the 1-q quantil of the conditional distribution
      f = function(B) {
        a    = gev2frech(B, loc = as.numeric(GEVparam[k,"loc"]), scale = as.numeric(GEVparam[k,"scale"]),
                         shape = as.numeric(GEVparam[k,"shape"]))
        help1 = log(b1[k]/a)/lh[k]
        help2 = log(b2[k]/a)/lh[k]
        P_sd_swe_1 = exp(-1/a*pnorm(lh[k]/2 + help1) - 1/b1[k]*pnorm(lh[k]/2 - help1))
        P_sd_swe_2 = exp(-1/a*pnorm(lh[k]/2 + help2) - 1/b2[k]*pnorm(lh[k]/2 - help2))
        ans = abs(1/q - (exp(-1/b2[k]) - exp(-1/b1[k]) - P_sd_swe_2 + P_sd_swe_1)/(exp(-1/b2[k]) - exp(-1/b1[k])))
        ans[which(ans == 1/q)] = Inf
        return(ans)
      }
      
      # calculate conditional return levels
      opt = suppressWarnings(nlminb(start_val[k], f))
      # if optimization was not successful (value of f >= 1e-03)), 
      # try to find better solution with different starting value
      opt1 = opt
      j = 0
      while (opt$objective >= 1e-03 && j <= 17) {
        l = 1.5 + (-1)^j*floor((j+1)/2)*0.1
        start_val_new = start_val[k]*l
        if (!is.infinite(f(start_val_new))) {
          opt = suppressWarnings(nlminb(start_val_new, f))
        }        
        j = j + 1
      }
      if (opt1$objective <= opt$objective) {
        opt = opt1
      }
      if (opt$objective >= 1e-03) {
        opt_seq = seq(0.8*start_val[k], 1.5*start_val[k], by = 0.001)
        opt_ind = which.min(f(opt_seq))
        if (all(is.na(f(opt_seq)))) {
          cond_rl[k]    = NA
          objectives[k] = NA
          warning("probability that conditioned variable is greater than 'cond_B' is zero -- change 'cond_B'")
        } else if (f(opt_seq[opt_ind]) <= opt$objective) {
          cond_rl[k]    = opt_seq[opt_ind]
          objectives[k] = f(opt_seq[opt_ind])
        } else {
          cond_rl[k]    = opt$par
          objectives[k] = opt$objective
        }
      } else {
        cond_rl[k]    = opt$par
        objectives[k] = opt$objective
      }
    }
  } else if (model == "ext-gauss") {
    # define the correlation function for the Extremal-Gaussian model
    matern_h = function(h,alpha,kappa) {
      2^(1-kappa)/gamma(kappa)*(h/alpha)^kappa*besselK(h/alpha, kappa)
    }

    # define the correlation parameters
    alpha     = cor_coeff["alpha"]
    sd_kappa  = cor_coeff["sd_kappa"]
    swe_kappa = cor_coeff["swe_kappa"]
    rho12     = cor_coeff["rho12"]

    # calculate correlations or cross-correlations, depending on 'same_var'
    if (same_var) {
      if (var == "sd") {
        matern = matern_h(h,alpha,sd_kappa)
      } else if (var == "swe") {
        matern = matern_h(h,alpha,swe_kappa)
      }
    } else{
      matern = rho12*matern_h(h,alpha,1/2*(sd_kappa + swe_kappa))
      matern[which(is.na(matern))] = rho12
    }

    # predefine conditional return level and objectives vector
    cond_rl    = rep(NA,times = n.site)
    objectives = rep(NA,times = n.site)

    # find the conditional return level for every location
    for (k in 1:n.site) {
      # define function f calculating abs(1/q - the conditional probability),
      #     which will be minimized in order to find the 1-q quantil of the conditional distribution
      f = function(B) {
        a    = gev2frech(B, loc = as.numeric(GEVparam[k,"loc"]), scale = as.numeric(GEVparam[k,"scale"]),
                         shape = as.numeric(GEVparam[k,"shape"]))
        help = (2/(1 - matern[k]^2))^(1/2)
        P_sd_swe_1 = exp(-1/a*pt(help*(b1[k]/a - matern[k]), df = 2) - 
                           1/b1[k]*pt(help*(a/b1[k] - matern[k]), df = 2))
        P_sd_swe_2 = exp(-1/a*pt(help*(b2[k]/a - matern[k]), df = 2) - 
                           1/b2[k]*pt(help*(a/b2[k] - matern[k]), df = 2))
        ans = abs(1/q - (exp(-1/b2[k]) - exp(-1/b1[k]) - P_sd_swe_2 + P_sd_swe_1)/(exp(-1/b2[k]) - exp(-1/b1[k])))
        ans[which(ans == 1/q)] = Inf
        return(ans)
      }

      # calculate conditional return levels
      opt = suppressWarnings(nlminb(start_val[k], f))
      # if optimization was not successful (value of f >= 1e-03)), 
      # try to find better solution with different starting value
      opt1 = opt
      j = 0
      while (opt$objective >= 1e-03 && j <= 17) {
        l = 1.5 + (-1)^j*floor((j+1)/2)*0.1
        start_val_new = start_val[k]*l
        if (!is.infinite(f(start_val_new))) {
          opt = suppressWarnings(nlminb(start_val_new, f))
        }        
        j = j + 1
      }
      if (opt1$objective <= opt$objective) {
        opt = opt1
      }
      if (opt$objective >= 1e-03) {
        opt_seq = seq(0.8*start_val[k], 1.5*start_val[k], by = 0.001)
        opt_ind = which.min(f(opt_seq))
        if (all(is.na(f(opt_seq)))) {
          cond_rl[k]    = NA
          objectives[k] = NA
          warning("probability that conditioned variable is greater than 'cond_B' is zero -- change 'cond_B'")
        } else if (f(opt_seq[opt_ind]) <= opt$objective) {
          cond_rl[k]    = opt_seq[opt_ind]
          objectives[k] = f(opt_seq[opt_ind])
        } else {
          cond_rl[k]    = opt$par
          objectives[k] = opt$objective
        }
      } else {
        cond_rl[k]    = opt$par
        objectives[k] = opt$objective
      }
    }
  } else if (model == "ext-t") {
    # define the correlation function for the Extremal-t model
    matern_h = function(h,alpha,kappa) {
      2^(1-kappa)/gamma(kappa)*(h/alpha)^kappa*besselK(h/alpha, kappa)
    }

    # define the correlation parameters
    alpha     = cor_coeff["alpha"]
    sd_kappa  = cor_coeff["sd_kappa"]
    swe_kappa = cor_coeff["swe_kappa"]
    rho12     = cor_coeff["rho12"]
    nu        = cor_coeff["nu"]

    # calculate correlations or cross-correlations, depending on 'same_var'
    if (same_var) {
      if (var == "sd") {
        matern = matern_h(h,alpha,sd_kappa)
      } else if (var == "swe") {
        matern = matern_h(h,alpha,swe_kappa)
      }
    } else{
      matern = rho12*matern_h(h,alpha,1/2*(sd_kappa + swe_kappa))
      matern[which(is.na(matern))] = rho12
    }

    # predefine conditional return level and objectives vector
    cond_rl    = rep(NA,times = n.site)
    objectives = rep(NA,times = n.site)

    # find the conditional return level for every location
    for (k in 1:n.site) {
      # define function f calculating abs(1/q - the conditional probability),
      #     which will be minimized in order to find the 1-q quantil of the conditional distribution
      f = function(B) {
        a    = gev2frech(B, loc = as.numeric(GEVparam[k,"loc"]), scale = as.numeric(GEVparam[k,"scale"]),
                         shape = as.numeric(GEVparam[k,"shape"]))
        help = ((nu + 1)/(1 - matern[k]^2))^(1/2)
        P_sd_swe_1 = exp(-1/a*pt(help*((b1[k]/a)^(1/nu) - matern[k]), df = nu + 1) -
                           1/b1[k]*pt(help*((a/b1[k])^(1/nu) - matern[k]), df = nu + 1))
        P_sd_swe_2 = exp(-1/a*pt(help*((b2[k]/a)^(1/nu) - matern[k]), df = nu + 1) -
                           1/b2[k]*pt(help*((a/b2[k])^(1/nu) - matern[k]), df = nu + 1))
        ans = abs(1/q - (exp(-1/b2[k]) - exp(-1/b1[k]) - P_sd_swe_2 + P_sd_swe_1)/(exp(-1/b2[k]) - exp(-1/b1[k])))
        ans[which(ans == 1/q)] = Inf
        return(ans)
      }

      # calculate conditional return levels
      opt = suppressWarnings(nlminb(start_val[k], f))
      # if optimization was not successful (value of f >= 1e-03)), 
      # try to find better solution with different starting value
      opt1 = opt
      j = 0
      while (opt$objective >= 1e-03 && j <= 17) {
        l = 1.5 + (-1)^j*floor((j+1)/2)*0.1
        start_val_new = start_val[k]*l
        if (!is.infinite(f(start_val_new))) {
          opt = suppressWarnings(nlminb(start_val_new, f))
        }        
        j = j + 1
      }
      if (opt1$objective <= opt$objective) {
        opt = opt1
      }
      if (opt$objective >= 1e-03) {
        opt_seq = seq(0.8*start_val[k], 1.5*start_val[k], by = 0.001)
        opt_ind = which.min(f(opt_seq))
        if (all(is.na(f(opt_seq)))) {
          cond_rl[k]    = NA
          objectives[k] = NA
          warning("probability that conditioned variable is greater than 'cond_B' is zero -- change 'cond_B'")
        } else if (f(opt_seq[opt_ind]) <= opt$objective) {
          cond_rl[k]    = opt_seq[opt_ind]
          objectives[k] = f(opt_seq[opt_ind])
        } else {
          cond_rl[k]    = opt$par
          objectives[k] = opt$objective
        }
      } else {
        cond_rl[k]    = opt$par
        objectives[k] = opt$objective
      }
    }
  }

  if (printObjectives) {
    print(summary(objectives))
  }

  # Define function output ----------------------------------------------------------------------------------------

  return(cond_rl)

  #---------------------------------------------------------------------------------------------------------------#
}
