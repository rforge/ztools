#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: cond_returnlevel_plot ----------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

cond_returnlevel_plot = function(location, GEVparam, cond_location, cond_GEVparam, same_var,
                                 cor_coeff, cond_B, cond_B_2 = Inf, obs = NULL, var = NULL, 
                                 period_range = c(1,128), model = "ext-t", printObjectives = FALSE,
                                 plottitle = NULL, save_name = NULL, save_dir = getwd(),
                                 printPlot = TRUE) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: cond_returnlevel_plot(location = arl, GEVparam = sd_GEVparam_arl, cond_location = ibk,
  #                                cond_GEVparam = swe_GEVparam_ibk, same_var = FALSE, cond_B = cond_B,
  #                                cor_coeff = cor_coeff, model = "hr",
  #                                plottitle = "sd cond rl at arlberg given swe > swe_rl50 in ibk",
  #                                printPlot = FALSE)
  # output: a conditional return level plot
  #
  # --- this function creates a conditional return level plot
  #     the q-year conditional return level of variable Z1 at location x1 given variable Z2 at location x2 is
  #           defined as the threshold B, such that the conditional probability that Z1(x1) exceeds this threshold,
  #           given that Z2(x2) is in the interval ('cond_B','cond_B_2'), is 1/q : 
  #           Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)] = 1/q
  #
  # --- input:
  #           1. 'location': a vector with certain location characteristics as entries
  #                 characteristics are for example: longitude, latitude and altitude
  #           2. 'GEVparam': a named vector with the GEV parameters of the variable for which
  #                   conditional return levels are calculated
  #                 names are 'loc', 'scale' and 'shape'
  #           3. 'cond_location': a vector with the characteristics of the conditioned location
  #                 characteristics should be the same as in 'location'
  #           4. 'cond_GEVparam': a named vector with the GEV parameters of the conditioned variable
  #                 names are 'loc', 'scale' and 'shape'
  #           5. 'same_var': logical value
  #                 if TRUE, the conditioned variable is the same as the unconditioned variable
  #                 this has to be known in order to use the right correlation function
  #           6. 'cor_coeff': a named vector with the correlation parameters
  #                 for the Huesler-Reiss model: 'alpha', 'kappa' and 'lambda12'
  #                 for the Extremal-Gaussian model: 'alpha', 'sd_kappa', 'swe_kappa' and 'rho12'
  #                 for the Extremal-t model: 'alpha', 'sd_kappa', 'swe_kappa', 'rho12' and 'nu'
  #           7. 'cond_B': a real number as the lower barrier of the conditioned variable
  #                 (e.g. the return level)
  # --- optional input:
  #           8. 'cond_B_2': a real number as the upper barrier of the conditioned variable
  #                 default (if this input is missing) is Inf
  #           9. 'obs': a vector with empirical observations
  #                 if provided, sample quantiles of this vector are added to the plot as points
  #          10. 'var': a character string being 'sd' or 'swe' if 'same_var' is TRUE and model is
  #                 'ext-t' or 'ext-gauss'
  #                   that is, if the conditioned variable is the same as the unconditioned variable,
  #                   for the Extremal-Gaussian and Extremal-t model it has to be known which variable it is
  #          11. 'period_range': the range of the return period to be plotted
  #                 a vector with start and end point, which should be numbers greater or equal than 1
  #                 default (if this input is missing) is 'c(1,128)'
  #          12. 'model': a character string
  #                 chose which bivariate max-stable model should be used to calculate the conditional
  #                   return levels; either 'hr' for the Huesler-Reiss, 'ext-gauss' for the Extremal-Gaussian or
  #                   'ext-t' for the Extremal-t model
  #                 default (if this input is missing) is 'ext-t'
  #          13. 'printObjectives': logical value
  #                 if TRUE, a summary of the values 
  #                       f(B) = abs(1/q - Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)])
  #                 is printed, where B is the found conditional return level and 
  #                 Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)] is the conditional probability that Z1(x1) 
  #                 (eg. sd or swe) exceeds B, given that Z2(x2) (eg. swe or sd) is in between cond_B and cond_B_2;
  #                 by the minimaization of the function f we find the 1/q quantile of this conditional 
  #                 distribution, thus we want small values
  #                 default (if this input is missing) is FALSE
  #          14. 'plottitle': a character string defining the title of the plot
  #                 default (if this input is missing) is 'conditional return level plot'
  #          15. 'save_name': a character string defining the saving name of the plot
  #                 default (if this input is missing) is 'cond_return_level_plot'
  #          16. 'save_dir': a character string defining the directory for the plot to be saved
  #                 default (if this input is missing) is the working directory
  #          17. 'print_plot': logical value
  #                 if TRUE, the plot is printed
  #                 default (if this input is missing) is TRUE
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
  # # load required package 'ggplot2'
  # if (inherits(try(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE), silent = TRUE), "try-error")) {
  #   message("required package 'ggplot2' is not installed yet -- trying to install package")
  #   install.packages("ggplot2", quiet = TRUE)
  #   if (inherits(try(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE), silent = TRUE), "try-error")) {
  #     stop("package 'ggplot2' couldn't be installed")
  #   } else {
  #     message("package successfully installed and loaded")
  #   }
  # }

  # check whether 'location' and 'cond_location' are vectors
  if (!is.vector(location) || !is.vector(cond_location)) {
    stop("'location' and 'cond_location' have to be vectors")
  }

  # check whether the length of 'location' and 'cond_location' coincide
  if (length(location) != length(cond_location)) {
    stop("'location' and 'cond_location' must have the same length")
  }

  # check whether 'GEVparam' and 'cond_GEVparam' are vectors
  if (!is.vector(GEVparam) || !is.vector(cond_GEVparam)) {
    stop("'GEVparam' and 'cond_GEVparam' have to be vectors")
  }

  # names of 'GEVparam' must be 'loc', 'scale' and 'shape'
  if (!(all(c("loc","scale","shape") %in% names(GEVparam)))) {
    stop("names of 'GEVparam' must be 'loc', 'scale' and 'shape'")
  }

  # names of 'cond_GEVparam' must be 'loc', 'scale' and 'shape'
  if (!(all(c("loc","scale","shape") %in% names(cond_GEVparam)))) {
    stop("names of 'cond_GEVparam' must be 'loc', 'scale' and 'shape'")
  }

  # 'cond_B' has to be a real number
  if (!is.vector(cond_B) || length(cond_B) != 1 || !is.numeric(cond_B)) {
    stop("'cond_B' has to be a real number")
  }
  
  # 'cond_B_2' has to be a real number
  if (!is.vector(cond_B_2) || length(cond_B_2) != 1 || !is.numeric(cond_B_2)) {
    stop("'cond_B_2' has to be a real number")
  }
  
  # 'cond_B_2' has to be greater than 'cond_B'
  if (cond_B_2 <= cond_B) {
    stop("'cond_B_2' has to be greater than 'cond_B'")
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

  # check whether 'obs' is a vector
  if (!missing(obs) && !is.vector(obs)) {
    stop("'obs' has to be a vector")
  }

  # 'printObjectives' has to be a logical value
  if (!(is.logical(printObjectives))) {
    stop(sprintf("'printObjectives' has to be TRUE or FALSE -- '%s' is not allowed",printObjectives))
  }

  # 'period_range' has to be a vector of length 2 with numbers greater or equal than 1
  if (length(period_range) != 2) {
    stop("please state start and end point of return level period: 'period_range = c(start,end)'")
  }
  if (any(period_range < 1)) {
    stop("numbers in 'period_range' have to be greater or equal than 1")
  }
  if (period_range[1] >= period_range[2]) {
    stop("endpoint of 'period_range' is less or equal its startpoint")
  }

  # 'printPlot' has to be a logical value
  if (!(is.logical(printPlot))) {
    stop(sprintf("'printPlot' has to be TRUE or FALSE -- '%s' is not allowed",printPlot))
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

  rl <- NULL
  
  #---------------------------------------------------------------------------------------------------------------#

  # Perform calculations ------------------------------------------------------------------------------------------

  # calculate the spatial lag between the two locations
  h = sqrt(sum(location - cond_location)^2)

  # if 'same_var' is TRUE, 'locations' and 'cond_locations' can't be the same (h != 0)
  if (same_var) {
    if (h == 0) {
      stop(c("same variable ('same_var = TRUE') can't be conditioned on same location",
             "\n  either condition on different location or use different variables"))
    }
  }

  # define non equidistant sequence of return level period
  a = period_range[1]:period_range[2]
  # define number of grid points with n
  n = period_range[2] - period_range[1]
  # partition the logarithmized intervall into equidistant points
  loga     = log(a)
  log_int  = max(loga) - min(loga)
  log_equi = loga[1] + log_int*(1:n)/n
  # transform log intervall back (q has to be greater than 1)
  if (period_range[1] == 1) {
    q = exp(log_equi)
  } else {
    q = c(period_range[1], exp(log_equi))
  }

  # transform conditional barriers 'cond_B' and 'cond_B_2' to unit Frechet
  b1 = gev2frech(cond_B, loc = as.numeric(cond_GEVparam["loc"]), scale = as.numeric(cond_GEVparam["scale"]),
                 shape = as.numeric(cond_GEVparam["shape"]))
  b2 = gev2frech(cond_B_2, loc = as.numeric(cond_GEVparam["loc"]), scale = as.numeric(cond_GEVparam["scale"]),
                 shape = as.numeric(cond_GEVparam["shape"]))

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
    cond_rl    = rep(NA,times = length(q))
    objectives = rep(NA,times = length(q))

    # find the conditional return level for every return period in q
    for (k in 1:length(q)) {
      # define function f calculating abs(1/q[k] - the conditional probability),
      #     which will be minimized in order to find the 1-q[k] quantil of the conditional distribution
      f = function(B) {
        a    = gev2frech(B, loc = as.numeric(GEVparam["loc"]), scale = as.numeric(GEVparam["scale"]),
                         shape = as.numeric(GEVparam["shape"]))
        help1 = log(b1/a)/lh
        help2 = log(b2/a)/lh
        P_sd_swe_1 = exp(-1/a*pnorm(lh/2 + help1) - 1/b1*pnorm(lh/2 - help1))
        P_sd_swe_2 = exp(-1/a*pnorm(lh/2 + help2) - 1/b2*pnorm(lh/2 - help2))
        ans = abs(1/q[k] - (exp(-1/b2) - exp(-1/b1) - P_sd_swe_2 + P_sd_swe_1)/(exp(-1/b2) - exp(-1/b1)))
        ans[which(ans == 1/q[k])] = Inf
        return(ans)
      }
      
      # define the starting value for the optimization to find the conditional return level
      start_val = returnlevels(GEVparam = GEVparam, q = q[k])

      # calculate conditional return levels
      opt = suppressWarnings(nlminb(start_val, f))
      # if optimization was not successful (value of f >= 1e-03)), 
      # try to find better solution with different starting value
      opt1 = opt
      j = 0
      while (opt$objective >= 1e-03 && j <= 17) {
        l = 1.5 + (-1)^j*floor((j+1)/2)*0.1
        start_val_new = start_val*l
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
    cond_rl    = rep(NA,times = length(q))
    objectives = rep(NA,times = length(q))

    # find the conditional return level for every return period in q
    for (k in 1:length(q)) {
      # define function f calculating abs(1/q[k] - the conditional probability),
      #     which will be minimized in order to find the 1-q[k] quantil of the conditional distribution
      f = function(B) {
        a    = gev2frech(B, loc = as.numeric(GEVparam["loc"]), scale = as.numeric(GEVparam["scale"]),
                         shape = as.numeric(GEVparam["shape"]))
        help = (2/(1 - matern^2))^(1/2)
        P_sd_swe_1 = exp(-1/a*pt(help*(b1/a - matern), df = 2) -
                           1/b1*pt(help*(a/b1 - matern), df = 2))
        P_sd_swe_2 = exp(-1/a*pt(help*(b2/a - matern), df = 2) -
                           1/b2*pt(help*(a/b2 - matern), df = 2))
        ans = abs(1/q[k] - (exp(-1/b2) - exp(-1/b1) - P_sd_swe_2 + P_sd_swe_1)/(exp(-1/b2) - exp(-1/b1)))
        ans[which(ans == 1/q[k])] = Inf
        return(ans)
      }

      # define the starting value for the optimization to find the conditional return level
      start_val = returnlevels(GEVparam = GEVparam, q = q[k])

      # calculate conditional return levels
      opt = suppressWarnings(nlminb(start_val, f))
      # if optimization was not successful (value of f >= 1e-03)), 
      # try to find better solution with different starting value
      opt1 = opt
      j = 0
      while (opt$objective >= 1e-03 && j <= 17) {
        l = 1.5 + (-1)^j*floor((j+1)/2)*0.1
        start_val_new = start_val*l
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
    cond_rl    = rep(NA,times = length(q))
    objectives = rep(NA,times = length(q))

    # find the conditional return level for every return period in q
    for (k in 1:length(q)) {
      # define function f calculating abs(1/q[k] - the conditional probability),
      #     which will be minimized in order to find the 1-q[k] quantil of the conditional distribution
      f = function(B) {
        a    = gev2frech(B, loc = as.numeric(GEVparam["loc"]), scale = as.numeric(GEVparam["scale"]),
                         shape = as.numeric(GEVparam["shape"]))
        help = ((nu + 1)/(1 - matern^2))^(1/2)
        P_sd_swe_1 = exp(-1/a*pt(help*((b1/a)^(1/nu) - matern), df = nu + 1) -
                           1/b1*pt(help*((a/b1)^(1/nu) - matern), df = nu + 1))
        P_sd_swe_2 = exp(-1/a*pt(help*((b2/a)^(1/nu) - matern), df = nu + 1) -
                           1/b2*pt(help*((a/b2)^(1/nu) - matern), df = nu + 1))
        ans = abs(1/q[k] - (exp(-1/b2) - exp(-1/b1) - P_sd_swe_2 + P_sd_swe_1)/(exp(-1/b2) - exp(-1/b1)))
        ans[which(ans == 1/q[k])] = Inf
        return(ans)
      }
      
      # define the starting value for the optimization to find the conditional return level
      start_val = returnlevels(GEVparam = GEVparam, q = q[k])

      # calculate conditional return levels
      opt = suppressWarnings(nlminb(start_val, f))
      # if optimization was not successful (value of f >= 1e-03)), 
      # try to find better solution with different starting value
      opt1 = opt
      j = 0
      while (opt$objective >= 1e-03 && j <= 17) {
        l = 1.5 + (-1)^j*floor((j+1)/2)*0.1
        start_val_new = start_val*l
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

  # define dataframe with 'cond_rl' and 'q'
  df = data.frame(cond_rl = cond_rl, q = q)

  # define plottitle
  if (missing(plottitle)) {
    plottitle = "conditional return level plot"
  }

  # define breaks and labels
  b = period_range[1]*2^(0:10)
  b = c(b[which(b < period_range[2])],period_range[2])
  l = b

  # prepare GEV parameters for the ggplot
  GEV = paste0("loc = ",round(GEVparam["loc"], digits = 2),
               ", scale = ",round(GEVparam["scale"], digits = 2),
               ", shape = ",round(GEVparam["shape"], digits = 2))

  # create plot
  if (missing(obs)) {
    plot = ggplot(df,aes(x = q, y = cond_rl)) +
      geom_line() +
      scale_x_log10(name = "return period",
                    breaks = b,
                    labels = l) +
      scale_y_continuous(name = "conditional return level") +
      geom_label(aes(label = GEV, x = b[1], y = max(cond_rl, na.rm = TRUE)), hjust = 0, vjust = 1, size = 3) +
      ggtitle(plottitle)
  } else {
    r = quantile(obs, 1-1/q, names = FALSE)
    df.obs = data.frame(rl = r, q = q)
    plot = ggplot(df,aes(x = q, y = cond_rl)) +
      geom_line() +
      geom_point(data = df.obs, aes(x = q, y = rl)) +
      scale_x_log10(name = "return period",
                    breaks = b,
                    labels = l) +
      scale_y_continuous(name = "conditional return level") +
      geom_label(aes(label = GEV, x = b[1], y = max(max(df$cond_rl, na.rm = TRUE), max(df.obs$rl, na.rm = TRUE))), 
                 hjust = 0, vjust = 1, size = 3) +
      ggtitle(plottitle)
  }

  # define save name
  if (!missing(save_name)) {
    plotname = paste0(save_name,".pdf",sep = "")
    # plot is saved to given directory
    path = save_dir
    message(sprintf("plot was saved as '%s' to directory: \n  '%s'",plotname,path))
    # save plot
    ggsave(filename = plotname, plot = plot, path = path, width = 7, height = 7)
    
  }
  
  # print plot if 'printPlot' is TRUE
  if (printPlot) {
    print(plot)
  }

  #---------------------------------------------------------------------------------------------------------------#

}
