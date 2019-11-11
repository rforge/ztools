returnlevel_plot = function(GEVparam, period_range = c(1,128), obs = NULL, plottitle = NULL,
                            save_name = NULL, save_dir = getwd(), printPlot = TRUE) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: returnlevel_plot(GEVparam = sd_GEVparam_arl, save_name = "return_level_plot_arlberg",
  #                           printPlot = FALSE)
  # output: a return level plot
  #
  # --- this function creates a return level plot for a given location
  #
  # --- input:
  #           1. 'GEVparam': a named vector with the GEV parameters of the location for which
  #                   a return level plot is wanted
  #                 names are 'loc', 'scale' and 'shape'
  # --- optional input:
  #           2. 'period_range': the range of the return period to be plotted
  #                 a vector with start and end point, which should be numbers greater or equal than 1
  #                 default (if this input is missing) is 'c(1,128)'
  #           3. 'obs': a vector with empirical observations
  #                 if provided, sample quantiles of this vector are added to the plot as points
  #           4. 'plottitle': a character string defining the title of the plot
  #                 default (if this input is missing) is 'return level plot'
  #           5. 'save_name': a character string defining the saving name of the plot
  #                 default (if this input is missing) is 'return_level_plot'
  #           6. 'save_dir': a character string defining the directory for the plot to be saved
  #                 default (if this input is missing) is the working directory
  #           7. 'print_plot': logical value
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

  # check whether 'GEVparam' is a vector
  if (!is.vector(GEVparam)) {
    stop("'GEVparam' has to be a vector")
  }

  # names of 'GEVparam' must be 'loc', 'scale' and 'shape'
  if (!(all(c("loc","scale","shape") %in% names(GEVparam)))) {
    stop("names of 'GEVparam' must be 'loc', 'scale' and 'shape'")
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

  # check whether 'obs' is a vector
  if (!missing(obs) && !is.vector(obs)) {
      stop("'obs' has to be a vector")
  }

  # 'printPlot' has to be a logical value
  if (!(is.logical(printPlot))) {
    stop(sprintf("'printPlot' has to be TRUE or FALSE -- '%s' is not allowed",printPlot))
  }

  #---------------------------------------------------------------------------------------------------------------#

  # Perform calculations ------------------------------------------------------------------------------------------

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

  # calculate return levels
  p = 1/q
  rl = qgev(1-p, loc = GEVparam["loc"], scale = GEVparam["scale"],
            shape = GEVparam["shape"])

  # define dataframe with 'rl' and 'q'
  df = data.frame(rl = rl, q = q)

  
  # define plottitle
  if (missing(plottitle)) {
    plottitle = "return level plot"
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
    plot = ggplot(df,aes(x = q, y = rl)) +
      geom_line() +
      scale_x_log10(name = "return period",
                    breaks = b,
                    labels = l) +
      scale_y_continuous(name = "return level") +
      geom_label(aes(label = GEV, x = b[1], y = max(rl, na.rm = TRUE)), hjust = 0, vjust = 1, size = 3) +
      ggtitle(plottitle)
  } else {
    r = quantile(obs, 1-p, names = FALSE)
    df.obs = data.frame(rl = r, q = q)
    plot = ggplot(df,aes(x = q, y = rl)) +
      geom_line() +
      geom_point(data = df.obs, aes(x = q, y = rl)) +
      scale_x_log10(name = "return period",
                    breaks = b,
                    labels = l) +
      scale_y_continuous(name = "return level") +
      geom_label(aes(label = GEV, x = b[1], y = max(max(df$rl, na.rm = TRUE), max(df.obs$rl, na.rm = TRUE))), hjust = 0, vjust = 1, size = 3) +
      ggtitle(plottitle)
  }

  # define save name
  if (!missing(save_name)) {
    plotname = paste0(save_name,".pdf",sep = "")
    # plot is saved to given directory
    path = save_dir
    message(sprintf("plot was saved as '%s' to directory: \n  '%s'",plotname,path))
    # save plot
    ggsave(filename = plotname, plot = plot, path = path, width = 10, height = 5.2)
    
  }

  # print plot if 'printPlot' is TRUE
  if (printPlot) {
    print(plot)
  }

  #---------------------------------------------------------------------------------------------------------------#

}
