#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: model_difference_map -----------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

model_difference_map = function(covariables, rl_model_1, rl_model_2, name_model_1 = "model_1",
                                name_model_2 = "model_2", plottitle = NULL,
                                save_name = NULL, save_dir = getwd(),
                                printPlot = TRUE) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: model_difference_map(covariables = covariables, rl_model_1 = sd_rl100_hr,
  #                               rl_model_2 = sd_rl100_smooth, name_model_1 = "hr", name_model_2 = "smooth")
  # output: a map of Austria with the relative difference of return levels for two given models
  #
  # --- this function creates a map of Austria with the relative difference of return levels for two given models
  #           the relative difference is calculated like: (rl_model_1 - rl_model_2)/rl_model_1
  #
  # --- input:
  #           1. 'covariables': a named matrix with the covariables
  #                 each row corresponds to one location,
  #                 columns should include at least 'lon' and 'lat'
  #           2. 'rl_model_1': a vector with the return level for every location and the first model
  #           3. 'rl_model_2': a vector with the return level for every location and the second model
  # --- optional input:
  #           4. 'name_model_1': a character string defining the name of the first model
  #                 default (if this input is missing) is 'model_1'
  #           5. 'name_model_2': a character string defining the name of the second model
  #                 default (if this input is missing) is 'model_2'
  #           6. 'plottitle': a character string defining the title of the plot
  #                 default (if this input is missing) is 'relative difference, name_model_1 - name_model_2',
  #                   depending on 'name_model_1' and 'name_model_2'
  #           7. 'save_name': a character string defining the saving name of the map
  #                 default (if this input is missing) is 'model_difference_map'
  #           8. 'save_dir': a character string defining the directory for the map to be saved
  #                 default (if this input is missing) is the working directory
  #           9. 'print_plot': logical value
  #                 if TRUE, the plot is printed
  #                 default (if this input is missing) is TRUE
  #---------------------------------------------------------------------------------------------------------------#

  # Check required packages and input parameters ------------------------------------------------------------------

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

  # check whether 'covariables' is a matrix
  if (!is.matrix(covariables)) {
    stop("'covariables' has to be a matrix")
  }

  # columns of 'covariables' have to be named and include at least 'lon' and 'lat'
  if (length(colnames(covariables)) != ncol(covariables)) {
    stop("columns of 'covariables' have to be named")
  }
  if (!(all(c("lon","lat") %in% colnames(covariables)))) {
    stop("colnames of 'covariables' have to include at least 'lon' and 'lat'")
  }

  # check whether the number of rows in 'covariables' and the length of 'rl_model_1' coincide
  if (nrow(covariables) != length(rl_model_1)) {
    stop(sprintf("number of rows (%i) in 'covariables' and length of 'rl_model_1' (%i) don't match up",
                 nrow(covariables), length(rl_model_1)))
  }

  # check whether the length of 'rl_model_1' and 'rl_model_2' coincide
  if (length(rl_model_1) != length(rl_model_2)) {
    stop(sprintf("length of 'rl_model_1' (%i) and length of 'rl_model_2' (%i) don't match up",
                 length(rl_model_1), length(rl_model_2)))
  }

  # 'printPlot' has to be a logical value
  if (!(is.logical(printPlot))) {
    stop(sprintf("'printPlot' has to be TRUE or FALSE -- '%s' is not allowed",printPlot))
  }

  lon <- lat <- border.at <- long <- group <- NULL
  
  #---------------------------------------------------------------------------------------------------------------#

  # Perform calculations ------------------------------------------------------------------------------------------

  # define dataframe with 'lon', 'lat' and 'd'
  df = data.frame(lon = covariables[,"lon"], lat = covariables[,"lat"], d = (rl_model_1 - rl_model_2)/rl_model_1)

  # define breaks, colors and plottitle
  if (missing(plottitle)) {
    plottitle = paste0("relative difference, ",name_model_1," - ",name_model_2)
  }
  b = c(-1,-0.6,-0.4,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,0.4,0.6)
  pal_returnlevels = c("#000000","#1B244B","#48548B","#7782B7", "#A4AAC9",
                       "#004512","#299443","#73DD8D","#93FCAD",
                       "#CCA0A8","#B8707F","#9A304E","#690129")

  # define labels and levels
  l         = b[-1]
  df$diff   = cut(df$d,breaks=b,labels=l)

  message("range and distribution of the relative difference:")
  print(summary(df$diff))

  # load country borders
  data("border.at",envir = environment())

  # create plot
  plot = ggplot() +
    geom_point(data = df,aes(lon,lat,col=diff), size = 0.5, shape = 15) +
    scale_color_manual(values = pal_returnlevels,
                       limits = b,
                       breaks = rev(l),
                       name = names(df)[4],
                       na.value = NA) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    geom_polygon(data = border.at, aes(x = long,y = lat, group = group), fill = NA, size = 0.7, col = "#585858") +
    ggtitle(plottitle)

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
