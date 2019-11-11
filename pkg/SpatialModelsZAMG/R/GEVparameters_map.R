#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: GEVparameters_map --------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

GEVparameters_map = function(covariables, GEVparam, sd_or_swe, parameter = "shape",
                             plottitle = NULL, save_name = NULL, save_dir = getwd(), 
                             printPlot = TRUE) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: GEVparameters_map(covariables = covariables, GEVparam = sd_GEVparam, sd_or_swe = "sd", 
  #                            parameter = "shape", printPlot = FALSE)
  # output: a map of Austria with the chosen GEV parameters
  #
  # --- this function creates a map of Austria with the chosen GEV parameter 
  #
  # --- input:
  #           1. 'covariables': a named matrix with the covariables
  #                 each row corresponds to one location, 
  #                 columns should include at least 'lon' and 'lat'
  #           2. 'GEVparam': a matrix with the GEV parameters
  #                 each row corresponds to one location
  #                 columns are 'loc', 'scale' and 'shape'
  #           3. 'sd_or_swe': a character string
  #                 chose snow depth ('sd') or snow water equivalent ('swe')
  # --- optional input:
  #           4. 'parameter': a character string
  #                 chose which GEV parameter should be plotted, either 'loc', 'scale' or 'shape'
  #                 default (if this input is missing) is 'shape'
  #           5. 'plottitle': a character string defining the title of the plot
  #                 default (if this input is missing) is 'sd_or_swe shape/scale/loc parameter', 
  #                   depending on 'sd_or_swe' and 'parameter'
  #           6. 'save_name': a character string defining the saving name of the map
  #                 default (if this input is missing) is 'sd_or_swe_parameter_map', 
  #                   depending on 'sd_or_swe' and 'parameter'
  #           7. 'save_dir': a character string defining the directory for the map to be saved
  #                 default (if this input is missing) is the working directory
  #           8. 'print_plot': logical value
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
  
  # check whether 'GEVparam' is a matrix
  if (!is.matrix(GEVparam)) {
    stop("'GEVparam' has to be a matrix")
  }
  
  # check whether the number of rows in 'covariables' and 'GEVparam' coincide
  if (nrow(covariables) != nrow(GEVparam)) {
    stop(sprintf("number of rows (%i) in 'covariables' and number of rows (%i) in 'GEVparam' don't match up", 
                 nrow(covariables), nrow(GEVparam)))
  }
  
  # columns of 'covariables' have to be named and include at least 'lon' and 'lat'
  if (length(colnames(covariables)) != ncol(covariables)) {
    stop("columns of 'covariables' have to be named")
  }
  if (!(all(c("lon","lat") %in% colnames(covariables)))) {
    stop("colnames of 'covariables' have to include at least 'lon' and 'lat'")
  }
  
  # columns of 'GEVparam' have to be 'loc', 'scale' and 'shape'
  if(!(all(c("loc","scale","shape") %in% colnames(GEVparam)))) {
    stop("columns of 'GEVparam' have to be 'loc', 'scale' and 'shape'")
  }
  
  # 'sd_or_swe' has to be either 'sd' or 'swe'
  if (!(sd_or_swe %in% c("sd","swe"))) {
    stop(sprintf("'sd_or_swe' has to be either 'sd' or 'swe' -- '%s' is not allowed",sd_or_swe))
  }
  
  # 'parameter' has to be either 'loc', 'scale' or 'shape'
  if (!(parameter %in% c("loc","scale","shape"))) {
    stop(sprintf("'parameter' has to be either 'loc', 'scale' or 'shape' -- '%s' is not allowed",parameter))
  }
  
  # 'printPlot' has to be a logical value
  if (!(is.logical(printPlot))) {
    stop(sprintf("'printPlot' has to be TRUE or FALSE -- '%s' is not allowed",printPlot))
  }
  
  
  lon <- lat <- border.at <- long <- group <- NULL
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Perform calculations ------------------------------------------------------------------------------------------
  
  # define dataframe with 'lon', 'lat' and the chosen GEV parameter
  df = data.frame(lon = covariables[,"lon"], lat = covariables[,"lat"], param = GEVparam[,parameter])
  
  # define different breaks for the different parameters and 'sd' or 'swe'
  if (parameter == "shape") {
    b = c(min(df$param,na.rm = TRUE) - 0.01, -0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)
  } else if (parameter == "loc" && sd_or_swe == "sd") {
    b = c(min(df$param,na.rm = TRUE) - 0.1,10,20,30,40,60,80,100,150,200,500,1000)
  } else if (parameter == "loc" && sd_or_swe == "swe") {
    b = c(min(df$param,na.rm = TRUE) - 0.1,15,30,50,75,100,150,200,500,1000,2000,3000)
  } else if (parameter == "scale" && sd_or_swe == "sd") {
    b = c(0,5,10,15,20,30,40,50,100,150,300,500)
  } else if (parameter == "scale" && sd_or_swe == "swe") {
    b = c(0,10,20,30,40,60,80,100,150,200,800,1500)
  }
  
  # define plottitle, colors, labels and levels
  if (missing(plottitle)) {
    plottitle = paste0(sd_or_swe," ",parameter," parameter",sep = "")
  }
  col_val   = c("#000000","#1B244B","#333B62","#48548B","#5767AC","#A4AAC9",
                "#CCA0A8","#B8707F","#A84E64","#9A304E","#8E063B","#690129")
  l         = b[-1]
  df$p      = cut(df$param,breaks=b,labels=l)
  names(df) = c(names(df)[-4],parameter)
  
  message(paste0("range and distribution of ",sd_or_swe," ",parameter," parameter:",sep = ""))
  print(summary(df[,4]))
  
  # load country borders
  data("border.at",envir = environment())
  
  # create plot
  plot = ggplot() + 
    geom_point(data = df,aes(lon,lat,col=df[,4]), size = 0.5, shape = 15) +
    scale_color_manual(values = col_val,
                       limits = b,
                       breaks = rev(l),
                       name = names(df)[4],
                       na.value = NA) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    geom_polygon(data = border.at, aes(x = long,y = lat, group = group), fill = NA, size = 0.5, col = "#585858") +
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