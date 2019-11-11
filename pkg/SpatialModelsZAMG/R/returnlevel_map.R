#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: returnlevel_map ----------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

returnlevel_map = function(covariables, rl, q, sd_or_swe, plottitle = NULL,
                           save_name = NULL, save_dir = getwd(), printPlot = TRUE) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: returnlevel_map(covariables = covariables, rl = sd_rl, q = 100, sd_or_swe = "sd", printPlot = FALSE)
  # output: a map of Austria with the return levels
  #
  # --- this function creates a map of Austria with return levels
  #
  # --- input:
  #           1. 'covariables': a named matrix with the covariables
  #                 each row corresponds to one location, 
  #                 columns should include at least 'lon' and 'lat'
  #           2. 'rl': a vector with the return level for every location
  #           3. 'q': the return period -- must be a number greater than 1
  #           4. 'sd_or_swe': a character string
  #                 chose snow depth ('sd'), snow water equivalent ('swe') or the quotient swe/sd ('quot')
  # --- optional input:
  #           5. 'plottitle': a character string defining the title of the plot
  #                 default (if this input is missing) is 'sd_or_swe return level map', depending on 'sd_or_swe'
  #           6. 'save_name': a character string defining the saving name of the map
  #                 default (if this input is missing) is 'sd_rl_map', 'swe_rl_map' or 'swe_per_sd_rl_map',
  #                   depending on 'sd_or_swe'
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
  
  # columns of 'covariables' have to be named and include at least 'lon' and 'lat'
  if (length(colnames(covariables)) != ncol(covariables)) {
    stop("columns of 'covariables' have to be named")
  }
  if (!(all(c("lon","lat") %in% colnames(covariables)))) {
    stop("colnames of 'covariables' have to include at least 'lon' and 'lat'")
  }
  
  # check whether the number of rows in 'covariables' and the length of 'rl' coincide
  if (nrow(covariables) != length(rl)) {
    stop(sprintf("number of rows (%i) in 'covariables' and length of 'rl' (%i) don't match up", 
                 nrow(covariables), length(rl)))
  }
  
  # 'q' has to be a number greater than 1
  if (q <= 1) {
    stop(sprintf("'q' has to be a number greater than 1 -- '%s' is not allowed",q))
  }
  
  # 'sd_or_swe' has to be either 'sd', 'swe' or 'quot'
  if (!(sd_or_swe %in% c("sd","swe","quot"))) {
    stop(sprintf("'sd_or_swe' has to be either 'sd', 'swe' or 'quot' -- '%s' is not allowed",sd_or_swe))
  }
  
  # 'printPlot' has to be a logical value
  if (!(is.logical(printPlot))) {
    stop(sprintf("'printPlot' has to be TRUE or FALSE -- '%s' is not allowed",printPlot))
  }
  
  lon <- lat <- border.at <- long <- group <- NULL
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Perform calculations ------------------------------------------------------------------------------------------
  
  # define dataframe with 'lon', 'lat' and 'rl'
  df = data.frame(lon = covariables[,"lon"], lat = covariables[,"lat"], rl = rl)
  
  # define different breaks, colors and plottitles for 'sd', 'swe' and 'quot', respectively
  if (sd_or_swe == "sd") {
    b = c(0,50,100,200,300,400,500,600,700,800,900,1000,1200,1400)
    pal_returnlevels =
      c("#000000","#1B244B","#333B62","#48548B","#5767AC","#7782B7","#A4AAC9",
        "#CCA0A8","#C88592","#B8707F","#A84E64","#9A304E","#8E063B","#690129")
    if (missing(plottitle)) {
      plottitle = "sd return level map"
    }
  } else if (sd_or_swe == "swe") {
    b = c(0,100,150,200,350,500,700,900,1200,1500,2000,2500,3000,5000)
    pal_returnlevels =
      c("#000000","#1B244B","#333B62","#48548B","#5767AC","#7782B7","#A4AAC9",
        "#CCA0A8","#C88592","#B8707F","#A84E64","#9A304E","#8E063B","#690129")
    if (missing(plottitle)) {
      plottitle = "swe return level map"
    }
  } else if (sd_or_swe == "quot") {
    b = c(min(df$rl) - 0.01,1.8,1.9,2,2.1,2.3,2.5,3.5) 
    pal_returnlevels =
      c("#000000","#CCA0A8","#C88592","#B8707F","#A84E64","#9A304E","#8E063B","#690129")
    if (missing(plottitle)) {
      plottitle = "swe/sd return level map"
    }
  }
  
  # define labels and levels
  l         = b[-1]
  df$rl_q   = cut(df$rl,breaks=b,labels=l)
  names(df) = c(names(df)[-4],paste0("rl", q, sep = ""))
  
  if (sd_or_swe == "quot") {
    message(paste0("range and distribution of swe per sd return levels with q  = ",q,":",sep = ""))
  } else {
    message(paste0("range and distribution of ",sd_or_swe," return levels with q = ",q,":",sep = ""))
  }
  print(summary(df[,4]))
  
 
  # load country borders
  data("border.at",envir = environment())
  
  # create plot
  plot = ggplot() +
    geom_point(data = df,aes(lon,lat,col=df[,4]), size = 0.5, shape = 15) +
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