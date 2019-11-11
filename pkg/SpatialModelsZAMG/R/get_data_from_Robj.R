#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#---------------------------------------- function: get_data_from_Robj -------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

get_data_from_Robj = function(dataname, measurements_barrier = 10) {
  # Information ---------------------------------------------------------------------------------------------------
  # example: get_data_from_Robj(dataname = "snow_data.robj", measurements_barrier = 3)
  # output: a list with 'snow_data', 'sd_max_data', 'swe_max_data' and 'covariables' (see output details below)
  #
  # --- this function extracts the snow data from a saved robj file 
  #           which consists of a list of stations (see below for an example)
  #
  # --- input:
  #           1. 'dataname': the dataname of the saved robj file as a character string
  # --- optional input:
  #           2. 'measurements_barrier': how many measurements does each station has to have at least
  #                 those stations are used later on in the optimization (see functions 'optimizer_..._model')
  #                 default (if this input is missing) is 10 
  # --- output: a list with
  #           1. the data of the saved file (a list): 'snow_data'
  #                 all stations with less than 'measurements_barrier' measurements have been deleted
  #           2. a matrix with the yearly maxima of snow depth (sd): 'sd_max_data'
  #                 each row corresponds to one station, 
  #                 columns are the corresponding years
  #                 matrix might contain NA's
  #           3. a matrix with the yearly maxima of snow water equivalent (swe): 'swe_max_data'
  #                 each row corresponds to one station, 
  #                 columns are the corresponding years
  #                 matrix might contain NA's
  #           4. a matrix with the covariables for each station: 'covariables'
  #                 each row corresponds to one station, 
  #                 columns are 'lon' (longitude), 'lat' (latitude), 'alt' (altitude), 
  #                     'mdday' (mean of the time difference in days of sd and swe maxima),
  #                     'sd_mmax' (mean maxima of sd over all years) and
  #                     'swe_mmax' (mean maxima of swe over all years)
  #
  # --- the saved robj file should include a list of stations with 
  #         dataframes/matrices 'meta' and 'data' according to this example:
  #
  # --- Example: Station UNKEN in Salzburg
  #
  # --- $'5045303'
  # --- $'5045303'$meta
  # ---      stnr stname   region provider   lon      lat alt length ybeg yend missing_years
  # --- 1 5045303  UNKEN SALZBURG     ZAMG 12.75 47.66667 545     42 1948 1989             0
  # 
  # --- $'5045303'$data
  # ---    year  sd    swe    date.sd   date.swe
  # --- 1  1948  41  82.57 1949-02-02 1949-02-02
  # --- 2  1949  40  52.31 1950-01-05 1950-02-09
  # --- .   ... ...    ...        ...        ...
  # ---
  # --- end example
  #---------------------------------------------------------------------------------------------------------------#
  
  # Check required packages and input parameters ------------------------------------------------------------------
  
  # # load required package 'lubridate'
  # if (inherits(try(library(lubridate, warn.conflicts = FALSE, quietly = TRUE), silent = TRUE), "try-error")) {
  #   message("required package 'lubridate' is not installed yet -- trying to install package")
  #   install.packages("lubridate", quiet = TRUE)
  #   if (inherits(try(library(lubridate, warn.conflicts = FALSE, quietly = TRUE), silent = TRUE), "try-error")) {
  #     stop("package 'lubridate' couldn't be installed")
  #   } else {
  #     message("package successfully installed and loaded")
  #   }
  # }
  
  # 'measurements_barrier' has to be a positive integer
  if (measurements_barrier < 0 || (measurements_barrier != as.integer(measurements_barrier))) {
    stop(sprintf(
      "'measurements_barrier' has to be a positive integer -- '%s' is not allowed",measurements_barrier))
  }
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Perform calculations ------------------------------------------------------------------------------------------
  
  # load data and define the number of stations with q
  snow_data = get(load(dataname))
  q = length(snow_data)
  
  # find the number of data for each station
  nrow_data = rep(NA,times = q)
  for (k in 1:q) {
    nrow_data[k] = nrow(snow_data[[k]]$data)
  }
  
  # find all stations with less than b measurements and delete them
  b = measurements_barrier
  delete = which(nrow_data < b)
  if (length(delete) > 0) {
    snow_data = snow_data[-delete]
    nrow_data = nrow_data[-delete]
    q = length(snow_data)
  }
  if (q == 0) {
    stop(sprintf(
      "all stations have less than %i measurements -- 'measurements_barrier' being %i is too high",b,b))
  }
  
  # define a vector with stationnumbers and 
  #     find the first and the last year for which a value exists
  statnum = rep(NA,times = q)
  ybegs   = rep(NA,times = q)
  yends   = rep(NA,times = q)
  for (k in 1:q) {
    statnum[k] = snow_data[[k]]$meta$stnr
    ybegs[k]   = min(as.numeric(as.vector(snow_data[[k]]$data[,"year"])), na.rm = TRUE)
    yends[k]   = max(as.numeric(as.vector(snow_data[[k]]$data[,"year"])), na.rm = TRUE)
  }
  ystart = min(ybegs, na.rm = TRUE)
  ystop  = max(yends, na.rm = TRUE)
  
  # define the extreme value matrices 'sd_max_data' and 'swe_max_data'
  nyears = length(ystart:ystop)
  sd_max_data  = matrix(NA, nrow = q, ncol = nyears, dimnames = list(statnum, ystart:ystop))
  swe_max_data = matrix(NA, nrow = q, ncol = nyears, dimnames = list(statnum, ystart:ystop))
  for (k in 1:q) {
    # define vector with all years of the station for which data exists
    statyears = as.character(snow_data[[k]]$data$year)
    # get the data for max_data matrices
    sd_max_data[k,statyears]  = as.numeric(as.vector(snow_data[[k]]$data$sd))
    swe_max_data[k,statyears] = as.numeric(as.vector(snow_data[[k]]$data$swe))
  }
  
  # define the covariable matrix 'covariables'
  covariables = matrix(NA, nrow = q, ncol = 6, dimnames = 
                         list(statnum, c("lon","lat","alt","mdday","sd_mmax","swe_mmax")))
  for (k in 1:q) {
    diff_days = abs(ymd(as.character(snow_data[[k]]$data$date.sd)) - 
                      ymd(as.character(snow_data[[k]]$data$date.swe)))
    mdday     = mean(diff_days, na.rm = TRUE)
    sd_mmax   = mean(as.numeric(as.vector(snow_data[[k]]$data$sd)),  na.rm = TRUE)
    swe_mmax  = mean(as.numeric(as.vector(snow_data[[k]]$data$swe)), na.rm = TRUE)
    covariables[k,] = c(as.numeric(as.vector(snow_data[[k]]$meta[c("lon","lat","alt")])), 
                        as.numeric(mdday), sd_mmax, swe_mmax)
  }
  
  #---------------------------------------------------------------------------------------------------------------#
  
  # Define function output ----------------------------------------------------------------------------------------
  
  ans = list(snow_data = snow_data, sd_max_data = sd_max_data, 
             swe_max_data = swe_max_data, covariables = covariables)
  
  return(ans)
  
  #---------------------------------------------------------------------------------------------------------------#
}