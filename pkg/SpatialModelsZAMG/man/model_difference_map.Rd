\name{model_difference_map}
\alias{model_difference_map}
\title{Create model difference maps for Austria}
\description{this function creates a map of Austria with the relative difference of return levels for two given models. the relative difference is calculated like:
\deqn{\frac{{\rm{rl\_model\_1}} - {\rm{rl\_model\_2}}}{{\rm{rl\_model\_1}}}\, .}{(rl_model_1 - rl_model_2)/rl_model_1 .}
}
\usage{
model_difference_map(covariables, rl_model_1, rl_model_2,
                     name_model_1 = "model_1",
                     name_model_2 = "model_2",
                     plottitle = NULL, save_name = NULL,
                     save_dir = getwd(), printPlot = TRUE)}
\arguments{
  \item{covariables}{a named matrix with the covariables. each row corresponds to one location, columns should include at least \var{lon} and \var{lat}}
  \item{rl_model_1}{a vector with the return level for every location and the first model}
  \item{rl_model_2}{a vector with the return level for every location and the second model}
  \item{name_model_1}{a character string defining the name of the first model. \cr default
  is \code{name_model_1 = "model_1"}}
  \item{name_model_2}{a character string defining the name of the second model. \cr default
  is \code{name_model_2 = "model_2"}}
    \item{plottitle}{a character string defining the title of the plot. default is \cr \code{plottitle = "relative difference, name_model_1 - name_model_2"}, depending on \code{name_model_1} and \code{name_model_2}. if no plottitle is wanted, use \code{plottitle = ""}}
  \item{save_name}{a character string defining the saving name of the map. }
  \item{save_dir}{a character string defining the directory for the map to be saved.\cr default is the working directory}
  \item{printPlot}{logical value; if \code{TRUE} (default), the plot is printed}

}
\value{a map of Austria with the relative difference of return levels for two given models}
\seealso{
\code{\link{returnlevels}}}
\examples{
# get covariables
lon = get(data(lon.at))
lon = as.vector(t(lon))
lat = get(data(lat.at))
lat = as.vector(t(lat))
alt = get(data(alt.at))
alt = as.vector(t(alt))

x = get(data(sample_grid_data))

mdday    = x$mdday
sd_mmax  = x$mmsd
swe_mmax = x$mmswe

# take only locations of the (Austrian) domain
lon = lon[which(!is.na(mdday))]
lat = lat[which(!is.na(mdday))]
alt = alt[which(!is.na(mdday))]
mdday    = mdday[which(!is.na(mdday))]
sd_mmax  = sd_mmax[which(!is.na(sd_mmax))]
swe_mmax = swe_mmax[which(!is.na(swe_mmax))]

# define matrix 'covariables'
covariables = cbind("lon" = lon, "lat" = lat, "alt" = alt,
                    "mdday" = mdday, "sd_mmax" = sd_mmax,
                    "swe_mmax" = swe_mmax)

# load function output from model_selection
sd_m_select  = get(data("sd_m_select"))
swe_m_select = get(data("swe_m_select"))

# perform optimization with Ext-Gaussian model
optim_gauss =
  optimizer_biv_ext_gauss_model(sd_m_select = sd_m_select,
                                swe_m_select = swe_m_select,
                                method = "bobyqa")
sd_coeff_gauss = optim_gauss$coefficients$sd_coeff

# load optimization results from HR model
data(optim_hr)
sd_coeff_hr = optim_hr$coefficients$sd_coeff

# calculate GEV parameters from linear models
sd_GEVparam_gauss = 
  GEVparameters_from_models(covariables = covariables, 
                            coefficients = sd_coeff_gauss)
sd_GEVparam_hr = 
  GEVparameters_from_models(covariables = covariables, 
                            coefficients = sd_coeff_hr)

# calculate return levels
q = 100
sd_rl_gauss  = 
  returnlevels(GEVparam = sd_GEVparam_gauss,  q = q)
sd_rl_hr      = 
  returnlevels(GEVparam = sd_GEVparam_hr,  q = q)

# create relative difference map
model_difference_map(covariables = covariables, 
                     rl_model_1 = sd_rl_hr, 
                     rl_model_2 = sd_rl_gauss,
                     printPlot = FALSE)

}
\keyword{functions}
\keyword{plots}
\keyword{modelselection}