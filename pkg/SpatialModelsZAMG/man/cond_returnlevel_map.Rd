\name{cond_returnlevel_map}
\alias{cond_returnlevel_map}
\title{Create conditional return level maps for Austria}
\description{this function creates a map of Austria with conditional return levels}
\usage{
cond_returnlevel_map(covariables, cond_rl, q, sd_or_swe,
                     plottitle = NULL, save_name = NULL,
                     save_dir = getwd(), printPlot = TRUE)}
\arguments{
  \item{covariables}{a named matrix with the covariables. each row corresponds to one location, columns should include at least \var{lon} and \var{lat}}
  \item{cond_rl}{a vector with the conditional return level for every location}
  \item{q}{the return period -- must be a number greater than 1}
  \item{sd_or_swe}{a character string; you can either chose snow depth (\code{"sd"}) or snow water equivalent (\code{"swe"})}
    \item{plottitle}{a character string defining the title of the plot. \cr default is \code{plottitle = "sd_or_swe conditional return level map"}, depending on \code{sd_or_swe}. if no plottitle is wanted, use \code{plottitle = ""}}
  \item{save_name}{a character string defining the saving name of the map.}
  \item{save_dir}{a character string defining the directory for the map to be saved. \cr default is the working directory}
  \item{printPlot}{logical value; if \code{TRUE} (default), the plot is printed}
}
\value{a map of Austria with the conditional return levels}
\seealso{
\code{\link{cond_returnlevels}}, \code{\link{cond_returnlevel_plot}}, \code{\link{returnlevel_map}}
}
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

# load function output from GEVparameters_from_models
sd_GEVparam  = get(data("sd_GEVparam"))

# load function output from optimizer_biv_hr_model
data("optim_hr")

# define coefficients
sd_coeff  = optim_hr$coefficients$sd_coeff
cor_coeff = optim_hr$coefficients$cor_coeff

# define return period (same as in sd_GEVparam)
q = 100

# calculate return levels
sd_rl  = returnlevels(GEVparam = sd_GEVparam,  q = q)

# sd given sd in ibk
ibk = c("lon" = 11.392778, "lat" = 47.267222, "alt" = 574,
        "mdday" = 6.12, "sd_mmax" = 16.8, "swe_mmax" = 29)
sd_GEVparam_ibk  = GEVparameters_from_models(ibk, sd_coeff)
sd_rl_ibk = returnlevels(GEVparam = sd_GEVparam_ibk, q = q)

cond_rl_ibk =
  cond_returnlevels(locations = covariables,
                    GEVparam = sd_GEVparam, q = q,
                    cond_locations = ibk,
                    cond_GEVparam = sd_GEVparam_ibk,
                    same_var = TRUE, cond_B = sd_rl_ibk,
                    cor_coeff = cor_coeff, model = "hr")

# create conditional return level map
cond_returnlevel_map(covariables = covariables,
                     cond_rl = cond_rl_ibk, q = q,
                     sd_or_swe = "sd", printPlot = TRUE)
}
\keyword{functions}
\keyword{conditionals}
\keyword{plots}
