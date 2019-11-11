\name{returnlevel_map}
\alias{returnlevel_map}
\title{Create return level maps for Austria}
\description{this function creates a map of Austria with return levels}
\usage{
returnlevel_map(covariables, rl, q, sd_or_swe,
                plottitle = NULL, save_name = NULL,
                save_dir = getwd(), printPlot = TRUE)}
\arguments{
  \item{covariables}{a named matrix with the covariables. each row corresponds to one location, columns should include at least \var{lon} and \var{lat}}
  \item{rl}{a vector with the return level for every location}
  \item{q}{the return period -- must be a number greater than 1}
  \item{sd_or_swe}{a character string; you can either chose snow depth (\code{"sd"}), snow water equivalent (\code{"swe"}) or the quotient swe/sd (\code{"quot"})}
  \item{plottitle}{a character string defining the title of the plot. \cr default is \code{plottitle = "sd_or_swe return level map"}, depending on \code{sd_or_swe}. if no plottitle is wanted, use \code{plottitle = ""}}
  \item{save_name}{a character string defining the saving name of the map. }
  \item{save_dir}{a character string defining the directory for the map to be saved. \cr default is the working directory}
  \item{printPlot}{logical value; if \code{TRUE} (default), the plot is printed}
}
\value{a map of Austria with the return levels}
\seealso{
\code{\link{returnlevels}}, \code{\link{returnlevel_plot}}, \code{\link{cond_returnlevel_map}}
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
swe_GEVparam = get(data("swe_GEVparam"))

# calculate returnlevels
q = 100
sd_rl  = returnlevels(GEVparam = sd_GEVparam,  q = q)
swe_rl = returnlevels(GEVparam = swe_GEVparam, q = q)

# create return level maps
returnlevel_map(covariables = covariables,
                rl = sd_rl, q = q,
                sd_or_swe = "sd",
                printPlot = FALSE)
returnlevel_map(covariables = covariables,
                rl = swe_rl, q = q,
                sd_or_swe = "swe",
                printPlot = FALSE)
}
\keyword{functions}
\keyword{plots}
