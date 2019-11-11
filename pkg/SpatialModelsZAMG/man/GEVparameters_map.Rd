\name{GEVparameters_map}
\alias{GEVparameters_map}
\title{Create GEV parameter maps for Austria}
\description{this function creates a map of Austria with the chosen GEV parameter}
\usage{
GEVparameters_map(covariables, GEVparam, sd_or_swe,
                  parameter = "shape", plottitle = NULL,
                  save_name = NULL, save_dir = getwd(),
                  printPlot = TRUE)}
\arguments{
  \item{covariables}{a named matrix with the covariables. each row corresponds to one location, columns should include at least \var{lon} and \var{lat}}
  \item{GEVparam}{a matrix with the GEV parameters. each row corresponds to one location, columns are \var{loc}, \var{scale} and \var{shape}}
  \item{sd_or_swe}{a character string; you can either chose snow depth (\code{"sd"}) or snow water equivalent (\code{"swe"})}
  \item{parameter}{a character string; chose which GEV parameter should be plotted; either \code{"loc"}, \code{"scale"} or \code{"shape"} (default)}
  \item{plottitle}{a character string defining the title of the plot. \cr default is \code{plottitle = "sd_or_swe shape/scale/loc parameter"}, depending on \code{sd_or_swe} and \code{parameter}. if no plottitle is wanted, use \code{plottitle = ""}}
  \item{save_name}{a character string defining the saving name of the map. }
  \item{save_dir}{a character string defining the directory for the map to be saved.\cr default is the working directory}
  \item{printPlot}{logical value; if \code{TRUE} (default), the plot is printed}

}
\value{a map of Austria with the chosen GEV parameters}
\seealso{
\code{\link{GEVparameters_from_models}}}
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

# create GEV parameter map
GEVparameters_map(covariables = covariables,
                  GEVparam = sd_GEVparam, sd_or_swe = "sd",
                  parameter = "shape", printPlot = FALSE)
}
\keyword{functions}
\keyword{plots}
\keyword{GEVparameters}
