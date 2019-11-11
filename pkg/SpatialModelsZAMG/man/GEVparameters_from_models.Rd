\name{GEVparameters_from_models}
\alias{GEVparameters_from_models}
\title{Calculate GEV parameters according to given linear models}
\description{this function calculates the GEV parameters according to the linear models \cr specified by \code{coefficients}}
\usage{
GEVparameters_from_models(covariables, coefficients)}
\arguments{
  \item{covariables}{a named matrix or vector with the covariables. each row corresponds to one location, columns are the corresponding covariables (names of the coefficients) from the linear models; the intercept is added seperately and thus doesn't have to be included}
  \item{coefficients}{a list with the named coefficients of the linear models including\cr
  \var{loccoeff} (coefficients of location parameter) -- must be named\cr
  \var{scalecoeff} (coefficients of scale parameter) -- must be named\cr
  \var{shapecoeff} (coefficients of shape parameter) -- must be named}
}
\value{
  \item{GEVparam }{a matrix or vector with the GEV parameters. each row corresponds to one location, columns are \var{loc}, \var{scale} and \var{shape}}
}
\seealso{
\code{\link{GEVparameters_map}}, \code{\link{returnlevels}}, \code{\link{cond_returnlevels}}, \code{\link{returnlevel_plot}}, \cr \code{\link{cond_returnlevel_plot}}
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

# load function output from optimizer_biv_hr_model
data("optim_hr")

# define coefficients
sd_coeff  = optim_hr$coefficients$sd_coeff
swe_coeff = optim_hr$coefficients$swe_coeff
cor_coeff = optim_hr$coefficients$cor_coeff

# Calculate the GEV parameters from linear models
sd_GEVparam =
  GEVparameters_from_models(covariables = covariables,
                            coefficients = sd_coeff)

swe_GEVparam =
  GEVparameters_from_models(covariables = covariables,
                            coefficients = swe_coeff)
}
\keyword{functions}
\keyword{GEVparameters}
