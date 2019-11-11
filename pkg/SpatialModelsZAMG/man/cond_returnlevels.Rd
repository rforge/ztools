\name{cond_returnlevels}
\alias{cond_returnlevels}
\title{Calculate conditional return levels}
\description{this function calculates conditional return levels. \cr

the \eqn{q}-year conditional return level of variable \eqn{Z_1}{Z1} at location \eqn{x_1}{x1} given variable \eqn{Z_2}{Z2} at location \eqn{x_2}{x2}, is defined as the threshold \eqn{B}, such that the conditional probability that \eqn{Z_1(x_1)}{Z1(x1)} exceeds this threshold, given that \eqn{Z_2(x_2)}{Z2(x2)} is in the interval \eqn{({\rm{cond\_B}}\,,\,{\rm{cond\_B\_2}})}{(cond_B,cond_B_2)}, is \eqn{1/q} :
\deqn{{\rm{Pr}}\bigl[Z_1(x_1) > B\, \big|\, Z_2(x_2) \in ({\rm{cond\_B}}\,,\,{\rm{cond\_B\_2}})\bigr] = \frac{1}{q} \, .}{Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)] = 1/q .}
}
\usage{
cond_returnlevels(locations, GEVparam, q,
                  cond_locations, cond_GEVparam, same_var,
                  cor_coeff, cond_B, cond_B_2 = Inf, 
                  var = NULL, model = "ext-t", 
                  printObjectives = FALSE)
}
\arguments{
  \item{locations}{
  a matrix or vector with certain location characteristics as columns/entries.
  each row corresponds to one location, columns are for example: longitude, latitude and altitude}
  \item{GEVparam}{a named matrix or vector with the GEV parameters of the variable for which the conditional return levels are calculated.
  each row corresponds to one location (same ones as in \code{locations}), columns are \var{loc}, \var{scale} and \var{shape}}
  \item{q}{the return period for the calculation of the conditional return levels -- must be a number greater than 1}
  \item{cond_locations}{a matrix or a vector with the location characteristics of the conditioned locations.
  each row corresponds to one location, columns should be the same as in \code{locations}}
  \item{cond_GEVparam}{a named matrix or vector with the GEV parameters of the conditioned variable.
  each row corresponds to one location (same ones as in \code{cond_locations}), columns are \var{loc}, \var{scale} and \var{shape}}
    \item{same_var}{logical value; if \code{TRUE}, the conditioned variable is the same as the unconditioned variable. this has to be known in order to use the right correlation function}
  \item{cor_coeff}{a named vector with the correlation parameters. \cr
  for the Huesler-Reiss model: \var{alpha}, \var{kappa} and \var{lambda12} ; \cr
  for the Extremal-Gaussian model: \var{alpha}, \var{sd_kappa}, \var{swe_kappa} and \var{rho12} ;\cr
  for the Extremal-t model: \var{alpha}, \var{sd_kappa}, \var{swe_kappa}, \var{rho12} and \var{nu}}
  \item{cond_B}{a vector of real numbers as the lower barriers of the conditioned variable (e.g. the return levels)}
  \item{cond_B_2}{a vector of real numbers as the upper barriers of the conditioned variable. \cr default is \code{cond_B_2 = Inf}}
  \item{var}{a character string being either \code{"sd"} or \code{"swe"} whenever \code{same_var = TRUE} and \code{model = "ext-t"} or \code{"ext-gauss"}. that is, if the conditioned variable is the same as the unconditioned variable, for the Extremal-Gaussian and Extremal-t model it has to be known which variable it is}
  \item{model}{a character string; chose which bivariate max-stable model should be used to calculate the conditional return levels. this should be either \code{"hr"} for the Huesler-Reiss, \code{"ext-gauss"} for the Extremal-Gaussian or \code{"ext-t"} (default) for the Extremal-t model}
  \item{printObjectives}{logical value; if \code{TRUE}, a summary of the values \deqn{f(B) = \Bigl|\,1/q - {\rm{Pr}}\bigl[Z_1(x_1) > B\, \big|\, Z_2(x_2) \in ({\rm{cond\_B}}\,,\,{\rm{cond\_B\_2}})\bigr]\,\Bigr|}{f(B) = abs(1/q - Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)])} is printed, where \eqn{B} is the found conditional return level and \deqn{{\rm{Pr}}\bigl[Z_1(x_1) > B\, \big|\, Z_2(x_2) \in ({\rm{cond\_B}}\,,\,{\rm{cond\_B\_2}})\bigr]}{Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)]} is the conditional probability that \eqn{Z_1(x_1)}{Z1(x1)} (e.g. sd or swe) exceeds \eqn{B}, given that \eqn{Z_2(x_2)}{Z2(x2)} (e.g. swe or sd) is in between \eqn{{\rm{cond\_B}}}{cond_B} and \eqn{{\rm{cond\_B\_2}}}{cond_B_2}; \cr by the minimaization of the function \eqn{f} we find the \eqn{1/q} quantile of this conditional distribution, thus we want small values. \cr default is \code{FALSE}
  }
}
\value{
a vector with the conditional return levels
}
\references{
Genton, M.G. & Padoan, S.A. & Sang, H. (2015): Multivariate max-stable spatial processes. Biometrika 102(1): 215-230.

\url{http://repository.kaust.edu.sa/kaust/bitstream/10754/552385/1/2013.GPS.Biometrika.Rev_14.pdf}
}
\seealso{
\code{\link{GEVparameters_from_models}},
\code{\link{cond_returnlevel_plot}},
\code{\link{cond_returnlevel_map}}, \cr
\code{\link{returnlevels}}
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

# load function output from optimizer_biv_hr_model
data("optim_hr")

# define coefficients
sd_coeff  = optim_hr$coefficients$sd_coeff
swe_coeff = optim_hr$coefficients$swe_coeff
cor_coeff = optim_hr$coefficients$cor_coeff

# define return period (same as in sd_GEVparam)
q = 100

# calculate return levels
sd_rl  = returnlevels(GEVparam = sd_GEVparam,  q = q)
swe_rl = returnlevels(GEVparam = swe_GEVparam, q = q)

# sd given swe at same locations
cond_rl_sd =
  cond_returnlevels(locations = covariables,
                    GEVparam = sd_GEVparam,
                    q = q, cond_locations = covariables,
                    cond_GEVparam = swe_GEVparam,
                    same_var = FALSE, cond_B = swe_rl,
                    cor_coeff = cor_coeff, model = "hr",
                    printObjectives = TRUE)

# swe given sd at same location
cond_rl_swe =
  cond_returnlevels(locations = covariables,
                    GEVparam = swe_GEVparam,
                    q = q, cond_locations = covariables,
                    cond_GEVparam = sd_GEVparam,
                    same_var = FALSE, cond_B = sd_rl,
                    cor_coeff = cor_coeff, model = "hr",
                    printObjectives = TRUE)

# sd given sd in ibk
ibk = c("lon" = 11.392778, "lat" = 47.267222, "alt" = 574,
        "mdday" = 6.12, "sd_mmax" = 16.8, "swe_mmax" = 29)
sd_GEVparam_ibk  = GEVparameters_from_models(ibk, sd_coeff)
sd_rl_ibk = returnlevels(GEVparam = sd_GEVparam_ibk, q = q)

cond_rl_ibk =
  cond_returnlevels(locations = covariables,
                    GEVparam = sd_GEVparam,
                    q = q, cond_locations = ibk,
                    cond_GEVparam = sd_GEVparam_ibk,
                    same_var = TRUE, cond_B = sd_rl_ibk,
                    cor_coeff = cor_coeff, model = "hr",
                    printObjectives = TRUE)
}
\keyword{functions}
\keyword{conditionals}
\keyword{returnlevels}
