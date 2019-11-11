\name{returnlevels}
\alias{returnlevels}
\title{Calculate return levels}
\description{this function calculates the return levels according to the given GEV parameters \code{GEVparam} and a return period of \code{q} years}
\usage{
returnlevels(GEVparam, q)
}
\arguments{
  \item{GEVparam}{a named matrix or vector with the GEV parameters. each row corresponds to one location, columns are \var{loc}, \var{scale} and \var{shape}}
  \item{q}{the return period for the calculation of return levels -- must be a number greater than 1}
}
\value{
a vector with the return levels for each location
}
\seealso{
\code{\link{GEVparameters_from_models}},
\code{\link{returnlevel_plot}},
\code{\link{returnlevel_map}},
\code{\link{cond_returnlevels}}
}
\examples{
# load function output from GEVparameters_from_models
sd_GEVparam  = get(data("sd_GEVparam"))
swe_GEVparam = get(data("swe_GEVparam"))

# calculate returnlevels
q = 100
sd_rl  = returnlevels(GEVparam = sd_GEVparam,  q = q)
swe_rl = returnlevels(GEVparam = swe_GEVparam, q = q)
}
\keyword{functions}
\keyword{returnlevels}
