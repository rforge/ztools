\name{MEVparameters_from_models}
\alias{MEVparameters_from_models}
\title{Calculate MEV parameters according to given linear models}
\description{this function calculates the MEV parameters according to the linear models \cr specified by \code{coefficients}}
\usage{
MEVparameters_from_models(covariables, coefficients)}
\arguments{
  \item{covariables}{a named matrix or vector with the covariables. each row corresponds to one location, columns are the corresponding covariables (names of the coefficients) from the linear models; the intercept is added seperately and thus doesn't have to be included}
  \item{coefficients}{a list with the named coefficients of the linear models including\cr
  \var{scalecoeff} (coefficients of scale parameter) -- must be named\cr
  \var{shapecoeff} (coefficients of shape parameter) -- must be named}
}
\value{
  \item{MEVparam }{a matrix or vector with the MEV parameters. each row corresponds to one location, columns are \var{scale} and \var{shape}}
}
\seealso{
\code{\link{GEVparameters_from_models}}, \code{\link{GEVparameters_map}}, \code{\link{returnlevels}}, \code{\link{cond_returnlevels}}, \code{\link{returnlevel_plot}}, \cr \code{\link{cond_returnlevel_plot}}
}

\keyword{functions}
\keyword{MEV,parameters}
