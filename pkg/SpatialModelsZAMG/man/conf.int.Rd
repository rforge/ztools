\name{conf.int}
\alias{conf.int}
\title{Confidence intervals for GEV and MEV distributions}
\description{Computes confidence intervals for GEV and MEV distributions using a parametric bootstrap method.}
\usage{conf.int(x, alpha = 0.05, return.periods = c(2,10,30,50,100),
         dist = c("GEV","MEV"), method = "boot", R = 502, verbose = FALSE)}
\arguments{
  \item{x}{data.frame with either columns "W", "C", "n" for the MEV distribution, or "loc","scale","shape" for the GEV distribution (see details).}
  \item{alpha}{Number between zero and one giving the 1 - alpha confidence level.}
  \item{return.periods}{Numeric vector giving the desired return periods.}
  \item{dist}{Character string naming which EVD to calculate confidence intervals from. }
  \item{method}{character naming which method for obtaining confidence intervals should be used. Option "boot" employs a parametric bootstrap that simulates data from the fitted model, and then fits the EVD to each simulated data set to obtain a sample of parameters or return levels. Currently only "boot" is available.}
  \item{R}{the number of bootstrap iterations.}
  \item{verbose}{Should additional information be printed during rundtime? Can be \code{TRUE} or \code{FALSE}, which is the default. }

}


\details{
\code{x} must be a data.frame, whose content depends on the distribution \code{dist}:
\cr\cr
\code{MEV}: If \code{dist='MEV'}, then the data.frame \code{x} must contain columns "W", "C", "n". They represent shape and scale parameter of the MEV, and mean number of wet days.
\cr\cr
\code{GEV}: If \code{dist='GEV'}, columns "loc","scale","shape" must be location, scale and shape parameter of the GEV .
}

\value{Either a numeric vector of length 3 (if only one return level is used) or a matrix.}

\author{Harald Schellander}

\references{
Schellander, H., Lieb, A. and Hell, T. (2019) 'Error Structure of Metastatistical and Generalized Extreme Value Distributions for Modeling Extreme Rainfall in Austria', Earth and Space Science, 6, pp. 1616-1632. doi: 10.1029/2019ea000557.
}

\keyword{extreme value distribution, confidence interval}
