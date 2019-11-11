\name{model_selection_MEV}
\alias{model_selection_MEV}
\title{Find best linear models for MEV distributed data according to AIC}
\description{this function calculates the pointwise MEV parameter estimates for each station and finds the best linear models according to Akaike Information Criterion (AIC)}
\usage{
model_selection_MEV(fitted_par, covariables, plot_station_distr = FALSE)
}

\arguments{
  \item{fitted_par}{ data.frame with MEV parameters w,c,n as columns. Each row corresponds to one station.}
  \item{covariables}{ data.frame with covariables for each station. Each row corresponds to one station, columns should include at least 'lon' (longitude), 'lat' (latitude) and 'alt' (altitude).}
  \item{plot_station_distr}{logical value; if \code{TRUE}, a plot with the distribution of the used stations for the model selection is generated. \cr default is \code{FALSE}}
}
\value{a list with
  \item{fitted_par}{the given \code{fitted_par} matrix}
  \item{covariables}{the given \code{covariables} matrix}
  \item{models}{a list of \code{lm}-class objects with the best fitted linear models for the MEV parameters:\cr
  scale_model: \eqn{{\rm{scale}} \sim \ldots}{scale ~ \ldots} , \cr
  shape_model: \eqn{{\rm{shape}} \sim \ldots}{shape ~ \ldots}}
}
\seealso{
\code{\link{optimizer_smooth_model_MEV_pwm}}, \code{\link{optimizer_smooth_model_MEV_mle}}
}

\keyword{functions}
\keyword{modelselection}
