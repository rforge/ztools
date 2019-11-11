\name{model_selection}
\alias{model_selection}
\title{Find best linear models according to AIC}
\description{this function calculates the pointwise GEV parameter estimates for each station and finds the best linear models according to Akaike Information Criterion (AIC)}
\usage{
model_selection(max_data, covariables, at_least_data = 30,
                plot_station_distr = FALSE)}
\arguments{
  \item{max_data}{a matrix with the yearly maxima of snow depth (sd) or snow water equivalent (swe). each row corresponds to one station, columns are the corresponding years. matrix might contain \code{NA}'s}
  \item{covariables}{a matrix with the covariables for each station. each row corresponds to one station, columns should include at least \var{lon} (longitude), \var{lat} (latitude) and \var{alt} (altitude)}
  \item{at_least_data}{how many measurements does each station has to have at least for the model selection. \cr default is \code{at_least_data = 30}}
  \item{plot_station_distr}{logical value; if \code{TRUE}, a plot with the distribution of the used stations for the model selection is generated. \cr default is \code{FALSE}}
}
\value{a list with
  \item{max_data}{the given \code{max_data} matrix}
  \item{covariables}{the given \code{covariables} matrix}
  \item{point_est}{a matrix with the pointwise GEV parameter estimates. each row corresponds to one station, columns are \var{loc} (location parameter), \var{scale} (scale parameter) and \var{shape} (shape parameter)}
  \item{models}{a list of \code{lm}-class objects with the best fitted linear models for the GEV parameters:\cr
  loc_model:   \eqn{{\rm{loc}}   \sim \ldots}{loc   ~ \ldots} , \cr
  scale_model: \eqn{{\rm{scale}} \sim \ldots}{scale ~ \ldots} , \cr
  shape_model: \eqn{{\rm{shape}} \sim \ldots}{shape ~ \ldots}}
  \item{used_for_model_selection}{the indices of the stations which were used for the model selection. you can use it for example as \code{max_data[used_for_model_selection,]} to get the maxima data of all stations used for the model selection}
}
\seealso{
\code{\link{get_data_from_Robj}}, \code{\link{optimizer_smooth_model}}, \code{\link{optimizer_biv_hr_model}}
}
\examples{
# load function output from get_data_from_Robj
get_data = get(data("get_data"))

# define covariables for sd and swe in order
# to perform model selection (drop swe_mmax
# for sd and sd_mmax for swe model selection)

sd_covariables  = get_data$covariables[,-6]
swe_covariables = get_data$covariables[,-5]

# perform model selection
sd_m_select =
  model_selection(max_data = get_data$sd_max_data,
                  covariables = sd_covariables,
                  at_least_data = 10)

swe_m_select =
  model_selection(max_data = get_data$swe_max_data,
                  covariables = swe_covariables,
                  at_least_data = 10)
}
\keyword{functions}
\keyword{modelselection}
