\name{returnlevel_plot}
\alias{returnlevel_plot}
\title{Create return level plots}
\description{this function creates a return level plot for a given location}
\usage{
returnlevel_plot(GEVparam, period_range = c(1,128),
                 obs = NULL, plottitle = NULL, 
                 save_name = NULL, save_dir = getwd(), 
                 printPlot = TRUE)}
\arguments{
  \item{GEVparam}{a named vector with the GEV parameters of the location for which a return level plot is wanted. names are \var{loc}, \var{scale} and \var{shape}}
  \item{period_range}{the range of the return period to be plotted. a vector with start and end point, which should be numbers greater or equal than 1. \cr default is \code{period_range = c(1,128)}}
  \item{obs}{a vector with empirical observations. if provided, sample quantiles of this vector are added to the plot as points}
  \item{plottitle}{a character string defining the title of the plot. \cr default is \code{plottitle = "return level plot"}. if no plottitle is wanted, use \code{plottitle = ""}}
  \item{save_name}{a character string defining the saving name of the plot. }
  \item{save_dir}{a character string defining the directory for the plot to be saved. \cr default is the working directory}
  \item{printPlot}{logical value; if \code{TRUE} (default), the plot is printed}
}
\value{a return level plot}
\seealso{
\code{\link{GEVparameters_from_models}}, \code{\link{returnlevel_map}}, \code{\link{cond_returnlevel_plot}}
}
\examples{
# load function output from optimizer_biv_hr_model
data("optim_hr")

# define coefficients
sd_coeff  = optim_hr$coefficients$sd_coeff
swe_coeff = optim_hr$coefficients$swe_coeff
cor_coeff = optim_hr$coefficients$cor_coeff

# arlberg
arl = c("lon" = 10.211, "lat" = 47.13, "alt" = 1979,
        "mdday" = 15.781, "sd_mmax" = 245.6,
        "swe_mmax" = 506.75)
sd_GEVparam_arl  = GEVparameters_from_models(arl, sd_coeff)
swe_GEVparam_arl = GEVparameters_from_models(arl, swe_coeff)

# create return level plots
returnlevel_plot(GEVparam = sd_GEVparam_arl,
                 save_name = "sd_return_level_plot_arlberg",
                 printPlot = FALSE)
returnlevel_plot(GEVparam = swe_GEVparam_arl,
                 save_name = "swe_return_level_plot_arlberg",
                 printPlot = FALSE)
}
\keyword{functions}
\keyword{plots}
