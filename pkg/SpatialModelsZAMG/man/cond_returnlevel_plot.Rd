\name{cond_returnlevel_plot}
\alias{cond_returnlevel_plot}
\title{Create conditional return level plots}
\description{this function creates a conditional return level plot. \cr

the \eqn{q}-year conditional return level of variable \eqn{Z_1}{Z1} at location \eqn{x_1}{x1} given variable \eqn{Z_2}{Z2} at location \eqn{x_2}{x2}, is defined as the threshold \eqn{B}, such that the conditional probability that \eqn{Z_1(x_1)}{Z1(x1)} exceeds this threshold, given that \eqn{Z_2(x_2)}{Z2(x2)} is in the interval \eqn{({\rm{cond\_B}}\,,\,{\rm{cond\_B\_2}})}{(cond_B,cond_B_2)}, is \eqn{1/q} :
\deqn{{\rm{Pr}}\bigl[Z_1(x_1) > B\, \big|\, Z_2(x_2) \in ({\rm{cond\_B}}\,,\,{\rm{cond\_B\_2}})\bigr] = \frac{1}{q} \, .}{Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)] = 1/q .}
}
\usage{
cond_returnlevel_plot(location, GEVparam, cond_location,
                      cond_GEVparam, same_var, cor_coeff,
                      cond_B, cond_B_2 = Inf, obs = NULL, 
                      var = NULL, period_range = c(1,128), 
                      model = "ext-t", 
                      printObjectives = FALSE, 
                      plottitle = NULL, save_name = NULL, 
                      save_dir = getwd(), printPlot = TRUE)
}
\arguments{
  \item{location}{
  a vector with certain location characteristics as entries. characteristics are for example: longitude, latitude and altitude}
  \item{GEVparam}{a named vector with the GEV parameters of the variable for which conditional return levels are calculated. names are \var{loc}, \var{scale} and \var{shape}}
  \item{cond_location}{a vector with the characteristics of the conditioned location. characteristics should be the same as in \code{location}}
  \item{cond_GEVparam}{a named vector with the GEV parameters of the conditioned variable. names are \var{loc}, \var{scale} and \var{shape}}
  \item{same_var}{logical value; if \code{TRUE}, the conditioned variable is the same as the unconditioned variable. this has to be known in order to use the right correlation function}
    \item{cor_coeff}{a named vector with the correlation parameters. \cr
  for the Huesler-Reiss model: \var{alpha}, \var{kappa} and \var{lambda12} ; \cr
  for the Extremal-Gaussian model: \var{alpha}, \var{sd_kappa}, \var{swe_kappa} and \var{rho12} ;\cr
  for the Extremal-t model: \var{alpha}, \var{sd_kappa}, \var{swe_kappa}, \var{rho12} and \var{nu}}
    \item{cond_B}{a real number as the lower barrier of the conditioned variable (e.g. the return level)}
    \item{cond_B_2}{a real number as the upper barrier of the conditioned variable. \cr default is \code{cond_B_2 = Inf}}
    \item{obs}{a vector with empirical observations. if provided, sample quantiles of this vector are added to the plot as points}
    \item{var}{a character string being either \code{"sd"} or \code{"swe"} whenever \code{same_var = TRUE} and \code{model = "ext-t"} or \code{"ext-gauss"}. that is, if the conditioned variable is the same as the unconditioned variable, for the Extremal-Gaussian and Extremal-t model it has to be known which variable it is}
  \item{period_range}{the range of the return period to be plotted. a vector with start and end point, which should be numbers greater or equal than 1. \cr default is \code{period_range = c(1,128)}}
  \item{model}{a character string; chose which bivariate max-stable model should be used to calculate the conditional return levels. this should be either \code{"hr"} for the Huesler-Reiss, \code{"ext-gauss"} for the Extremal-Gaussian or \code{"ext-t"} (default) for the Extremal-t model}
  \item{printObjectives}{logical value; if \code{TRUE}, a summary of the values \deqn{f(B) = \Bigl|\,1/q - {\rm{Pr}}\bigl[Z_1(x_1) > B\, \big|\, Z_2(x_2) \in ({\rm{cond\_B}}\,,\,{\rm{cond\_B\_2}})\bigr]\,\Bigr|}{f(B) = abs(1/q - Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)])} is printed, where \eqn{B} is the found conditional return level and \deqn{{\rm{Pr}}\bigl[Z_1(x_1) > B\, \big|\, Z_2(x_2) \in ({\rm{cond\_B}}\,,\,{\rm{cond\_B\_2}})\bigr]}{Pr[Z1(x1) > B | Z2(x2) in (cond_B,cond_B_2)]} is the conditional probability that \eqn{Z_1(x_1)}{Z1(x1)} (e.g. sd or swe) exceeds \eqn{B}, given that \eqn{Z_2(x_2)}{Z2(x2)} (e.g. swe or sd) is in between \eqn{{\rm{cond\_B}}}{cond_B} and \eqn{{\rm{cond\_B\_2}}}{cond_B_2}; \cr by the minimaization of the function \eqn{f} we find the \eqn{1/q} quantile of this conditional distribution, thus we want small values. \cr default is \code{FALSE}
  }
    \item{plottitle}{a character string defining the title of the plot. \cr default is \code{plottitle = "conditional return level plot"}. if no plottitle is wanted, use \code{plottitle = ""}}
  \item{save_name}{a character string defining the saving name of the plot.}
  \item{save_dir}{a character string defining the directory for the plot to be saved. \cr default is the working directory}

  \item{printPlot}{logical value; if \code{TRUE} (default), the plot is printed}
}
\value{
a conditional return level plot
}
\references{
Genton, M.G. & Padoan, S.A. & Sang, H. (2015): Multivariate max-stable spatial processes. Biometrika 102(1): 215-230.

\url{http://repository.kaust.edu.sa/kaust/bitstream/10754/552385/1/2013.GPS.Biometrika.Rev_14.pdf}
}
\seealso{
\code{\link{GEVparameters_from_models}}, \code{\link{cond_returnlevel_map}}, \code{\link{returnlevel_plot}}
}
\examples{
# load function output from optimizer_biv_hr_model
data("optim_hr")

# define coefficients
sd_coeff  = optim_hr$coefficients$sd_coeff
swe_coeff = optim_hr$coefficients$swe_coeff
cor_coeff = optim_hr$coefficients$cor_coeff

# conditional sd return level at arlberg
# given swe > swe_rl50 in ibk

# arlberg
arl = c("lon" = 10.211, "lat" = 47.13, "alt" = 1979,
        "mdday" = 15.781, "sd_mmax" = 245.6,
        "swe_mmax" = 506.75)
sd_GEVparam_arl  = GEVparameters_from_models(arl, sd_coeff)
swe_GEVparam_arl = GEVparameters_from_models(arl, swe_coeff)

# innsbruck
ibk = c("lon" = 11.392778, "lat" = 47.267222, "alt" = 574,
        "mdday" = 6.12, "sd_mmax" = 16.8, "swe_mmax" = 29)
sd_GEVparam_ibk  = GEVparameters_from_models(ibk, sd_coeff)
swe_GEVparam_ibk = GEVparameters_from_models(ibk, swe_coeff)

# create plot
cond_B = returnlevels(GEVparam = swe_GEVparam_ibk, q = 50)
cond_returnlevel_plot(location = arl, GEVparam = sd_GEVparam_arl,
                      cond_location = ibk,
                      cond_GEVparam = swe_GEVparam_ibk,
                      same_var = FALSE, cond_B = cond_B,
                      cor_coeff = cor_coeff, model = "hr",
                      save_name = "sd_cond_return_level_plot_arl_given_ibk",
                      plottitle = "sd cond rl at arlberg given swe > swe_rl50 in ibk",
                      printPlot = FALSE)
}
\keyword{functions}
\keyword{conditionals}
\keyword{plots}
