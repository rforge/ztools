\name{1 SpatialModelsZAMG}
\alias{SpatialModelsZAMG}
\docType{package}

\title{Package \pkg{SpatialModelsZAMG}}

\description{The package \pkg{SpatialModelsZAMG} was developed for spatial modeling of extreme snow depth (sd) and snow water equivalent (swe) in Austria, but can in principle be used for any kind of extreme events.

There is a function for model selection, and several functions for optimization of the model coefficients, like via smooth modeling, bivariate Huesler-Reiss, bivariate Extremal-Gaussian or bivariate Extremal-t model.

Furthermore, return levels can be calculated and maps can be plotted.

For the bivariate models, there is also the possibility to calculate and plot conditional return levels.}

\details{The package provides the following functions:
\enumerate{
    \item \code{\link{get_data_from_Robj}}: extract the snow data from a saved \code{robj} file
    \item \code{\link{model_selection}}: find best linear models according to AIC
    \item \code{\link{optimizer_smooth_model}}: optimize the coefficients of the best fitted linear models via smooth modeling
    \item \code{\link{optimizer_biv_hr_model}}: optimize the coefficients of the best fitted linear models via bivariate Huesler-Reiss model
    \item \code{\link{optimizer_biv_ext_gauss_model}}: optimize the coefficients of the best fitted linear models via bivariate Extremal-Gaussian model
    \item \code{\link{optimizer_biv_ext_t_model}}: optimize the coefficients of the best fitted linear models via bivariate Extremal-t model
    \item \code{\link{GEVparameters_from_models}}: calculate GEV parameters according to given linear models
    \item \code{\link{GEVparameters_map}}: create GEV parameter maps for Austria
    \item \code{\link{returnlevels}}: calculate return levels
    \item \code{\link{returnlevel_plot}}: create return level plots
    \item \code{\link{returnlevel_map}}: create return level maps for Austria
    \item \code{\link{cond_returnlevels}}: calculate conditional return levels
    \item \code{\link{cond_returnlevel_plot}}: create conditional return level plots
    \item \code{\link{cond_returnlevel_map}}: create conditional return level maps for Austria
    \item \code{\link{model_difference_map}}: create model difference maps for Austria \cr
    \cr
  }

  The package also includes some data files, which are used in the example codes:
  \enumerate{
  \item \file{\link{lon.at}}: longitude data from Austria
  \item \file{\link{lat.at}}: latitude data from Austria
  \item \file{\link{alt.at}}: altitude data from Austria
  \item \file{\link{border.at}}: border data from Austria
  \item \file{\link{sample_grid_data}}: sample grid data for sd and swe 
  \item \file{\link{sample_snow_data}}: sample snow data from 36 examplary stations distributed over Austria
  \item \file{\link{get_data}}: exemplary function output from \code{\link{get_data_from_Robj}}
  \item \file{\link{sd_m_select}}: exemplary function output from \code{\link{model_selection}} for sd
  \item \file{\link{swe_m_select}}: exemplary function output from \code{\link{model_selection}} for swe
  \item \file{\link{optim_hr}}: exemplary function output from \code{\link{optimizer_biv_hr_model}}
  \item \file{\link{sd_GEVparam}}: exemplary function output from \code{\link{GEVparameters_from_models}} for sd
  \item \file{\link{swe_GEVparam}}: exemplary function output from \code{\link{GEVparameters_from_models}} for swe\cr
  \cr
  
  }
}
\author{
Simon Gstoehl, \email{simon.gstoehl@vol.at} \cr
\cr

\bold{Contact:}

\emph{University of Innsbruck:} Tobias Hell, \email{tobias.hell@uibk.ac.at}

\emph{ZAMG Innsbruck:} Harald Schellander, \email{harald.schellander@zamg.ac.at}\cr
\cr

}
\keyword{packageinfo}