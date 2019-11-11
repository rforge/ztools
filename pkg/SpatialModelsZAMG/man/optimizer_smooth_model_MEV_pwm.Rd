\name{optimizer_smooth_model_MEV_pwm}
\alias{optimizer_smooth_model_MEV_pwm}
\title{Optimize the coefficients of the best fitted linear MEV models via smooth modeling and probability weighted moments}
\description{this function optimizes the coefficients of the best fitted linear MEV models \cr (from the function \code{\link{model_selection_MEV}}) via probability weighted moments optimization}
\usage{
optimizer_smooth_model_MEV_pwm(m_select, data, method = c("nlminb","BFGS",
                       "ucminf","Nelder-Mead"),
                       follow.on = FALSE, itnmax = NULL,
                       printParam = FALSE)} 
\arguments{
  \item{m_select}{this input should be a list including \var{max_data}, \var{covariables} and \var{models} as in the output of the function \code{\link{model_selection}}}
  \item{data}{list whose elements are vectors including all observed daily values at one station, the stations have to be the same and used in the same order as the stations used for model selection.}
  \item{method}{optimization method(s) for external function \code{\link[optimx]{optimx}}, this can also be a vector. possible methods are: \var{Nelder-Mead}, \var{BFGS}, \var{CG}, \var{L-BFGS-B}, \var{nlm}, \var{nlminb}, \var{spg}, \var{ucminf}, \var{newuoa}, \var{bobyqa}, \var{nmkb}, \var{hjkb}, \var{Rcgmin}, \var{Rvmmin}\cr
  default is \code{method = c("nlminb","BFGS","ucminf","Nelder-Mead")}}
  \item{follow.on}{logical value; if \code{TRUE}, and there are multiple methods, then the last set of coefficients from one method is used as the starting set for the next \cr default is \code{FALSE}}
  \item{itnmax}{if provided as a vector of the same length as the length of  \code{method}, this gives the maximum number of iterations or function values for the corresponding method. if a single number is provided, this will be used for all methods}
  \item{printParam}{logical value; if \code{TRUE}, the GEV parameters during the optimization are printed. this might be useful to check the proper functioning of the optimization (shape parameter should be approximately between -0.5 and 0.5) \cr default is \code{FALSE}}
  
}
\value{a list with
  \item{summary}{a summary of the optimization results, including an information message whether the optimization was successful or not and which method delivered the best coefficients}
  \item{coefficients}{a list with the optimized coefficients. \cr
  containing: \cr
  \var{scalecoeff}, \var{shapecoeff} and \var{all_coeff}}
}
\references{
Blanchet, J. & Lehning, M. (2010): Mapping snow depth return levels: smooth spatial modeling versus station interpolation. Hydrology and Earth System Sciences 14(12): 2527-2544.
\cr\cr
Schellander, H., Lieb, A. and Hell, T. (2019) 'Error Structure of Metastatistical and Generalized Extreme Value Distributions for Modeling Extreme Rainfall in Austria', Earth and Space Science, 6, pp. 1616-1632. doi: 10.1029/2019ea000557.

}
\seealso{
\code{\link{model_selection_MEV}}, \code{\link{optimizer_smooth_model_MEV_pwm}}
}

\keyword{functions}
\keyword{optimizers}
