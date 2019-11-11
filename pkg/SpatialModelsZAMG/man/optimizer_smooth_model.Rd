\name{optimizer_smooth_model}
\alias{optimizer_smooth_model}
\title{Optimize the coefficients of the best fitted linear models via smooth modeling}
\description{this function optimizes the coefficients of the best fitted linear models \cr (from the function \code{\link{model_selection}}) via smooth modeling \cr and with maximum likelihood estimation}
\usage{
optimizer_smooth_model(m_select, method = c("nlminb","BFGS",
                       "ucminf","Nelder-Mead"),
                       follow.on = FALSE, itnmax = NULL,
                       printParam = FALSE)}
\arguments{
  \item{m_select}{this input should be a list including \var{max_data}, \var{covariables} and \var{models} as in the output of the function \code{\link{model_selection}}}
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
  \var{loccoeff}, \var{scalecoeff}, \var{shapecoeff} and \var{all_coeff}}
}
\references{
Blanchet, J. & Lehning, M. (2010): Mapping snow depth return levels: smooth spatial modeling versus station interpolation. Hydrology and Earth System Sciences 14(12): 2527-2544.

\url{https://www.hydrol-earth-syst-sci.net/14/2527/2010/hess-14-2527-2010.pdf}
}
\seealso{
\code{\link{model_selection}}, \code{\link{optimizer_biv_hr_model}}
}
\examples{
# load function output from model_selection
sd_m_select  = get(data("sd_m_select"))
swe_m_select = get(data("swe_m_select"))

# perform optimization
sd_optim =
  optimizer_smooth_model(m_select = sd_m_select,
                         method = c("nlminb","ucminf"))

swe_optim =
  optimizer_smooth_model(m_select = swe_m_select,
                         method = "nlminb",
                         itnmax = 500)
}
\keyword{functions}
\keyword{optimizers}
