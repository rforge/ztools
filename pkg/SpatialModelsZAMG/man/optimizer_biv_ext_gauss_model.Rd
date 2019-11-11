\name{optimizer_biv_ext_gauss_model}
\alias{optimizer_biv_ext_gauss_model}
\title{Optimize the coefficients of the best fitted linear models via bivariate Extremal-Gaussian model}
\description{this function optimizes the coefficients of the best fitted linear models \cr (from the function \code{\link{model_selection}}) via bivariate Extremal-Gaussian model \cr (= Extremal-t model with fixed \eqn{{\rm{nu}} = 1}{nu = 1}) and with composite likelihood inference}
\usage{
optimizer_biv_ext_gauss_model(sd_m_select, swe_m_select,
                              method = "ucminf",
                              follow.on = FALSE,
                              itnmax = NULL,
                              printParam = FALSE)}
\arguments{
  \item{sd_m_select}{this input should be a list including \var{max_data}, \var{covariables} and \var{models} as in the output of the function \code{\link{model_selection}}}
  \item{swe_m_select}{this input should be a list including \var{max_data}, \var{covariables} and \var{models} as in the output of the function \code{\link{model_selection}}}
  \item{method}{optimization method(s) for external function \code{\link[optimx]{optimx}}, this can also be a vector, but optimization might take several hours for one method! possible methods are: \var{Nelder-Mead}, \var{BFGS}, \var{CG}, \var{L-BFGS-B}, \var{nlm}, \var{nlminb}, \var{spg}, \var{ucminf} (default), \var{newuoa}, \var{bobyqa}, \var{nmkb}, \var{hjkb}, \var{Rcgmin}, \var{Rvmmin}}
  \item{follow.on}{logical value; if \code{TRUE}, and there are multiple methods, then the last set of coefficients from one method is used as the starting set for the next \cr default is \code{FALSE}}
  \item{itnmax}{if provided as a vector of the same length as the length of \code{method}, this gives the maximum number of iterations or function values for the corresponding method. if a single number is provided, this will be used for all methods}
  \item{printParam}{logical value; if \code{TRUE}, the GEV parameters during the optimization are printed. this might be useful to check the proper functioning of the optimization (shape parameter should be approximately between -0.5 and 0.5) \cr default is \code{FALSE}}
}
\value{a list with
  \item{complete_sd_max_data}{a completed submatrix of the \var{max_data} matrix for sd, where the completation was performed via function \code{\link[mice]{complete}} of the package \code{\link[mice]{mice}}. each row corresponds to one station (unchanged), columns (years) were chosen such that less than 50\% of the original column-entries were \code{NA}'s}
  \item{complete_swe_max_data}{completed submatrix of the \var{max_data} matrix for swe, where the completation was performed via function \code{\link[mice]{complete}} of the package \code{\link[mice]{mice}}. each row corresponds to one station (unchanged), columns (years) were chosen such that less than 50\% of the original column-entries were \code{NA}'s}
  \item{summary}{a summary of the optimization results, including an information message whether the optimization was successful or not and which method delivered the best coefficients}
  \item{coefficients}{a list with the optimized coefficients. \cr
  containing: \cr
   \var{sd_coeff} (\var{sd_loccoeff}, \var{sd_scalecoeff}, \var{sd_shapecoeff}), \cr
    \var{swe_coeff} (\var{swe_loccoeff}, \var{swe_scalecoeff}, \var{swe_shapecoeff}), \cr
    \var{cor_coeff} (\var{alpha}, \var{sd_kappa}, \var{swe_kappa}, \var{rho}) and \cr
    \var{all_coeff}}
}
\references{
Genton, M.G. & Padoan, S.A. & Sang, H. (2015): Multivariate max-stable spatial processes. Biometrika 102(1): 215-230.

\url{http://repository.kaust.edu.sa/kaust/bitstream/10754/552385/1/2013.GPS.Biometrika.Rev_14.pdf}
}
\seealso{
\code{\link{model_selection}}, \code{\link{optimizer_biv_ext_t_model}}, \code{\link{optimizer_biv_hr_model}}
}

\keyword{functions}
\keyword{optimizers}
