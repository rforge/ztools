\name{optim_hr}
\alias{optim_hr}
\docType{data}
\title{
Exemplary function output from \code{\link{optimizer_biv_hr_model}}
}
\description{
this data file consists of the function output from \code{\link{optimizer_biv_hr_model}} stemming from the data of \file{\link{sample_snow_data}}
}
\usage{data("optim_hr")}
\format{
The format is: \cr

\code{List of 4} \cr
\code{ $ complete_sd_max_data : num [1:27, 1:41] 33 69 50 9 40 117 342 175 18 22 ...}\cr
\code{  ..- attr(*, "dimnames")=List of 2}\cr
\code{  .. ..$ : chr [1:27] "101" "102" "103" "104" ...}\cr
\code{  .. ..$ : chr [1:41] "1972" "1973" "1974" "1975" ...}\cr
\cr
\code{ $ complete_swe_max_data: num [1:27, 1:41] 48.2 122.5 94.5 11.2 47.6 ...}\cr
\code{  ..- attr(*, "dimnames")=List of 2}\cr
\code{  .. ..$ : chr [1:27] "101" "102" "103" "104" ...}\cr
\code{  .. ..$ : chr [1:41] "1972" "1973" "1974" "1975" ...}\cr
\cr
\code{ $ summary      :List of 3}\cr
\code{  ..$ message: chr} \cr
\code{      "optimization was successful, coefficients were chosen from method 'ucminf'"}\cr
\code{  ..$ method : chr "ucminf"}\cr
\code{  ..$ result :Classes 'optimx' and 'data.frame':	1 obs. of  33 variables ...}\cr
\cr
\code{ $ coefficients :List of 4}\cr
\code{  ..$ sd_coeff :List of 3}\cr
\code{  .. ..$ loccoeff  : Named num [1:3] 3.328 0.522 0.157}\cr
\code{  .. .. ..- attr(*, "names")= chr [1:3] "(Intercept)" "sd_mmax" "mdday"}\cr
\code{  .. ..$ scalecoeff: Named num [1:3] 4.232 0.391 0.111}\cr
\code{  .. .. ..- attr(*, "names")= chr [1:3] "(Intercept)" "sd_mmax" "mdday"}\cr
\code{  .. ..$ shapecoeff: Named num [1:5] 2.0718 -0.0588 -0.6298 0.0215 0.7808}\cr
\code{  .. .. ..- attr(*, "names")= chr [1:5] "(Intercept)" "lat" "loc" "sd_mmax" ...}\cr
\code{  ..$ swe_coeff:List of 3}\cr
\code{  .. ..$ loccoeff  : Named num [1:3] 8.5681 0.4757 -0.0205}\cr
\code{  .. .. ..- attr(*, "names")= chr [1:3] "(Intercept)" "swe_mmax" "mdday"}\cr
\code{  .. ..$ scalecoeff: Named num [1:3] 9.9517 0.0022 0.388}\cr
\code{  .. .. ..- attr(*, "names")= chr [1:3] "(Intercept)" "alt" "swe_mmax"}\cr
\code{  .. ..$ shapecoeff: Named num [1:5] 3.2771 -0.0657 -0.0739 0.3381 -0.2246}\cr
\code{  .. .. ..- attr(*, "names")= chr [1:5] "(Intercept)" "lat" "swe_mmax" "loc" ...}\cr
\code{  ..$ cor_coeff: Named num [1:3] 0.000126 0.152239 1.477005}\cr
\code{  .. ..- attr(*, "names")= chr [1:3] "alpha" "kappa" "lambda12"}\cr
\code{  ..$ all_coeff: Named num [1:25] 0.000126 0.152239 1.477005 3.3279 0.521775 ...}\cr
\code{  .. ..- attr(*, "names")= chr [1:25] "alpha" "kappa" "lambda12" "loc1_Int" ...}\cr
}
\seealso{
\file{\link{sd_m_select}}, \file{\link{swe_m_select}}
}
\examples{
data("optim_hr")

# define coefficients
sd_coeff  = optim_hr$coefficients$sd_coeff
swe_coeff = optim_hr$coefficients$swe_coeff
cor_coeff = optim_hr$coefficients$cor_coeff
}
\keyword{datasets}

