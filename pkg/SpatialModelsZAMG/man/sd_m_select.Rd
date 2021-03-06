\name{sd_m_select}
\alias{sd_m_select}
\docType{data}
\title{
Exemplary function output from \code{\link{model_selection}} for sd
}
\description{
this data file consists of the function output from \code{\link{model_selection}} stemming from the data of \file{\link{sample_snow_data}} for the variable snow depth
}
\usage{data("sd_m_select")}
\format{
The format is:\cr

\code{List of 5} \cr
\code{ $ max_data     : num [1:36, 1:68] NA NA NA 24 NA NA NA NA NA 40 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : chr [1:36] "101" "102" "103" "104" ...} \cr
\code{  .. ..$ : chr [1:68] "1948" "1949" "1950" "1951" ...} \cr
\cr
\code{ $ covariables  : num [1:36, 1:5] 14.3 14.9 14.8 13.9 14.4 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : chr [1:36] "101" "102" "103" "104" ...}\cr
\code{  .. ..$ : chr [1:5] "lon" "lat" "alt" "mdday" ...} \cr
\cr
\code{ $ point_est    : num [1:36, 1:3] 27.9 40.1 21.9 19.3 15.8 ...}\cr
\code{  ..- attr(*, "dimnames")=List of 2}\cr
\code{  .. ..$ : chr [1:36] "101" "102" "103" "104" ...}\cr
\code{  .. ..$ : chr [1:3] "loc" "scale" "shape"}\cr
\cr
\code{ $ models       :List of 3}\cr
\code{  ..$ loc_model  :List of 30}\cr
\code{  ..$ scale_model:List of 30} \cr
\code{  ..$ shape_model:List of 30} \cr
\cr
\code{ $ used_for_model_selection: Named int [1:32] 1 2 3 4 5 6 7 8 9 10 ...} \cr
\code{  ..- attr(*, "names")= chr [1:32] "101" "102" "103" "104" ... }\cr
}
\seealso{
\file{\link{swe_m_select}}, \file{\link{optim_hr}}
}
\examples{
sd_m_select = get(data("sd_m_select"))
}
\keyword{datasets}
