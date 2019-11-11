\name{get_data}
\alias{get_data}
\docType{data}
\title{
Exemplary function output from \code{\link{get_data_from_Robj}}
}
\description{
this data file consists of the function output from \code{\link{get_data_from_Robj}} stemming from the data of \file{\link{sample_snow_data}}
}
\usage{data("get_data")}
\format{
  The format is: \cr

\code{List of 4} \cr
\code{
 $ snow_data   :List of 36} \cr
  \code{..$ 101:List of 2} \cr
  \code{.. ..$ meta:'data.frame':	1 obs. of  4 variables:} \cr
  \code{.. .. ..$ stnr: num 101 }\cr
  \code{.. .. ..$ lon : num 14.3} \cr
  \code{.. .. ..$ lat : num 46.6} \cr
  \code{.. .. ..$ alt : num 436 }\cr
  \code{.. ..$ data:'data.frame':	63 obs. of  5 variables: }\cr
  \code{.. .. ..$ year    : Factor w/ 63 levels "1950","1951", ... }\cr
  \code{.. .. ..$ sd      : num [1:63] 34 108 29 50 43 19 27 49 16 23 ...} \cr
  \code{.. .. ..$ swe     : num [1:63] 80.8 154.9 60.7 78.9 66.7 ... }\cr
  \code{.. .. ..$ date.sd : Factor w/ 63 levels "1951-02-05","1952-02-15", ... }\cr
  \code{.. .. ..$ date.swe: Factor w/ 63 levels "1951-02-05","1952-02-15", ...} \cr
  \code{.. .. ..- attr(*, "na.action")=Class 'omit'  Named int 64 }\cr
  \code{.. .. .. .. ..- attr(*, "names")= chr "64"} \cr
 \code{[... ...] }\cr
  \cr
 \code{$ sd_max_data : num [1:36, 1:68] NA NA NA 24 NA NA NA NA NA 40 ...} \cr
  \code{..- attr(*, "dimnames")=List of 2} \cr
  \code{.. ..$ : chr [1:36] "101" "102" "103" "104" ... }\cr
 \code{ .. ..$ : chr [1:68] "1948" "1949" "1950" "1951" ... }\cr
  \cr
 \code{$ swe_max_data: num [1:36, 1:68] NA NA NA 33 NA ...} \cr
  \code{..- attr(*, "dimnames")=List of 2 }\cr
  \code{.. ..$ : chr [1:36] "101" "102" "103" "104" ...} \cr
  \code{.. ..$ : chr [1:68] "1948" "1949" "1950" "1951" ... }\cr
  \cr
 \code{$ covariables : num [1:36, 1:6] 14.3 14.9 14.8 13.9 14.4 ... }\cr
  \code{..- attr(*, "dimnames")=List of 2} \cr
  \code{.. ..$ : chr [1:36] "101" "102" "103" "104" ...} \cr
  \code{.. ..$ : chr [1:6] "lon" "lat" "alt" "mdday" ...} \cr
  }
  \seealso{
\file{\link{sd_m_select}}, \file{\link{swe_m_select}}
}
\examples{
get_data = get(data("get_data"))
}
\keyword{datasets}
