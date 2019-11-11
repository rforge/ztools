\name{sample_snow_data}
\alias{sample_snow_data}
\docType{data}
\title{
Sample snow data from 36 examplary stations distributed over Austria
}
\description{
this data file consists of a list with 36 examplary stations including randomized snow data, like the maximal snow depth and maximal snow water equivalent with the respective date of occurance as well as geographical data like longitude, latituide and altitude.
}
\details{
stations are ordered according to a station number \var{stnr}. in the sample station numbers are from 101 to 136. each station contains meta data (data frame with \var{stnr}, \var{lon}, \var{lat} and \var{alt}) and the snow data (data frame with columns \var{year}, \var{sd}, \var{swe}, \var{date.sd} and \var{date.swe})
}
\usage{data("sample_snow_data")}
\format{
  The format is: \cr

  \code{List of 36} \cr
  \code{$ 101:List of 2} \cr
  \code{..$ meta:'data.frame':	1 obs. of  4 variables:} \cr
  \code{.. ..$ stnr: num 101}\cr
  \code{.. ..$ lon : num 14.3}\cr
  \code{.. ..$ lat : num 46.6}\cr
  \code{.. ..$ alt : num 436}\cr
  \code{..$ data:'data.frame':	63 obs. of  5 variables:}\cr
  \code{.. ..$ year    : Factor w/ 63 levels "1950","1951", ... }\cr
  \code{.. ..$ sd      : num [1:63] 34 108 29 50 43 19 27 49 16 23 ...}\cr
  \code{.. ..$ swe     : num [1:63] 80.8 154.9 60.7 78.9 66.7 ...}\cr
  \code{.. ..$ date.sd : Factor w/ 63 levels "1951-02-05","1952-02-15", ...}\cr
  \code{.. ..$ date.swe: Factor w/ 63 levels "1951-02-05","1952-02-15", ...}\cr
  \code{.. ..- attr(*, "na.action")=Class 'omit'  Named int 64}\cr
  \code{.. .. .. ..- attr(*, "names")= chr "64"}\cr
 \code{[...]}
}
\seealso{
\file{\link{sample_grid_data}}
}
\examples{
snow_data = get(data("sample_snow_data"))
}
\keyword{datasets}


