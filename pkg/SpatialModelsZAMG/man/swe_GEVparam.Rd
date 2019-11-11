\name{swe_GEVparam}
\alias{swe_GEVparam}
\docType{data}
\title{
Exemplary function output from \code{\link{GEVparameters_from_models}} for swe
}
\description{
this data file consists of the function output from \code{\link{GEVparameters_from_models}} stemming from the data of \file{\link{sample_snow_data}} for the variable snow water equivalent
}
\usage{data("swe_GEVparam")}
\format{
  The format is: \cr

\code{ num [1:92166, 1:3] 22.7 17.6 58.4 65.8 114.6 ...}\cr
\code{ - attr(*, "dimnames")=List of 2}\cr
\code{  ..$ : NULL}\cr
\code{  ..$ : chr [1:3] "loc" "scale" "shape"}
}
\seealso{
\file{\link{sd_GEVparam}}, \file{\link{optim_hr}}
}
\examples{
swe_GEVparam = get(data("swe_GEVparam"))
}
\keyword{datasets}
