\name{sd_GEVparam}
\alias{sd_GEVparam}
\docType{data}
\title{
Exemplary function output from \code{\link{GEVparameters_from_models}} for sd
}
\description{
this data file consists of the function output from \code{\link{GEVparameters_from_models}} stemming from the data of \file{\link{sample_snow_data}} for the variable snow depth
}
\usage{data("sd_GEVparam")}
\format{
  The format is:\cr

\code{ num [1:92166, 1:3] 17.6 14.2 30.7 38.7 71.4 ...}\cr
\code{ - attr(*, "dimnames")=List of 2}\cr
\code{  ..$ : NULL}\cr
\code{  ..$ : chr [1:3] "loc" "scale" "shape"}\cr
}
\seealso{
\file{\link{swe_GEVparam}}, \file{\link{optim_hr}}
}
\examples{
sd_GEVparam = get(data("sd_GEVparam"))
}
\keyword{datasets}

