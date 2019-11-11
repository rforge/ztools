\name{sample_grid_data}
\alias{sample_grid_data}
\docType{data}
\title{
Sample grid data for snow depth and snow water equivalent
}
\description{
this data file includes a \eqn{183690 \times 3}{183690 x 3} data frame where the three columns are \var{mdday} (mean of the time difference in days of sd and swe maxima), \var{mmsd} (mean maxima of sd) and \var{mmswe} (mean maxima of swe). the rows correspond to grid points covering Austria. the grid points are stemming from the same \eqn{585 \times 314}{585 x 314} grid as in \file{\link{lon.at}}, \file{\link{lat.at}} and \file{\link{alt.at}} when transforming them to a vector of length \eqn{585 \cdot 314 = 183690}{585*314 = 183690} , e.g.: \code{as.vector(t(lon.at))}.
}
\usage{data("sample_grid_data")}
\format{
  A data frame with 183690 observations of the following 3 variables.
  \describe{
    \item{\code{mmsd}}{a numeric vector: mean maxima of sd}
    \item{\code{mmswe}}{a numeric vector: mean maxima of swe}
    \item{\code{mdday}}{a numeric vector: mean of the time difference in days of sd and swe maxima}
  }
}
\seealso{
\file{\link{lon.at}}, \file{\link{lat.at}}, \file{\link{alt.at}}, \file{\link{sample_snow_data}}
}
\examples{
x = get(data("sample_grid_data"))
}
\keyword{datasets}

