\name{border.at}
\alias{border.at}
\docType{data}
\title{
Border data from Austria
}
\description{
a \eqn{15134 \times 7}{15134 x 7} data frame with border data from Austria. this data file is used within the functions that create maps, like \code{\link{returnlevel_map}}.
}
\usage{data("border.at")}
\format{
  A data frame with 15134 observations on the following 7 variables.
  \describe{
    \item{\code{long}}{a numeric vector: longitude data}
    \item{\code{lat}}{a numeric vector: latitude data}
    \item{\code{order}}{a numeric vector}
    \item{\code{hole}}{a logical vector}
    \item{\code{piece}}{a factor with levels \code{1} and \code{2}}
    \item{\code{id}}{a character vector}
    \item{\code{group}}{a factor with levels \code{1.1} and \code{1.2}}
  }
}
\seealso{
\file{\link{lon.at}}, \file{\link{lat.at}}, \file{\link{alt.at}}
}
\examples{
border = get(data("border.at"))
}
\keyword{datasets}

