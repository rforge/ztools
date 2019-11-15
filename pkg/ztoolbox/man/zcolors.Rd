\name{zcolors}
\alias{zcolors}
\alias{show_palettes}

\title{Convenient selection of colors based on HCL color space}

\description{
 Provides a wrapper for simple selection of colors for typical meteorological parameters.
}
\usage{
zcolors(n, parameter = c("std","rain","snow","temperature","terrain"), 
  type=c("qualitative","sequential","diverge"), plot = FALSE, ...)
  
show_palettes(n = 10)  
 

}

\arguments{
     \item{n}{ number of colors.}
     \item{parameter}{ one of \code{std}, \code{rain}, \code{snow}, \code{temperature}, \code{terrain}. The default is \code{std}.}
     \item{type}{ type of color palette. One of \code{sequential}, \code{diverge}, \code{qualitative}. The default is \code{qualitative} (see details).}
     \item{plot}{ Shall an image be plotted? Can be \code{TRUE} or \code{FALSE}, which is the default.}
     \item{\dots}{ Other parameters that are transferred to the HCL routines (see details).}
}

\details{
Colors are provided based on the HCL color space, which is used as reference. See package \code{\link{colorspace}}.
Function \code{show_palettes} shows an overview of currently available colorsets.
}

\value{
  A character string with HEX representation of the \code{n} colors. 
  If \code{plot = TRUE} an image with the chosen colorset is shown. Function \code{show_palettes} only produces an overview image. 

}

\keyword{zamg, color, colorspace, hcl}

\examples{

zcolors(8,plot=TRUE)
show_palettes(12)

}
