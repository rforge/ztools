\name{nearest_gridpoint}
\alias{nearest_gridpoint}

\title{Calculates the nearest point on a grid to a reference point}

\description{
 Calculates nearest point on a grid to a reference point in longitude, latitude degree decimal.
}
\usage{
nearest_gridpoint(lon, lat, grid, measure = c("euclidian","haversine"))

}

\arguments{
     \item{lon, lat}{ longitude and latitude of the reference point in degrees decimal.}
     \item{grid}{ matrix or data.frame with at leats two columns named \code{lon} and \code{lat}. The function looks for the nearest grid point to the reference point.}
     \item{measure}{ method for distance calculation. Currently \code{euclidian}, which is the default and \code{haversine} are possible (see details).}
}

\details{
The distance of each point on \code{grid} is calculated either as a euclidian distance, or with the Haversine Formula, which computes the great-circle distance between two points on a sphere.

}

\value{
  An list with components:
  \item{idx}{ index of the point nearest to coordinates \code{lon, lat} on \code{grid}.}
  \item{lon}{ longitude of reference point on \code{grid}.}
  \item{lat}{ latitude of reference point on \code{grid}.}

}

