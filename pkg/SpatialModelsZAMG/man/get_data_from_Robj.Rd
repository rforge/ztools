\name{get_data_from_Robj}
\alias{get_data_from_Robj}
\title{Extract the snow data from a saved \code{robj} file}
\description{this function extracts the snow data from a saved \code{robj} file which consists of a list of stations like described in the \emph{Details} section}
\usage{
get_data_from_Robj(dataname, measurements_barrier = 10)
}
\arguments{
  \item{dataname}{
  the dataname of the saved \code{robj} file as a character string}
  \item{measurements_barrier}{how many measurements does each station has to have at least. those stations are used later on in the optimization (see functions \code{optimizer_\ldots_model}). \cr default is \code{measurements_barrier = 10}}
}
\value{
a list with
  \item{snow_data }{the data of the saved file (as a list). \cr all stations with less than \code{measurements_barrier} measurements have been deleted}
  \item{sd_max_data }{a matrix with the yearly maxima of snow depth (sd). each row corresponds to one station, columns are the corresponding years. matrix might contain \code{NA}'s}
  \item{swe_max_data }{a matrix with the yearly maxima of snow water equivalent (swe). each row corresponds to one station, columns are the corresponding years. matrix might contain \code{NA}'s}
  \item{covariables }{a matrix with the covariables for each station. each row corresponds to one station, columns are \var{lon} (longitude), \var{lat} (latitude), \var{alt} (altitude), \var{mdday} (mean of the time difference in days of sd and swe maxima), \var{sd_mmax} (mean maxima of sd over all years) and \var{swe_mmax} (mean maxima of swe over all years)}
}
\details{
the saved \code{robj} file should include a list of stations with dataframes \var{meta} and \var{data}. stations are ordered via station number.

\var{meta} should include \var{stnr}, \var{lon}, \var{lat} and \var{alt}

\var{data} should include \var{year}, \var{sd}, \var{swe}, \var{date.sd} and \var{date.swe}

see also the sample snow data from 36 examplary stations distributed over Austria: \cr \file{\link{sample_snow_data}}
}
\seealso{
\code{\link{model_selection}}, \code{\link{optimizer_smooth_model}}, \code{\link{optimizer_biv_hr_model}}
}
\examples{
# load sample_snow_data and save it to working directory
data(sample_snow_data)
save(sample_snow_data, file = "sample_snow_data.robj")

get_data =
  get_data_from_Robj(dataname = "sample_snow_data.robj",
                     measurements_barrier = 3)
}
\keyword{functions}
