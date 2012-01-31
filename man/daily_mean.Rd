\name{daily_mean}
\alias{daily_mean}
\title{Production of daily means from hourly temperature series}
\usage{
  daily_mean(hourly_list, series_names = NULL)
}
\arguments{
  \item{hourly_list}{the list of hourly temperatures}

  \item{series_names}{names of the serie to be averaged (if
  \code{NULL} (default): all series)}
}
\value{
  A list of daily averaged series
}
\description{
  The function works on a list of hourly temperature
  series. The hourly list is the output of the
  interpolation function \code{Th_int}, called iteratively
  to produce a list where each component represents one
  interpolated series.
}
\note{
  The first element of \code{hourly_list} must be a data
  frame named "Date" and its columns "year", "month", "day"
  (a fourth column ("hours") is not used in this function)
}
\examples{
data(Trentino_hourly_T)
# generates daily means for series T0001 and T0129:
Tm_list <- daily_mean(hourly_list = Th_int_list, series_names = c("T0001", "T0129"))
}
\author{
  Emanuele Eccel, Emanuele Cordano
  \email{emanuele.eccel@iasma.it}
}
\references{
  Eccel, E., 2010: What we can ask to hourly temperature
  recording. Part II: hourly interpolation of temperatures
  for climatology and modelling. Italian Journal of
  Agrometeorology XV(2):45-50

  Eccel, E., 2010: What we can ask to hourly temperature
  recording. Part I: statistical vs. meteorological meaning
  of minimum temperature. Italian Journal of
  Agrometeorology XV(2):41-43.
}
\seealso{
  \code{\link{Th_int_list}}
}

