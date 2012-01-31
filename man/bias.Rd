\name{bias}
\alias{bias}
\title{Calculates mean bias (difference between (max+min)/2 and 24-hour averages) in mean daily temperature series}
\usage{
  bias(TMIN, TMAX, TMEAN, min_valid = 21)
}
\arguments{
  \item{TMIN}{data frame with daily minimum temperatures in
  columns. The first 3 columns are skipped (dates as year,
  month and day are supposed to be stored in these
  columns)}

  \item{TMAX}{same for TMAX}

  \item{TMEAN}{same for TMEAN. Should come from 24-hour
  daily means.}

  \item{min_valid}{min nr. of valid days in a month for
  retaining its average value (if valid days are fewer,
  monthly value is \code{NA}). Default is 21.}
}
\value{
  A vector of means of daily biases, where the \code{TMEAN}
  is considered the "true" (reference) value
}
\description{
  Calculates the average difference between the series of
  mean daily temperatures calculated by (Tmax + Tmin) / 2
  and the average calculated by 24 hourly values a day, as
  resulting from the interpolation (or from measurements).
  The function works on data tables with series on columns.
}
\note{
  Biases are calculated only on columns that are present in
  both \code{TMIN}/\code{TMAX} and \code{TMEAN}
}
\examples{
data(Trentino_hourly_T)
mo_bias <- bias(TMIN = Tn, TMAX = Tx, TMEAN = Tm_list, min_valid = 20)
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
  \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf},\url{www.agrometeorologia.it}

  See also: Eccel, E., 2010: What we can ask to hourly
  temperature recording. Part I: statistical vs.
  meteorological meaning of minimum temperature. Italian
  Journal of Agrometeorology XV(2):41-43.
  \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf},\url{www.agrometeorologia.it}
}
\seealso{
  \code{\link{daily_mean}}
}

