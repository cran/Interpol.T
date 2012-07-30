\name{Th_int_series}
\alias{Th_int_series}
\title{Hourly interpolation of multiple daily temperature series}
\usage{
  Th_int_series(cal_times, TMIN, TMAX, start_year,
    end_year, cal_shape = NULL, active_IDs = NULL,
    min_mo.length = 21, full.24.hrs.span_min = TRUE,
    silent = FALSE)
}
\arguments{
  \item{cal_times}{calibration table of "time" parameters
  for the specific series}

  \item{cal_shape}{calibration table for "shape" parameter
  for the specific series; default is \code{NULL}(don't
  use)}

  \item{TMIN}{minimum temperature daily table}

  \item{TMAX}{maximum temperature daily table}

  \item{start_year}{year of simulation start}

  \item{end_year}{year of simulation end}

  \item{active_IDs}{a set of series IDs to be interpolated.
  If \code{NULL} (default), all series are interpolated}

  \item{min_mo.length}{minimum number of days necessary to
  calculate monthly dtr values}

  \item{full.24.hrs.span_min}{logical, if set to
  \code{FALSE} does not allow to shift minimum time to the
  late hours of the day}

  \item{silent}{logical: if set to \code{TRUE} suppresses
  any warning issue}
}
\value{
  A list of interpolated hourly temperatures for the
  \code{active_IDs} series
}
\description{
  The function creates sets of hourly temperature series,
  in a specified period, from minimum and maximum daily
  series. The calibration files, created by the functions
  \code{par_calibration} and \code{shape_calibration}, are
  loaded and passed to function \code{Th_interp}.
}
\note{
  \code{TMIN} and \code{TMAX} are data frames with the
  first three columns devoted to the date. The names of
  these columns must be "year", "month", and "day",
  irrespective of their order. Data series range from
  column 4 to the last.

  If the series ID coincides with one with non-null results
  of the \code{par_calibration} function (enough data for
  calibration) its table is passed to the interpolation
  function, otherwise the average (\code{cal_table}) is
  used.

  \code{full.24.hrs.span_min} is \code{TRUE} as default. If
  must be set to \code{FALSE} only if minimum values of the
  daily series have been calculated on a restricted time
  band, which is included in \code{band_min} (see function
  \code{par_calibration}). If this is the case, the minimum
  of the interpolated curve will always fall within
  \code{band_min} (early hours of the day). If this option
  is erroneously chosen, errors as large as 0.6 deg C can
  arise in the average of mean daily T.

  Tmin of the day before the first is set = to Tmin of the
  first day and Tmin of the day after the last = Tmin of
  the last day.

  Since the first value of T at sunset (of the day before)
  is \code{NULL}, the first hourly values produced till
  \code{time_min} are = Tmin.
}
\examples{
library(Interpol.T)
data(Trentino_hourly_T)
stations <- c("T0001","T0010","T0129")
# interpolation of temperature for series T0001 and T0129, from 2004 to 2005
Th_int_list <- Th_int_series(cal_times = calibration_l, cal_shape = calibration_shape, TMIN=Tn, TMAX=Tx, start_year = 2004, end_year = 2005, active_IDs = stations)
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
  \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf}

  Original algorithm from: Cesaraccio, C., Spano, D., Duce,
  P., Snyder, R.L., 2001. An improved model for determining
  degree-day values from daily temperature data. Int. J.
  Biometeorol. 45: 161-169.
  \url{http://www.springerlink.com/content/qwctkmlq3tebthek/}

  See also: Eccel, E., 2010: What we can ask to hourly
  temperature recording. Part I: statistical vs.
  meteorological meaning of minimum temperature. Italian
  Journal of Agrometeorology XV(2):41-43.
  \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf}
}
\seealso{
  \code{\link{Th_interp}}, \code{\link{Mo.Th.Ra.}},
  \code{\link{date.time}}
}

