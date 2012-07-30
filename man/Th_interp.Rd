\name{Th_interp}
\alias{Th_interp}
\title{24-hourly interpolation of temperature}
\usage{
  Th_interp(Tmin, Tmax, Tsuns = NULL, Th_24_before = NULL,
    day, tab_calibr, dtr_month = NULL, ratio_dtr = NULL,
    late_min = TRUE)
}
\arguments{
  \item{Tmin}{a daily table of 4 named columns, the first 3
  being year, month, day, the 4th minimum temperature. The
  column names "month" and "T" are mandatory}

  \item{Tmax}{same for Tmax}

  \item{Tsuns}{temperature at sunset time}

  \item{Th_24_before}{temperature at time 24 of the
  previous day (time 00 of the present day)}

  \item{day}{progressive number of the day (row of both
  \code{Tmin} and \code{Tmax}), corresponding to a day}

  \item{tab_calibr}{"hour" parameter calibration table for
  the specific series. See \code{\link{par_calibration}}}

  \item{dtr_month}{monthly daily thermal range table (see
  function \code{Mo.Th.Ra.})}

  \item{ratio_dtr}{parameter for the choice of the night
  curve shape; it is \code{NULL} if no calibration_shape is
  passed to the function by \code{Th_int_series}}

  \item{late_min}{logical; allows to shift the time of
  occurrence of minima to the late hours of the day
  (assumes the value of \code{full.24.hrs.span_min}, passed
  by functions \code{Th_int_series)} and
  \code{shape_calibration})}
}
\value{
  A vector containing the values from hour = 0 (element 1)
  to hour = 23 (element 24)
}
\description{
  The function creates 24 values of hourly temperature from
  minimum and maximum daily values. This function applies
  to single series and to single day couples of minimum and
  maximum temperature. It is called by functions
  \code{Th_int_series} and \code{shape_calibration}. The
  function uses four different curves: from time 00 to the
  minimum time: a horizontal-axis parabola (a line, if this
  choice is enabled and according to the daily thermal
  range of the day); from minimum to maximum time: an
  increasing sinusoidal curve; from maximum time to sunset:
  a decreasing sinusoidal curve; from sunset to time = 23:
  a horizontal-axis parabola (a line, if this choice is
  enabled and according to the daily thermal range of the
  day). Calibration parameters are series- and
  monthly-specific. This function is operationally called
  by \code{Th_int_series}, which requires the daily series
  and the calibration table as input (plus other
  parameters). A general user will conveniently use the
  latter function.
}
\note{
  The function is called by \code{Th_int_series}.

  If the series ID coincides with one with non-null results
  of the \code{par_calibration} function (enough data for
  calibration) its table is passed to the interpolation
  function, otherwise the average (\code{cal_table}) is
  used.

  A non-NULL value for \code{ratio_dtr} enables the
  function to interpolate night values with a line, if the
  conditions on Daily Thermal Range occur. This may give
  rise to a sharp change for the hours following the min.

  Tmin of the day before the first is set = to Tmin of the
  first day and Tmin of the day after the last = Tmin of
  the last day.

  If T from sunset falls below the minimum of the day
  (temptatively attributed to \code{time_min}), an
  adjustement is done: the early hours assume a constant T
  = T[00] and the minimum is shifted so that T[23] =
  \code{Tmin} for that day

  Since the very first value of T series at sunset (of the
  day before) is \code{NULL}, the first hourly values
  produced till \code{time_min} are = \code{Tmin} of the
  day.
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

  Original algorithm from: Cesaraccio, C., Spano, D., Duce,
  P., Snyder, R.L., 2001. An improved model for determining
  degree-day values from daily temperature data. Int. J.
  Biometeorol. 45: 161-169.
  \url{http://www.springerlink.com/content/qwctkmlq3tebthek/}

  See also: Eccel, E., 2010: What we can ask to hourly
  temperature recording. Part I: statistical vs.
  meteorological meaning of minimum temperature. Italian
  Journal of Agrometeorology XV(2):41-43.
  \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf},\url{www.agrometeorologia.it}
}
\seealso{
  \code{\link{Th_int_list}}
}

