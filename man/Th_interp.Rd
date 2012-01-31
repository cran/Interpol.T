\name{Th_interp}
\alias{Th_interp}
\title{24-hourly interpolation of temperature}
\usage{
  Th_interp(Tmin, Tmax, Tsuns = NULL, day, tab_calibr,
    dtr_month, ratio_dtr, delta.T_night = NULL)
}
\arguments{
  \item{Tmin}{a daily table of 4 named columns, the first
  are year, month, day, the 4th is minimum temperature. The
  column names "month" and "T" are mandatory}

  \item{Tmax}{same for Tmax}

  \item{Tsuns}{temperature at sunset time}

  \item{day}{progressive number of the day (row of both
  \code{Tmin} and \code{Tmax}), corresponding to a day}

  \item{tab_calibr}{"hour" parameter calibration table for
  the specific series. See \code{\link{par_calibration}}}

  \item{dtr_month}{monthly daily thermal range table (see
  function \code{Mo.Th.Ra.})}

  \item{ratio_dtr}{parameter for the choice of the night
  curve shape}

  \item{delta.T_night}{parameter (optional) for the
  detection and correction of minima occurring at late
  hours of the day}
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
  function uses three different curves: from minimum to
  maximum times: an increasing sinusoidal curve; from
  maximum time to sunset: a decreasing sinusoidal curve;
  from sunset to the minimum of the following day: either a
  horizontal-axis parabola or a line, according to the
  daily thermal range of the day. Calibration parameters
  are series- and monthly-specific. This function is
  operationally called by \code{Th_int_series}, which
  requires the daily series and the calibration files as
  input (plus other parameters). A general user will
  conveniently use the latter function.
}
\note{
  The function is called by \code{Th_int_series}.

  If the series is one with non-null results of the
  \code{par_calibration} function (enough data for
  calibration) its table is passed to the interpolation
  function, otherwise the average (\code{cal_table}) is
  used

  \code{delta.T_night}: used only if night_adjust is
  \code{TRUE}

  \code{delta.T_cl.nig.} is used (passed to function
  \code{Th_interp}) only if \code{night_adjust} is
  \code{TRUE}

  Tmin of the day before the first is set = to Tmin of the
  first day and Tmin of the day after the last = Tmin of
  the last day.

  Since the first value of T at sunset (of the day before)
  is \code{NULL}, the first hourly values produced till
  \code{time_min} are \code{NA}. # interpolates for each
  series, as called from function \code{Th_int_series}:
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

