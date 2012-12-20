\name{shape_calibration}
\alias{shape_calibration}
\title{Calibrates the shape of the night interpolating curve}
\usage{
  shape_calibration(meas, date.format = "ymd",
    cal_times_list, band_min = 0:23, band_max = 0:23,
    ratio_dtr_range = c(0, 6), nr_cycles = 10,
    min_mo.length = 21, full.24.hrs.span_min = TRUE,
    silent = FALSE)
}
\arguments{
  \item{meas}{measured hourly values file (table), where
  the first column is the series' ID}

  \item{date.format}{input date format (formats for
  function \code{chron})}

  \item{cal_times_list}{calibration list of "time"
  parameters (output of \code{\link{par_calibration}})}

  \item{band_min}{band of hours of occurrence of day
  minimum in the daily series (continuous). See
  \code{Note}}

  \item{band_max}{same for maximum time}

  \item{ratio_dtr_range}{range for seeking the optimal
  value of \code{ratio_dtr}}

  \item{nr_cycles}{number of calibration trials within the
  calibration ranges (all)}

  \item{min_mo.length}{minimum number of days to calculate
  any monthly values of dtr (is passed to function
  \code{\link{Mo.Th.Ra.}})}

  \item{full.24.hrs.span_min}{logical, if set to
  \code{FALSE} does not allow to shift minimum time to the
  late hours of the day}

  \item{silent}{logical, if set to \code{TRUE} suppresses
  any warning issue}
}
\value{
  A list containing the optimum values of \code{ratio_dtr}
}
\description{
  Calibrates the shape of the night interpolating curve,
  either horizontal-axe parabola or line, by changing the
  exponent z (see reference). It functions according to the
  comparison of the daily thermal range and the climate
  (reference) monthly one.
}
\note{
  \code{meas} must be organized as 4-field records, all
  series in the same file, no headers. Column order:
  station ID, date, time (hour), T, [others fields, if
  any...] separated by spaces. This field order is
  mandatory.

  Default date format is "ymd" (yyyy/mm/dd). Different
  combinations can be passed to function with
  \code{date.format}, but separator must be "/"

  \code{band_min} and \code{band_max} are the time bands
  according to which the minimum and maximum temperature
  were calculated in the daily series to be interpolated.
  In general, they range from 0 to 23, unless the series
  has had some restriction in the calculation of minimum
  and maximum values. Hence, these bands can be different
  from those used to calibrate the most frequent occurrence
  of min and max.

  The optimal value of \code{ratio_dtr} (\code{k}, eq. 7,
  in the quoted reference Eccel (2010a)) is chosen as the
  one with the (absolute) minimum value of the bias
  (irrespective of its sign). \code{ratio_dtr} is the ratio
  between the daily thermal range of the day to be
  interpolated and the mean monthly value for that series.
  The corresponding values of mean absolute error and RMSE
  can be checked in the resulting list.

  \code{min_mo.length} (passed to function
  \code{Mo.Th.Ra.}) refers to the sum of days along all the
  series for each specific month, not for any single month
  in one year.

  \code{full.24.hrs.span_min} is \code{TRUE} as default. If
  must be set to \code{FALSE} only if minimum values of the
  daily series have been calculated on a restricted time
  band, which is included in \code{band_min} (see function
  \code{par_calibration}). If this is the case, the minimum
  of the interpolated curve will always fall within
  \code{band_min} (early hours of the day). If this option
  is erroneously chosen, errors as large as 0.6 deg C can
  arise in the average of mean daily T.
}
\examples{
library(Interpol.T)
data(Trentino_hourly_T)

stations <- c("T0001","T0010","T0129")

calibration_shape <- shape_calibration(meas = h_d_t[h_d_t$V1 \%in\% stations,],
						cal_times_list = calibration_l[stations],
                     band_min = 0:23, band_max = 0:23, ratio_dtr_range = c(0,4),
						min_mo.length=21)
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
  \code{\link{par_calibration}}, \code{\link{Th_interp}}
}

