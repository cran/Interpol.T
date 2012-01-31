\name{par_calibration}
\alias{par_calibration}
\title{Calibration of "hour" parameters for temperature interpolation}
\usage{
  par_calibration(meas, date.format = "ymd",
    cal_period = NULL, missing_value_code = NA,
    min_valid_yrs = 1, band_min = 4:9, band_max = 12:16,
    band_suns = 14:20, silent = FALSE, aver_series = NULL)
}
\arguments{
  \item{meas}{measured hourly values file (table), where
  the first column is the series' ID}

  \item{date.format}{date format in \code{meas}. Default is
  \code{"ymd"} (year month day)}

  \item{cal_period}{calibration period for the series.
  Default is \code{NULL} (use whole period)}

  \item{missing_value_code}{allows to deal missing data if
  a code (either numeric or character) is passed. Default
  is \code{NA}}

  \item{min_valid_yrs}{minimum fraction of year(s) on any
  series ID necessary for carrying out calibration}

  \item{band_min}{(continuous) band of hours to seek day
  minimum}

  \item{band_max}{same for maximum time}

  \item{band_suns}{same for sunset time}

  \item{silent}{if set to \code{TRUE} removes notice of
  insufficient length for the calibration of single short
  series}

  \item{aver_series}{set of series IDs (chr) used to
  calculate the average calibration (if \code{NULL}: all
  stations with valid data are included)}
}
\value{
  a list of calibration tables, one for each series, plus
  one average table at the bottom (named "Average")
}
\description{
  Calibrates the monthly parameters for the definition of
  the interpolation algorithm values of: time of minimum
  temperature; time of maximum temperature; time of sunset;
  parameter 'c' (see references for details). It works on
  more than one series at once. It also calculates an
  average calibration table, for simulations of series with
  no calibration.
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

  See also: Eccel, E., 2010: What we can ask to hourly
  temperature recording. Part I: statistical vs.
  meteorological meaning of minimum temperature. Italian
  Journal of Agrometeorology XV(2):41-43.
}
\examples{
data(Trentino_hourly_T)
stations <- c("T0001","T0010","T0129")

calibration_l <- par_calibration(meas = h_d_t[h_d_t$V1 \%in\% stations,], missing_value_code = -999.9, band_min = 4:8, band_max = 13:16, band_suns = 16:20, cal_period = NULL)
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
}
\seealso{
  \code{\link{shape_calibration}}
}

