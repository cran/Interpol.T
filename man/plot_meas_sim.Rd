\name{plot_meas_sim}
\alias{plot_meas_sim}
\title{Plots hourly curves}
\usage{
  plot_meas_sim(meas, sim, series = NULL, chart.start,
    chart.end, date.format = "ymd", missing_code = NA,
    wait = 1, plot.leg = TRUE, leg.pos = "bottomright")
}
\arguments{
  \item{meas}{measured hourly values file (table), where
  the first column is the series' ID}

  \item{sim}{simulated hourly list (output of
  Th_int_series)}

  \item{series}{names of the seris to plot. If \code{NULL},
  plots all the series}

  \item{chart.start}{start date for the plotting. Format
  example: "1Jan2000"}

  \item{chart.end}{end date for the plotting. Format
  example: "1Jan2000"}

  \item{date.format}{input date format for measurements
  (formats for function \code{chron}). Default is "ymd"}

  \item{missing_code}{code (either real or character) for
  missing values in measurements. Default is \code{NA}}

  \item{wait}{lag time (seconds) between plot appearance on
  the screen (default is 1 second)}

  \item{plot.leg}{logical: if \code{TRUE} (default) legends
  are plotted}

  \item{leg.pos}{position of legends (only if
  \code{plot.leg} = \code{TRUE}). Default is "bottomright".
  Values for \code{par} can be passed to function}
}
\value{
  A plot with two curves: measured values and hourly
  interpolations
}
\description{
  Plots curves for comparison between measures and
  simulations (hourly interpolations)
}
\note{
  If daily minimum and maximum values are the absolute ones
  for each day, the interpolated curve will be generally
  either higher (maximum) or lower (minimum) than hourly
  measurements, which are mean hourly values (hence,
  slightly lower and higher than the absolute daily ones,
  respectively). This may explain an apparent mismatch
  between the two curves.
}
\examples{
data(Trentino_hourly_T)
stations <- c("T0001","T0010","T0129")
plot_meas_sim(meas = h_d_t, sim = Th_int_list, series=stations,
              missing_code = -999.9, chart.start = "1Feb2004",
             chart.end = "29Feb2004", leg.pos = "top")
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
  \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf}
}

