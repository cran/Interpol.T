\name{Trentino_hourly_T}
\docType{data}
\alias{Trentino_hourly_T}
\alias{Tn}
\alias{Tx}
\alias{h_d_t}
\alias{mo_bias}
\alias{Th_int_list}
\alias{Tm_list}
\alias{calibration_l}
\alias{calibration_shape}

\title{
Dataset of hourly temperature for Trentino, Italy
}

\usage{
data(Trentino_hourly_T)
}

\description{
Contains the following objects: 
\describe{
\item{\code{Tn}}{Data frame with \code{year},\code{month}, \code{day} and daily minimum temperature in 39 stations in Trentino, Italy. This series has been directly calculated from the hourly one (\code{h_d_t}), so it contains the daily minima and maxima of hourly measurements, not daily absolute minima and maxima.}
\cr
\item{\code{Tx}}{Data frame containing \code{year},\code{month} , \code{day} and daily maximum temperature in 39 stations in Trentino, Italy. This series has been directly calculated from the hourly one (\code{h_d_t}), so it contains the daily minima of hourly measurements.}
\cr
\item{\code{h_d_t}}{Data frame containing station id, date (yyyy/mm/dd), hour (integer 0,...23), hourly temperature, and a quality flag (not used) in 39 stations in Trentino.}
\cr
\item{\code{mo_bias}}{Data frame containing the monthly (lines 1-12) and annual (line 13) average bias between mean daily temperatures calculated as (Tmin + Tmax)/2 and the corresponding 24-value mean, for 39 stations in Trentino.}
\cr
\item{\code{Th_int_list}}{List containing the simulated (interpolated) hourly values. The first element (\code{Date}) is a data frame of \code{year}, \code{month}, \code{day}, and \code{hour}. All other elements, each having the name of one station id, are numeric vectors of temperature.}
\cr
\item{\code{Tm_list}}{List containing the daily means, calculated from the hourly interpolations. The first element (\code{Date}) is a data frame of \code{year}, \code{month}, and \code{day}. All other elements, each having the name of one station id, are numeric data frames of daily temperature.}
\cr
\item{\code{calibration_l}}{List containing the "hour" calibration parameters. Each element (name: station id) reports the following parameters, for each month: \code{time_min}: mode value of the time of occurrence of minimum temperature;  \code{time_max}: mode value of the time of occurrence of maximum temperature; \code{time_suns}: mode value of the time of occurrence of sunset; \code{C_m}: value of "c" (see \code{\link{par_calibration}} and the quoted reference). An unreported station id means insufficient data for calibration.}
\cr
\item{\code{calibration_shape}}{List containing the "shape" calibration parameters for the night portion of the curve. The list has two data frames. The first (name: \code{ratio}) has in every line (one line per station id) the following: 1. the value of \code{ratio_dtr} that minimizes the mean error. 2., 3., 4. the mean error, mean absolute error, root mean square error, respectively, corresponding to this value (all in "deg C"). The second data frame (name: \code{delta.T_night}) has the same data organization of the first one, but refers to the delta of temperature betwen time 0 and minimum of the same day to apply the night correction algorithm (minimum time is shifted form the calibrated time of minimum to h=23), if this option is required. For details see \code{\link{shape_calibration}}.}
\cr

}
}
\details{
Dataset from Trentino, Italy, with examples of hourly temperature series, the daily series (min and max) obtained from these, and the results of the calibration and application of the interpolation algorithm: the calibration lists (parameters), the interpolated list of series, and the bias between interpolated 24-hour daily means and the means obtained by (Tmin + Tmax)/2.
The user can easily use the package with his/her own data after replacing the values of such variables.
}
\format{
Data frames and lists.
}
\source{
Original data are provided by Provincia Autonoma di Trento, Italy (\url{http://www.meteotrentino.it/}).

This dataset is intended for research purposes only, being distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
}

\references{
Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50
\url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf},\url{www.agrometeorologia.it}
}

\keyword{dataset}
