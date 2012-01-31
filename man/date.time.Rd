\name{date.time}
\alias{date.time}
\title{Generation of dates and hours between start and end terms}
\usage{
  date.time(day.begin, day.end, h.begin = 0, h.end = 23,
    date.format = "y/m/d")
}
\arguments{
  \item{day.begin}{begin date (day) - format specified in
  \code{date.format}}

  \item{day.end}{end date (day) - format specified in
  \code{date.format}}

  \item{h.begin}{begin time (hour - integer)}

  \item{h.end}{end time (hour - integer)}

  \item{date.format}{input date format (formats for
  function \code{chron})}
}
\value{
  A 4-column table having the folloging fields (integer):
  "year", "month", "day", "hour"
}
\description{
  Generates dates and hours between start and end terms.
  The result is a data frame with years, months, days and
  hours. Called in function \code{\link{Th_int_series}}
}
\note{
  Input start and end dates as \code{character} (default
  format "yyyy/mm/dd"), hours as \code{integers} (0 to 23).

  Date format can be changed according to package
  \code{chron}'s standard, e.g.  "y/m/d" (default) or
  "m/d/y"
}
\examples{
data(Trentino_hourly_T)
date<-date.time(day.begin="01/01/2004", day.end="31/12/2005", date.format= "d/m/y")
}
\author{
  Emanuele Eccel, Emanuele Cordano
  \email{emanuele.eccel@iasma.it}
}
\seealso{
  \code{\link{Th_int_series}}
}

