\name{Mo.Th.Ra.}
\alias{Mo.Th.Ra.}
\title{Creates a table of mean monthly "daily thermal range" (dtr)}
\usage{
  Mo.Th.Ra.(Tmin, Tmax, name, min_mo.length = 21,
    silent = FALSE)
}
\arguments{
  \item{Tmin}{a table of 4 named columns, the first are
  year, month, day, the 4th is minimum temperature. The
  column names "month" and "T" are mandatory}

  \item{Tmax}{same for Tmax}

  \item{name}{name of the series dtr is calculated for}

  \item{min_mo.length}{minimum number of days necessary to
  calculate monthly dtr values}

  \item{silent}{logical. If \code{TRUE} no warning is
  issued}
}
\value{
  A vector of 12 dtr values (1 for January, ... 12 for
  December)
}
\description{
  Creates a table of mean monthly DTR (difference between
  maximum and minimum temperature) over the series' period,
  to use a reference for establishing the exponent "z" of
  the night stretch of the interpolation curve.
}
\examples{
data(Trentino_hourly_T)
id<-"T0001"
Tmin<-data.frame(Tn[,1:3], T=Tn[,id])
Tmax<-data.frame(Tx[,1:3], T=Tx[,id])
dtr <- Mo.Th.Ra.(Tmin = Tmin, Tmax = Tmax, name = id)
}
\author{
  Emanuele Eccel, Emanuele Cordano
  \email{emanuele.eccel@iasma.it}
}
\seealso{
  \code{\link{shape_calibration}},
  \code{\link{Th_int_series}}
}

