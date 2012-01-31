NULL
#' 
#' The function creates sets of hourly temperature series, in a specified period, from minimum and maximum daily series. The calibration files, created by the functions \code{par_calibration} and \code{shape_calibration}, are loaded and passed to function \code{Th_interp}.
#'
#' @title Hourly interpolation of multiple daily temperature series
#' 
#' 
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param cal_times  calibration table of "time" parameters for the specific series
#' @param cal_shape  calibration table for "shape" parameter(s) for the specific series
#' @param TMIN  minimum temperature daily table
#' @param TMAX  maximum temperature daily table
#' @param start_year  year of simulation start
#' @param end_year  year of simulation end
#' @param night_adjust  logical: if set to \code{TRUE} applies the correction for minima occurring at late hours of the day
#' @param active_IDs  a set of series IDs to be interpolated. If \code{NULL} (default), all series are interpolated
#' @param min_mo.length minimum number of days necessary to calculate monthly dtr values
#' @param silent logical, if set to \code{TRUE} suppresses any warning issue

#' @export 

#' @return A list of interpolated hourly temperatures for the \code{active_IDs} series

#' @references 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50 
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf}
#' 
#' 
#' Original algorithm from: Cesaraccio, C., Spano, D., Duce, P., Snyder, R.L., 2001. An improved model for determining degree-day values from daily temperature data. Int. J. Biometeorol. 45: 161-169.
#' \url{http://www.springerlink.com/content/qwctkmlq3tebthek/}
#' 
#'  
#' See also: Eccel, E., 2010: What we can ask to hourly temperature recording. Part I: statistical vs. meteorological meaning of minimum temperature. Italian Journal of Agrometeorology XV(2):41-43.
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf}
#' 
#' @note
#' \code{TMIN} and \code{TMAX} are data frames with the first three columns devoted to the date. The names of these columns must be "year", "month", and "day", irrespective of their order. Data series range from column 4 to the last.
#'
#' If the series is one with non-null results of the \code{par_calibration} function (enough data for calibration) its table is passed to the interpolation function, otherwise the average (\code{cal_table}) is used.
#'
#' \code{delta.T_night} is used only if night_adjust is \code{TRUE}.
#'
#' Tmin of the day before the first is set = to Tmin of the first day and Tmin of the day after the last = Tmin of the last day.
#'
#' Since the first value of T at sunset (of the day before) is \code{NULL}, the first hourly values produced till \code{time_min} are \code{NA}.

#' @examples
#' data(Trentino_hourly_T)
#' stations <- c("T0001","T0010","T0129")
#' # interpolation of temperature for series T0001 and T0129, from 1970 to 1972
#' Th_int_list <- Th_int_series(cal_times = calibration_l, cal_shape = calibration_shape, TMIN=Tn, TMAX=Tx, start_year = 2004, end_year = 2005, active_IDs = stations)

#' @seealso \code{\link{Th_interp}},  \code{\link{Mo.Th.Ra.}},  \code{\link{date.time}}


###################################################################
# INTERPOLATES DAILY TEMPERATURE **SERIES**, MORE SERIES AT ONCE 
###################################################################

Th_int_series<-function(cal_times, cal_shape, TMIN, TMAX, start_year, end_year, night_adjust=FALSE, active_IDs=NULL, min_mo.length=21, silent=FALSE)

{


if(length(cal_shape) == 2)   # delta.T_night included
 delta.T_night_aver<-mean(cal_shape$delta.T_night$delta.T_night.min, na.rm=T)

if(is.null(active_IDs))
  IDs<-names(TMIN)[-c(1:3)]  else IDs <- active_IDs

IDs_OK<-names(cal_times)[-length(names(cal_times))]  # all except "Average"

month.year.day<-TMIN[,c(1:3)]
d<-paste("0",month.year.day$day[1],sep=""); m<-paste("0",month.year.day$month[1],sep="")
start<-paste(month.year.day$year[1],"/", substr(m,nchar(m)-1, nchar(m)), "/", substr(d,nchar(d)-1, nchar(d)), sep="")
d<-paste("0",month.year.day$day[nrow(month.year.day)],sep=""); m<-paste("0",month.year.day$month[nrow(month.year.day)],sep="")
end<-paste(month.year.day$year[nrow(month.year.day)],"/", substr(m,nchar(m)-1, nchar(m)), "/", substr(d,nchar(d)-1, nchar(d)), sep="")
date<-date.time(day.begin=start, day.end=end)
if(date$year[1] > start_year | date$year[nrow(date)] < end_year) print("Warning: interpolation period required exceeds daily series limits!", quote=F) 
date<-date[date$year>=start_year & date$year<=end_year,]

interp_list<-list(date)

for(id in IDs)  

{
if(silent==FALSE)
 print(paste("Now processing series", id),quote=F)

Tmin<-data.frame(TMIN[,1:3], T=TMIN[,id])
Tmax<-data.frame(TMAX[,1:3], T=TMAX[,id])
  
dtr<-Mo.Th.Ra.(Tmin=Tmin, Tmax=Tmax, name=id, min_mo.length=min_mo.length, silent=silent)  

Tmin<-Tmin[Tmin$year>=start_year & Tmin$year<=end_year,]
Tmax<-Tmax[Tmax$year>=start_year & Tmax$year<=end_year,]

# if series OK loads the corresponding calibr. table and substit. missing parameters with the average ones; if series not OK chooses the average table
if(id %in% IDs_OK) {
calibr<-cal_times[[id]] 
calibr[is.na(calibr)]<- cal_times$Average[is.na(calibr)] } else calibr<-cal_times$Average

ratio_dtr_aver<-mean(cal_shape$ratio$ratio.min, na.rm=T)
ratio_dtr<-cal_shape$ratio[id,1]; if(is.na(ratio_dtr)) ratio_dtr<-ratio_dtr_aver

delta.T_night<-NULL
if(night_adjust)
 { delta.T_night<-cal_shape$delta.T_night[id,1]; if(is.na(delta.T_night)) delta.T_night<-delta.T_night_aver }

# interpolates for each series
Th_int<-NULL
Ts<-NULL
for(d in 1:nrow(Tmin)) 
{  
 if(exists("Th_l")) Ts<-Th_l$Tsuns 
 Th_l<-Th_interp(Tmin=Tmin, Tmax=Tmax, day=d, Tsuns=Ts, tab_calibr=calibr, dtr_month=dtr, ratio_dtr=ratio_dtr, delta.T_night=delta.T_night)
 Th<-Th_l$Th
 Th_int<-append(Th_int, round(Th,1))
}
interp_list<-append(interp_list, list(Th_int))

} # stations


names(interp_list)<-c("Date", IDs)  

return(interp_list)

}

