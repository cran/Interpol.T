NULL
#' 
#' The function creates 24 values of hourly temperature from minimum and maximum daily values. This function applies to single series and to single day couples of minimum and maximum temperature. It is called by functions \code{Th_int_series} and \code{shape_calibration}.
#' The function uses three different curves: from minimum to maximum times: an increasing sinusoidal curve; from maximum time to sunset: a decreasing sinusoidal curve; from sunset to the minimum of the following day: either a horizontal-axis parabola or a line, according to the daily thermal range of the day.
#' Calibration parameters are series- and monthly-specific.
#' This function is operationally called by \code{Th_int_series}, which requires the daily series and the calibration files as input (plus other parameters). A general user will conveniently use the latter function. 
#' @title 24-hourly interpolation of temperature
#'
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param Tmin  a daily table of 4 named columns, the first are year, month, day, the 4th is minimum temperature. The column names "month" and "T" are mandatory
#' @param Tmax  same for Tmax
#' @param Tsuns  temperature at sunset time
#' @param day  progressive number of the day (row of both \code{Tmin} and \code{Tmax}), corresponding to a day
#' @param tab_calibr  "hour" parameter calibration table for the specific series. See \code{\link{par_calibration}}
#' @param dtr_month  monthly daily thermal range table (see function \code{Mo.Th.Ra.})
#' @param ratio_dtr  parameter for the choice of the night curve shape
#' @param delta.T_night  parameter (optional) for the detection and correction of minima occurring at late hours of the day


#' @export 

#' @return A vector containing the values from hour = 0 (element 1) to hour = 23 (element 24)

#' @references 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50 
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf},\url{www.agrometeorologia.it}
#' 
#' 
#' Original algorithm from: Cesaraccio, C., Spano, D., Duce, P., Snyder, R.L., 2001. An improved model for determining degree-day values from daily temperature data. Int. J. Biometeorol. 45: 161-169.
#' \url{http://www.springerlink.com/content/qwctkmlq3tebthek/}
#'  	
#' See also: Eccel, E., 2010: What we can ask to hourly temperature recording. Part I: statistical vs. meteorological meaning of minimum temperature. Italian Journal of Agrometeorology XV(2):41-43.
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf},\url{www.agrometeorologia.it}
#' 
#' 
#' @note
#' The function is called by \code{Th_int_series}.
#'
#' If the series is one with non-null results of the \code{par_calibration} function (enough data for calibration) its table is passed to the interpolation function, otherwise the average (\code{cal_table}) is used
#'
#' \code{delta.T_night}: used only if night_adjust is \code{TRUE}
#'
#' \code{delta.T_cl.nig.} is used (passed to function \code{Th_interp}) only if \code{night_adjust} is \code{TRUE}
#'
#' Tmin of the day before the first is set = to Tmin of the first day and Tmin of the day after the last = Tmin of the last day.
#'
#' Since the first value of T at sunset (of the day before) is \code{NULL}, the first hourly values produced till \code{time_min} are \code{NA}.



#' # interpolates for each series, as called from function \code{Th_int_series}:
#' @seealso \code{\link{Th_int_list}}



#########################################################
# INTERPOLATES DAILY TEMPERATURE (ONE DAY, ONE SERIES)
# CALLED BY Th_int_series AND BY shape_calibration
#########################################################


Th_interp<-function(Tmin, Tmax, Tsuns=NULL, day, tab_calibr, dtr_month, ratio_dtr, delta.T_night=NULL)

{
Th<-NULL
mm<-Tmin$month[day]
Tn<-Tmin$T[day] 
Tx<-Tmax$T[day] 
Tsuns_d.before<-Tsuns
if(!is.na(Tmin[day+1,1])) Tn_after<-Tmin$T[day+1] else Tn_after<-Tmin$T[day]
cloudy_morning<-!is.na(Tx-Tn) & (Tx-Tn)/dtr_month[Tmin$month[day]]<=ratio_dtr  
cloudy_night<-!is.na(Tx-Tn_after) & (Tx-Tn_after)/dtr_month[Tmin$month[day]]<=ratio_dtr 
Tsuns<-Tx - tab_calibr$C_m[mm]*(Tx-Tn_after)  
time_min<-tab_calibr$time_min[mm]
time_max<-tab_calibr$time_max[mm]
time_suns<-tab_calibr$time_suns[mm]
Th_24_before <- Tn
Tn_after <- Tn


# h = 0 to min
Tn_original<-Tn
if(cloudy_morning==TRUE){z<-1} else z<-0.5 
if(!is.null(delta.T_night) )  
 if(Tn < Th_24_before - delta.T_night & !is.na(Th_24_before))
  Tn <- Th_24_before - delta.T_night
b_d.before<-(Tn-Tsuns_d.before)/((time_min + 24-time_suns)^z)
for(h in 0:(time_min))  #  Th[1] = Th at time h = 0 and so on
 Th[h+1]<-Tsuns_d.before + b_d.before*((h+24-time_suns)^z)

# h = min to max
for(h in time_min:time_max) 
 Th[h+1]<-Tn+(Tx-Tn)/2 * ( 1 + sin((h-time_min)/(time_max-time_min)*pi - pi/2)) 

# h = max to sunset
for(h in (time_max):(time_suns)) 
 Th[h+1]<- Tsuns + (Tx-Tsuns)*sin(pi/2*(1+(h-time_max)/(time_suns-time_max)))     

# h = sunset to 23
if(cloudy_night==TRUE) {z<-1} else z<-0.5  
if(!is.null(delta.T_night) )  # apply the correction
 if(Tn_after  > Tn_original - delta.T_night & !is.na(Tn_original > Tn_after) )
   Tn_after <- Tn_original - delta.T_night
b<-(Tn_after - Tsuns) / ((time_min + 24-time_suns)^z)
for(h in time_suns:23)
 Th[h+1]<-Tsuns + b*((h-time_suns)^z)
Th_24_before<-Th[24]
 
Th_day.before<-Th 

Th_list<-list(Th,Tsuns); names(Th_list)<-c("Th", "Tsuns")

return(Th_list)

}
