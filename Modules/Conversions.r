#***************************** General conversion routines ******************************

#Latitude <- 48.3230 #Latitude in ddd.mmsssss format (degress, minutes, decimal seconds). Negative for southern hemisphere. Moved to preferences file.

func.DMS <- function(Latitude){ #It's a bit of overkill to put this in as a function, but it means that we can ask the user for the latitude, or call it from a file. 
	Ldegrees <- trunc(Latitude) #Loses the decimals
	Lminutes <- trunc((Latitude - Ldegrees)*100) #extracts minutes
	Lseconds <- 10000 * (Latitude - (Lminutes/100 + Ldegrees)) #extracts seconds
	Latitude <- Ldegrees + (Lminutes / 60) + (Lseconds / 3600) #This section converts dms to decimal degrees.
	Latitude <- Latitude * pi/180 #and convert to radians
	return(Latitude)
}#end function

#Converts month name to a number and adds the current proportion of the month. eg: Jan and halfway through monthly timesteps returns 1.5.
#note, this has change from Chimble 1.
#Month <- Month + (Timestep)%%1 is now Month + Rows[Time]%%1 * Timestep * 10 (grab decimals, and switch to fraction of month passed).
func.Month <- function (Time){
	if (class(Met[Time,"Month"]) == "character"){
	Month <- sapply((substr (MonthName,1,3)),function(x) grep(paste("(?i)",x,sep=""),month.abb)) 
	} else {
	Month <- Met[Time,"Month"]
	}
	if (Datarate == "Monthly") {
	Rows <- as.numeric(rownames(Met))
	Month <- Month + Rows[Time]%%1 * Timestep * 10
	} else {
	Month <- Month + (as.numeric(Met[Time,"Day"])/30) #assumes 30 days per month. Should be close enough. 
	}
	return(Month)
	#The grep component removes the case sensitivity from the match. This function assigns months by their name with a number. Jan = 1. This is used for Ra calculation and possibly elsewhere.
}# end function 

func.MonthNum <- function (Month){
	if (class(Month) == "character"){
	Month <- sapply((substr (Month,1,3)),function(x) grep(paste("(?i)",x,sep=""),month.abb)) 
	}
	return(Month)
	#The grep component removes the case sensitivity from the match. This function assigns months by their name with a number. Jan = 1. This is used for Ra calculation and possibly elsewhere.
}# end function 



func.RaVal <- function(Latitude){ #Ra from Valiantzas' paper. Checked against example on P697.
	#There is an error here. All latitudes should be in Radians - now fixed.
	#However, there is also a second error in equation 2 (118 * N * etc). This yields unrealistic values for tropical regions. Therefore func.Ra has been rewritten (below) to use the FAO calculations - Allen 1998 FAO Chapter 3. eq 28.
	N <- 4 * (func.DMS(Latitude)) * sin(0.53 * func.Month(Time) - 1.65) + 12 #Number of hours of sunlight. Required by Ra calc following.
	if (Mod(func.DMS(Latitude)) > 23.5*pi/180){ #Tropic check
		Ra <- 3 * N * sin((0.131 * N) - (0.95 * Mod(func.DMS(Latitude))))
	} else if (Mod(func.DMS(Latitude)) <= 23.5 * pi/180){ #end if
		Ra <- 118*(N^0.2)*sin(0.131*N - (0.2 * Mod(func.DMS(Latitude))))
	}#end else
	return (Ra)
}#end function


func.Ra <- function(Latitude){ #Ra from Valiantzas' paper. Checked against example on P697.
	#There is an error here. All latitudes should be in Radians.
	JulianDayApprox = (Met[Time,"Month"]-1)*30 + Met[Time,"Day"]
	Dr <- 1+0.033*cos(2*pi*JulianDayApprox/365)
	SolarDec <- 0.409*sin(2*pi*JulianDayApprox/365-1.39)
	Sunset <- acos(-tan(func.DMS(Latitude))*tan(SolarDec))
	Ra <- 24*60/pi * 0.082 * Dr * (Sunset * sin(func.DMS(Latitude))*sin(SolarDec)+cos(func.DMS(Latitude))*cos(SolarDec)*sin(Sunset))
	return (Ra)
}#end function


#Returns number of days in a month.
func.monthlength <- function (month){
if (Datarate == "Monthly") {
		#month <- func.Month(Time)
		if (month == 1|month == 3|month == 5|month == 7|month == 8|month == 10|month == 12){ #This section (and same in func.PET) changes the daily rate to a per month rate, based on the number of days in each month. This is different from Steinman who uses 30 days per month.
			Runtime <- 31
		}#end if
		else if (month ==  4|month == 6|month == 9|month == 11){
			Runtime <- 30
		} else {#end else if
			Runtime <- 28
		}#end else
} else {
Runtime = 1
}
return (Runtime)
}

#Better version of above
DaysInMonth <- function(month, year){
  leap = 0
  if ((year %% 4 == 0 & year %% 100 != 0) | year %% 400 == 0)
    leap = 1
  # Return the number of days in the month
  return(switch(month,
                '01' = 31,
                '02' = 28 + leap,  # adds 1 if leap year
                '03' = 31,
                '04' = 30,
                '05' = 31,
                '06' = 30,
                '07' = 31,
                '08' = 31,
                '09' = 30,
                '10' = 31,
                '11' = 30,
                '12' = 31))
} #Thanks to Jacob Amos. Nice tip - https://stackoverflow.com/questions/6243088/find-out-the-number-of-days-of-a-month-in-r?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa

#Convert month names (Jan or January) to number
MonthNum <- function(Month) {grep(paste("(?i)",substr(Month,1,3),sep=""),month.abb)}


#And we'll also apply a few of the conversions so that we don't have to run the functions regularly. A few of these conversions only have to be done once.
LatitudeRadians <- func.DMS(Latitude) #Lakes don't move much. :D


#Salinity converters
ECToTDS <- function(EC, ConversionFactor){
TDS <- EC*ConversionFactor
return (TDS)
}

PracToRefSalinity <- function (PracSalinity){
RefSalinity <- PracSalinity * 35.16504/35
return (RefSalinity)
}

PracToAbsSalinity <- function (PracSalinity){
}