#read meteorological data
#Adding a check here to make sure that the Timestep is setup the right way.

message ("Loading meteorological timeseries.")
if (MetDataPath == ""){
	message("\n*****************\nThe second file required is meteorological values for the time period of interest")
	message("This file requires a header line, followed by a tab separated rows with the following values: \nMonth\nPrecipitation (mm)\nTemperature (ºC)\nRelative Humidity (%)\nSolar radiation (MJ m-2d-1)\nWind speed (m/s)\nd18O in precipitation (0/00)\ndDeuterium in precipitation (0/00)\nStratification depth (m)\n Lake air temperature offset (ºC)\n")
	readline(prompt = "To select the meteorological file, hit <Enter>")
	MetDataPath <- tryCatch(file.choose(), error = function(e) "")
	message("The file ",basename(MetDataPath) ," has been loaded")

	#Now to assign some variables from these files.
	#One point to remember is that growing a dataframe in R is very slow. It would be better to precalculate the size of the data frame, and leaving it empty. We can do this simply by looking at the length of the Metfile dataframe.
	#Just an aside, I'll be using data frames rather than lists to allow for the use of other data types such as strings used for notes, etc in the future.
} else {#end if
#	message ("*** Reading meteorological data file specified in Preferences")
}#end else
MetDf <- LoadFile(MetDataPath, "Timeseries")
if (Prefs["Datarate","Value"] == "Monthly"){
message ("*** ", nrow(MetDf), " months identified in the Meteorological data.")
} else {
message ("*** ",nrow(MetDf), " days identified in the Meteorological data.")
}

#This is a check to avoid any situations where RH = 100%, as this condition breaks the CG model.
if (length (MetDf[,"RH"][MetDf[,"RH"] > 95])){
message ("*** ", length (MetDf[,"RH"][MetDf[,"RH"] > 95]), " days/months had relative humidity values above 95%. Values set to 95% to avoid breaking the isotope fractionation equations.")
}
MetDf[,"RH"][MetDf[,"RH"] > 95] = 95

#read stratification data
if (Prefs["StratificationModule","Value"] == "No"){
message ("Loading stratification timeseries.")
StratData <- LoadFile(StratDataPath, "Timeseries")
} else {
message ("Loading stratification module.")
}

#read GWdata
if (Prefs["GroundwaterModule","Value"] == "No"){
message ("Loading groundwater timeseries.")
GWData <- LoadFile(GWDataPath, "Timeseries")
} else {
message ("Loading groundwater module.")
}

#Combine timeseries data into the Met File

message ("Combining timeseries data into single dataframe (Met dataframe).")
if (Prefs["StratificationModule","Value"] == "No"){
	if (nrow(StratData) == nrow (MetDf)){
	MetDf$SL <- StratData$SL
	} else {
	message ("Stratification and Metfiles are of different length.")
	}
} else {
	MetDf$SL <- 0
	MetDf$Toff <-0
}

if (Prefs["GroundwaterModule","Value"] == "No"){
	if (nrow(GWData) == nrow (MetDf)){
	MetDf$Groundwater <- GWData$Groundwater
	MetDf$d18Ogw <- GWData$d18Ogw
	MetDf$dDgw <- GWData$dDgw
	} else {
	message ("Groundwater and Metfiles are of different length.")
	}
} else {
	MetDf$Groundwater <- 0
	MetDf$d18Ogw <- 0
	MetDf$dDgw <- 0
}

if (Prefs["StratificationModule","Value"] == "No"){
source("Modules/Stratification Smoother.r")
	if (Timestep < 1){
	SLsmooth <- func.SLsmooth()$SL #Replace the stepped SL with the smoothed one.
	} else {
	SLsmooth <- MetDf$SL
	}#SLsmooth takes a while to run. This prevents it running if not required (Timestep = 1).
}

#This next step modifies the timeseries data to account for multiple iterations during a month/day.
MetDf <- MetDf[rep(1:nrow(MetDf),each=(1/Timestep)),] #This duplicates the metrows by the timestep factor, so that we don't require multiple indexes. 
Metrows <- nrow(MetDf)
if (Prefs["StratificationModule","Value"] == "No"){
MetDf$SL <- SLsmooth 
} else {
MetDf$SL <- 0
}
#The above is messy (with regard to SL calculation and should be tidied)

#This next line adds a single column to the Met file for incoming events
MetDf$Events <- 0

#Convert the Met Dataframe to a Matrix for faster running of the main loop. 
Met <- MetDf
#Months to integers (If not already numeric)
if (!is.numeric(Met[1,"Month"])){
Met$Month <- sapply((substr (Met$Month,1,3)),function(x) grep(paste("(?i)",x,sep=""),month.abb)) 
}
#And days to the total number of days recorded for that month

if (Prefs["Datarate","Value"] == "Monthly"){
Met$Day <- mapply (DaysInMonth,Met$Month,Met$Year)
#Note: This currently represents the number of days in each month. We could easily modify it with the Timestep variable to represent days per timestep. 
#Which we now do below, as it makes it easier to manage events.
#Method is Ceiling(Total month days * Timestep * decimal of rowname*10 + (half the timestep difference))
Met$Day <- ceiling((Met$Day)*Timestep * (as.numeric(rownames(Met)) %% 1)*10 + (((Met$Day)*Timestep)/2))
}
#and this gives us a day in the approximate middle of each timestep

Met <- data.matrix (Met, rownames.force=TRUE)