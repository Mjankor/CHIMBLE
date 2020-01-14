#What do we need to do here.
#Event manager function part 1
#Parse events list
#If daily then match day, month, year with Met and assign a 1 to that Met["Row","Events"]
#If monthly then match day to nearest timestep, month, year with Met and assign a 1 to that Met["Row","Events"]
#May be worth having an error catch here in case someone uses numerical.


#Convert the dates to numerical month
if (exists("ScenarioPath")){
if (is.character(ScenarioEvents$Month)){
ScenarioEvents$Month <- mapply (MonthNum,ScenarioEvents$Month)
}
}

#Event planner, places a flag in the metfile for each event.
EventFlagger <- function(mute){
if (nrow(ScenarioEvents) > 0){
	if (Datarate == "Daily"){
		for (i in 1:nrow(ScenarioEvents)){
		EventDate <- which(Met[,"Year"] == ScenarioEvents$Year[i] & Met[,"Month"] == ScenarioEvents$Month[i] & Met[,"Day"] == ScenarioEvents$Day[i])
		Met[EventDate,"Events"] <<- 1
		if (mute != "Mute"){
		message ("Event: ", ScenarioEvents$Name[i], " set to change to ",ScenarioEvents$Value[i]," on " ,Met[EventDate,"Day"],"/",Met[EventDate,"Month"],"/",Met[EventDate,"Year"])
		}
		}
		}

	if (Datarate == "Monthly"){
		for (i in 1:nrow(ScenarioEvents)){
		PossibleDays <- which(Met[,"Year"] == ScenarioEvents$Year[i] & Met[,"Month"] == ScenarioEvents$Month[i])
		PossibleDays <- c(Met[PossibleDays[1]:PossibleDays[length(PossibleDays)],"Day"])
		NearestDay <- stepfun(PossibleDays, 0:length(PossibleDays)) #Cool little sort function
		NearestDay <- PossibleDays[NearestDay(ScenarioEvents$Day[i])]
		ScenarioEvents$Day[i] <<- NearestDay
		#Two options here. We can just replace the scenario events day with the nearest one, or we can write a search routine. We'll use the former for simplicity and speed.
		EventDate <- which(Met[,"Year"] == ScenarioEvents$Year[i] & Met[,"Month"] == ScenarioEvents$Month[i] & Met[,"Day"] == ScenarioEvents$Day[i])
		Met[EventDate,"Events"] <<- 1
			if (mute != "Mute"){
		message ("Event: ", ScenarioEvents$Name[i], " set to change to ",ScenarioEvents$Value[i]," on " ,Met[EventDate,"Day"],"/",Met[EventDate,"Month"],"/",Met[EventDate,"Year"])
		}
		}
	}
}
}

#Event compare can be called using VarName and VarValue from the ScenarioEvents file. 
#This has been significantly changed from the previous version. Previous version updated variables in all datafiles. Path names, etc. As this new version will just be loading new variables and files into the mass balance these changes will be treated as non-permanent and only for that model run.
EventCompare <- function (VarName, VarValue, i){
if (nrow(ScenarioEvents) > 0){
VariableToggle <- 0
#FileToggle <- 0
FilesChanged <- c()
RerunTimeSeries <- 0
#If name matches variable then variable is updated.
	if (VarName %in% VarList$Variable){
		#Then we need to find it in the list
		#VariableToggle <- 1
		VarRow <- which(VarList$Variable == VarName)
		TempFile<-get(VarList$File[VarRow])
		PreValue <- TempFile[VarName,"Value"]
		#TempFile[VarName,"Value"] <- VarValue
		#assign(VarList$File[VarRow], TempFile, envir= .GlobalEnv)
		#Change ^^^ to place into eventslist
		if (grepl("^[A-z]{1,}$", VarValue)){
		ScenarioEventsList[[i]] <<- VarValue
		} else {
		ScenarioEventsList[[i]] <<- as.numeric(VarValue)
		}
		#message ("Event: ",VarName," from file ",VarList$File[VarRow]," updated with value ",VarValue,". Previous value was ",PreValue)
	}

	if (VarName %in% row.names(Paths)){
		#FileToggle <- 1
		#FilesChanged <- c(FilesChanged,VarName) No longer need a list of files. Instead we process them one at a time.
		#TempFile<-Paths
		OldFile <- Paths[VarName,"Filename"]
			if (Paths[VarName,"Type"] == "Timeseries"){
			RerunTimeSeries <- 1
			MetBackup <- Met #Backup the Met File
			TempFile<-Paths
			TempFile[VarName,"Filename"] <- VarValue
			assign("Paths", TempFile, envir= .GlobalEnv)
			RebuildPaths()
			}
		#TempFile[VarName,"Filename"] <- VarValue	
		#^^^ This updates the filename in the paths object. Not sure if still relevant if we switch to paired lists and RCPP code.
		#assign("Paths", TempFile, envir= .GlobalEnv)
		#message ("Event: ",VarName, " updated with file ",VarValue,", replacing file ",OldFile)
		#paths file is no longer updated. Now we need to process the file into something useful.
		UpdatedFile <- LoadFile(paste(Paths[VarName,"Subfolder"],"/",Paths[VarName,"Filename"],sep=""),Paths[VarName,"Type"])
			if (RerunTimeSeries == 1){
				message ("Rebuilding timeseries datafiles")
				source ("Modules/Read timeseries.r", local = FALSE)
				#Note. this won't have modifiers applied after the reload.
				
				EventFlagger("Mute")
				ScenarioEventsList[[i]] <<- Met
				#Met <<- MetBackup #removed so that the Met file is continually updated.
			} else {
			ScenarioEventsList[[i]] <<- UpdatedFile
			}
	}
	
	#And special treatment for the modifiers. These are applied to the entire metfile on load, so to apply them as an event, requires the creation of a new metfile, with the values applied.
	if (VarName == "PMod"){
	MetBackup <- Met
	Met[,"P"] <- PrimaryMetBackup[,"P"] * VarValue
	ScenarioEventsList[[i]] <<- Met
	#Met <<- MetBackup
	}
	if (VarName == "TaMod"){
	MetBackup <- Met
	Met[,"Ta"] <- PrimaryMetBackup[,"Ta"] + VarValue
	#And validate lake temperatures to make sure it doesn't freeze
	source("Modules/Lake temperature validation.r")
	ScenarioEventsList[[i]] <<- Met
	#Met <<- MetBackup
	}	
	if (VarName == "TwMod"){
	MetBackup <- Met
	Met[,"Toff"] <- PrimaryMetBackup[,"Toff"] + VarValue
	#And validate lake temperatures to make sure it doesn't freeze
	source("Modules/Lake temperature validation.r")
	ScenarioEventsList[[i]] <<- Met
	#Met <<- MetBackup
	}	
	if (VarName == "WSMod"){
	MetBackup <- Met
	Met[,"WS"] <- PrimaryMetBackup[,"WS"] * VarValue
	ScenarioEventsList[[i]] <<- Met
	#Met <<- MetBackup
	}	
	if (VarName == "RHMod"){
	MetBackup <- Met
	Met[,"RH"] <- PrimaryMetBackup[,"RH"] * VarValue
	#And a few fixes to resolve data that's out of range. eg: humidity can't be <0 or >100%
	Met[,"RH"][Met[,"RH"] > 100] = 100
	Met[,"RH"][Met[,"RH"] < 0] = 0
	ScenarioEventsList[[i]] <<- Met
	#Met <<- MetBackup
	}	
	if (VarName == "RsMod"){
	MetBackup <- Met
	Met[,"Rs"] <- PrimaryMetBackup[,"Rs"] * VarValue
	ScenarioEventsList[[i]] <<- Met
	#Met <<- MetBackup
	}
	if (VarName == "GWMod"){
	MetBackup <- Met
	Met[,"Groundwater"] <- PrimaryMetBackup[,"Groundwater"] + VarValue
	ScenarioEventsList[[i]] <<- Met
	#Met <<- MetBackup
	}
	if (VarName == "Stratmod"){
	MetBackup <- Met
	Met[,"SL"] <- PrimaryMetBackup[,"SL"] * VarValue
	ScenarioEventsList[[i]] <<- Met
	#Met <<- MetBackup
	}
	#Update the Met file with the revised values
	Met <<- Met

}
}
#if (VariableToggle == 1){
	#Update all active variables. 
#	LoadAllVariables()
#	}
	
#if (FileToggle == 1){
	#Reload file.
#	message (FilesChanged)
#	for (i in 1:length(FilesChanged)){
	#update the file. Need to assign this outside the function, so we'll use assign again.
#UpdatedFile <- LoadFile(paste(Paths[FilesChanged[i],"Subfolder"],"/",Paths[FilesChanged[i],"Filename"],sep=""),Paths[FilesChanged[i],"Type"])
#assign (VarName, UpdatedFile, envir= .GlobalEnv)

	#Update the new file paths. We can also call these direct from the Paths dataframe, but having the short FilenamePath variables is quite convenient
#	RebuildPaths()
	#The above is probably no longer needed, Change the assign line to instead add the file to a list.
#	if (Paths[FilesChanged[i],"Type"] == "Timeseries"){
#	RerunTimeSeries <- 1
#	}
	#If name matches filename then file reloaded and rerun
	#Timeseries data may be a problem as we'll need to rerun the timeseries module.
	#Note: This does not include changing data directories. I don't forsee any examples where that would be a good, or useful idea.
	#Need to update this to preserve the original Met File. We should do this in the calling function, rather than here, so it only has to be done once, rather than every time a timeseries file is changed.
#	}

#	}
#}

#This function parses the ScenarioEvents list (which should be prepped with dates the same as the Met file from running this source). For all events with the same date, it runs the EventCompare() routine. 

FindEvents <- function (){
	#backup the Met file before setting up the events.
	PrimaryMetBackup <<- Met #Backup the Met File
	for (i in 1:nrow(ScenarioEvents)){
	EventCompare(ScenarioEvents$Name[i], ScenarioEvents$Value[i], i)
	}
	#restore the original Met file
	Met <<- PrimaryMetBackup
}

FindEventsOrig <- function (line){
	message ("Time is ",line)
	EventDate <- which(ScenarioEvents$Year == Met[line,"Year"] & ScenarioEvents$Month == Met[line,"Month"] & ScenarioEvents$Day == Met[line,"Day"])
	for (i in 1:length(EventDate)){
	EventLine <- EventDate[i]
	EventCompare(ScenarioEvents$Name[EventLine], ScenarioEvents$Value[EventLine], i)
	}
}
#If called, then parse the ScenarioEvents list
#For daily check day, month, year
#Check the date against all rows (Which?). If date matches a row, then check variable name. 
#For monthly as above but use same calculation as Event manager part 1
#If if matches something in the 'paths', then reload the file. Also, perhaps, reform the metfile.
#If it doesn't, then just redo the variable in question. write message.
