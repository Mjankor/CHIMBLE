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


#Event planner, places a flag in the metfile for each event.
if (Datarate == "Daily"){
	for (i in 1:nrow(ScenarioEvents)){
	EventDate <- which(Met[,"Year"] == ScenarioEvents$Year[i] & Met[,"Month"] == ScenarioEvents$Month[i] & Met[,"Day"] == ScenarioEvents$Day[i])
	Met[EventDate,"Events"] <- 1
	message ("Event: ", ScenarioEvents$Name[i], " set to change to ",ScenarioEvents$Value[i]," on " ,Met[EventDate,"Day"],"/",Met[EventDate,"Month"],"/",Met[EventDate,"Year"])
	}
	}

if (Datarate == "Monthly"){
	for (i in 1:nrow(ScenarioEvents)){
	PossibleDays <- which(Met[,"Year"] == ScenarioEvents$Year[i] & Met[,"Month"] == ScenarioEvents$Month[i])
	PossibleDays <- c(Met[PossibleDays[1]:PossibleDays[length(PossibleDays)],"Day"])
	NearestDay <- stepfun(PossibleDays, 0:length(PossibleDays)) #Cool little sort function
	NearestDay <- PossibleDays[NearestDay(ScenarioEvents$Day[i])]
	ScenarioEvents$Day[i] <- NearestDay
	#Two options here. We can just replace the scenario events day with the nearest one, or we can write a search routine. We'll use the former for simplicity and speed.
	EventDate <- which(Met[,"Year"] == ScenarioEvents$Year[i] & Met[,"Month"] == ScenarioEvents$Month[i] & Met[,"Day"] == ScenarioEvents$Day[i])
	Met[EventDate,"Events"] <- 1
	message ("Event: ", ScenarioEvents$Name[i], " set to change to ",ScenarioEvents$Value[i]," on " ,Met[EventDate,"Day"],"/",Met[EventDate,"Month"],"/",Met[EventDate,"Year"])
	}
}
}

#Event compare can be called using VarName and VarValue from the ScenarioEvents file. 
EventCompare <- function (VarName, VarValue){
VariableToggle <- 0
FileToggle <- 0
FilesChanged <- c()
RerunTimeSeries <- 0
#If name matches variable then variable is updated.
if (VarName %in% VarList$Variable){
	#Then we need to find it in the list
	VariableToggle <- 1
	VarRow <- which(VarList$Variable == VarName)
	TempFile<-get(VarList$File[VarRow])
	PreValue <- TempFile[VarName,"Value"]
	TempFile[VarName,"Value"] <- VarValue
	assign(VarList$File[VarRow], TempFile, envir= .GlobalEnv)
	message ("Event: ",VarName," from file ",VarList$File[VarRow]," updated with value ",VarValue,". Previous value was ",PreValue)
	#This will update the prefs table, so make sure that any variables that are commonly called without referring to the prefs table, eg: Timestep, DataRate are updated too.
	#To do this I'll update the variable definition section in CHIMBLE Setup.R to be a function.
}

if (VarName %in% row.names(Paths)){
	FileToggle <- 1
	FilesChanged <- c(FilesChanged,VarName)
	TempFile<-Paths
	OldFile <- TempFile[VarName,"Value"]
	TempFile[VarName,"Filename"] <- VarValue	
	assign("Paths", TempFile, envir= .GlobalEnv)
	message ("Event: ",VarName, " updated with file ",VarValue,", replacing file ",OldFile)
	#paths file is updated. Now we need to process the file into something useful.
}

if (VariableToggle == 1){
	#Update all active variables. 
	LoadAllVariables()
	}
	
if (FileToggle == 1){
	#Reload file.
	message (FilesChanged)
	for (i in 1:length(FilesChanged)){
	#update the file. Need to assign this outside the function, so we'll use assign again.
UpdatedFile <- LoadFile(paste(Paths[FilesChanged[i],"Subfolder"],"/",Paths[FilesChanged[i],"Filename"],sep=""),Paths[FilesChanged[i],"Type"])
assign (VarName, UpdatedFile, envir= .GlobalEnv)
	#Update the new file paths. We can also call these direct from the Paths dataframe, but having the short FilenamePath variables is quite convenient
	RebuildPaths()
	if (Paths[FilesChanged[i],"Type"] == "Timeseries"){
	RerunTimeSeries <- 1
	}
	#If name matches filename then file reloaded and rerun
	#Timeseries data may be a problem as we'll need to rerun the timeseries module.
	#Note: This does not include changing data directories. I don't forsee any examples where that would be a good, or useful idea.
	}
		if (RerunTimeSeries == 1){
			message ("Rebuilding timeseries datafiles")
		source ("Modules/Read timeseries.r", local = FALSE)
		}
	}
}

#This function parses the ScenarioEvents list (which should be prepped with dates the same as the Met file from running this source). For all events with the same date, it runs the EventCompare() routine. 

FindEvents <- function (line){
	message ("Time is ",line)
	EventDate <- which(ScenarioEvents$Year == Met[line,"Year"] & ScenarioEvents$Month == Met[line,"Month"] & ScenarioEvents$Day == Met[line,"Day"])
	for (i in 1:length(EventDate)){
	EventLine <- EventDate[i]
	EventCompare(ScenarioEvents$Name[EventLine], ScenarioEvents$Value[EventLine])
	}
}


#If called, then parse the ScenarioEvents list
#For daily check day, month, year
#Check the date against all rows (Which?). If date matches a row, then check variable name. 
#For monthly as above but use same calculation as Event manager part 1
#If if matches something in the 'paths', then reload the file. Also, perhaps, reform the metfile.
#If it doesn't, then just redo the variable in question. write message.
