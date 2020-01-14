#This module sets up the file paths, reads the preferences and parameters files and applies any variable modifications.

#Read data from scenario file (We use this routine to avoid the code and notes at the top of the file).
if (exists("ScenarioPath")){
message ("loading Scenario data file: ",substr(ScenarioFile,1,nchar(ScenarioFile)-2)," data.txt", sep = "")
ScenarioData <- read.table (paste(ScenarioDir,"/", substr(ScenarioFile,1,nchar(ScenarioFile)-2)," data.txt", sep = ""),header=TRUE, sep="\t", row.names = 1, stringsAsFactors = FALSE)
if (nrow(ScenarioData) < 1)
	message ("No data in Scenario Events file")
} else {
message ("No scenario data file loaded.")
ScenarioData <- data.frame(Name="",Value="",Notes="")
}

if (exists("ScenarioPath")){
message ("loading Scenario events file: ",substr(ScenarioFile,1,nchar(ScenarioFile)-2)," events.txt", sep = "")
ScenarioEvents <- read.table (paste(ScenarioDir,"/", substr(ScenarioFile,1,nchar(ScenarioFile)-2)," events.txt", sep = ""),header=TRUE, sep="\t", row.names = NULL, stringsAsFactors = FALSE)
ScenarioEventsList <- list()
	if (nrow(ScenarioEvents) < 1){
		message ("No data in Scenario Events file")
	}
} else {
	message ("No scenario events file loaded.")
	ScenarioEvents <- data.frame(Year=numeric(0),Month=numeric(0),Day=numeric(0),Name=character(0),Value=character(0),Notes=character(0))
	ScenarioEventsList <- list()
}

#Function to compare variables against the ones found in the scenario file. If a variable is active in the scenario file, then it is used in preference to the default. 
VarCompare <- function (defvar, default){
if (defvar %in% row.names(ScenarioData)){
var <- ScenarioData[defvar,"Value"]
message ("Scenario data used for: ",defvar, ", using: ", var)
} else {
var <- default
}
return (var)
}

#This function loads files based on a filepath.
LoadFile <- function(FileName, Type){
	if (Type == "Static"){
	NewFile <- read.table(FileName,header=TRUE,sep="\t", row.names = 1, stringsAsFactors = FALSE)
	} else if (Type == "Timeseries"){
	NewFile <- read.table(FileName,header=TRUE,sep="\t", stringsAsFactors = FALSE)	
	} else if (Type == "Grid"){
	NewFile <- read.table(FileName,header=FALSE,sep="\t", stringsAsFactors = FALSE)	
	}
	return(NewFile)
}

#Define details for all filepaths in the filepaths file
Paths <- read.table("Preferences/Filepaths.txt",header=TRUE, sep="\t", row.names = 1, stringsAsFactors = FALSE)
#Define all filepaths
RebuildPaths <- function(){
for (i in 1:nrow(Paths)){
assign(paste(rownames(Paths)[i],"Path", sep="") ,paste(Paths$Subfolder[i],"/",VarCompare(rownames(Paths)[i],Paths[i,"Filename"]),sep = "")
,envir= .GlobalEnv)
}
}

#And build the paths
RebuildPaths()

#Load preference files
Prefs <- LoadFile(ChimblePrefsPath, "Static")
GWPrefs <- LoadFile(GWPrefsPath, "Static")
ChemPrefs <- LoadFile(ChemPrefsPath, "Static")
HistPrefs <- LoadFile(HistPrefsPath, "Static")
StratPrefs <-LoadFile(StratPrefsPath, "Static")
MetPrefs <- LoadFile(MetPrefsPath, "Static")
HypsoPrefs <- LoadFile(HypsoPrefsPath, "Static")
IsoPrefs <- LoadFile(IsoPrefsPath, "Static")

message ("Input files loaded.")

#Compare variables in preferences to the scenario file and update any changes.
#We do this here to make it easy to chase down any preferences missed by the scenario VarCompare tester
VarCompareAll <- function (variabletable){
if (nrow(variabletable) > 0){
for (i in 1: nrow(variabletable)){
if (row.names(variabletable)[i] %in% row.names(ScenarioData)){
variabletable[i,"Value"] <- ScenarioData[row.names(variabletable)[i], "Value"]
message ("Scenario data used for: ",row.names(variabletable)[i], ", using: ", ScenarioData[row.names(variabletable)[i], "Value"])
}
}
}
return (variabletable)
}

Prefs <- VarCompareAll(Prefs)
GWPrefs <- VarCompareAll(GWPrefs)
ChemPrefs <- VarCompareAll(ChemPrefs)
HistPrefs <- VarCompareAll(HistPrefs)
StratPrefs <- VarCompareAll(StratPrefs)
MetPrefs <- VarCompareAll(MetPrefs)
IsoPrefs <- VarCompareAll(IsoPrefs)
HypsoPrefs <- VarCompareAll(HypsoPrefs)


#Load Parameters
Params <- LoadFile(ParamsPath, "Static")
#And the parameter file too.
Params <- VarCompareAll(Params)

#Load initial conditions
#Load isotope initial values
IsoData <- LoadFile(IsoDataPath, "Static")

#Load initial chemistry
ChemData <- LoadFile(ChemDataPath, "Static")

#Make a list of all variables

CreateVarList <- function (){
VarList <- ""
FilesList <- ""
VarList <- if (nrow(Prefs) > 0) c(row.names(Prefs)) else ""
FilesList <- if (nrow(Prefs) > 0) (replicate(nrow(Prefs),"Prefs")) else ""
VarList <- if (nrow(GWPrefs) > 0) c(VarList,(row.names(GWPrefs))) else VarList
FilesList <- if (nrow(GWPrefs) > 0) c(FilesList, (replicate(nrow(GWPrefs),"GWPrefs"))) else FilesList 
VarList <- if (nrow(ChemPrefs) > 0) c(VarList,(row.names(ChemPrefs))) else VarList
FilesList <- if (nrow(ChemPrefs) > 0) c(FilesList, (replicate(nrow(ChemPrefs),"ChemPrefs"))) else FilesList
VarList <- if (nrow(HistPrefs) > 0) c(VarList,(row.names(HistPrefs))) else VarList
FilesList <- if (nrow(HistPrefs) > 0) c(FilesList, (replicate(nrow(HistPrefs),"HistPrefs"))) else FilesList
VarList <- if (nrow(StratPrefs) > 0) c(VarList,(row.names(StratPrefs))) else VarList
FilesList <- if (nrow(StratPrefs) > 0) c(FilesList, (replicate(nrow(StratPrefs),"StratPrefs"))) else FilesList
VarList <- if (nrow(MetPrefs) > 0) c(VarList,(row.names(MetPrefs))) else VarList
FilesList <- if (nrow(MetPrefs) > 0) c(FilesList, (replicate(nrow(MetPrefs),"MetPrefs"))) else FilesList
VarList <- if (nrow(HypsoPrefs) > 0) c(VarList,(row.names(HypsoPrefs))) else VarList
FilesList <- if (nrow(HypsoPrefs) > 0) c(FilesList, (replicate(nrow(HypsoPrefs),"HypsoPrefs"))) else FilesList
VarList <- if (nrow(IsoPrefs) > 0) c(VarList,(row.names(IsoPrefs))) else VarList
FilesList <- if (nrow(IsoPrefs) > 0) c(FilesList, (replicate(nrow(IsoPrefs),"IsoPrefs"))) else FilesList
VarList <- if (nrow(Params) > 0) c(VarList,(row.names(Params))) else VarList
FilesList <- if (nrow(Params) > 0) c(FilesList, (replicate(nrow(Params),"Params"))) else FilesList
#And make them a dataframe
VarList <- data.frame(Variable = VarList, File = FilesList, stringsAsFactors = FALSE)
return (VarList)
}
VarList <- CreateVarList()

CheckNames <- function(){
	for (i in 1:nrow(VarList)){
		if (VarList$Variable[i] %in% row.names(Paths)){
		message ("Shared name for variable and file identifier for ",VarList$Variable[i],". Please fix this in the Filepaths or ",VarList$File[i]," preference file.")
		}
	}
	for (i in 1:nrow(VarList)){
		if (VarList$Variable[i] %in% VarList$Variable[-i]){
		message ("Shared name for variable: ",VarList$Variable[i],". Please fix this in the ", VarList$File[i]," preference file.")
		}
	}
}

CheckNames()



#Current project - Variables list
#Then use that to event manage

#And we'll do another for files.
#CreateFilesList <- function(){
##Name <- c(
#Folder <- c(
#Type <- c(
#FileList <- data.frame(Name = "", Folder = "", Type = "", stringsAsFactors = FALSE)
#return(FileList)
#}
#FileList <- CreateFilesList()