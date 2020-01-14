#Now for the parameters file. We want to create a data table with the same number of rows as the Metfile +1. The plus one is because most of these values are the result of the observations (meteoroogical or otherwise) for the month prior. eg: if we start with a lake level of 2m for January 2017, and have 0.1m of rainfall observed for Jan 2017, then we need 2 months in the lake and other files, one for the 2m of January, and a second for the 2.1m of February 2017. 

#This function (currently overkill) does that. It's also there to provide a framework later for filling in the rest of the blank data in the final row of the running matrixes.
RunningMatrixLength <- function (type){
if (type == "Full"){
	RowLength <- length(Met[,1])+1
} else if (type == "Clean") {
FirstYearMonths <- 13 - Met[1,"Month"]
FinalYearMonths <- Met[nrow(Met),"Month"]
RowLength <- FirstYearMonths + FinalYearMonths + 12*(length(unique(MetDf$Year))-2)+1
}
return (RowLength)
}

#This function just generates a final date and timestep for the end of the running matrixes. The running matrixes are always 1 row longer than the Met timeseries file. This is easily explainable if we consider just 1 data row. The lake begins with a volume 100 in timestep 1. Timestep 1 has 10 units of rainfall, which means that the lake begins Timestep 2 with a volume of 110. The other option would be to avoid filling in the parameters to timestep 1, and treat them separately (as Timestep 0). However, I find it easier to follow if the initial parameters and timeseries data all starts at the same time.
FinalLine <- function(){
if (Datarate == "Daily"){
MetFinalDate <- as.Date(paste(Met[nrow(Met),"Year"],Met[nrow(Met),"Month"],Met[nrow(Met),"Day"],sep = "-"))
ClosingDate <- MetFinalDate+1
ClosingDateList <- c(as.numeric(row.names(Met)[nrow(Met)]) + (as.numeric(row.names(Met)[nrow(Met)]) -  as.numeric(row.names(Met)[nrow(Met)-1])))
ClosingDateList <- c(ClosingDateList, as.numeric(format(as.Date(ClosingDate, format="%d/%m/%Y"),"%Y")))
ClosingDateList <- c(ClosingDateList, as.numeric(format(as.Date(ClosingDate, format="%d/%m/%Y"),"%m")))
ClosingDateList <- c(ClosingDateList, as.numeric(format(as.Date(ClosingDate, format="%d/%m/%Y"),"%d")))
return (ClosingDateList)
}
if (Datarate == "Monthly"){
MetFinalDate <- as.Date(paste(Met[nrow(Met),"Year"],Met[nrow(Met),"Month"],Met[nrow(Met),"Day"],sep = "-"))
ClosingDate <- MetFinalDate+ceiling(30*Timestep)
ClosingDateList <- c(ceiling(as.numeric(row.names(Met)[nrow(Met)]) + (as.numeric(row.names(Met)[nrow(Met)]) -  as.numeric(row.names(Met)[nrow(Met)-1]))))
ClosingDateList <- c(ClosingDateList, as.numeric(format(as.Date(ClosingDate, format="%d/%m/%Y"),"%Y")))
ClosingDateList <- c(ClosingDateList, as.numeric(format(as.Date(ClosingDate, format="%d/%m/%Y"),"%m")))
ClosingDateList <- c(ClosingDateList, as.numeric(format(as.Date(ClosingDate, format="%d/%m/%Y"),"%d")))
return (ClosingDateList)
}
}

ReservoirList <- list("RESsl","RESdl","RESss","RESds","RESsp","RESin","GW")

Lake <- matrix(data = as.numeric(0),nrow = RunningMatrixLength("Full"), ncol = length(ReservoirList) + 7, dimnames = list(c(),c("Timestep","Year","Month","Day","Depth","Volume","Area", ReservoirList)))
message ("Lake running matrix created")
#need to deal with odd months. Such a function will be useful in the future too.
#length(unique(Met$Year)) gives us all unique years. However, we also need to know whether the start or end years are a full year (unlikely).

if (Datarate == "Monthly") {
Clean <- matrix(data = as.numeric(0),nrow = RunningMatrixLength("Clean"), ncol = length(ReservoirList) + 11, dimnames = list(c(),c("Timestep","Year","Month","Day","Depth","Volume","Area",ReservoirList, "Strat_Depth", "d18O_sampledepth", "dD_sampledepth","Cl_sampledepth")))
}else {
Clean <- matrix(data = as.numeric(0),nrow = RunningMatrixLength("Clean"), ncol = length(ReservoirList) + 11, dimnames = list(c(),c("Timestep","Year","Month","Day","Depth","Volume","Area",ReservoirList, "Strat_Depth", "d18O_sampledepth", "dD_sampledepth","Cl_sampledepth")))
}
#Both of these are the same currently. The second (daily datarate) may be modified in the future as required.

#Now let's populate the first line. 
#Lake$CA[1] <- ("CA")
#Lake$AWCss[1] <- ("AWCss")
#Lake$AWCds[1] <- ("AWCds")
#Lake$Cin[1] <- ("Cin")
#Lake$Csr[1] <- ("Csr")
#Lake$RESsl[1] <- ("RESsl")
#Lake$RESdl[1] <- ("RESdl")
#Lake$RESss[1] <- ("RESss")
#Lake$RESds[1] <- ("RESds")
#Lake$RESin[1] <- ("RESss")
#Lake$RESsp[1] <- ("RESsp")
#Lake$Volume[1] <- ("Volume")
#Lake$Area[1] <- ("NA")
#Lake$Depth[1] <- ("NA")
#I've removed this temporarily to avoid having to convert strings to characters. We'll put the header data in just before we write out the file. 

#Now for the flux data frame. We'll make this the same structure as the Lake file, just in case we have to do any integration to fluxes over months. Later we may drop it to a single row that gets overwritten each time. We don't need to have a header line, but to keep the structure easy, we'll keep it.

Flux <- matrix(data = as.numeric(0),nrow = RunningMatrixLength("Full"), ncol = 29, dimnames = list(c(),c("Timestep","Year","Month","Day","FpC","FpL","Fsm","Fe","Fsos","Fdos","Fdlm","Fssi","Fsse","Fssd","Fdse","Fdsd","Fro","Fin","Fsf","Fgw","FExt","Fof","GW_FHead", "GW_Ext", "GW_Recharge", "GW_ET", "E","PET", "GW_Err")))
message ("Flux running matrix created")
Flux18O <- matrix(data = as.numeric(0),nrow = RunningMatrixLength("Full"), ncol = 29, dimnames = list(c(),c("Timestep","Year","Month","Day","FpC","FpL","Fsm","Fe","Fsos","Fdos","Fdlm","Fssi","Fsse","Fssd","Fdse","Fdsd","Fro","Fin","Fsf","Fgw","FExt","Fof","GW_FHead", "GW_Ext", "GW_Recharge", "GW_ET", "E","PET", "GW_Err")))
FluxD <- matrix(data = as.numeric(0),nrow = RunningMatrixLength("Full"), ncol = 29, dimnames = list(c(),c("Timestep","Year","Month","Day","FpC","FpL","Fsm","Fe","Fsos","Fdos","Fdlm","Fssi","Fsse","Fssd","Fdse","Fdsd","Fro","Fin","Fsf","Fgw","FExt","Fof","GW_FHead", "GW_Ext", "GW_Recharge", "GW_ET", "E","PET", "GW_Err")))

#These two routines sort out the initial data (and column names) for isotopes and chemistry. 
IsotopeInitial <- function (){
ReservoirList <- c(ReservoirList, "E","PET")
counter <- 1
IsoPrep <- data.frame(Name = character(length(ReservoirList)*nrow(IsoData)), Value = numeric(length(ReservoirList)*nrow(IsoData)), stringsAsFactors = FALSE)
for (i in 1:length(ReservoirList)){
	for (j in 1:nrow(IsoData)){
	IsoPrep$Name[counter] <- paste(ReservoirList[i],row.names(IsoData)[j], sep = "_")
	counter <- counter + 1
	}
	}
	counter <- 1
	for (i in 1:ncol(IsoData)){
	for (j in 1:nrow(IsoData)){
		IsoPrep$Value[counter] <- IsoData[j,ReservoirList[[i]]]
	counter <- counter + 1
	}
	}
	return (IsoPrep)
	}

ChemistryInitial <- function (){
counter <- 1
ReservoirList <- c(ReservoirList, "P","Ext")
ChemPrep <- data.frame(Name = character(length(ReservoirList)*nrow(ChemData)), Value = numeric(length(ReservoirList)*nrow(ChemData)), stringsAsFactors = FALSE)
	for (j in 1:nrow(ChemData)){
		for (i in 1:length(ReservoirList)){
	ChemPrep$Name[counter] <- paste(ReservoirList[i],row.names(ChemData)[j], sep = "_")
	ChemPrep$Value[counter] <- ChemData[j,ReservoirList[[i]]]
	counter <- counter + 1
	}
	}
	return (ChemPrep)
	}

#Create isotope running file
Isotopes <- matrix(data = as.numeric(0),nrow = RunningMatrixLength("Full"), ncol = nrow(IsotopeInitial()) + 4, dimnames = list(c(),c("Timestep","Year","Month","Day",IsotopeInitial()[,"Name"])))
Isotopes[1,5:(ncol(Isotopes))] <- IsotopeInitial()[,"Value"]

message ("Isotopes running matrix created")

#Create chemistry running file
Chemistry <- matrix(data = as.numeric(0),nrow = RunningMatrixLength("Full"), ncol = nrow(ChemistryInitial()) + 4, dimnames = list(c(),c("Timestep","Year","Month","Day",ChemistryInitial()[,"Name"])))
Chemistry[1,5:ncol(Chemistry)] <- ChemistryInitial()[,"Value"]
message ("Chemistry running matrix created")

#Fill out timestep data for each matrix
Lake[1:nrow(Lake)-1,"Timestep"] <- as.numeric(rownames(Met))
Lake[1:nrow(Lake)-1,"Year"] <- as.numeric(Met[,"Year"])
Lake[1:nrow(Lake)-1,"Month"] <- as.numeric(Met[,"Month"])
Lake[1:nrow(Lake)-1,"Day"] <- as.numeric(Met[,"Day"])
Lake[nrow(Lake),c(1,2,3,4)] <- FinalLine()

Flux[1:nrow(Flux)-1,"Timestep"] <- as.numeric(rownames(Met))
Flux[1:nrow(Flux)-1,"Year"] <- as.numeric(Met[,"Year"])
Flux[1:nrow(Flux)-1,"Month"] <- as.numeric(Met[,"Month"])
Flux[1:nrow(Flux)-1,"Day"] <- as.numeric(Met[,"Day"])
Flux[nrow(Lake),c(1,2,3,4)] <- FinalLine()


Isotopes[1:nrow(Isotopes)-1,"Timestep"] <- as.numeric(rownames(Met))
Isotopes[1:nrow(Isotopes)-1,"Year"] <- as.numeric(Met[,"Year"])
Isotopes[1:nrow(Isotopes)-1,"Month"] <- as.numeric(Met[,"Month"])
Isotopes[1:nrow(Isotopes)-1,"Day"] <- as.numeric(Met[,"Day"])
Isotopes[nrow(Lake),c(1,2,3,4)] <- FinalLine()


Chemistry[1:nrow(Chemistry)-1,"Timestep"] <- as.numeric(rownames(Met))
Chemistry[1:nrow(Chemistry)-1,"Year"] <- as.numeric(Met[,"Year"])
Chemistry[1:nrow(Chemistry)-1,"Month"] <- as.numeric(Met[,"Month"])
Chemistry[1:nrow(Chemistry)-1,"Day"] <- as.numeric(Met[,"Day"])
Chemistry[nrow(Lake),c(1,2,3,4)] <- FinalLine()

#Define initial conditions for each dataset.
for (i in 1:length(ReservoirList)){
Lake[1,ReservoirList[[i]]] <- as.numeric(Params[ReservoirList[[i]],"Value"])
Lake[1,"Volume"] <- Lake[1,"RESsl"] + Lake[1,"RESdl"]
Lake[1,"Area"] <- func.Lakearea(Lake[1,"Volume"])
Lake[1,"Depth"] <- func.Lakedepth(Lake[1,"Volume"])
}

