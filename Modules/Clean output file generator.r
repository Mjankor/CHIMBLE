#This little section just grabs some useful information from the various dataframes and compiles it into a nice simple dataframe for ease of reading.
#Will have to remove/replace CL- section depending on chemistry options.
func.Sample <- function(MetTime, Time){
if (Sampledepth >= Met[MetTime,"SL"]){
	SampleD <- Isotopes[Time, "RESdl_D"]
	Sample18O <- Isotopes[Time, "RESdl_18O"]
	SampleCl <- Chemistry[Time, "RESdl_Cl"]
} else {
	SampleD <- Isotopes[Time, "RESsl_D"]
	Sample18O <- Isotopes[Time, "RESsl_18O"]
	SampleCl <- Chemistry[Time, "RESsl_Cl"]
}
Sample <- list(SampleD, Sample18O, SampleCl)
return(Sample)
}

#Used for old chimble
#Lake$Timestep[nrow(Lake)] <- floor(as.numeric(Lake$Timestep[nrow(Lake)-1]))+1
##Flux$Timestep[nrow(Flux)] <- floor(as.numeric(Lake$Timestep[nrow(Flux)]))+1
#ISOD$Timestep[nrow(ISOD)] <- floor(as.numeric(Lake$Timestep[nrow(ISOD)-1]))+1
#ISOD$Timestep[nrow(ISOD)] <- floor(as.numeric(Lake$Timestep[nrow(ISOD)-1]))+1

if (Datarate == "Monthly") {
for (Time in 1:nrow(Lake)){ #Parallelizable
	if (Time > nrow(Met)){
	MetTime <- Time-1 #Met is one row shorter than lake matrix. Need to double up on the last Met value. 
	}else {
	MetTime <- Time
	}
	if (as.numeric(Lake[Time,"Timestep"])%%1 == 0){
	temptime <- as.numeric(Lake[Time,"Timestep"])
	Clean[temptime,"Timestep"] <- as.numeric(Lake[Time,"Timestep"])
	Clean[temptime,"Year"] <- as.numeric(Met[MetTime,"Year"])
	Clean[temptime,"Month"] <- as.numeric((Met[MetTime,"Month"]))
	Clean[temptime,"Day"] <- as.numeric(Met[MetTime,"Day"])
	Clean[temptime,"RESsl"] <- Lake[Time,"RESsl"]
	Clean[temptime,"RESdl"] <- Lake[Time,"RESdl"]
	Clean[temptime,"RESss"] <- Lake[Time,"RESss"]
	Clean[temptime,"RESds"] <- Lake[Time,"RESds"]
	Clean[temptime,"RESin"] <- Lake[Time,"RESin"]
	Clean[temptime,"GW"] <- Lake[Time,"GW"]
	Clean[temptime,"Volume"] <- Lake[Time,"Volume"]
	Clean[temptime,"Area"] <- Lake[Time,"Area"]
	Clean[temptime,"Depth"] <- Lake[Time,"Depth"]
	Clean[temptime,"Strat_Depth"] <- Met[MetTime,"SL"]
	Clean[temptime,"dD_sampledepth"] <- as.numeric(func.Sample(MetTime, Time)[1])
	Clean[temptime,"d18O_sampledepth"] <- as.numeric(func.Sample(MetTime, Time)[2])
	Clean[temptime,"Cl_sampledepth"] <- as.numeric(func.Sample(MetTime, Time)[3])
	}
	}
} else {
Cleanindex <- 1
for (Time in 1:nrow(Lake)){ #Parallelizable
	if (Time > nrow(Met)){
	MetTime <- Time-1  #Met is one row shorter than lake matrix. Need to double up on the last Met value. 
	}else {
	MetTime <- Time
	}
	if (as.numeric(Lake[Time,"Day"]) == 1 || as.numeric(Lake[Time,"Day"]) == 0){
	Clean[Cleanindex,"Timestep"] <- as.numeric(Lake[Time,"Timestep"])
	Clean[Cleanindex,"Year"] <- Met[MetTime,"Year"]
	Clean[Cleanindex,"Month"] <- as.numeric((Met[MetTime,"Month"]))
	Clean[Cleanindex,"Day"] <- Met[MetTime,"Day"]
	Clean[Cleanindex,"RESsl"] <- Lake[Time,"RESsl"]
	Clean[Cleanindex,"RESdl"] <- Lake[Time,"RESdl"]
	Clean[Cleanindex,"RESss"] <- Lake[Time,"RESss"]
	Clean[Cleanindex,"RESds"] <- Lake[Time,"RESds"]
	Clean[Cleanindex,"RESin"] <- Lake[Time,"RESin"]
	Clean[Cleanindex,"GW"] <- Lake[Time,"GW"]
	Clean[Cleanindex,"Volume"] <- Lake[Time,"Volume"]
	Clean[Cleanindex,"Area"] <- Lake[Time,"Area"]
	Clean[Cleanindex,"Depth"] <- Lake[Time,"Depth"]
	Clean[Cleanindex,"Strat_Depth"] <- Met[MetTime,"SL"]
	Clean[Cleanindex,"dD_sampledepth"] <- as.numeric(func.Sample(MetTime, Time)[1])
	Clean[Cleanindex,"d18O_sampledepth"] <- as.numeric(func.Sample(MetTime, Time)[2])
	Clean[Cleanindex,"Cl_sampledepth"] <- as.numeric(func.Sample(MetTime, Time)[3])
	Cleanindex <- Cleanindex +1
	}
	}
}

message ("Model run complete\n****************")
message ("Final lake depth = ", round(Clean[nrow(Clean), "Depth"], digits=2))
