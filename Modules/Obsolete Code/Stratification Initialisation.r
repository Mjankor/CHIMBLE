Stratification <- matrix(data = as.numeric(0), ncol = 15, nrow = (ceiling(max(Lakevolumes[,1])) - floor(min(Lakevolumes[,1])))/LayerThickness, dimnames = list(c(), c("LowerDepth", "UpperDepth", "Thickness", "PreviousTemp", "PreviousD", "Previous18O", "PreviousSalinity", "Temp", "D", "18O", "Salinity", "Density", "SpecHeat", "Diffusion", "Area")))

MaxLayers <- nrow(Stratification)
Stratification[1,"LowerDepth"] <- floor(max(Lakevolumes[,1]))
#Stratification[1,"UpperDepth"] <- ceiling(max(Lakevolumes[,1]))
Stratification[1,"UpperDepth"] <- (max(Lakevolumes[,1])) #For When we get dynamic Layers
Stratification[MaxLayers,"LowerDepth"] <- floor(min(Lakevolumes[,1]))
Stratification[MaxLayers,"UpperDepth"] <- ceiling(min(Lakevolumes[,1]))
#Stratification[MaxLayers,"UpperDepth"] <- (min(Lakevolumes[,1])) #For When we get dynamic Layers


Stratification[1:(MaxLayers-1),"LowerDepth"] = floor(max(Lakevolumes[,1])):ceiling(min(Lakevolumes[,1]))
Stratification[2:(MaxLayers-1),"UpperDepth"] = Stratification[2:(MaxLayers-1),"LowerDepth"] + LayerThickness

Stratification[,"Thickness"] = Stratification[,"UpperDepth"] - Stratification[,"LowerDepth"]

InitThermoclineDepth <- Lake[1,"Depth"] - Met[1,"SL"]
InitThermoclineLayer <- which(Stratification[,"LowerDepth"] == floor(InitThermoclineDepth))

Stratification[,"Temp"] = InitialLakeTemp
Stratification[1:InitThermoclineLayer,"D"] = Isotopes[1,"RESsl_D"]
Stratification[InitThermoclineLayer:MaxLayers,"D"] = Isotopes[1,"RESdl_D"]
Stratification[1:InitThermoclineLayer,"18O"] = Isotopes[1,"RESsl_18O"]
Stratification[InitThermoclineLayer:MaxLayers,"18O"] = Isotopes[1,"RESdl_18O"]

for (i in 1:nrow(Stratification)){
	indexupper <- which(Lakevolumes[,1] == Stratification[i,"UpperDepth"])
	#indexlower <- which(Lakevolumes[,1] == Stratification[i,"LowerDepth"])
	#Stratification[i,"Area"] <- (Lakevolumes[indexlower,"Area"] + Lakevolumes[indexupper,"Area"])/2
	Stratification[i,"Area"] <- Lakevolumes[indexupper,"Area"]
	#also do temps
	if (Stratification[i,"UpperDepth"]+Stratification[i,"Thickness"]/2 <= Lake[1,"Depth"]){
		Stratification[i,"Temp"] = InitialLakeTemp 
	}
}


#Stratification[1:InitThermoclineLayer,"PreviousSalinity"] = 
# Stratification[InitThermoclineLayer:MaxLayers,"PreviousSalinity"] = 
