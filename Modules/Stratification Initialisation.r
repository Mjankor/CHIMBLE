Stratification <- matrix(data = as.numeric(0), ncol = 15, nrow = (ceiling(max(Lakevolumes[,1])) - floor(min(Lakevolumes[,1])))/LayerThickness, dimnames = list(c(), c("LowerDepth", "UpperDepth", "Thickness", "PreviousTemp", "PreviousD", "Previous18O", "PreviousSalinity", "Temp", "D", "18O", "Salinity", "Density", "SpecHeat", "Diffusion", "Area")))
#need an additional check to make sure there's enough lake depth.

MaxLayers <- nrow(Stratification)
if (max(Lakevolumes[,1])%%1 == 0){
Stratification[1,"LowerDepth"] <- (max(Lakevolumes[,1]))-1
Stratification[1,"UpperDepth"] <- (max(Lakevolumes[,1])) #For When we get dynamic Layers
FirstFullLayerUpper <- (max(Lakevolumes[,1]))-1
} else {
Stratification[1,"LowerDepth"] <- floor(max(Lakevolumes[,1]))
Stratification[1,"UpperDepth"] <- (max(Lakevolumes[,1])) #For When we get dynamic Layers
FirstFullLayerUpper <- floor(max(Lakevolumes[,1]))
}

if (min(Lakevolumes[,1])%%1 == 0){
Stratification[MaxLayers,"LowerDepth"] <- (min(Lakevolumes[,1]))
Stratification[MaxLayers,"UpperDepth"] <- (min(Lakevolumes[,1]))+1
LastFullLayer <- (min(Lakevolumes[,1]))+1
} else{
Stratification[MaxLayers,"LowerDepth"] <- floor(min(Lakevolumes[,1]))
Stratification[MaxLayers,"UpperDepth"] <- ceiling(min(Lakevolumes[,1]))
LastFullLayer <- floor(min(Lakevolumes[,1]))+1
}

Stratification[1:(MaxLayers-1),"LowerDepth"] = FirstFullLayerUpper:LastFullLayer
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
	LakeVolume <- func_Lakevolume(Stratification[i,"UpperDepth"], Lakevolumes)
	Stratification[i,"Area"] <- func_Lakearea(LakeVolume, Lakevolumes)
	#also do temps
	if (Stratification[i,"UpperDepth"]+Stratification[i,"Thickness"]/2 <= Lake[1,"Depth"]){
		Stratification[i,"Temp"] = InitialLakeTemp 
	}
}


#Stratification[1:InitThermoclineLayer,"PreviousSalinity"] = 
# Stratification[InitThermoclineLayer:MaxLayers,"PreviousSalinity"] = 
