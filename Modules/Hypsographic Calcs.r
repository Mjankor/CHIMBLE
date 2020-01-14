
message("Loading hypsographic data")
#This first section is used to convert a table of lake data (Height, Volumes, Area) to a function.
#The file chosen should be tab separated text, with headers, and columns of Depth(height), Volume and Area.

#Determine interpolation type
Interpolationtype <- HypsoPrefs["Interpolationtype","Value"]

if (HypsoDataPath == ""){
message("*****************\nThis program is designed to convert lake data (Height/Depth, Volumes and Areas to a Loess function.\nIt requires a tab separated text file with three columns for height, volume and area. \nPS: Make sure your lake data extends well above the height of the current lake level (preferably to include the entire catchment or to overflow level. The highest lake level/maximum volume should be at the top of the file.)\n")
readline(prompt = "To select your lake data file, please hit <Enter>")
#LakeAVH <- function(){
HypsoDataPath <- tryCatch(file.choose(), error = function(e) "")
} else { # end if
#message ("*****************\n*** Reading lake volume file specified in Preferences")
}#end else
Lakevolumes <- LoadFile(HypsoDataPath,"Timeseries")
Lakevolumes <- as.matrix(Lakevolumes) # change to matrix so that we can use c++
message("The file ",basename(HypsoDataPath) ," has been loaded. CHIMBLE will now fit ", if (HypsoPrefs["Interpolationtype","Value"] == "Loess") "a Loess function" else "line segments", " to the hypsographic data.")
message ("*** To display the hypsographic curves, use the 'ShowHypsoCurves()' function")

options(digits=7) # Just to keep things legible.
options(scipen=999) #This bit of code pushes R to display things in decimal rather than scientific format.

#********************************* Interpolation **************************************

func.Lakearea <- function (volume){
if (Interpolationtype == "Loess"){
Area <- predict(LakeVA, volume)
return (Area)
} else if (Interpolationtype == "Segments"){
Mindex <-  which(Lakevolumes[,"Volume"] == min(Lakevolumes[,"Volume"][(Lakevolumes[,"Volume"] > volume)]))
Slope <- (Lakevolumes[Mindex,"Area"] - Lakevolumes[Mindex+1,"Area"])/(Lakevolumes[Mindex,"Volume"] - Lakevolumes[Mindex+1,"Volume"])
Intercept <- Lakevolumes[Mindex,"Area"] - Lakevolumes[Mindex,"Volume"]*Slope
Area <- volume*Slope + Intercept
return (Area)
} else {
message ("Interpolation type unrecognised (Lakearea Function)")
}
} #end func

func.Lakedepth <- function (volume){
if (Interpolationtype == "Loess"){
Depth <- predict(LakeVH, volume)
return (Depth)
} else if (Interpolationtype == "Segments"){
Mindex <- which(Lakevolumes[,"Volume"] == min(Lakevolumes[,"Volume"][(Lakevolumes[,"Volume"] > volume)]))
Slope <- (Lakevolumes[Mindex,"Depth"]  - Lakevolumes[Mindex+1,"Depth"])/(Lakevolumes[Mindex,"Volume"] - Lakevolumes[Mindex+1,"Volume"])
Intercept <- Lakevolumes[Mindex,"Depth"] - Lakevolumes[Mindex,"Volume"]*Slope
Depth <- volume*Slope + Intercept
return (Depth)
} else {
message ("Interpolation type unrecognised (Lakedepth Function)")
}
} #end func

func.Lakevolume <- function (depth){
if (Interpolationtype == "Loess"){
Volume <- predict(LakeHV, depth)
return (Volume)
} else if (Interpolationtype == "Segments"){
Mindex <- Mindex <- which(Lakevolumes[,"Depth"] == min(Lakevolumes[,"Depth"][(Lakevolumes[,"Depth"] > depth)]))
Slope <- (Lakevolumes[Mindex,"Volume"] - Lakevolumes[Mindex+1,"Volume"])/(Lakevolumes[Mindex,"Depth"] - Lakevolumes[Mindex +1,"Depth"])
Intercept <- Lakevolumes[Mindex,"Volume"] - Lakevolumes[Mindex,"Depth"]*Slope
Volume <- depth*Slope + Intercept
return (Volume)
} else {
message ("Interpolation type unrecognised (Lakevolume Function)")
}
} #end func

#********************************* Loess Interpolation **********************************
#These three lines are just to enable them as global variables
LakeVA <- 0
LakeVH <- 0
LakeHV <- 0
if (Interpolationtype == "Loess"){
LakeVA <- loess (Lakevolumes[,3] ~ Lakevolumes[,2], span = Fitspan, surface = "direct")
message ("\nVolume to area best fit completed with a loess function. Standard error = ",summary(LakeVA)[5])

LakeVH <- loess (Lakevolumes[,1] ~ Lakevolumes[,2], span = Fitspan, surface = "direct")
message ("Volume to lake depth best fit completed with a loess function. Standard error = ",summary(LakeVH)[5])

#This function is required for the SVC calculation. 
LakeHV <- loess (Lakevolumes[,2] ~ Lakevolumes[,1], span = Fitspan, surface = "direct")
}
#********************************* Interpolation Plots **********************************

ShowHypsoCurves <- function(){
	if (Prefs["HeadlessMode","Value"] == "No"){
	quartz()
	if (Interpolationtype == "Loess"){
	par(mfrow=c(2,2))
	} else {
	par(mfrow=c(2,1))
	}


	#The following section is commented out for headless mode
	plot (Lakevolumes[,2], Lakevolumes[,1], main="Volumes & Lake Depth", sub="", xlab="Volume", ylab="Depth")
	if (Interpolationtype == "Loess"){
	lines (Lakevolumes[,2], (predict(LakeVH)), col="red")
	plot (Lakevolumes[,2], resid(LakeVH), main="Volumes & Lake Depth Residuals", sub="", xlab="Volume", ylab="Depth")
	} else {
	lines (Lakevolumes[,2], Lakevolumes[,1], col="red")
	}

	plot(Lakevolumes[,2], Lakevolumes[,3], main="Volumes & Surface Area", sub="", xlab="Volume", ylab="Surface Area")
	if (Interpolationtype == "Loess"){
	lines (Lakevolumes[,2], (predict(LakeVA)), col="red")
	plot (Lakevolumes[,2], resid(LakeVA), main="Volumes & Surface Area Residuals", sub="", xlab="Volume", ylab="Surface Area")
	} else {
	lines (Lakevolumes[,2], Lakevolumes[,3], col="red")
	}
	}
}