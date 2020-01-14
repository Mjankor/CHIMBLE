require(compiler)
enableJIT(3)
options(digits=7)
options(scipen=999)
#library(rbenchmark)
library(Rcpp)
#sourceCpp("Modules/ChimbleGW.cpp") Removed. Now part of CHIMBLERcpp package
#sourceCpp("Assignments.cpp")
#sourceCpp("Vectorcheck.cpp")
#Read in preferences. Might as well use the same structure as CHIMBLE.
#GWPrefs <- read.table("BM Groundwater 50M/Preferences.txt",header=TRUE, sep="\t", row.names = 1, stringsAsFactors = FALSE)
#^^^Moved to CHIMBLE setup.r
#Prefs <- read.table("GWTest/Preferences.txt",header=TRUE, sep="\t", row.names = 1, stringsAsFactors = FALSE)

#************************** Read Preferences *************************************
#Things to do.
#Double check files and preferences against those used in Chimble. Remove any unused.
#Reassign preferences to new locations as per Chimble 2
#File reading location can be removed from here and left in CHIMBLE Setup.r
#May potentially move here, if it turns out to be annoying, but having NA for values should be fine.

#Read in grids. Each grid will simply be a tab separated file denoting a square celled grid. This means that they can be generated in Excel.


K1 <- as.numeric(GWPrefs ["K1", "Value"])
S1 <- as.numeric(GWPrefs ["S1", "Value"])
K2 <- as.numeric(GWPrefs ["K2", "Value"])
S2 <- as.numeric(GWPrefs ["S2", "Value"])
K3 <- as.numeric(GWPrefs ["K3", "Value"])
S3 <- as.numeric(GWPrefs ["S3", "Value"])
K4 <- as.numeric(GWPrefs ["K4", "Value"])
S4 <- as.numeric(GWPrefs ["S4", "Value"])
K5 <- as.numeric(GWPrefs ["K5", "Value"])
S5 <- as.numeric(GWPrefs ["S5", "Value"])
K6 <- as.numeric(GWPrefs ["K6", "Value"])
S6 <- as.numeric(GWPrefs ["S6", "Value"])
K7 <- as.numeric(GWPrefs ["K7", "Value"])
S7 <- as.numeric(GWPrefs ["S8", "Value"])
K8 <- as.numeric(GWPrefs ["K8", "Value"])
S8 <- as.numeric(GWPrefs ["S8", "Value"])
K9 <- as.numeric(GWPrefs ["K9", "Value"])
S9 <- as.numeric(GWPrefs ["S9", "Value"])

#GWTimestep <- as.numeric(GWPrefs ["Timestep", "Value"]) Now using CHIMBLE timestep
#GWRuntime <- as.numeric(GWPrefs ["Runtime", "Value"])
Cellxdim <- as.numeric(GWPrefs ["Cellxdim", "Value"])
Cellydim <- as.numeric(GWPrefs ["Cellydim", "Value"])
Accuracy <- as.numeric(GWPrefs ["Accuracy", "Value"])
Maxiterations <- as.numeric(GWPrefs ["Maxiterations", "Value"])


#*********************** Convert to dataframes/Matrixes ********************************

#convert the grids to dataframes
LakeDepths <- data.matrix(LoadFile(GWLakesPath,"Static"))
Topogrid <- data.matrix(LoadFile(TopogridPath,"Grid"))
GWPumpinggrid <- data.matrix(LoadFile(GWPumpinggridPath,"Grid"))
GWgeogrid <- data.matrix(LoadFile(GWgeogridPath,"Grid"))
GWcatchmentgrid <- data.matrix(LoadFile(GWcatchmentgridPath,"Grid"))
GWbasegrid <- data.matrix(LoadFile(GWbasegridPath,"Grid"))
GWboundarygrid <- data.matrix(LoadFile(GWboundarygridPath,"Grid"))
GWheadgrid <- data.matrix(LoadFile(GWheadgridPath,"Grid"))
GWlakesedimentgrid <- data.matrix(LoadFile(GWlakesedimentgridPath,"Grid"))
GWlakekgrid <- data.matrix(LoadFile(GWlakekgridPath,"Grid"))
# and generate a classification grid too.
Classgrid <- matrix(data = 0, nrow (Topogrid), ncol (Topogrid)) 
Lakecondgrid <- matrix(data = 0, nrow (Topogrid), ncol (Topogrid)) 
#and backup the boundary grid, as it will be modified regularly.
GWboundarybackup <- GWboundarygrid

#We should also combine the K and S values into a nice little matrix.
Geomatrix <- matrix(c(K1,S1,K2,S2,K3,S3,K4,S4,K5,S5,K6,S6,K7,S7,K8,S8,K9,S9), ncol = 2, byrow = TRUE)


#Make an ETgrid to apply Evapotranspiration to each cell.
#GWETgrid <- as.matrix(read.table(GWETgrid,header = FALSE))
GWETgrid <- matrix(data = 0, nrow = nrow(Topogrid), ncol = ncol(Topogrid))

#And a recharge grid to for recharge from the catchment
GWRechargegrid <- matrix(data = 0, nrow = nrow(Topogrid), ncol = ncol(Topogrid))

message ("Groundwater Module loaded.")

#*********************** Define lake depths ********************************

func.Lakedepths <- function(){#and the same for lake depths
CurrentLakeDepth <- (func.Lakedepth(Lake[Time,"Volume"]))
	for (i in 1:nrow(LakeDepths)){
		if (CurrentLakeDepth + LakeDepths[i,"LakeDelta"] > LakeDepths[i,"MinDepth"]){
		LakeDepths[i,"CurrentDepth"] <- CurrentLakeDepth + LakeDepths[i,"LakeDelta"]
		} else {
		LakeDepths[i,"CurrentDepth"] <- MinDepth
		}
	}
	return (LakeDepths)
}


#*********************** Groundwater stats matrix ********************************
# These give the C functions something to write to. In particular, GWNumbers is the volumes of flows in and out for each GW model run. 
GWNumbers <- matrix(data = 0, nrow = 9, ncol = 5, dimnames = list(c(),c("LakeID", "SeepageIn", "SeepageOut", "SurfaceFlows", "NetFlows")))
GWStats <- list("MassInit" = 0, "MassFinal" = 0, "MassError(%)" = 0, "FixedHeadFlux" = 0, "Recharge Flux" = 0, "ExtFlows(Abstractions)" = 0, "LakeSeepage_In" = 0, "LakeSeepage_Out" = 0, "SurfaceFlows" = 0, "NonCatchmentOutflows" = 0, "ET" = 0, "Timestep" = 0, "Runtime" = 0, "Recharge" = 0, "Conductivity" = 0)


#*********************** Create K and S Matrixes ********************************

#and from that, generate Kgrid and Sgrid matrixes for conductivity and specific storage. It's easier to do that in R, rather than in C and will only have to be done once.
GWKgrid <- matrix(data = 0, nrow (Topogrid), ncol (Topogrid)) 
for (i in 1:nrow(Geomatrix)){
GWKgrid[GWgeogrid == i] <- Geomatrix[i,1]
} #end for

GWSgrid <- matrix(data = 0, nrow (Topogrid), ncol (Topogrid)) 
for (i in 1:nrow(Geomatrix)){
GWSgrid[GWgeogrid == i] <- Geomatrix[i,2]
} #end for

message("The groundwater grids have been loaded.")

#now we should calculate the grid size and perhaps area covered.
GWgridx <- ncol(Topogrid)
GWgridy <- nrow(Topogrid)
GWarea <- (GWgridx*Cellxdim) * (GWgridy*Cellydim)

message("Total size = x:", (GWgridx*Cellxdim), "m, y:",(GWgridy*Cellydim), "m, for a total area of ",GWarea, "m^2")

message("Catchment area: ", (Cellxdim*Cellydim * length(GWcatchmentgrid[GWcatchmentgrid == 1])))


#*********************** Calculate initial volume ********************************
Lake[1,"GW"] <- sum((GWheadgrid - GWbasegrid)*GWSgrid*Cellxdim *Cellydim)

#**************************** CHIMBLE specific parameters*** *************************
#Setup things for Chimble. Parameters that have to be aquired from CHIMBLE rather than from the prefs file.
#Timesteps moved to GW C++ file.
#GWRuntime <- func.monthlength(Time)
#GWRuntime <- Met[Time,"Day"]
#If Met day field is modified, change this to func.monthlength(Met[1,"Month"])
#GWTimestep <- GWRuntime*Timestep #need to update this section to account for multiple timesteps. 
#*********************** Create a "depth" grid (head - base) *************************
#No longer needed. All sorted in C++ file
#func.updatedepth <- function (){
#GWdepth <- GWheadgrid - GWbasegrid
#}
#GWdepth <- func.updatedepth()

#************************** Check grid dimensions **************************************
#12D sometimes generates grids that are a column and row short. This just checks that all the groundwater grids are the same dimensions.
grids <- c("Topogrid", "GWRechargegrid", "GWETgrid", "GWPumpinggrid", "GWgeogrid", "GWcatchmentgrid", "GWboundarygrid", "GWbasegrid", "GWheadgrid", "GWlakesedimentgrid", "GWlakekgrid")
for (i in 1:length(grids)){
	if (nrow(Topogrid) - nrow(get(grids[i])) != 0){
		message ("Grid ", grids[i], " does not have the same number of rows as the Topogrid. Topogrid rows = ", nrow(Topogrid), ", ", grids[1], " rows = ", nrow(get(grids[i])))
	}
	if (ncol(Topogrid) - ncol(get(grids[i])) != 0){
		message ("Grid ", grids[i], " does not have the same number of columns as the Topogrid. Topogrid columns = ", nrow(Topogrid), ", ", grids[1], " columns = ", nrow(get(grids[i])))
	}
}

#************************** Check aquifer thickness **************************************
#This is a check to ensure that GWbasegrid lies below Topogrid - GWlakesedimentgrid. If not, then we get negative numbers, and the base layer might end up above the head height. 
#If there are areas where the base lies above, then the base is set as Topogrid - GWlakesedimentgrid - 1 (so 1 m thickness).
if (min(Topogrid - GWlakesedimentgrid - GWbasegrid) < 0){
message ("There is a problem with GWbasegrid. It has elevations above the top of the aquifer, resulting in a negative aquifer thickness. Changing elevations to sit 1m below top of aquifer.")
GWbasegrid[(Topogrid -GWlakesedimentgrid - GWbasegrid) < 0] <- (Topogrid[(Topogrid -GWlakesedimentgrid - GWbasegrid) < 0] - 1)
}
#************************** Determine Catchment **************************************

#This is done using Catchment divide.R It generates a binary grid of 1 and 0 with 0 being exterior and 1 being the catchment. In the case of dual lake system, the second catchment is defined with value 2.

#************************** Cell recharge **************************************

#This function adds the recharge to each grid cell
#Acquire directly from func.soil Fdsd.
func.GWrecharge <- function(){ 
#need to use the gridded area estimate, or we'll end up with a discrepancy between Fro & Fdsd and the recharge volume.
GWRechargegrid <- matrix(data = func.Soil(Time)$Fdsd/((sum(func.class() == 11))*Cellxdim*Cellydim*GWTimestep), nrow (Topogrid), ncol (Topogrid))
##substitute .03 for the mm recharge value from Chimble.
return (GWRechargegrid)
}

func.GWETgrid <- function(){ 
#need to use the gridded area estimate, or we'll end up with a discrepancy between Fro & Fdsd and the recharge volume.
GWETgrid <- matrix(data = func.PET(Time)/func.monthlength(Time), nrow (Topogrid), ncol (Topogrid))
##substitute .03 for the mm recharge value from Chimble.
return (GWETgrid)
}	

#************************** Determine Area of lake *******************************


#First, a function to define the area that is lake. Need two grids, topo and classification.
#This needs to be updated to also constrain the lake area inside the catchment. Otherwise, lowpoints outside the catchment, which may be critical as a limit on groundwater, will be modelled as though at lake height. - done

#Major update. Lake classes 0-9. Catchments 10-19. Outside catchments 20.
func.class <- function (){
Classgrid[GWcatchmentgrid == 1] <- Topogrid[GWcatchmentgrid == 1] - (func.Lakedepth(Lake[Time,"Volume"]))
#hack for gnotuk
Classgrid[GWcatchmentgrid == 2] <- Topogrid[GWcatchmentgrid == 2] - (func.Lakedepth(Lake[Time,"Volume"]) - 40)
#define lake (Has to be catchment first, or catchment class overrides lake)
Classgrid[Classgrid > 0 & GWcatchmentgrid == 1] <- 11
Classgrid[Classgrid <= 0 & GWcatchmentgrid == 1] <- 1
#define catchment
#hack for gnotuk
Classgrid[Classgrid > 0 & GWcatchmentgrid == 2] <- 12 #catchment
Classgrid[Classgrid <= 0 & GWcatchmentgrid == 2] <- 2 #lake
#Define outside of catchment
Classgrid[is.na(Classgrid)] <- 20
return (Classgrid)
}

#Commented out lines - original class function. No Gnotuk option.
#func.class <- function (){
#Classgrid <- Topogrid - Lakedepth
#Classgrid[Classgrid <= 0] <- 0
#Classgrid[Classgrid > 0] <- 1
#Classgrid[GWcatchmentgrid == 0] <- 2
##Convert to Boolean
#return (Classgrid)
#}

#****************************** Define lake conductance *********************************
#This can be done once so it's easiest to do it in R at the start of the routine.
func.Lakeconductance <- function (){
Lakecondgrid <- (GWlakekgrid * Cellxdim* Cellydim)/GWlakesedimentgrid
Lakecondgrid[!is.finite(Lakecondgrid)] <- 0
return (Lakecondgrid)
}


#************************ Set lake cell heads to lake heights ***************************

#this function will set the lake area to have a head height the same as the lake.
#not needed anymore as lakes are now head dependent.
#func.GWlake <- function (){
#Lakedepth <- (func.Lakedepth(Lake$Volume[Time]))
#GWheadgrid[func.class() == 1] <- Lakedepth
##hack for gnotuk
#Gnotukdelta <- 40
#Gnotukdepth <- if (Lakedepth - Gnotukdelta > 86) Lakedepth -Gnotukdelta else 86
#GWheadgrid[func.class() == 2] <- Gnotukdepth
##Hack to allow for Gnotuk
##GWheadgrid[func.class() == 2]
#return (GWheadgrid)
#}
 # Setup fixed heads for lake
#Now we should have updated head values based on recharge and an updated lake level. 

#boundary conditions. Tricky as the lake forms a boundary condition. However, it should also be self regulating to some degree. 
#Let's describe a few boundary conditions first.
#0 is no boundary condition, default cell.
#1 defines a dirichlet condition. (fixed head)
#2 defines a noflow boundary (Noflow also considered as extent of grid). Probably not needed.
#Specified flux conditions (Neumann) are implemented already through recharge. We can expand on that to include wells if required.
#Fixed head/Dirichlet condition. Applicable to the lake. Requires updating of the Boundary grid with each new lake level.
#Boundaries around the model are difficult to define. I want them to rise and fall as part of the regional watertable, but they still need to form the boundary. One tricky solution is to define their head heights based on the slope of the lower aquifer surface (thereby maintaining the regional groundwater flow), and with a head height defined as some sort of weighted average of nearby cells. We also need to quantify how much water balance changes due to groundwater movement at the boundaries.  
#What we'll do instead is simply define the area outside the model as no-flow, and rely on flow + P + surface seepage to constrain the water.

#************************ Define lake cells as dirichlet condition ***********************
#No longer needed as lake is now defined as a head dependent (Cauchy) seepage condition.
#func.boundary <- function (){
#GWboundarygrid <- GWboundarybackup #restore original boundary grid
#GWboundarygrid[func.class() == 0] <- 1 #update boundary grid with lake cells
##hack for gnotuk
#GWboundarygrid[func.class() == 3] <- 1
#return (GWboundarygrid)
#}

#********************************* Define lake depths**********************************

#*************************** Run R functions and C code *******************************
#Run R functions
GWheadgrid <- GWheadgrid * 1 #dodgy hack to stop Rcpp losing the address of GWheadgrid and not updating it.
#Now here we want the C function.
func.Lakeflow <- function (){
Lakeflow <- Headcalc(Topogrid, GWRechargegrid, GWPumpinggrid, GWKgrid, GWSgrid, GWbasegrid, GWheadgrid, GWboundarygrid, GWlakesedimentgrid, Lakecondgrid, GWETgrid, Classgrid, GWTimestep, GWRuntime, ETexdepth, Lakedepth, Cellxdim, Cellydim, Accuracy, Maxiterations)
#Now we need to strip out the flows to the lake
#[func.class() == 0]) Leaving out class, picks up overland flows
return (Lakeflow)
#Totalruntime <- Runtime + Totalruntime
}

#microbenchmark(func.Lakeflow(), times = 1000)
#Test 100,000 timesteps takes 18.3 seconds, using system.time
#By moving the entire loop system into C++, it now takes 1.76 seconds to do 100,000 iterations.
	

	#******************************************************************
	#===================== End Groundwater Setup ======================
	#******************************************************************
