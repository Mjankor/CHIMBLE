#Chimble Mk2
#Version 2.0
#Date 17/4/18
#Things to do
#1) Daily/Monthly check for Strat and GW data. 

#Initial Setup
require(compiler)
enableJIT(3)
library(Rcpp)
library(CHIMBLERcpp)
options(digits=2) # Just to keep things legible.
options(scipen=999) #This bit of code pushes R to display things in decimal rather than scientific format.
options(max.print=2000) #Just to make it easy not to be overwhelmed by large matrixes and dataframes.



#Find folders, preferences and setup CHIMBLE
source ("Modules/CHIMBLE Setup.r")

#Load up some preferences. Particularly those that define time and program function.
LoadVariables <- function (File){
	for (i in 1:nrow(File)){
		if (grepl("^[A-z]{1,}$", File[i,"Value"])){ #check if alpha or numeric variable
		assign(paste(rownames(File)[i]), File[i,"Value"], envir = .GlobalEnv)
		} else {
		 assign(paste(rownames(File)[i]), as.numeric(File[i,"Value"]), envir = .GlobalEnv)
		}
	}
}

LoadAllVariables  <- function(){
	Files <- c(unique(VarList$File))
	for (i in 1:length(Files)){
		LoadVariables (get(Files[i]))
	}

#add a few checks on formatting.
	if (Timestep > 1){
	message ("Timestep format incorrect. Using reciprocal (", 1/Timestep, ") instead.")
	Timestep <<- 1/Timestep
	}
	
	message("Variables Loaded")
}

LoadAllVariables()

#Load some useful conversion routines
source ("Modules/Conversions.r")

#Load historical data
source ("Modules/Load historical data.r")

#Define Hypsographic Curves. This module is loaded for historical purposes but is not required. It has been replaced with CHIMBLERcpp.
source ("Modules/Hypsographic Calcs.r")
#ShowHypsoCurves() #uncomment to show the graphs

#And load general conversion routines
source("Modules/Conversions.r")

#Set time to 1 to start the model run.
#Time is initially 1. This refers to the first line of the Met data frame.
#We can iterate through the Met data frame as long as Time <= Metrows
Time <- 1

#Chimble V1 stores isotope, groundwater data, and stratification data in the Metfile. This makes it easy to keep tabs on what aligns with what. However, it also means that every modification of these parameters requires a new MetFile. By splitting these datasets out, it makes it clearer as to which parameters are being modified for a differing scenario.

#One thing we need to do to speed things up is coerce our dataframes into matrixes. This will make it easy to take them into C++, and should make the R stuff a little bit faster.
#Any of the time series files should use Year, Month, Day. (These are Met, Strat, GW) We could simply split them off and use their row number to process the files. Then, once the main mass balance is finished, add the year, month, day data back.
#The above is a future possibility. Currently the Met file (and others) keeps the date fields, with the day field having the days in the month. 
#All we should need to do is convert the Met data into a matrix, ditching the date fields.
source ("Modules/Read timeseries.r")

#Next step is to setup the running files. These are matrixes used to track reservoir and flux values.
#Tables from Chimble V1 are Lake, Clean (daily and monthly), Flux, Flux18O, FluxD, ISOD, ISO18O, 
#We need Lake, Clean (daily and monthly), Flux?, Chemistry (~70 fields), Isotopes (~6 fields). 
#Problem - When we switch to thermal lake model, we may end up with more reservoirs.
#Chemistry should be recorded monthly. #Isotopes as per timestep.
#We should use matrixes instead. Then we can use them in C++ without problems.

#Structure
#Lake: Timestep	Year	Month	Day	depth	volume	areas	RESsl	RESdl	RESss	RESds	RESsp	RESin
#Clean:  Timestep	Year	Month	Day	depth	volume	areas	RESsl	RESdl	RESss	RESds	RESsp	RESin	d18O	dD
#Flux: 	Timestep	Year	Month	Day	Fr	Fp	Fsm	Fe	Fsos	Fdos	Fdlm	Fssi	Fsse	Fssd	Fdse	Fdsd	Fro	Fin	Fsf	Fgw	Fof	E	PET
#Chemistry: 
#Isotopes: Timestep	Year	Month	Day	18O_RESsl	D_RESsl	18O_RESdl	D_RESdl	18O_RESss	D_RESss	18O_RESds	D_RESds	18O_RESsp	D_RESsp	18O_RESin	18O_GW	D_GW	D_RESin	18O_E	D_E	18O_PET	D_PET
source("Modules/Create running matrixes.r")

#Load the salinity parameters
source("Modules/Salinity_prep.r")

#Activate groundwater module.
if (Prefs["GroundwaterModule","Value"] == "Yes" || Prefs["GroundwaterModule","Value"] == "yes"){
source ("Modules/Groundwater Module Initialisation.r")
} else {
source ("Modules/Groundwater Placeholder.r")
}

#Activate stratification moddule
source ("Modules/Stratification Initialisation.r")
#This is run regardless as the file is needed as an argument for the mass balance. However, it's not used if Stratification is turned off.

#Load the event manager. This runs once, to flag all the events. Then, once the main loop is running, if Met$Events is TRUE, then MassBalance will parse the EventScenarios for any new events.
source ("Modules/Event manager.r")
EventFlagger("Loud")
FindEvents()

#Next, Apply modifiers to the various datasets. eg: PMod, WSMod, etc.
#Event driven modifiers are applied to the basic metfile, but any initial modifiers set in prefs aren't included (as they aren't an event). Therefore this function has to be applied again here, after the event manager.
source("Modules/Apply Modifiers.r")

#Load hydrological functions. This module is loaded for historical purposes but is not required. It has been replaced with CHIMBLERcpp.
source ("Modules/Hydrology functions.r")

#Load isotope calcs. This module is loaded for historical purposes but is not required. It has been replaced with CHIMBLERcpp.
source("Modules/Isotope functions.r")

#Load some plotting functions
source ("Modules/MapPlotter.r")

#And a mass balance check function

#load the mass balance
library(CHIMBLERcpp)

if (GroundwaterModule == "No"){
	if (nrow(Met)/20000 < 60){
	message ("There are ",nrow(Met), " integrations. This will take around ", round (nrow(Met)/20000, digits = 2), " seconds to run.")
	}else {
	message ("There are ",nrow(Met), " integrations. This will take around ", round (nrow(Met)/1200000, digits = 0), " minutes to run.")
	}
} else {
	if (nrow(Met)/9 < 60){
	message ("There are ",nrow(Met), " integrations. This will take around ", round (nrow(Met)/9, digits = 2), " seconds to run.")
	#Minor problem. This is dependent upon how much of a "jump" the groundwater has to cover each timestep. For monthly it's more like nrow(Met)/2.
	}else {
	message ("There are ",nrow(Met), " integrations. This will take around ", round (nrow(Met)/540, digits = 0), " minutes to run.")
	}
}
readline(prompt = "To run the lake model please hit enter.")

MassBalance (KcSS, ALBearth, ALBlake, CA, Cin, Csr, KcDS, KcLsum, KcLwin, LatitudeRadians, RunoffRatio, PWF, Timestep, AWCds, AWCss, Time, Flux, Isotopes, Lake, Met, Lakevolumes, Datarate, Interpolationtype, GroundwaterModule, Topogrid, GWRechargegrid, GWPumpinggrid, GWKgrid, GWSgrid, GWbasegrid, GWheadgrid, GWboundarygrid, GWlakesedimentgrid, Lakecondgrid, GWETgrid, Classgrid, ETexdepth, Cellxdim, Cellydim, Accuracy, Maxiterations, GWcatchmentgrid, GWlakekgrid, LakeDepths, Chemistry, DirectExternalFlows, AtmosphericShift, Turbulence, Feedback, FluxD, Flux18O, GWStats, GWNumbers, ScenarioEventsList, ScenarioEvents, ChemConcentrations, ChemParams, NeutralDragCoeff, SWExtCoef, ThermoclineTargetThickness, ThermoclineMaximum, StratificationModule, Stratification, GWLakeID, Theta)

#Do some chacks and generate clean file.
source ("Modules/Clean output file generator.r")
source ("Modules/MassBal Checks.r")
source ("Modules/Plotting routines.r")
Dates <- SyncDates()
LoadGraphs()
message ("Use 'LoadGraphs()', 'func.plan()' and 'profplot(0,0)' to display graphs")

#Final step. Clear the scenario path so that scenarios are only loaded when specified.
if (exists("ScenarioPath")){
rm(ScenarioPath)
}

