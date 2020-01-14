#We need a function to calculate the actual catchment area as the initial value in the Lake data frame is the total catchment which consists of lake surface + earth catchment. This will vary depending on lake level. This function should be used instead of Lake$CA[Time] for most equations. 
func.CAe <- function(Time){ #Checked
	CAe <- as.numeric(Params["CA","Value"]) - Lake[Time,"Area"]#}#end if
	return (CAe) #end else
}#end function

#Equation 10: To call this function Fsm <- func.Fsm(time) and you'll get the value for Fsm for the current month. 
#Monthly data requirement
func.Fsm <- function(Time){ #Checked
	if (Datarate == "Monthly") {
	snowmelt <- 0.021 #mm per month
	} else {
	snowmelt <- 0.00069 #mm per day
	}
	#This section has been replaced. Now, if Ta > -2 and RESsp > 0, then snowmelt is calculated. The issue may be that Stella may not be able to do negative fluxes and that therefore having a negative flux isn't an issue. In R it's a problem.
	if (Lake[Time,"RESsp"] > (snowmelt * (Met[Time,"Ta"] + 2) * func.CAe(Time) * Timestep) && Met[Time,"Ta"] > -2) {
		Fsm <- snowmelt * (Met[Time,"Ta"] + 2) * func.CAe(Time) * Timestep
		Fsm18O <- Fsm * Isotopes[Time,"RESsp_18O"]
		FsmD <- Fsm * Isotopes[Time,"RESsp_D"]
	} else if (Lake[Time,"RESsp"] < (snowmelt * (Met[Time,"Ta"] + 2) * func.CAe(Time) * Timestep) && Met[Time,"Ta"] > -2){
		Fsm <- Lake[Time,"RESsp"]
		Fsm18O <- Fsm * Isotopes[Time,"RESsp_18O"]
		FsmD <- Fsm * Isotopes[Time,"RESsp_D"]
	} else {
		Fsm <- 0
		Fsm18O <- 0
		FsmD <- 0
}
Fsm <- list(Fsm = Fsm,Fsm18O = Fsm18O, FsmD = FsmD)
return (Fsm)
}#end func

#groundwater function - groundwater is run directly into the RESin at the moment. Volume and isotopes come from meteorological data.
func.Fgw <- function(Time){ #checked
	Fgw <- Met[Time,"Groundwater"] * Timestep
	Fgw18O <- Fgw*Met[Time,"d18Ogw"]
	FgwD <- Fgw*Met[Time,"dDgw"]
	Fgw <- list(Fgw = Fgw,Fgw18O = Fgw18O, FgwD = FgwD)
return (Fgw)
}#end function
#This function can probably be updated to call the groundwater module if groundwater module is activated.


#equation 11-16
func.Soil <- function (Time){ #Checked
SSmax <- func.CAe(Time) * as.numeric(Params["AWCss","Value"])#Max volume of soil
DSmax <- func.CAe(Time) * as.numeric(Params["AWCds","Value"])
Influx <- as.numeric(func.Fr(Time)["Fr"]) + as.numeric(func.Fsm(Time)["Fsm"]) #this has to balance
SSflux <- as.numeric(Lake[Time,"RESss"]) + Influx - as.numeric(func.Evap(Time)["Fsse"]) #RESss + fluxes
SSexcess <- SSflux - SSmax #excess from RESss
DSflux <- as.numeric(Lake[Time,"RESds"]) - as.numeric(func.Evap(Time)["Fdse"])



#Easy solution for this is that the incoming flux is split up amongst all the outgoing fluxes with no fractionation. All we need to do for the isotope concentrations is break up the incoming flux (Fr+Fsm) by the percentages of each outgoing flux.
Influx18O <- as.numeric(func.Fr(Time)["Fr18O"]) + as.numeric(func.Fsm(Time)["Fsm18O"]) #this has to balance
InfluxD <- as.numeric(func.Fr(Time)["FrD"]) + as.numeric(func.Fsm(Time)["FsmD"]) #this has to balance
if (Influx == 0) {
	Fssi <- 0
	Fssi18O <- 0
	FssiD <- 0
	Fssd <- 0
	Fssd18O <- 0
	FssdD <- 0 
	Fdsd <- 0
	Fdsd18O <- 0
	FdsdD <- 0
	Fro <- 0
	Fro18O <- 0
	FroD <- 0
} else if (SSmax >= SSflux){#RESss can hold all flux in and out
	Fssi <- Influx
	Fssi18O <- Influx18O
	FssiD <- InfluxD
	Fssd <- 0
	Fssd18O <- 0
	FssdD <- 0 
	Fdsd <- 0
	Fdsd18O <- 0
	FdsdD <- 0
	Fro <- 0
	Fro18O <- 0
	FroD <- 0
	} else if (SSmax < SSflux && DSmax >= DSflux + SSexcess*RunoffRatio/100) {
	#RESss cannot hold all flux, so half is diverted to RESds, half to Fro. RESds is able to hold all additional flux.
	Fssi <- SSmax - as.numeric(Lake[Time,"RESss"]) + as.numeric(func.Evap(Time)["Fsse"])
	Fssd <- SSexcess*RunoffRatio
	Fdsd <- 0
	Fro <- SSexcess*(1-RunoffRatio)
	Fssi18O <- Influx18O * Fssi/(Fssi+Fssd+Fdsd+Fro)
	Fssd18O <- Influx18O * Fssd/(Fssi+Fssd+Fdsd+Fro)
	Fdsd18O <- Influx18O * Fdsd/(Fssi+Fssd+Fdsd+Fro)
	Fro18O <- Influx18O * Fro/(Fssi+Fssd+Fdsd+Fro)
	FssiD <- InfluxD * Fssi/(Fssi+Fssd+Fdsd+Fro)
	FssdD <- InfluxD * Fssd/(Fssi+Fssd+Fdsd+Fro)
	FdsdD <- InfluxD * Fdsd/(Fssi+Fssd+Fdsd+Fro)
	FroD <- InfluxD * Fro/(Fssi+Fssd+Fdsd+Fro)
	} else if ((SSmax < SSflux && DSmax < DSflux + SSexcess*RunoffRatio)){
	#RESss cannot hold all flux, so half is diverted to RESds, half to Fro. RESds is also unable to hold all additional flux. 
	Fssi <- SSmax - as.numeric(Lake[Time,"RESss"]) + as.numeric(func.Evap(Time)["Fsse"])
	Fssd <- DSmax - as.numeric(Lake[Time,"RESds"]) + as.numeric(func.Evap(Time)["Fdse"])
	Fdsd <- SSexcess*RunoffRatio - Fssd
	Fro <- SSexcess*(1-RunoffRatio)
	Fssi18O <- Influx18O * Fssi/(Fssi+Fssd+Fdsd+Fro)
	Fssd18O <- Influx18O * Fssd/(Fssi+Fssd+Fdsd+Fro)
	Fdsd18O <- Influx18O * Fdsd/(Fssi+Fssd+Fdsd+Fro)
	Fro18O <- Influx18O * Fro/(Fssi+Fssd+Fdsd+Fro)
	FssiD <- InfluxD * Fssi/(Fssi+Fssd+Fdsd+Fro)
	FssdD <- InfluxD * Fssd/(Fssi+Fssd+Fdsd+Fro)
	FdsdD <- InfluxD * Fdsd/(Fssi+Fssd+Fdsd+Fro)
	FroD <- InfluxD * Fro/(Fssi+Fssd+Fdsd+Fro)
	} else {
	message ("You missed an option in the func.Soil routine")
	}
	Soil <- list(Fssi = Fssi,Fssd = Fssd, Fdsd = Fdsd, Fro= Fro, Fssi18O = Fssi18O, FssiD = FssiD, Fssd18O = Fssd18O, FssdD = FssdD, Fdsd18O = Fdsd18O, FdsdD = FdsdD, Fro18O = Fro18O, FroD = FroD)
	return (Soil)
}

func.Evap <- function(Time){ #Checked
#equations 12 and 14. These need to add up to the PET. As RESss is drawn down and unable to fulfil the PET, then evap occurs from RESds.
#We need to estimate the actual evapotranspiration ET.  We'll do this based on Van Boxel's model, rather than Palmer. 1965 which is approximately what Steinman does.
#First we will use two coefficients to create a linear relationship between ET and Soil water content.
#Additional future work would be to include saturation, as well as field capacity for soil, so that the RESss can go above field capacity and be withdrawn at the full rate.

SSETC <- Lake[Time,"RESss"]/(func.CAe(Time) * as.numeric(Params["AWCss","Value"]))
DSETC <- Lake[Time,"RESds"]/(func.CAe(Time) * as.numeric(Params["AWCds","Value"])) # + Lake$AWCss[Time]

SSETa <- (func.CAe(Time) * func.PET(Time) * as.numeric(MetPrefs["KcSS","Value"]) * SSETC)
DSETa <- (func.CAe(Time) * func.PET(Time) * as.numeric(MetPrefs["KcDS","Value"]) * DSETC)
#These two coefficients will result in a two step evaporation rate with most evaporation coming from upper soil and a lower rate from the deeper soil. Of interest is that this method also removes the dependency of evaporation in stages (first layer 1, then layer 2), as the lower layer now simply represents the evapotranspiration from deep rooted plants, thereby simplifying the code. 

if (Lake[Time,"RESss"] > (SSETa * Timestep)) { 
#so if RESss is greater than PET, then all evap comes from RESss.
	Fsse <- SSETa * Timestep
	} else {
	Fsse <- Lake[Time,"RESss"]
	}#endif
Fsse18O <- Fsse * Isotopes[Time,"RESss_18O"]
FsseD <- Fsse *  Isotopes[Time,"RESss_D"]
	
if (Lake[Time,"RESds"] > (DSETa * Timestep)) { 
#so if RESds is greater than ETa, then all evap comes from RESds.
	Fdse <- DSETa * Timestep
	} else {
	Fdse <- Lake[Time,"RESds"]
	}#endif
Fdse18O <- Fdse*Isotopes[Time,"RESds_18O"]
FdseD <- Fdse*Isotopes[Time,"RESds_D"]
#as evapotranspiration is not supposed to result in fractionation (steinman 2010), we simply set the isotopic flux to be the same as the soil reservoir.

		Evap <- list(Fsse = Fsse, Fdse = Fdse, Fsse18O = Fsse18O, 
		FsseD=FsseD,Fdse18O=Fdse18O,FdseD=FdseD)
		#Evap <-	lapply(Evap, function(x) x*0.5) # fudge here for testing the effect if some water remains in soil.
		return (Evap)
}


#Function for Fr over the catchment (catchment area - lake surface area)
func.Fr <- function(Time){ #Checked
	if (Met[Time,"Ta"] > 0){
		Fr <- ((Met[Time,"P"]*0.001) * func.CAe(Time) * Timestep)
		Fr18O <- Fr * Met[Time,"d18Op"]
		FrD <- Fr * Met[Time,"dDp"]
	} else {
		Fr <- 0
		Fr18O <- 0
		FrD <- 0
	}
	Fr <- list (Fr=Fr, Fr18O = Fr18O, FrD=FrD)
	return (Fr)
}#end function

#Equation 9.
func.Fsf <- function(Time){ #Checked
	if (Met[Time,"Ta"] <= 0){
		Fsf <- ((Met[Time,"P"]*0.001) * func.CAe(Time) * Timestep)
		Fsf18O <- Fsf * Met[Time,"d18Op"]
		FsfD <- Fsf * Met[Time,"dDp"]
	} else {
		Fsf <- 0
		Fsf18O <- 0
		FsfD <- 0
	} #end else
	Fsf <- list (Fsf = Fsf, Fsf18O = Fsf18O, FsfD = FsfD)
	return (Fsf)
}#end func

#equation 17
#Monthly data requirement
#Cin has been adapted to be 12/365 times the rate if daily data is used. 
#half distance problem here. RESin will never drain.
func.Fin <- function(Time){ #Checked
	if (Lake[Time,"RESin"] > (func.CAe(Time) * 0.001)) { #This just gives us a value that assumes that RESin will drain until AWC in the runoff and deeper drainage drops to 0.001mm over the catchment. 
		Fin <- Lake[Time,"RESin"] * as.numeric(Params["Cin","Value"]) * Timestep
		Fin18O <- Fin * Isotopes[Time,"RESin_18O"]
		FinD <- Fin * Isotopes[Time,"RESds_D"]
	} else {
		Fin <- 0
		Fin18O <- 0
		FinD <-	0	
	}
	Fin <- list (Fin = Fin, Fin18O = Fin18O, FinD = FinD)
	return (Fin)
}#end function

#equation 18
###In S2010 equations Deep_lake_depth is derived from the stratification profiles. From that a volume is calculated (Deep_lake_volume_m3). This volume is then subtracted from the total lake volume to give SVC. 
###In our version, SVC is volume determined by Total volume - Volume(Lake Height - Stratification depth). When the depth is 0, SVC = 0. 
###If SVC > 0, then SVC - RESsl is used. If SVC > RESsl, then Fdlm is positive, water flows into the upper layer. If SVC < RESsl, then Fdlm is negative, and water moves into the lower layer. It should be possible to have the Fslm and Fdlm fluxes as a single positive/negative flux. 
###turns out that inverse predictions are quite tricky. Probably best just to create a new fit with height vs volume. This is found at the top of the file with the other loess functions.
###We cannot rely on the predicted LakeHV function when the lake level gets very low, as the slope can go negative. Therefore, we will assume that if the lake gets to within 0.5m of the minimum level, that the entire lake is mixed into RESsl.	

func.SVC <- function(Time){ #Checked
	Basemix <- ifelse (Interpolationtype == "Loess", 0.6, 0)
	if (Met[Time,"SL"] == 0){
		SVC <- 0
	} else if (func.Lakedepth(Lake[Time,"Volume"]) - min(Lakevolumes[,1]) - Met[Time,"SL"] <= Basemix) {
		SVC <- Lake[Time,"Volume"]
	} else {
		SVC <- (Lake[Time,"Volume"] - func.Lakevolume(Lake[Time,"Depth"] - Met[Time,"SL"]))
	}
	return (SVC)
}#end function

#func.Fslm <- function(Time) This flux should no longer be needed. It's function will be taken up by the Fdlm as a positive/negative function.
Fslm <- 0 

#Moves water from one Epilimnion to Hypolimnion and vice versa
#Half distance problem here if we use the timestep.
##SVCfunc=0 This step is moved into the main program loops, as it's required to determine whether the fluxes go into RESsl or RESsl
func.Fdlm <- function(Time){ #Checked
	Fdlm <- (func.SVC(Time) - Lake[Time,"RESsl"])# * Timestep
		if (Fdlm < 0){
		Fdlm18O <- Fdlm* Isotopes[Time,"RESsl_18O"]
		FdlmD <- Fdlm*Isotopes[Time,"RESsl_D"]
		} else {
		Fdlm18O <- Fdlm*Isotopes[Time,"RESdl_18O"]
		FdlmD <- Fdlm*Isotopes[Time,"RESdl_D"]
		}
		Fdlm <- list (Fdlm = Fdlm, Fdlm18O = Fdlm18O, FdlmD = FdlmD)
		return(Fdlm)
	#so if SVC is larger than RESsl, Fdlm is positive showing water moving from deep to shallow, else, negative
}#end function

#Equation 20 and 21. These two will have to be updated to include the updated outseepage model used in Steinman 2012. 
func.Fsos <- function (time){ #Checked
	Fsos <- as.numeric(Params["Csr","Value"]) * Lake[Time,"RESsl"] * Timestep
	Fsos18O <- Fsos * Isotopes[Time,"RESsl_18O"]
	FsosD <- Fsos * Isotopes[Time,"RESsl_D"] 
	Fsos <- list(Fsos = Fsos, Fsos18O = Fsos18O, FsosD = FsosD)
	return (Fsos)}
	
#equation 21 
func.Fdos <- function (time){ #Checked
	Fdos <- as.numeric(Params["Csr","Value"]) * Lake[Time,"RESdl"] * Timestep
	Fdos18O <- Fdos * Isotopes[Time,"RESdl_18O"]
	FdosD <- Fdos * Isotopes[Time,"RESdl_D"]
	Fdos <- list(Fdos = Fdos, Fdos18O = Fdos18O, FdosD = FdosD)
	return (Fdos)
}

#Evaporation section. Note that evaporation is calculated on a daily basis, and then multiplied by 30 in Steinman's model. In our case, we'll probably multiply it by the number of days in the month, as that's relatively easy to code for (feb = 28 days).
#F_E = E*Lake_SA__m2/1000 <- and here's the equation from Steinman's Stella model. :D
#Ra is probably defined from the simplified expressions on P696 in Valiantzas' paper.
#This means that Ra can be derived entirely from the site latitude. (pretty cool huh)

#Equation 22
#Monthly data requirement
#ALBlake <- 0.08 #albedo for lake (from Steinman 2010) These three variables have been moved into the preferences file.
#ALBearth <- 0.25 #albedo for surface (from Steinman 2010)
#PWF <- 1 # Penman Wind Constant. 1 for original Penman wind function. 0.5 for reduced Penman wind function. 0 for Linacre wind function. 
func.E <- function(Time){ #Checked
	if (Met[Time,"Ta"] <=0) {
	E <- 0
	} else {
	E <- 0.051 * (1-as.numeric(MetPrefs["ALBlake","Value"])) * Met[Time,"Rs"] * ((Met[Time,"Ta"] + 9.5)^0.5) - (2.4 * ((Met[Time,"Rs"]/func.Ra(as.numeric(Prefs["Latitude","Value"])))^2)) + 0.052 * (Met[Time,"Ta"] + 20) * (1 - Met[Time,"RH"]/100) * (as.numeric(MetPrefs["PWF","Value"]) - 0.38 + 0.54 * Met[Time,"WS"])
	}
	E <- E *0.001 #convert to metres #Daily evap rate.
	if (Datarate == "Monthly"){
	E <- DaysInMonth(Met[Time,"Month"],Met[Time,"Year"]) * E
	}
	return (E)
}#end function

#Equation 23
#Major change from Chimble 1. PWF/2 was used in Chimble 1, changed to 0.5 here.
func.PET <- function(Time){ #Checked
	if (Met[Time,"Ta"] <=0) {
	E <- 0
	} else {
	E <- 0.051 * (1-as.numeric(MetPrefs["ALBearth","Value"])) * Met[Time,"Rs"] * ((Met[Time,"Ta"] + 9.5)^0.5) - (2.4 * ((Met[Time,"Rs"]/func.Ra(as.numeric(Prefs["Latitude","Value"])))^2)) + 0.048 * (Met[Time,"Ta"] + 20) * (1 - Met[Time,"RH"]/100) * (0.5 + 0.536 * Met[Time,"WS"]) #Daily evap rate.
	}
	E <- E *0.001 #convert to metres
	if (Datarate == "Monthly"){
	E <- DaysInMonth(Met[Time,"Month"],Met[Time,"Year"]) * E
	}
	#Major change from CHIMBLE 1. Error in previous version may have limited all months to 28 days. Time used instead of Month. 
	return (E)
}#end function

#Fp is precip over the lake area
func.Fp <- function(Time){#Checked
	Fp <- Met[Time,"P"] * 0.001 * Lake[Time,"Area"] * Timestep
	Fp18O <- Fp * Met[Time,"d18Op"]
	FpD <- Fp * Met[Time,"dDp"]
	Fp <- list(Fp = Fp, Fp18O = Fp18O, FpD = FpD)
	return (Fp)
}#end function

#This routines sets a seasonal variable Kc, based on month, with a slight offset to account for the lag in average lake and air temperatures. 
func.KcLake <- function (Time){ #Checked
	Kcamplitude <- (as.numeric(MetPrefs["KcLwin","Value"]) - as.numeric(MetPrefs["KcLsum","Value"]))/2
	Kcseasonoffset <- 0.8
	Kcshift <- (as.numeric(MetPrefs["KcLwin","Value"])+as.numeric(MetPrefs["KcLsum","Value"]))/2
	KcLake <- Kcshift + (sin(((func.Month(Time)/12)+Kcseasonoffset)*pi*2)*Kcamplitude)
	# was - Kcshift + (sin(((Met[Time,"Month"] + (Met[Time,"Day"]/DaysInMonth(Met[Time,"Month"],Met[Time,"Year"])))/12 +Kcseasonoffset)*pi*2)*Kcamplitude)
return (KcLake)
}

func.Fe <- function(Time) { #Checked
	Fe <- func.E(Time) * Lake[Time,"Area"] * Timestep * func.KcLake(Time)
	Fe18O <- Fe * Isotopes[Time,"E_18O"]
	FeD <- Fe *  Isotopes[Time,"E_D"]
	Fe <- list(Fe = Fe, Fe18O = Fe18O, FeD = FeD)
	return (Fe)
}#end function
#This is the end of the hydrology functions. The next step is to build a routine to tie them all together. 