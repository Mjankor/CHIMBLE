#************************** Main isotopic functions *********************************
#This is the isotopic functions. 
#equation 24.
func.dE <- function(Time, Isotype){
Tlake <- Met[Time,"Ta"]+Met[Time,"Toff"]
#We'll start with equations 31 and 32 and work back from there.
if (Isotype == "18O"){
	C <- 14.25
	ALPHA <- exp ((350410 /(273.15 + Tlake)^3) - (1666.4/(273.15 + Tlake)^2) + (6.7123 * 1/(273.15 + Tlake)) - 0.007685) 
	#equation 28 #Checked against graphs from Horita et al 2008
} else if (Isotype == "D"){
	C <- 12.55
	ALPHA <- exp ((1.1588 * ((273.15 + Tlake)^3)/10^9) - (1.6201 * ((273.15 + Tlake)^2)/10^6) + (0.79484 * (273.15 + Tlake)/1000) + (2.9992 * (10^6)/(273.15 + Tlake)^3) - 0.16104) 
	# equation 29 #Checked against graphs from Horita et al 2008
} #end else if
Esa <- 6.108 * exp((17.27*Met[Time,"Ta"])/(Met[Time,"Ta"] + 237.7)) #Equation 26
Esw <- 6.108 * exp((17.27*Tlake)/(Tlake + 237.7)) 
RALPHA <- 1/ALPHA # equation 30
hn = Met[Time,"RH"] * (Esa/Esw) *.01 # equation 25 
#Modifiers to Evaporation fractionation calculation
theta <- (1 - ((Met[Time,"RH"]/100)*(1-Feedback) + Feedback))/(1 - (Met[Time,"RH"])/100)

Ek <- C * (1 - hn)*theta #equation 32 
#Eeq <- 1000 * (ALPHA - 1)  #Gonfianti convention
Eeq <- 1000 * (1 - RALPHA) #equation 31 
Etot <- Ek + Eeq
func.df <- function(){if(Isotype == "D")"dDp" else "d18Op"}


dA <- Met[Time,func.df()]*RALPHA - Eeq*AtmosphericShift#Steinman fixed
#func.df allows us to grab data from a particular dataframe. eg: ISOD[Time,column]
#need to set this to take from RESdl if SL = 0
#This has been updated with isotopic enrichment version from Gibson, 2002 also used by Steinman.

#Needs work. Have to modify the following to deal with both isotopes.


limit <- (hn*dA + Etot)/(hn - 0.001*Etot) #Not needed in the model as back diffusion should be taken care of by the evaporation calculation.
if (func.SVC(Time) == 0){
	dE <- ((RALPHA * Isotopes[Time,"RESdl_18O"] - (hn * dA) - Etot)/(1 - hn + 0.001*Ek))
	} else if (func.SVC(Time)){
	dE <- ((RALPHA * Isotopes[Time,"RESsl_18O"] - (hn * dA) - Etot)/(1 - hn + 0.001*Ek))
}

#add a percentage of evap back into the atmospheric vapour.
dA <- (Met[Time,func.df()]*RALPHA - Eeq*AtmosphericShift)*(1-Feedback) + Feedback*dE
#and process again
if (func.SVC(Time) == 0){
	dE <- ((RALPHA * Isotopes[Time,"RESdl_18O"] - (hn * dA) - Etot)/(1 - hn + 0.001*Ek))
	} else if (func.SVC(Time)){
	dE <- ((RALPHA * Isotopes[Time,"RESsl_18O"] - (hn * dA) - Etot)/(1 - hn + 0.001*Ek))
}

return (dE)
} #end function. The values for O18 appear comparable with those of Horita 2008 (-30pmil)

