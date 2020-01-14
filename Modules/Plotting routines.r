LoadGraphs <- function(){
Dates <- SyncDates() #Calculate dates for the first of each month to avoid data overload.
LoadBasicGraphs()
if (HistoricalLevels == "Yes"){
	HistoricalLevelData <- HistoricalModelPrep()
	Calcheck(HistoricalLevelData, "nograph")
	Histplot(HistoricalLevelData)
}
if (HistoricalTemps == "Yes"){
	HistoricalTempsData <- HistoricalTempsPrep()
	Tempsplot(HistoricalTempsData)
}
}

HistPlots <- function(){
if (HistoricalLevels == "Yes"){
	HistoricalLevelData <- HistoricalModelPrep()
	Calcheck(HistoricalLevelData, "nograph")
	Histplot(HistoricalLevelData)
}
if (HistoricalTemps == "Yes"){
	HistoricalTempsData <- HistoricalTempsPrep()
	Tempsplot(HistoricalTempsData)
}

if (HistoricalChem == "Yes"){
	HistoricalChemData <- HistoricalChemPrep()
	HistChemplot(HistoricalChemData)
}
if (HistoricalIsotopes == "Yes"){
	HistoricalIsoData <- HistoricalIsoPrep()
	HistIsoplot(HistoricalIsoData)
}

}

SyncDates <- function (){
if (Datarate == "Monthly"){
Dates <- Lake[Lake[,"Month"] %in% c(1) & Lake[,"Day"] %in% c(1:ceiling (30*Timestep)),]
} else {
Dates <- Lake[Lake[,"Month"] %in% c(1) & Lake[,"Day"] %in% c(1),]
}
return (Dates)
}

LoadBasicGraphs <- function(){
quartz(width = 12, height = 10)
par(mfrow=c(2,2))

plot (1:nrow(Lake),Lake[,"Depth"],main="Lake depth over time",xlab="Time", ylab="Depth", type = "l", col = "blue", xaxt = 'n')
axis (side = 1, at = (Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), labels = Dates[,"Year"])
abline (v=(Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), col="grey", lwd = 0.5)
#abline (h=seq(ceiling(min(Lake$Depth)),floor(max(Lake$Depth)), 1) , col="grey83")
abline (h=seq(min(Lake[,"Depth"]),max(Lake[,"Depth"]), 1) , col="grey83")
textloc1 <- 8000
textloc2 <- (4 * max(Lake[,"Depth"]) + 1 * min(Lake[,"Depth"]))/5
text(textloc1, textloc2, paste("GroundwaterModule = ", GroundwaterModule, "\nKcSS = ", KcSS, "\nKcDS = ", KcDS, "\nCsr =",  Csr, "\nCin =", Cin, "\nPWF =", PWF, sep = ''), pos = 4)
plot (1:nrow(Isotopes),Isotopes[,"RESsl_18O"],main="RESsl & RESdl O18 delta",xlab="Time", type = "l", ylab="dO18", col = "blue", xaxt = 'n')
axis (side = 1, at = (Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), labels = Dates[,"Year"] )
abline (v=(Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), col="grey", lwd = 0.5)
#abline (h=seq(ceiling(min(Lake$Depth)),floor(max(Lake$Depth)), 1) , col="grey83")
lines (1:nrow(Isotopes),Isotopes[,"RESdl_18O"])

plot (1:nrow(Clean),Clean[,"Depth"],main="Lake depth over time (Cleaned)",xlab="Time", type = "l", ylab="Depth", col = "blue", xaxt = 'n')
axis (side = 1, at = Dates[,"Timestep"], labels = Dates[,"Year"])
abline (v=Dates[,"Timestep"], col="grey", lwd = 0.5)
abline (h=seq(min(Lake[,"Depth"]), max(Lake[,"Depth"]), 1) , col="grey83")
#abline (h=seq(ceiling(min(Lake$Depth)),floor(max(Lake$Depth)), 1) , col="grey83")
#lines(1:nrow(Clean),Clean[,"Depth"])
textloc1 <- 8000
textloc2 <- (4 * max(Lake[,"Depth"]) + 1 * min(Lake[,"Depth"]))/5
text(textloc1, textloc2, paste("GroundwaterModule = ", GroundwaterModule, "\nKcSS = ", KcSS, "\nKcDS = ", KcDS, "\nCsr =",  Csr, "\nCin =", Cin, "\nPWF =", PWF, sep = ''), pos = 4)
plot (1:nrow(Isotopes),Isotopes[,"RESsl_D"],main="RESsl & RESdl D delta",xlab="Time", type = "l", ylab="dD", col = "blue", xaxt = 'n')
axis (side = 1, at = Dates[,"Timestep"]/Timestep, labels = Dates[,"Year"] )
abline (v=(Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), col="grey", lwd = 0.5)
#abline (h=seq(ceiling(min(Lake$Depth)),floor(max(Lake$Depth)), 1) , col="grey83")
lines (1:nrow(Isotopes),Isotopes[,"RESdl_D"])

quartz(width = 12, height = 10)
par(mfrow=c(2,1))
plot (1:nrow(Clean),Clean[,"d18O_sampledepth"],main="RESsl & RESdl d18O (Clean)", sub = paste("Sample Depth: ", Sampledepth),xlab="Time", type = "l", ylab="dO18", col = "blue",xaxt = 'n')
axis (side = 1, at = Dates[,"Timestep"], labels = Dates[,"Year"] )
abline (v=Dates[,"Timestep"], col="grey", lwd = 0.5)
#abline (h=seq(ceiling(min(Lake$Depth)),floor(max(Lake$Depth)), 1) , col="grey83")
plot (1:nrow(Clean),Clean[,"dD_sampledepth"],main="RESsl & RESdl dD (Clean)", sub = paste("Sample Depth: ", Sampledepth),xlab="Time", type = "l", ylab="dD", col = "blue",xaxt = 'n')
axis (side = 1, at = Dates[,"Timestep"], labels = Dates[,"Year"]  )
abline (v=Dates[,"Timestep"], col="grey", lwd = 0.5)
#abline (h=seq(ceiling(min(Lake$Depth)),floor(max(Lake$Depth)), 1) , col="grey83")

quartz(width = 12, height = 10)
par(mfrow=c(2,1))
LELfull <- lm(Clean[,"dD_sampledepth"]~ Clean[,"d18O_sampledepth"]) # Create a LEL from the D and D
plot (Clean[,"d18O_sampledepth"], Clean[,"dD_sampledepth"],main="Complete LEL", sub = paste("Sample Depth: ", Sampledepth),xlab="dD", ylab="dD", col = "blue")
abline(LELfull, col="red")
textloc1 <- min(Clean[,"d18O_sampledepth"])
textloc2 <- max(Clean[,"dD_sampledepth"])
text(textloc1, textloc2, paste("LEL slope = ", round(LELfull$coefficients[2], digits=2)), pos = 4)

LEL <- lm(Clean[36:nrow(Clean),"dD_sampledepth"] ~ Clean[36:nrow(Clean),"d18O_sampledepth"]) # Create a LEL from the D and D
plot (Clean[36:nrow(Clean),"d18O_sampledepth"],Clean[36:nrow(Clean),"dD_sampledepth"], main=paste("LEL 36:", nrow(Clean), "months."), sub = paste("Sample Depth: ", Sampledepth),xlab="dD", ylab="dD", col = "blue")
abline(LEL, col="red")
textloc1 <- min(Clean[36:nrow(Clean),"d18O_sampledepth"])
textloc2 <- max(Clean[36:nrow(Clean),"dD_sampledepth"])
text(textloc1, textloc2, paste("LEL slope = ", round(LEL$coefficients[2], digits=2)), pos = 4)
#add text
}
#************************** Check against historical data ******************************

#There are two separate routines for monthly or daily data. These can probably be condensed, with just a single "if" on the subset selection, but for the moment I'll keep the separate in case they require different graphs, etc.

Calcheck <- function(HistoricalLevelData, graph = "graph"){
HistoricalLevelData$Residuals <- HistoricalLevelData$Height - HistoricalLevelData$Modelledlevel
if (graph == "graph"){
quartz()
par(mfrow=c(2,1))
plot (HistoricalLevelData$Height, HistoricalLevelData$Modelledlevel,xlab="Historical Level", ylab="Modelled Level")
#plot (1:nrow(HistoricalLevelData), HistoricalLevelData$Residuals)
abline (a=0, b=1, col="black", lwd = 1)
#message ("D=",round(ISOD[6149,10],digits= 2),", D=",round(ISOD[6149,10], digits =2))
plot (HistoricalLevelData$Height, HistoricalLevelData$Residuals, main="Modelled Lake Level Residuals", sub="", xlab="Depth", ylab="Difference (m)")
}
#message ("Average of residuals = ", round(mean(HistoricalLevelData$Residuals), digits = 3), ". Standard deviation = ", round(sd(HistoricalLevelData$Residuals), digits = 3))
message ("Ave=", round(mean(HistoricalLevelData$Residuals), digits = 3), ". SD=", round(sd(HistoricalLevelData$Residuals), digits = 3))
return (c(round(mean(HistoricalLevelData$Residuals), digits = 3),round(sd(HistoricalLevelData$Residuals), digits = 3)))
} #end function


HistoricalModelPrep <- function () {
HistoricalLevelData$Modelledlevel <- 0
HistoricalLevelData$RecordNumber <- 0
	for (i in 1:nrow(HistoricalLevelData)){
		if (Datarate == "Monthly") {
		Daycheck <- c((HistoricalLevelData$Day[i]- ceiling(30 * Timestep)):(HistoricalLevelData$Day[i] + ceiling(30 * Timestep)))
		Record <- Lake[Lake[,"Year"] %in% HistoricalLevelData$Year[i] & Lake[,"Month"] %in% HistoricalLevelData$Month[i] & Lake[,"Day"] %in% Daycheck,"Depth"]
		RecordNumber <- which(Lake[,"Year"] %in% HistoricalLevelData$Year[i] & Lake[,"Month"] %in% HistoricalLevelData$Month[i] & Lake[,"Day"] %in% Daycheck)
		} else if (Datarate == "Daily"){
		Record <- Lake[Lake[,"Year"] %in% HistoricalLevelData$Year[i] & Lake[,"Month"] %in% HistoricalLevelData$Month[i] & Lake[,"Day"] %in% HistoricalLevelData$Day[i],"Depth"]
		RecordNumber <- which(Lake[,"Year"] %in% HistoricalLevelData$Year[i] & Lake[,"Month"] %in% HistoricalLevelData$Month[i] & Lake[,"Day"] %in% HistoricalLevelData$Day[i])
		}
		if (length(Record) == 0) {
		Record <- NA
		RecordNumber <- NA
		}
	HistoricalLevelData$Modelledlevel[i] <- Record[1]
	HistoricalLevelData$RecordNumber[i] <- RecordNumber[1]
	}
return (HistoricalLevelData)
} # end for

Histplot <- function (HistoricalLevelData) {
quartz(width = 10, height = 6)
plot (1:nrow(Lake),Lake[,"Depth"], type = "l", main="Lake depth over time",xlab="Time", ylab="Depth", col = "blue", xaxt = 'n', lwd = 1.5)
axis (side = 1, at = (Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), labels = Dates[,"Year"])
abline (v=(Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), col="grey", lwd = 0.5)
abline (h=seq(ceiling(min(Lake[,"Depth"])),floor(max(Lake[,"Depth"])), 1) , col="grey83")
z = 3
lines(y = HistoricalLevelData$Height[z:nrow(HistoricalLevelData)], x = HistoricalLevelData$RecordNumber[z:nrow(HistoricalLevelData)], lwd = 1.25)
points(y = HistoricalLevelData$Height[1:z], x = HistoricalLevelData$RecordNumber[1:z])
} #end if


#lines(x = 1:nrow(Lake), y = NoGWlake[,"Depth"], col = chartreuse4, lwd = 1.25)


#And for the temps
HistoricalTempsPrep <- function () {
HistoricalTempsData$ModelledTemps <- 0
HistoricalTempsData$RecordNumber <- 0
	for (i in 1:nrow(HistoricalTempsData)){
		if (Datarate == "Monthly") {
		Daycheck <- c((HistoricalTempsData$Day[i]- ceiling(30 * Timestep)):(HistoricalTempsData$Day[i] + ceiling(30 * Timestep)))
		Record <- Met[Met[,"Year"] %in% HistoricalTempsData$Year[i] & Met[,"Month"] %in% HistoricalTempsData$Month[i] & Met[,"Day"] %in% Daycheck,"Ta"] + Met[Met[,"Year"] %in% HistoricalTempsData$Year[i] & Met[,"Month"] %in% HistoricalTempsData$Month[i] & Met[,"Day"] %in% Daycheck,"Toff"]
		RecordNumber <- which(Met[,"Year"] %in% HistoricalTempsData$Year[i] & Met[,"Month"] %in% HistoricalTempsData$Month[i] & Met[,"Day"] %in% Daycheck)
		} else if (Datarate == "Daily"){
		Record <- Met[Met[,"Year"] %in% HistoricalTempsData$Year[i] & Met[,"Month"] %in% HistoricalTempsData$Month[i] & Met[,"Day"] %in% HistoricalTempsData$Day[i],"Ta"] + Met[Met[,"Year"] %in% HistoricalTempsData$Year[i] & Met[,"Month"] %in% HistoricalTempsData$Month[i] & Met[,"Day"] %in% HistoricalTempsData$Day[i],"Toff"]
		RecordNumber <- which(Met[,"Year"] %in% HistoricalTempsData$Year[i] & Met[,"Month"] %in% HistoricalTempsData$Month[i] & Met[,"Day"] %in% HistoricalTempsData$Day[i])
		}
		if (length(Record) == 0) {
		Record <- NA
		RecordNumber <- NA
		}
	HistoricalTempsData$ModelledTemps[i] <- Record[1]
	HistoricalTempsData$RecordNumber[i] <- RecordNumber[1]
	}
return (HistoricalTempsData)
} # end for

Tempsplot <- function (HistoricalTempsData) {
quartz()
minplot <- min(HistoricalTempsData$RecordNumber)
maxplot <- max(HistoricalTempsData$RecordNumber)

plot (1:nrow(Met),Met[,"Ta"]+Met[,"Toff"], type = "l", main="Lake surface temperature over time",xlab="Time", ylab="Depth", col = "blue", xaxt = 'n', xlim=c(minplot,maxplot), ylim = c(min(Met[,"Ta"]+Met[,"Toff"])-3,max(Met[,"Ta"]+Met[,"Toff"]) + 3))
axis (side = 1, at = Dates[,"Timestep"]/Timestep, labels = Dates[,"Year"])
abline (v=(Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), col="grey", lwd = 0.5)
abline (h=seq(ceiling(min(Met[,"Ta"]+Met[,"Toff"])),floor(max(Met[,"Ta"]+Met[,"Toff"])), 1) , col="grey83")
lines(y = HistoricalTempsData$Temp, x = HistoricalTempsData$RecordNumber)
text(minplot + 1, min((Met[,"Ta"]+Met[,"Toff"])-2), paste("Observed temps average, max, min", round(mean(HistoricalTempsData$Temp),2), ", ", min(HistoricalTempsData$Temp), ", ", max(HistoricalTempsData$Temp)), pos = 4)
text(minplot + 1, min((Met[,"Ta"]+Met[,"Toff"])-3), paste("Modelled temps average, max, min", round(mean(HistoricalTempsData$ModelledTemps),2), ", ", round(min(HistoricalTempsData$ModelledTemps),1), ", ", round(max(HistoricalTempsData$ModelledTemps)),1), pos = 4)

} #end if

HistoricalChemPrep <- function () {
HistoricalChemData$ModelledChem <- 0
HistoricalChemData$RecordNumber <- 0
	for (i in 1:nrow(HistoricalChemData)){
		if (Datarate == "Monthly") {
		Daycheck <- c((HistoricalChemData$Day[i]- ceiling(30 * Timestep)):(HistoricalChemData$Day[i] + ceiling(30 * Timestep)))
		#Record <- Met[Met[,"Year"] %in% HistoricalChemData$Year[i] & Met[,"Month"] %in% HistoricalChemData$Month[i] & Met[,"Day"] %in% Daycheck,"Ta"] + Met[Met[,"Year"] %in% HistoricalChemData$Year[i] & Met[,"Month"] %in% HistoricalChemData$Month[i] & Met[,"Day"] %in% Daycheck,"Toff"]
		RecordNumber <- which(Met[,"Year"] %in% HistoricalChemData$Year[i] & Met[,"Month"] %in% HistoricalChemData$Month[i] & Met[,"Day"] %in% Daycheck)
		} else if (Datarate == "Daily"){
		#Record <- Met[Met[,"Year"] %in% HistoricalChemData$Year[i] & Met[,"Month"] %in% HistoricalChemData$Month[i] & Met[,"Day"] %in% HistoricalChemData$Day[i],"Ta"] + Met[Met[,"Year"] %in% HistoricalChemData$Year[i] & Met[,"Month"] %in% HistoricalChemData$Month[i] & Met[,"Day"] %in% HistoricalChemData$Day[i],"Toff"]
		RecordNumber <- which(Met[,"Year"] %in% HistoricalChemData$Year[i] & Met[,"Month"] %in% HistoricalChemData$Month[i] & Met[,"Day"] %in% HistoricalChemData$Day[i])
		}
		#if (length(Record) == 0) {
		#Record <- NA
		#RecordNumber <- NA
		#}
	#HistoricalChemData$ModelledChem[i] <- Record[1]
	HistoricalChemData$RecordNumber[i] <- RecordNumber[1]
	}
return (HistoricalChemData)
} # end for

HistChemplot <- function (HistoricalChemData) {
quartz()
ChemTemp <- matrix(ncol = 1,nrow = nrow(Met))
#Need to redo this to use the available ions/tds, rather than hardcoded.
for (i in 1:(nrow(Chemistry)-1)){
	if (Sampledepth >= Met[i,"SL"]){
		ChemTemp[i] <- Chemistry[i, "RESdl_Na"]+Chemistry[i, "RESdl_Mg"]+Chemistry[i, "RESdl_Ca"]+Chemistry[i, "RESdl_K"]+Chemistry[i, "RESdl_Sr"]+Chemistry[i, "RESdl_Cl"]+Chemistry[i, "RESdl_SO4"]+Chemistry[i, "RESdl_HCO3"]+Chemistry[i, "RESdl_NO3"]+Chemistry[i, "RESdl_Br"]+Chemistry[i, "RESdl_CO3"]
	} else {
		ChemTemp[i] <- Chemistry[i, "RESsl_Na"]+Chemistry[i, "RESsl_Mg"]+Chemistry[i, "RESsl_Ca"]+Chemistry[i, "RESsl_K"]+Chemistry[i, "RESsl_Sr"]+Chemistry[i, "RESsl_Cl"]+Chemistry[i, "RESsl_SO4"]+Chemistry[i, "RESsl_HCO3"]+Chemistry[i, "RESsl_NO3"]+Chemistry[i, "RESsl_Br"]+Chemistry[i, "RESsl_CO3"]
	}
}

minplot <- min(HistoricalChemData$RecordNumber)
maxplot <- max(HistoricalChemData$RecordNumber)
minyplot <- min(HistoricalChemData$TDS)-2000
maxyplot <- max(HistoricalChemData$TDS)+2000
plot (1:nrow(ChemTemp),ChemTemp, type = "l", main="Lake TDS over time",xlab="Time", ylab="TDS", col = "blue", xaxt = 'n',xlim=c(minplot,maxplot), ylim =c(minyplot,maxyplot))
axis (side = 1, at = Dates[,"Timestep"]/Timestep, labels = Dates[,"Year"])
abline (v=(Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), col="grey", lwd = 0.5)

if (floor(max(ChemTemp)) > 10000){
MaxChem <- signif(floor(max(ChemTemp)), digits = 1)+2000
MinChem <- signif(ceiling(min(ChemTemp)), digits = 1)-2000
} else {
MaxChem <- signif(floor(max(ChemTemp)), digits = 2)+2000
MinChem <- signif(ceiling(min(ChemTemp)), digit = 2)-2000
}
message (MinChem, ", ",MaxChem)
abline (h=seq(MinChem,MaxChem, 1000), col="grey83")
z <- nrow(HistoricalChemData)#14
lines(y = HistoricalChemData$TDS[z:nrow(HistoricalChemData)], x = HistoricalChemData$RecordNumber[z:nrow(HistoricalChemData)])
points(y = HistoricalChemData$TDS[1:z], x = HistoricalChemData$RecordNumber[1:z], cex = 0.5)
lines (1:nrow(ChemTemp), ChemTemp, lwd = 1, col = "blue")
#lines (1:nrow(ChemTempNoGW), ChemTempNoGW, lwd = 1.5, col = "green")
} #end if


HistoricalIsoPrep <- function () {
HistoricalIsoData$ModelledIso <- 0
HistoricalIsoData$RecordNumber <- 0
	for (i in 1:nrow(HistoricalIsoData)){
		if (Datarate == "Monthly") {
		Daycheck <- c((HistoricalIsoData$Day[i]- ceiling(30 * Timestep)):(HistoricalIsoData$Day[i] + ceiling(30 * Timestep)))
		#Record <- Met[Met[,"Year"] %in% HistoricalIsoData$Year[i] & Met[,"Month"] %in% HistoricalIsoData$Month[i] & Met[,"Day"] %in% Daycheck,"Ta"] + Met[Met[,"Year"] %in% HistoricalIsoData$Year[i] & Met[,"Month"] %in% HistoricalIsoData$Month[i] & Met[,"Day"] %in% Daycheck,"Toff"]
		RecordNumber <- which(Met[,"Year"] %in% HistoricalIsoData$Year[i] & Met[,"Month"] %in% HistoricalIsoData$Month[i] & Met[,"Day"] %in% Daycheck)
		} else if (Datarate == "Daily"){
		#Record <- Met[Met[,"Year"] %in% HistoricalIsoData$Year[i] & Met[,"Month"] %in% HistoricalIsoData$Month[i] & Met[,"Day"] %in% HistoricalIsoData$Day[i],"Ta"] + Met[Met[,"Year"] %in% HistoricalIsoData$Year[i] & Met[,"Month"] %in% HistoricalIsoData$Month[i] & Met[,"Day"] %in% HistoricalIsoData$Day[i],"Toff"]
		RecordNumber <- which(Met[,"Year"] %in% HistoricalIsoData$Year[i] & Met[,"Month"] %in% HistoricalIsoData$Month[i] & Met[,"Day"] %in% HistoricalIsoData$Day[i])
		}
		#if (length(Record) == 0) {
		#Record <- NA
		#RecordNumber <- NA
		#}
	#HistoricalIsoData$Modelled18O[i] <- Record[1]
	HistoricalIsoData$RecordNumber[i] <- RecordNumber[1]
	}
return (HistoricalIsoData)
} # end for

HistIsoplot <- function (HistoricalIsoData) {
quartz(width = 8, height = 7)
IsoTemp <- matrix(ncol = 2,nrow = nrow(Met))
for (i in 1:(nrow(Isotopes)-1)){
	if (Sampledepth >= Met[i,"SL"]){
		IsoTemp[i,1] <- Isotopes[i, "RESdl_18O"]
		IsoTemp[i,2] <- Isotopes[i, "RESdl_D"]
	} else {
		IsoTemp[i,1] <- Isotopes[i, "RESsl_18O"]
		IsoTemp[i,2] <- Isotopes[i, "RESsl_D"]
	}
}

minplot <- min(HistoricalIsoData$RecordNumber)
maxplot <- max(HistoricalIsoData$RecordNumber)
minyplot <- min(HistoricalIsoData$Iso18O)-0.5
maxyplot <- max(HistoricalIsoData$Iso18O)+0.5
plot (1:nrow(IsoTemp),IsoTemp[,1], type = "l", main="Lake 18O over time",xlab="Time", ylab="18O", col = "blue", xaxt = 'n',xlim=c(minplot,maxplot), ylim =c(minyplot,maxyplot), lwd = 1.5)
axis (side = 1, at = Dates[,"Timestep"]/Timestep-((1/Timestep)-1), labels = Dates[,"Year"])
abline (v=(Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), col="grey", lwd = 0.5)
abline (h=seq(ceiling(min(IsoTemp)),floor(max(IsoTemp)), 100) , col="grey83")
lines(y = HistoricalIsoData$Iso18O, x = HistoricalIsoData$RecordNumber)

#and for D
quartz()
IsoTemp <- matrix(ncol = 2,nrow = nrow(Met))
for (i in 1:(nrow(Isotopes)-1)){
	if (Sampledepth >= Met[i,"SL"]){
		IsoTemp[i,1] <- Isotopes[i, "RESdl_18O"]
		IsoTemp[i,2] <- Isotopes[i, "RESdl_D"]
	} else {
		IsoTemp[i,1] <- Isotopes[i, "RESsl_18O"]
		IsoTemp[i,2] <- Isotopes[i, "RESsl_D"]
	}
}

minplot <- min(HistoricalIsoData$RecordNumber)
maxplot <- max(HistoricalIsoData$RecordNumber)
minyplot <- min(HistoricalIsoData$IsoD)-2
maxyplot <- max(HistoricalIsoData$IsoD)+2
plot (1:nrow(IsoTemp),IsoTemp[,2], type = "l", main="Lake D over time",xlab="Time", ylab="D", col = "blue", xaxt = 'n',xlim=c(minplot,maxplot), ylim =c(minyplot,maxyplot), lwd = 1.5)
axis (side = 1, at = Dates[,"Timestep"]/Timestep-((1/Timestep)-1), labels = Dates[,"Year"])
abline (v=(Dates[,"Timestep"]/Timestep)-((1/Timestep)-1), col="grey", lwd = 0.5)
abline (h=seq(ceiling(min(IsoTemp)),floor(max(IsoTemp)), 100) , col="grey83")
lines(y = HistoricalIsoData$IsoD, x = HistoricalIsoData$RecordNumber)
} #end if


#Profile plotter

profplot <- function(coffset,roffset){
HistoricalLevelData <- HistoricalModelPrep()
quartz (width = 10, height = 10)
#jpeg(width = 1200, height = 1200, units = "px", filename=paste ("Images/Scenario_", Scenario, "_Profile.jpg"))#Disable for headless mode
par(mfrow=c(2,1))
topo <- Topogrid[,floor(ncol(Topogrid)/2)-coffset]
Lakes <- Classgrid[,floor(ncol(Topogrid)/2)-coffset]
Base <- GWbasegrid[,floor(ncol(GWbasegrid)/2)-coffset]
Lakes[Lakes >= 10] <- NA
for (i in 1:nrow(LakeDepths)){
Lakes[Lakes == i] <- LakeDepths[i,4]
}
heads <- GWheadgrid[,floor(ncol(Topogrid)/2)-coffset]
plot (topo, type = "l", col = "black",main="Profile Plan North-South", xlab="Meters", ylab="Height",xaxt = 'n', ylim = c(min(topo-50),max(topo)))
lines (Lakes, col = "blue")
lines (Base, col = "grey")
lines (heads, col = "blue", lty = 2)
axis (side = 1, at=pretty(1:nrow(Topogrid), n = 10), labels = Cellydim * pretty(1:nrow(Topogrid), n = 10))

text(0, max(topo)-25, paste("K = ", GWKgrid[1,1], "\nLake K = ", GWlakekgrid[1,1], "\nSoil Layer =",  AWCss[1], "\nGeology Values = ", K1,", ", K2,", ",K3,", ",K4,", ",K5,", ",K6,", ",K7,", ",K8,", ",K9,", ", "\nCalcheck = ", Calcheck(HistoricalLevelData, "nograph")[1], ", SD = ", Calcheck(HistoricalLevelData, "nograph")[2], "\nBase = ", Paths["GWbasegrid","Filename"],"\nHeads = ", Paths["GWheadgrid","Filename"], "\nMet = ", Paths["MetData","Filename"], sep = ''), pos = 4, cex = 0.5	)

topo <- Topogrid[floor(nrow(Topogrid)/2)- roffset,]
Lakes <- Classgrid[floor(nrow(Topogrid)/2)- roffset,]
Base <- GWbasegrid[floor(nrow(GWbasegrid)/2)-roffset,]
Lakes[Lakes >= 10] <- NA
for (i in 1:nrow(LakeDepths)){
Lakes[Lakes == i] <- LakeDepths[i,4]
}
heads <- GWheadgrid[floor(nrow(Topogrid)/2)- roffset,]
plot (topo, type = "l", col = "black",main="Profile Plan East-West", xlab="Meters", ylab="Height",xaxt = 'n', ylim = c(min(topo-50),max(topo)))
lines (Lakes, col = "blue")
lines (heads, col = "blue", lty = 2)
lines (Base, col = "grey")
axis (side = 1, at=pretty(1:ncol(Topogrid), n = 10), labels = Cellydim * pretty(1:ncol(Topogrid), n = 10))

#And for the topo plan
func.rotatemat <- function(m) t(m)[,nrow(m):1]
func.profplan <- function(GWMatrix){
quartz (width = 8, height = 10)
#jpeg(width = 1200, height = 1200, units = "px", filename=paste ("Images/Scenario_", Scenario, "_Contour.jpg"))#Disable for headless mode
GWheadrotated <- func.rotatemat(GWMatrix)
#zlim <- seq(50,300, by = 5)
zlim <- pretty(seq(floor(min(Topogrid)), ceiling(max(Topogrid)), length.out = 50), n = 20)
colours <- colorRampPalette(c("black", "white"))( length(zlim) )

y <- Cellxdim*1:ncol(GWheadrotated)
x <- Cellydim*1:nrow(GWheadrotated)

#roffsets are inverted here as it filled contour uses 0,0 = bottom left. 
filled.contour(x,y,GWheadrotated, axes=TRUE, frame.plot=TRUE,
			plot.axes = { axis(1, seq(0, Cellxdim*ncol(Topogrid), by = 500), cex.axis = 1.5)
			axis(2, seq(0, Cellydim*nrow(Topogrid), by = 500), cex.axis = 1.5); rect (Cellxdim * (floor(ncol(Topogrid)/2-coffset)),0,Cellxdim * (floor(ncol(Topogrid)/2)+1-coffset),Cellydim * floor(nrow(Topogrid)), col = "blue") ;  rect (0,Cellxdim * (floor(nrow(Topogrid)/2 + roffset)),Cellxdim * (nrow(Topogrid)),Cellxdim * (floor(nrow(Topogrid)/2)+1 + roffset), col = "blue")},
            key.title = title(main = "Height\n(meters)", cex.main = 1.5,   font.main= 4),
			key.axes = axis(4, zlim, cex.axis = 1),
			levels = (zlim), col = colours,
			plot.title = title(main = "Topography", cex.main = 1.5, font.main= 4))
#dev.off()			
}
func.profplan(Topogrid)
}