#MSE for levels
Levels <- HistoricalModelPrep()
LevelsMSE <- mean((Levels[,5]-Levels[,4])^2)
LevelsRMSE <- sqrt(LevelsMSE)

#MSE for Temps
Temps <- HistoricalTempsPrep()
TempsMSE <- mean((Temps[,5]-Temps[,4])^2)
TempsRMSE <- sqrt(TempsMSE)

#MSE for chemistry
Chem <- HistoricalChemPrep()
ChemTemp <- matrix(ncol = 1,nrow = nrow(Met))
#Need to redo this to use the available ions/tds, rather than hardcoded.
for (i in 1:(nrow(Chemistry)-1)){
	if (Sampledepth >= Met[i,"SL"]){
		ChemTemp[i] <- Chemistry[i, "RESdl_Na"]+Chemistry[i, "RESdl_Mg"]+Chemistry[i, "RESdl_Ca"]+Chemistry[i, "RESdl_K"]+Chemistry[i, "RESdl_Sr"]+Chemistry[i, "RESdl_Cl"]+Chemistry[i, "RESdl_SO4"]+Chemistry[i, "RESdl_HCO3"]+Chemistry[i, "RESdl_NO3"]+Chemistry[i, "RESdl_Br"]+Chemistry[i, "RESdl_CO3"]
	} else {
		ChemTemp[i] <- Chemistry[i, "RESsl_Na"]+Chemistry[i, "RESsl_Mg"]+Chemistry[i, "RESsl_Ca"]+Chemistry[i, "RESsl_K"]+Chemistry[i, "RESsl_Sr"]+Chemistry[i, "RESsl_Cl"]+Chemistry[i, "RESsl_SO4"]+Chemistry[i, "RESsl_HCO3"]+Chemistry[i, "RESsl_NO3"]+Chemistry[i, "RESsl_Br"]+Chemistry[i, "RESsl_CO3"]
	}
}

for (i in 1:nrow(Chem)){
Chem[i,5] <- ChemTemp[Chem[i,6],1]
}
ChemMSE <- mean((Chem[180:nrow(Chem),5]-Chem[180:nrow(Chem),4])^2)
#Only using from record 180, to avoid using earlier highly variable observations.
ChemRMSE <- sqrt(ChemMSE)

message ("Levels RMSE: ", round(LevelsRMSE, digits = 2))
message ("Temps RMSE: ", round(TempsRMSE, digits = 1))
message ("Chem RMSE: ", round(ChemRMSE, digits = 0))