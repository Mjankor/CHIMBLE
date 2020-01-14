if (HistoricalLevels == "Yes"){
HistoricalLevelData <- LoadFile(HistoricalLevelDataPath, "Timeseries")
HistoricalLevelData$Month <- func.MonthNum(HistoricalLevelData$Month)
} else {
} 

#Convert months to numbers

if (HistoricalChem == "Yes"){
HistoricalChemData <- LoadFile(HistoricalChemDataPath, "Timeseries")
HistoricalChemData$Month <- func.MonthNum(HistoricalChemData$Month)
} else {
message ("***You have not specified a historical chemistry file in preferences.")
}

if (HistoricalIsotopes == "Yes"){
HistoricalIsoData <- LoadFile(HistoricalIsoDataPath, "Timeseries")
HistoricalIsoData$Month <- func.MonthNum(HistoricalIsoData$Month)
} else {
message ("***You have not specified a historical isotope file in preferences.")
}


if (HistoricalTemps == "Yes"){
HistoricalTempsData <- LoadFile(HistoricalTempsDataPath, "Timeseries")
HistoricalTempsData$Month <- func.MonthNum(HistoricalTempsData$Month)
} else {
} 
#Convert months to numbers

