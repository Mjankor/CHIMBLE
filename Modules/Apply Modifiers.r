#Apply modifiers to Met data
Met[,"P"] <- Met[,"P"] * as.numeric(Prefs["PMod","Value"])
Met[,"Ta"] <- Met[,"Ta"] + as.numeric(Prefs["TaMod","Value"])
Met[,"Toff"] <- Met[,"Toff"] + as.numeric(Prefs["TwMod","Value"])
Met[,"WS"] <- Met[,"WS"] * as.numeric(Prefs["WSMod","Value"])
Met[,"RH"] <- Met[,"RH"] * as.numeric(Prefs["RHMod","Value"])
Met[,"Rs"] <- Met[,"Rs"] * as.numeric(Prefs["RsMod","Value"])
Met[,"Groundwater"] <- Met[,"Groundwater"] + as.numeric(Prefs["GWMod","Value"])
Met[,"SL"] <- Met[,"SL"] * as.numeric(Prefs["StratMod","Value"])

#And a few fixes to resolve data that's out of range. eg: humidity can't be <0 or >100%
Met[,"RH"][Met[,"RH"] > 100] = 100
Met[,"RH"][Met[,"RH"] < 0] = 0

#And validate lake temperatures to make sure it doesn't freeze
source("Modules/Lake temperature validation.r")
