#************************** Checks and Reports *********************************
Func.LakeOnlyMassBalCheck <- function(){
FluxSum <- sum (Flux[,"Fsse"]) + sum (Flux[,"Fdse"])+sum (Flux[,"Fe"]) + sum (Flux[,"Fof"]) + sum (Flux[,"Fdos"]) + sum (Flux[,"Fsos"]) - (sum(Flux[,"FpL"]) + sum(Flux[,"Fsf"]) + sum(Flux[,"FpC"]) + sum(Flux[,"Fgw"]) + sum(Flux[,"FExt"]))

VolSum <- (Lake[1,"RESsl"] + Lake[1,"RESdl"] + Lake[1,"RESss"] + Lake[1,"RESds"] + Lake[1,"RESin"] + Lake[1,"RESsp"] + Lake[1,"GW"]) - (Lake[nrow(Lake),"RESsl"] + Lake[nrow(Lake),"RESdl"] + Lake[nrow(Lake),"RESss"] + Lake[nrow(Lake),"RESds"] + Lake[nrow(Lake),"RESin"] + Lake[nrow(Lake),"RESsp"] + Lake[nrow(Lake),"GW"])

message ("System Balance Check.") 

message ("(Lake) Sum of evaporation and outseepage - sum of rainfall, snowfall, external and groundwater = ", round(FluxSum, digits = 1))
message ("Initial system volume - final system volume = ", round(VolSum, digits = 1))
message ("Difference = ",round((FluxSum - VolSum), digits = 2))

SumD <- sum (FluxD[,"Fsse"]) + sum (FluxD[,"Fdse"])+sum (FluxD[,"Fe"]) + sum (FluxD[,"Fof"]) + sum (FluxD[,"Fdos"]) + sum (FluxD[,"Fsos"]) - (sum(FluxD[,"FpL"]) + sum(FluxD[,"Fsf"]) + sum(FluxD[,"FpC"]) + sum(FluxD[,"Fgw"]) + sum(FluxD[,"FExt"]))

Sum18O <- sum (Flux18O[,"Fsse"]) + sum (Flux18O[,"Fdse"])+sum (Flux18O[,"Fe"]) + sum (Flux18O[,"Fof"]) + sum (Flux18O[,"Fdos"]) + sum (Flux18O[,"Fsos"]) - (sum(Flux18O[,"FpL"]) + sum(Flux18O[,"Fsf"]) + sum(Flux18O[,"FpC"]) + sum(Flux18O[,"Fgw"]) + sum(Flux18O[,"FExt"]))

VolD <- (Lake[1,"RESsl"]*Isotopes[1,"RESsl_D"] + Lake[1,"RESdl"]*Isotopes[1,"RESdl_D"] + Lake[1,"RESss"]*Isotopes[1,"RESss_D"] + Lake[1,"RESds"]*Isotopes[1,"RESds_D"] + Lake[1,"RESin"]*Isotopes[1,"RESin_D"] + Lake[1,"RESsp"]*Isotopes[1,"RESsp_D"] + Lake[1,"GW"]*Isotopes[1,"GW_D"])- (Lake[nrow(Lake),"RESsl"]*Isotopes[nrow(Isotopes),"RESsl_D"] + Lake[nrow(Lake),"RESdl"]*Isotopes[nrow(Isotopes),"RESdl_D"] + Lake[nrow(Lake),"RESss"]*Isotopes[nrow(Isotopes),"RESss_D"] + Lake[nrow(Lake),"RESds"]*Isotopes[nrow(Isotopes),"RESds_D"] + Lake[nrow(Lake),"RESin"]*Isotopes[nrow(Isotopes),"RESin_D"] + Lake[nrow(Lake),"RESsp"]*Isotopes[nrow(Isotopes),"RESsp_D"] + Lake[nrow(Lake),"GW"]*Isotopes[nrow(Isotopes),"GW_D"])

Vol18O <- (Lake[1,"RESsl"]*Isotopes[1,"RESsl_18O"] + Lake[1,"RESdl"]*Isotopes[1,"RESdl_18O"] + Lake[1,"RESss"]*Isotopes[1,"RESss_18O"] + Lake[1,"RESds"]*Isotopes[1,"RESds_18O"] + Lake[1,"RESin"]*Isotopes[1,"RESin_18O"] + Lake[1,"RESsp"]*Isotopes[1,"RESsp_18O"] + Lake[1,"GW"]*Isotopes[1,"GW_18O"]) - (Lake[nrow(Lake),"RESsl"]*Isotopes[nrow(Isotopes),"RESsl_18O"] + Lake[nrow(Lake),"RESdl"]*Isotopes[nrow(Isotopes),"RESdl_18O"] + Lake[nrow(Lake),"RESss"]*Isotopes[nrow(Isotopes),"RESss_18O"] + Lake[nrow(Lake),"RESds"]*Isotopes[nrow(Isotopes),"RESds_18O"] + Lake[nrow(Lake),"RESin"]*Isotopes[nrow(Isotopes),"RESin_18O"] + Lake[nrow(Lake),"RESsp"]*Isotopes[nrow(Isotopes),"RESsp_18O"] + Lake[nrow(Lake),"GW"]*Isotopes[nrow(Isotopes),"GW_18O"])

message ("Sum of dD flux out - sum of d18O flux in = ", round(Sum18O,digits = 1))
message ("Sum of d18O volume change = ", round(Vol18O,digits = 1))
message ("Sum of D flux out - sum of dD flux in = ", round(SumD,digits = 1))
message ("Sum of dD volume change = ", round(VolD,digits = 1))
message ("Difference of d18O = ", round(Sum18O - Vol18O,digits = 1))
message ("Difference of dD = ", round(SumD - VolD,digits = 1))
}

Func.GWLakeMassBalCheck <- function(){
FluxSum <- (sum (Flux[,"FpL"]) + sum (Flux[,"Fsf"])+sum (Flux[,"FpC"]) + sum (Flux[,"Fgw"]) + sum (Flux[,"FExt"])) - (sum (Flux[,"Fe"]) + sum(Flux[,"Fof"]) + sum(Flux[,"Fdos"]) + sum(Flux[,"Fsos"]) + sum(Flux[,"Fsse"]) + sum(Flux[,"Fdse"]) + sum(Flux[,"Fdsd"]))

VolSum <- (Lake[nrow(Lake),"RESsl"] + Lake[nrow(Lake),"RESdl"] + Lake[nrow(Lake),"RESss"] + Lake[nrow(Lake),"RESds"] + Lake[nrow(Lake),"RESin"] + Lake[nrow(Lake),"RESsp"]) - (Lake[1,"RESsl"] + Lake[1,"RESdl"] + Lake[1,"RESss"] + Lake[1,"RESds"] + Lake[1,"RESin"] + Lake[1,"RESsp"])

SumD <- (sum (FluxD[,"FpL"]) + sum (FluxD[,"Fsf"])+sum (FluxD[,"FpC"]) + sum (FluxD[,"Fgw"]) + sum (FluxD[,"FExt"])) - (sum (FluxD[,"Fe"]) + sum(FluxD[,"Fof"]) + sum(FluxD[,"Fdos"]) + sum(FluxD[,"Fsos"]) + sum(FluxD[,"Fsse"]) + sum(FluxD[,"Fdse"]) + sum(FluxD[,"Fdsd"]))

Sum18O <- (sum (Flux18O[,"FpL"]) + sum (Flux18O[,"Fsf"])+sum (Flux18O[,"FpC"]) + sum (Flux18O[,"Fgw"]) + sum (Flux18O[,"FExt"])) - (sum (Flux18O[,"Fe"]) + sum(Flux18O[,"Fof"]) + sum(Flux18O[,"Fdos"]) + sum(Flux18O[,"Fsos"]) + sum(Flux18O[,"Fsse"]) + sum(Flux18O[,"Fdse"]) + sum(Flux18O[,"Fdsd"]))

VolD <- (Lake[nrow(Lake),"RESsl"]*Isotopes[nrow(Isotopes),"RESsl_D"] + Lake[nrow(Lake),"RESdl"]*Isotopes[nrow(Isotopes),"RESdl_D"] + Lake[nrow(Lake),"RESss"]*Isotopes[nrow(Isotopes),"RESss_D"] + Lake[nrow(Lake),"RESds"]*Isotopes[nrow(Isotopes),"RESds_D"] + Lake[nrow(Lake),"RESin"]*Isotopes[nrow(Isotopes),"RESin_D"] + Lake[nrow(Lake),"RESsp"]*Isotopes[nrow(Isotopes),"RESsp_D"]) - (Lake[1,"RESsl"]*Isotopes[1,"RESsl_D"] + Lake[1,"RESdl"]*Isotopes[1,"RESdl_D"] + Lake[1,"RESss"]*Isotopes[1,"RESss_D"] + Lake[1,"RESds"]*Isotopes[1,"RESds_D"] + Lake[1,"RESin"]*Isotopes[1,"RESin_D"] + Lake[1,"RESsp"]*Isotopes[1,"RESsp_D"])

Vol18O <- (Lake[nrow(Lake),"RESsl"]*Isotopes[nrow(Isotopes),"RESsl_18O"] + Lake[nrow(Lake),"RESdl"]*Isotopes[nrow(Isotopes),"RESdl_18O"] + Lake[nrow(Lake),"RESss"]*Isotopes[nrow(Isotopes),"RESss_18O"] + Lake[nrow(Lake),"RESds"]*Isotopes[nrow(Isotopes),"RESds_18O"] + Lake[nrow(Lake),"RESin"]*Isotopes[nrow(Isotopes),"RESin_18O"] + Lake[nrow(Lake),"RESsp"]*Isotopes[nrow(Isotopes),"RESsp_18O"]) - (Lake[1,"RESsl"]*Isotopes[1,"RESsl_18O"] + Lake[1,"RESdl"]*Isotopes[1,"RESdl_18O"] + Lake[1,"RESss"]*Isotopes[1,"RESss_18O"] + Lake[1,"RESds"]*Isotopes[1,"RESds_18O"] + Lake[1,"RESin"]*Isotopes[1,"RESin_18O"] + Lake[1,"RESsp"]*Isotopes[1,"RESsp_18O"])

FluxGW <- sum(Flux[,"GW_Recharge"]) + sum(Flux[,"Fsos"]) + sum(Flux[,"Fdos"]) - sum(Flux[,"GW_ET"]) + sum(Flux[,"GW_Ext"]) + sum(Flux[,"GW_FHead"]) 

FluxGWD <- sum(FluxD[,"GW_Recharge"]) + sum(FluxD[,"Fsos"]) + sum(FluxD[,"Fdos"]) - sum(FluxD[,"GW_ET"]) + sum(FluxD[,"GW_Ext"]) + sum(FluxD[,"GW_FHead"]) 
FluxGW18O <- sum(Flux18O[,"GW_Recharge"]) + sum(Flux18O[,"Fsos"]) + sum(Flux18O[,"Fdos"]) - sum(Flux18O[,"GW_ET"]) + sum(Flux18O[,"GW_Ext"]) + sum(Flux18O[,"GW_FHead"]) 

MaxGWError <- max(Flux[,"GW_Err"])
MinGWError <- min(Flux[,"GW_Err"])

VolGW <- Lake[1,"GW"] - Lake[nrow(Lake),"GW"]
VolGWD <- Lake[1,"GW"]*Isotopes[1,"GW_D"] - Lake[nrow(Lake),"GW"]*Isotopes[nrow(Isotopes),"GW_D"]
VolGWD <- Lake[1,"GW"]*Isotopes[1,"GW_18O"] - Lake[nrow(Lake),"GW"]*Isotopes[nrow(Isotopes),"GW_18O"]

message ("Model run complete\n****************\nSystem Balance Check.") 
message ("(Lake) Sum of evaporation and outseepage - sum of rainfall, snowfall, external and groundwater = ", round(FluxSum, digits = 1))
message ("Initial system volume - final system volume = ", round(VolSum, digits = 1))
message ("Difference = ",round((FluxSum - VolSum), digits = 2))
message ("Sum of d18O flux out - sum of d18O flux in = ", round(Sum18O,digits = 1))
message ("Sum of d18O volume change = ", round(Vol18O,digits = 1))
message ("Difference of d18O = ", round(Sum18O - Vol18O,digits = 1))
message ("Sum of D flux out - sum of dD flux in = ", round(SumD,digits = 1))
message ("Sum of dD volume change = ", round(VolD,digits = 1))
message ("Difference of dD = ", round(SumD - VolD,digits = 1))
message ("Max GW error (%) ", round(MaxGWError,digits = 4))
message ("Min GW error (%) ", round(MinGWError,digits = 4))
}
if (GroundwaterModule == "No"){
Func.LakeOnlyMassBalCheck()
} else {
Func.GWLakeMassBalCheck()
}