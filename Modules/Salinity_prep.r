#First row = Column number for the ions in Chemistry
#Second row = ResSL
#Third row = Column number for the ions in Chemistry
#Fourth row = ResDL
#Fifth row = Molar Masses.
ChemConcentrations <- matrix(data = as.numeric(0),nrow = 5, ncol = 14, dimnames = list(c(),c("Na","Mg","Ca","K","Cl","Br","NaCl", "MgCL2", "CaCl2", "KCl", "NaBr", "Mass", "Salinity", "Act_Water")))

ChemNames <-colnames(Chemistry)
for (i in 1:ncol(Chemistry)){
	if (ChemNames[i] == "RESsl_Na"){
	ChemConcentrations[1,"Na"] = i
	}
	if (ChemNames[i] == "RESsl_Mg"){
	ChemConcentrations[1,"Mg"] = i
	}
	if (ChemNames[i] == "RESsl_Ca"){
	ChemConcentrations[1,"Ca"] = i
	}
	if (ChemNames[i] == "RESsl_K"){
	ChemConcentrations[1,"K"] = i
	}
	if (ChemNames[i] == "RESsl_Cl"){
	ChemConcentrations[1,"Cl"] = i
	}
	if (ChemNames[i] == "RESsl_Br"){
	ChemConcentrations[1,"Br"] = i
	}

	if (ChemNames[i] == "RESdl_Na"){
	ChemConcentrations[3,"Na"] = i
	}
	if (ChemNames[i] == "RESdl_Mg"){
	ChemConcentrations[3,"Mg"] = i
	}
	if (ChemNames[i] == "RESdl_Ca"){
	ChemConcentrations[3,"Ca"] = i
	}
	if (ChemNames[i] == "RESdl_K"){
	ChemConcentrations[3,"K"] = i
	}
	if (ChemNames[i] == "RESdl_Cl"){
	ChemConcentrations[3,"Cl"] = i
	}
	if (ChemNames[i] == "RESdl_Br"){
	ChemConcentrations[3,"Br"] = i
	}

ChemConcentrations[5,"Na"] <- 22.989770
ChemConcentrations[5,"Mg"] <- 24.305
ChemConcentrations[5,"Ca"] <- 40.078
ChemConcentrations[5,"K"] <- 39.0983
ChemConcentrations[5,"Cl"] <- 35.453
ChemConcentrations[5,"Br"] <- 79.904
ChemConcentrations[5,"NaCl"] <- 58.44277
ChemConcentrations[5,"MgCL2"] <- 95.211
ChemConcentrations[5,"CaCl2"] <- 110.984
ChemConcentrations[5,"KCl"] <- 74.5513
ChemConcentrations[5,"NaBr"] <- 102.89377

ChemConcentrations[2,"Act_Water"] <- 1;
ChemConcentrations[4,"Act_Water"] <- 1;
}

ChemParams <- matrix(data = as.numeric(0),nrow = 13, ncol = 7, dimnames = list(c(),c("Salinity", "Density", "NaClAct", "MgCl2Act", "CaCl2Act", "KClAct", "NaBrAct")))
#Salinity densities based on seawater, from https://www.mt-oceanography.info/Utilities/density.html
#using 20ºC
#Activities derived using AIOMFAC (http://www.aiomfac.caltech.edu/model.html)
ChemParams[,1] = c(0,25,50,75,100,125,150,175,200,225,250,275,300)
ChemParams[,2] = c(1000, 1017.154, 1036.256, 1055.669, 1075.45, 1095.63, 1116.232, 1137.271, 1158.76, 1180.709, 1203.127, 1226.019, 1249.392)
ChemParams[,3] = c(1, 0.986, 0.970, 0.953, 0.934, 0.914, 0.891, 0.866, 0.839, 0.809, 0.776, 0.741, 0.704)
ChemParams[,4] = c(1, 0.987, 0.972, 0.952, 0.926, 0.896, 0.859, 0.816, 0.766, 0.711, 0.652, 0.588, 0.521)
ChemParams[,5] = c(1, 0.989, 0.977, 0.962, 0.945, 0.924, 0.899, 0.870, 0.836, 0.796, 0.752, 0.702, 0.648)
ChemParams[,6] = c(1, 0.989, 0.977, 0.965, 0.953, 0.939, 0.924, 0.909, 0.892, 0.874, 0.855, 0.834, 0.813)
ChemParams[,7] = c(1, 0.992, 0.983, 0.974, 0.964, 0.953, 0.940, 0.926, 0.911, 0.894, 0.875, 0.855, 0.832)

