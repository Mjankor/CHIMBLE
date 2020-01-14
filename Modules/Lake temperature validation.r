#If air temp + lake temp offset would set the lake to less than 4º, then lake goes to 4º (Same as Steinman 2010)

func.Toff <- function(temp){
if (temp < 4) {
temp <- 4 - Met[i,"Ta"]
} else {
temp <- Met[i,"Toff"]
}
return (temp)
}

for (i in 1:nrow(Met)){ #Parallelizable
Met[i,"Toff"] <- func.Toff(Met[i,"Toff"]+Met[i,"Ta"])
}

