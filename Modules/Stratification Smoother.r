# Smooth stratification data
func.SLsmooth <- function(){
SLsmooth <- data.frame(Timestep = numeric (nrow(MetDf) / Timestep) , SL=numeric (nrow(MetDf) / Timestep)) #(Met$SL[1])
for (Time in 1:(nrow(MetDf))){ #Possibly parallelizable
	for (step in 0:((1/Timestep)-1)) {
	Index <- ((Time-1)/Timestep)+step+1
	if (Time < nrow(MetDf)){
	SLsmooth$SL[Index] <- MetDf$SL[Time] + step*(MetDf$SL[Time+1] - MetDf$SL[Time])*Timestep
	} else {
	SLsmooth$SL[Index] <- MetDf$SL[Time]
	}
	}
	}
	return (SLsmooth)
}
