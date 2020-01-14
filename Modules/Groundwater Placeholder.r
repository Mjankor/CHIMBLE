#Rather than having separate functions for with/without the groundwater module, we'll simply create some placeholder data for when the groundwater module is not loaded. This data won't be used, but can be passed as an argument to the groundwater functions. 
Cellxdim <- 50
Cellydim <- 50
Accuracy <- .0001
Maxiterations <- 1000
LakeDepths <- matrix(data = 0, 1, 5) 
Topogrid <- matrix(data = 0, 1, 1) 
GWPumpinggrid <- matrix(data = 0, 1, 1) 
GWgeogrid <- matrix(data = 0, 1, 1) 
GWcatchmentgrid <- matrix(data = 0, 1, 1) 
GWbasegrid <- matrix(data = 0, 1, 1) 
GWboundarygrid <- matrix(data = 0, 1, 1) 
GWheadgrid <- matrix(data = 0, 1, 1) 
GWlakesedimentgrid <- matrix(data = 0, 1, 1) 
GWlakekgrid <- matrix(data = 0, 1, 1) 
# and generate a classification grid too.
Classgrid <- matrix(data = 0, 1, 1) 
Lakecondgrid <- matrix(data = 0, 1, 1) 
#and backup the boundary grid, as it will be modified regularly.
GWboundarybackup <- GWboundarygrid
GWRechargegrid<- matrix(data = 0, 1, 1) 
GWKgrid<- matrix(data = 0, 1, 1) 
GWSgrid<- matrix(data = 0, 1, 1) 
GWETgrid<- matrix(data = 0, 1, 1) 
ETexdepth <- 2
GWNumbers <- matrix(data = 0, 1, 1) 
GWStats <- c(1)