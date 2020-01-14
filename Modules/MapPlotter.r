#rotate the matrix
func.rotatemat <- function(m) t(m)[,nrow(m):1]

func.plan <- function(GWMatrix,lowest=NULL, highest=NULL){
GWheadrotated <- func.rotatemat(GWMatrix)

#zlim <- seq(50,300, by = 5)
if (is.null(lowest)) lowest <- floor(min(GWMatrix))
if (is.null(highest)) highest <- ceiling(max(GWMatrix))

message(lowest, ", ", highest)

zlim <- seq(lowest, highest, length.out = 50)
colours <- colorRampPalette(c("black", "white"))( length(zlim) )

#Need to expand this to include all scenarios
#if (floor (min(GWMatrix)) < 0){
#zlim1 <- seq(floor(min(GWMatrix)), 0, length.out = 25)
#zlim2 <- seq(0, ceiling(max(GWMatrix)), length.out = 25)
#}
#colours <- c(colorRampPalette(c("black", "blue"))(length(zlim1)),colorRampPalette(c("black", "red"))(length(zlim2)))



y <- 10*1:ncol(GWheadrotated)
x <- 10*1:nrow(GWheadrotated)

filled.contour(x,y,GWheadrotated, axes=TRUE, frame.plot=TRUE,
			plot.axes = { axis(1, seq(0, 10*ncol(Topogrid), by = 500), cex.axis = 1.5)
			axis(2, seq(0, 10*nrow(Topogrid), by = 500), cex.axis = 1.5) },
            key.title = title(main = "Height\n(meters)", cex.main = 1.5,   font.main= 4),
			key.axes = axis(4, seq(floor(min(GWMatrix)), ceiling(max(GWMatrix)), length.out = 50), cex.axis = 1),
			levels = (zlim), col = colours,
			plot.title = title(main = deparse(substitute(GWMatrix)), cex.main = 1.5,   font.main= 4))
}