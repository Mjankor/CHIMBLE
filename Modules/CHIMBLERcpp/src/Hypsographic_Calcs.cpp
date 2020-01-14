#include <Rcpp.h>
#include "Hypsographic_Calcs.hpp"

using namespace Rcpp;

//These functions interpolate between hypsographic data.

// [[Rcpp::export]]
double func_Lakearea(double Volume, NumericMatrix Lakevolumes){
double Mindex = 0;
double intercept;
double slope;
double Area;
// Update Aug2019 - we use 1 for the starting value in the for loop so that if the lake is above the maximum value, then mindex = 0, rather than -1 (which is very bad). 
for (int i = 1; i < Lakevolumes.nrow(); i++){
	if (Volume >= Lakevolumes(i,1)){
		Mindex = i-1;
		break;
	}
}
slope = (Lakevolumes(Mindex,2)- Lakevolumes(Mindex+1,2)) / (Lakevolumes(Mindex,1)- Lakevolumes(Mindex+1,1));
intercept = Lakevolumes(Mindex,2) - Lakevolumes(Mindex,1)*slope;
Area = Volume*slope + intercept;

return Area;
}

// [[Rcpp::export]]
double func_Lakedepth(double Volume, NumericMatrix Lakevolumes){
double Mindex = 0;
double intercept;
double slope;
double Depth;
for (int i = 1; i < Lakevolumes.nrow(); i++){
	if (Volume >= Lakevolumes(i,1)){
		Mindex = i-1;
		break;
	}
}
slope = (Lakevolumes(Mindex,0)- Lakevolumes(Mindex+1,0)) / (Lakevolumes(Mindex,1)- Lakevolumes(Mindex+1,1));
intercept = Lakevolumes(Mindex,0) - Lakevolumes(Mindex,1)*slope;
Depth = Volume*slope + intercept;

return Depth;
}

// [[Rcpp::export]]
double func_Lakevolume(double Depth, NumericMatrix Lakevolumes){
double Mindex = 0;
double intercept;
double slope;
double Volume;
for (int i = 1; i < Lakevolumes.nrow(); i++){
	if (Depth >= Lakevolumes(i,0)){
		Mindex = i-1;
		break;
	}
}
slope = (Lakevolumes(Mindex,1)- Lakevolumes(Mindex+1,1)) / (Lakevolumes(Mindex,0)- Lakevolumes(Mindex+1,0));
intercept = Lakevolumes(Mindex,1) - Lakevolumes(Mindex,0)*slope;
Volume = Depth*slope + intercept;

return Volume;
}