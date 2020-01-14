#ifndef GWCALCS
#define GWCALCS
using namespace Rcpp;

List Headcalc(NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakebedgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double Timestep, double Runtime, double ETexdepth, NumericMatrix Lakedepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, NumericMatrix GWNumbers);
void func_ClassGrid(long Time, NumericMatrix Topogrid, NumericMatrix LakeDepths, NumericMatrix GWcatchmentgrid, NumericMatrix Classgrid, NumericMatrix Lakevolumes);
void func_LakeConductance(NumericMatrix Lakecondgrid, int Cellxdim, int Cellydim, NumericMatrix GWlakebedgrid, NumericMatrix GWlakekgrid);
void func_RechargeGrid (long Time, int Cellxdim, int Cellydim, NumericMatrix GWRechargegrid, NumericMatrix Classgrid, NumericMatrix Flux, double Timestep);
void func_ETGrid (long Time, NumericMatrix GWETgrid, NumericMatrix Flux);
void func_LakeDepths(NumericMatrix Lake, long Time, NumericMatrix LakeDepths);


#endif // GWCALCS