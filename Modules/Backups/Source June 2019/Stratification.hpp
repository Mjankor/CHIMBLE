#ifndef STRATIFICATION_H
#define STRATIFICATION_H
using namespace Rcpp;

double func_Density (double LakeTemperature, double Salinity);
double func_LakeDrag (double LakeTemperature, double AirTemperature, double WindSpeed, double NeutralDragCoeff);
double func_LakeAlbedo (int Julian, double LatitudeRadians);
double func_IceAlbedo (double AirTemperature, double FreezingTemp);
double func_SnowAlbedo (bool MeltFlag);
double func_Teten(double Temperature);
NumericVector func_LatSensHeat (double LakeTemperature, double AirTemperature, double WindSpeed, double LakeIceThickness, double FreezingTemp, double SurfacePressure, double RH, double NeutralDragCoeff);
NumericVector TridiagSolver (int NumSystems, int Columns, int NumUnknowns, NumericMatrix SubdiagMatrix, NumericMatrix MainDiagMatrix, NumericMatrix SuperDiagMatrix, NumericMatrix RHS);
NumericVector func_IceSWRad (double SW, double LakeIceThickness, double SnowThickness);
double func_FreezingPoint(double SurfacePressure, double Salinity);
double func_OldDownwardLongwave(double Ra, double Rs, double AirTemperature, double WaterVaporPressure);
double func_DownwardLongwave(long Time, double LatitudeRadians, NumericMatrix Met);
double func_LongWave (double LakeIceThickness, double WaterTemperature, double TempIce, double Ra, double Rs, double AirTemperature, double RH, NumericMatrix Met, long Time, double LatitudeRadians);
void func_TempProfile (double LakeIceThickness, double LakeDepth, double LakeSurfaceArea, double LakeTemperature, double AirTemperature, double WindSpeed, double SurfacePressure, double RH, double Salinity, double TempIce, double Ra, double Rs, double LatitudeRadians, int Julian, NumericMatrix Stratification, double NeutralDragCoeff, double SWExtCoef, NumericMatrix Met, long Time);
void func_ConvectionMixing (NumericMatrix Stratification, double LakeDepth);
double func_Interpolate(double StartValue, double EndValue, int Days, int CurrentDay);
double func_CleanAirAct (NumericMatrix Stratification, double LakeDepth);
void func_RunStratification (long Time, NumericMatrix Met, NumericMatrix Lake, NumericMatrix ChemConcentrations, NumericMatrix Stratification, double LatitudeRadians, double NeutralDragCoeff, double SWExtCoef, double ThermoclineTargetThickness, double ThermoclineMaximum);



#endif // STRATIFICATION_H