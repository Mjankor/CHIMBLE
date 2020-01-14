#ifndef MASSBAL_H
#define MASSBAL_H
using namespace Rcpp;

std::string StrSplitFront (std::string StringToSplit, std::string Delim);
std::string StrSplitBack (std::string StringToSplit, std::string Delim);
void func_FluxFpCatchment(NumericMatrix Flux, NumericMatrix Met, long Time, double CA, NumericMatrix Lake, double Timestep);
void func_FluxFpLake(NumericMatrix Flux, long Time, NumericMatrix Met, NumericMatrix Lake, double Timestep);
void func_FluxFe (NumericMatrix Flux, long Time, double Timestep, NumericMatrix Lake, NumericMatrix Met, double ALBlake, NumericMatrix Isotopes, double PWF, double LatitudeRadians, std::string Datarate, double KcLwin, double KcLsum, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, double Feedback);
void func_FluxFsm (NumericMatrix Flux, double CA, long Time, NumericMatrix Lake,NumericMatrix Met, NumericMatrix Isotopes, std::string Datarate, double Timestep);
void func_FluxFext (NumericMatrix Flux, long Time, NumericMatrix Met, NumericMatrix Isotopes, double Timestep);
void func_FluxFdlm (NumericMatrix Flux, long Time, NumericMatrix Isotopes, NumericMatrix Met, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype);
void func_FluxSoil (NumericMatrix Flux, long Time, float AWCss, float AWCds, double CA, NumericMatrix Lake, NumericMatrix Met, NumericMatrix Isotopes, double KcDS, double  KcSS, double LatitudeRadians,float RunoffRatio, double ALBearth, double ALBlake, double PWF, std::string Datarate, double Timestep);
void func_FluxSoilEvap (NumericMatrix Flux, double CA, long Time, NumericMatrix Met, NumericMatrix Lake,NumericMatrix Isotopes, float AWCss, float AWCds, double KcDS, double LatitudeRadians, double KcSS, double ALBearth, double ALBlake, std::string Datarate, double PWF, double Timestep);
void func_FluxFin (NumericMatrix Flux, double CA, NumericMatrix Lake, long Time, NumericMatrix Isotopes, double Cin, double Timestep);
void func_FluxFsf (NumericMatrix Flux, double CA, NumericMatrix Lake, long Time, NumericMatrix Met, double Timestep);
void func_FluxFgw (NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakesedimentgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double Timestep, double ETexdepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, long Time, NumericMatrix Lake, NumericMatrix GWcatchmentgrid, NumericMatrix Lakevolumes, NumericMatrix GWlakekgrid, NumericMatrix Flux, NumericMatrix LakeDepths, double Csr, NumericMatrix Met, NumericMatrix Isotopes, std::string GroundwaterModule, List GWStats, NumericMatrix GWNumbers, int GWLakeID);
void func_FluxFof (NumericMatrix Flux, NumericMatrix Lake, NumericMatrix Lakevolumes, long Time);
void func_FluxE (NumericMatrix Flux, long Time, NumericMatrix Met, double ALBlake, double PWF, double LatitudeRadians, std::string Datarate, double Timestep, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, double Feedback);
void func_FluxPET (NumericMatrix Flux, long Time, NumericMatrix Met, double ALBearth, double PWF, double LatitudeRadians, std::string Datarate, double Timestep);
void func_Fluxes(double KcSS, double ALBearth, double ALBlake, double CA, double Cin, double Csr, double KcDS, double KcLsum, double KcLwin, double LatitudeRadians, float RunoffRatio, double PWF, double Timestep, float AWCds, float AWCss, long Time, NumericMatrix Flux, NumericMatrix Isotopes, NumericMatrix Lake, NumericMatrix Met, NumericMatrix Lakevolumes, std::string Datarate, std::string Interpolationtype, std::string GroundwaterModule, NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakesedimentgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double ETexdepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, NumericMatrix GWcatchmentgrid, NumericMatrix GWlakekgrid, NumericMatrix LakeDepths, List GWStats, NumericMatrix GWNumbers, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, int GWLakeID, double Feedback);
void func_Hydrology(NumericMatrix Flux, NumericMatrix Lake, NumericMatrix Isotopes, long Time, NumericMatrix Met, NumericMatrix Lakevolumes, std::string Interpolationtype, std::string DirectExternalFlows, std::string GroundwaterModule, List GWStats);
void MixChems(NumericMatrix Lake, long Time, NumericMatrix Flux, NumericMatrix Chemistry, int StratType);
void func_IsotopeFluxes (NumericMatrix Flux, long Time, std::string Isotype, NumericMatrix Isotopes, NumericMatrix Met, double AtmosphericShift, double Turbulence, double Feedback, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype, std::string Delim, std::string DirectExternalFlows, NumericMatrix FluxD, NumericMatrix Flux18O, std::string GroundwaterModule, List GWStats, NumericMatrix GWNumbers, NumericMatrix ChemConcentrations, double Theta);
void func_ChemistryFluxes (NumericMatrix Chemistry, long Time, NumericMatrix Flux, NumericMatrix Lake, NumericMatrix Met, std::string DirectExternalFlows, std::string GroundwaterModule, List GWStats, NumericMatrix GWNumbers);
void func_LakeMixer(NumericMatrix Flux, NumericMatrix Lake, NumericMatrix Chemistry, NumericMatrix Isotopes, long Time, NumericMatrix Met, NumericMatrix Lakevolumes, std::string Interpolationtype, NumericMatrix FluxD, NumericMatrix Flux18O);
void MassBalance (double KcSS, double ALBearth, double ALBlake, double CA, double Cin, double Csr, double KcDS, double KcLsum, double KcLwin, double LatitudeRadians, float RunoffRatio, double PWF, double Timestep, float AWCds, float AWCss, long Time, NumericMatrix Flux, NumericMatrix Isotopes, NumericMatrix Lake, NumericMatrix Met, NumericMatrix Lakevolumes, std::string Datarate, std::string Interpolationtype, std::string GroundwaterModule, NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakesedimentgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double ETexdepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, NumericMatrix GWcatchmentgrid, NumericMatrix GWlakekgrid, NumericMatrix LakeDepths, NumericMatrix Chemistry, std::string DirectExternalFlows, double AtmosphericShift, double Turbulence, double Feedback, NumericMatrix FluxD, NumericMatrix Flux18O, List GWStats, NumericMatrix GWNumbers, List ScenarioEventsList, DataFrame ScenarioEvents, NumericMatrix ChemConcentrations, NumericMatrix ChemParams, double NeutralDragCoeff, double SWExtCoef, double ThermoclineTargetThickness, double ThermoclineMaximum, std::string StratificationModule, NumericMatrix Stratification, int GWLakeID, double Theta);
void EventChecker(long Time, Function FindEvents);
void PointTest(int PointVar, Function RChangePointer);

#endif // MASSBAL_H