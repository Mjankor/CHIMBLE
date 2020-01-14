#ifndef HYDROCALCS
#define HYDROCALCS

using namespace Rcpp;

double func_CAe (double CA, NumericMatrix Lake, long Time) ;
NumericVector func_FpLake(long Time, NumericMatrix Met, NumericMatrix Lake, double Timestep);
NumericVector func_Fsm (double CA, long Time, NumericMatrix Lake,NumericMatrix Met, NumericMatrix Isotopes, std::string Datarate, double Timestep);
NumericVector func_Fext (long Time, NumericMatrix Met, NumericMatrix Isotopes, double Timestep);
NumericVector func_Fgw (NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakesedimentgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double Timestep, double ETexdepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, long Time, NumericMatrix Lake, NumericMatrix GWcatchmentgrid, NumericMatrix Lakevolumes, NumericMatrix GWlakekgrid, NumericMatrix Flux, NumericMatrix LakeDepths, double Csr, NumericMatrix Met, NumericMatrix Isotopes, std::string GroundwaterModule, List GWStats, NumericMatrix GWNumbers, int GWLakeID);
NumericVector func_FpCatchment(long Time, NumericMatrix Met, double CA, NumericMatrix Lake, double Timestep);
int CDaysInMonth(int Month, int Year);
int func_DayofYear(int Month, int Day);
double func_MonthFraction(long Time, NumericMatrix Met, std::string Datarate, double Timestep);
double func_Ra(long Time, double LatitudeRadians, NumericMatrix Met);
double func_E(long Time, NumericMatrix Met, double ALBlake, double PWF, double LatitudeRadians, std::string Datarate, double Timestep, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, double Feedback);
double func_PET(long Time, NumericMatrix Met, double ALBearth, double PWF, double LatitudeRadians, std::string Datarate, double Timestep);
NumericVector func_Evap(double CA, long Time, NumericMatrix Met, NumericMatrix Lake,NumericMatrix Isotopes, float AWCss, float AWCds, double KcDS, double LatitudeRadians, double KcSS, double ALBearth, double ALBlake, std::string Datarate, double PWF, double Timestep);
NumericVector func_Soil(long Time, float AWCss, float AWCds, double CA, NumericMatrix Lake, NumericMatrix Met, NumericMatrix Isotopes, double KcDS, double  KcSS, double LatitudeRadians,float RunoffRatio, double ALBearth, double ALBlake, double PWF, std::string Datarate, double Timestep);
NumericVector func_Fsf(double CA, NumericMatrix Lake, long Time, NumericMatrix Met, double Timestep);
NumericVector func_Fin(double CA, NumericMatrix Lake, long Time, NumericMatrix Isotopes, double Cin, double Timestep);
double func_SVC(long Time, NumericMatrix Met, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype);
NumericVector func_Fdlm(long Time, NumericMatrix Isotopes, NumericMatrix Met, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype);
double func_KcLake(long Time, double KcLwin, double KcLsum, NumericMatrix Met,std::string Datarate, double Timestep);
NumericVector func_Fe(long Time, double Timestep, NumericMatrix Lake, NumericMatrix Met, double ALBlake, NumericMatrix Isotopes, double PWF, double LatitudeRadians, std::string Datarate, double KcLwin, double KcLsum, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, double Feedback);
double func_FeedbackRH(NumericMatrix Met, long Time, double Feedback);
double func_Fof(NumericMatrix Lake, NumericMatrix Lakevolumes, long Time);
//NumericVector func_Fsos(long Time, double Csr, NumericMatrix Lake, double Timestep, NumericMatrix Isotopes);
//NumericVector func_Fdos(long Time, double Csr, NumericMatrix Lake, double Timestep, NumericMatrix Isotopes);


#endif // HYDROCALCS