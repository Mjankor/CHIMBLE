#ifndef ISOCALCS
#define ISOCALCS

using namespace Rcpp;

double func_dE(long Time, std::string Isotype, NumericMatrix Isotopes, NumericMatrix Met, double AtmosphericShift, double Turbulence, double Feedback, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype, NumericMatrix ChemConcentrations, double Theta);
double RH_Normalised(NumericMatrix Met, long Time, double Feedback, double Act_Water, double tLake);
#endif // ISOCALCS