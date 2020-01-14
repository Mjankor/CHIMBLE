#ifndef ISOCALCS
#define ISOCALCS

using namespace Rcpp;

double func_dE(long Time, std::string Isotype, NumericMatrix Isotopes, NumericMatrix Met, double AtmosphericShift, double Turbulence, double Feedback, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype, NumericMatrix ChemConcentrations);

#endif // ISOCALCS