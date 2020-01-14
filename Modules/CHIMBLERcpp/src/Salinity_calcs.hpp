#ifndef SALINITYCALCS_H
#define SALINITYCALCS_H
using namespace Rcpp;

double CalcMolarities(double Time, NumericMatrix Chemistry, NumericMatrix ChemConcentrations);
void EstimateSaltMolarities(NumericMatrix ChemConcentrations);
void CalcSalts(NumericMatrix ChemConcentrations);
void CalcSaltDensity (NumericMatrix ChemConcentrations, NumericMatrix ChemParams);
void CalcWaterActivity (NumericMatrix ChemConcentrations, NumericMatrix ChemParams);
void SalinityCalcs (double Time, NumericMatrix ChemConcentrations, NumericMatrix Chemistry, NumericMatrix ChemParams);


#endif // SALINITYCALCS_H
