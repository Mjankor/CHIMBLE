#ifndef HYPOCALCS_H
#define HYPOCALCS_H
using namespace Rcpp;

double func_Lakearea(double Volume, NumericMatrix Lakevolumes);
double func_Lakedepth(double Volume, NumericMatrix Lakevolumes);
double func_Lakevolume(double Depth, NumericMatrix Lakevolumes);

#endif // HYPOCALCS_H
