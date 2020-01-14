#include <Rcpp.h>
#include "Hydrology_functions.hpp"
#include "Hypsographic_Calcs.hpp"
#include "Isotope_Calcs.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
double func_dE(long Time, std::string Isotype, NumericMatrix Isotopes, NumericMatrix Met, double AtmosphericShift, double Turbulence, double Feedback, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype, NumericMatrix ChemConcentrations, double Theta){
double tLake;
float CKinetic=0;
double Alpha_Eq = 0;
double RAlpha_Eq;
//double Esa; Moved to RH_Normalisation function
//double Esw;
double ThetaVal;
double RH_Norm;
double EKin;
double EEq;
double ETotal;
//double Limit;
double dAtmos;
double dE;
double Act_Water = 1;
double Act_Isotopes = 0;


tLake = Met(Time,4) + Met(Time,10);
if (Isotype == "18O"){
	CKinetic = 28.5;
	Alpha_Eq = exp(350410 / pow((273.15 + tLake),3) - 1666.4 / pow((273.15 + tLake),2) + 6.7123 /(273.15 + tLake) - 0.007685);
	
	// determine salt effect	
	//check upper/lower reservoir.
	//Additive method of determining isotopes comes from Gat 1995 (P156), Values from Gat 2010.
	if (ChemConcentrations(1,12) != 0){
		Act_Water = ChemConcentrations(1,13);
		//Calc delta values for each salt based on salt molality.
		//Combine for a total delta value.
		Act_Isotopes = -1.11 * ChemConcentrations(1,7) - 0.47 * ChemConcentrations(1,8) + 0.16 * ChemConcentrations(1,9);
		//Apply to lake water in the isotope equation as per Gat 2010.
		//and apply activity to humidity.
	} else {
		//and get values from lower reservoir if upper is 0
		Act_Water = ChemConcentrations(3,13);
		Act_Isotopes = -1.11 * ChemConcentrations(3,7) - 0.47 * ChemConcentrations(3,8) + 0.16 * ChemConcentrations(3,9);
	}	
	
	
} else if (Isotype == "D"){
	CKinetic = 25.1;
	Alpha_Eq = exp(1.1588 * pow((273.15 + tLake),3)/pow(10,9) - 1.6201 * pow((273.15 + tLake),2)/pow(10,6) + 0.79484 * (273.15+tLake)/1000 + 2.9992 * pow(10,6)/pow((273.15 + tLake),3) - 0.16104);
	
	//calculate salt effect.
	if (ChemConcentrations(1,12) != 0){
		Act_Water = ChemConcentrations(1,13);
		//Calc delta values for each salt based on salt molality.
		//Combine for a total delta value.
		Act_Isotopes = 2.2 * ChemConcentrations(1,6) + 5.1 * ChemConcentrations(1,7) + 6.1 * ChemConcentrations(1,8) + 2.5 * ChemConcentrations(1,9);
		//Apply to lake water in the isotope equation as per Gat 2010.
		//and apply activity to humidity.
	} else {
		//and get values from lower reservoir if upper is 0
		Act_Water = ChemConcentrations(3,13);
		Act_Isotopes = 2.2 * ChemConcentrations(3,6) + 5.1 * ChemConcentrations(3,7) + 6.1 * ChemConcentrations(3,8) + 2.5 * ChemConcentrations(3,9);
	}	
	
}

RAlpha_Eq = 1/Alpha_Eq;
RH_Norm = RH_Normalised(Met, Time, Feedback, Act_Water, tLake);
//Rcpp::Rcout << "RHNorm" << RH_Norm << std::endl;
//Theta is defined by method from Gat 2010, or by a predefined parameter.
if (Theta == 999){
	if (Met(Time,5)==100){
		ThetaVal = 1;
		} else {
		ThetaVal = (1- func_FeedbackRH(Met, Time, Feedback))/(1- func_FeedbackRH(Met, Time, 0));
		} 
	} else {
	ThetaVal = Theta;
}

EKin = CKinetic * (1-RH_Norm)*ThetaVal * Turbulence;
EEq = 1000*(1-RAlpha_Eq);
ETotal = EKin + EEq;

if (Isotype == "18O"){
	dAtmos = Met(Time,8) * RAlpha_Eq - EEq * AtmosphericShift;
	if (func_SVC (Time, Met, Lakevolumes, Lake, Interpolationtype) == 0){
		dE = (RAlpha_Eq * Isotopes(Time,7) + Act_Isotopes - RH_Norm*dAtmos - ETotal)/(1 - RH_Norm + 0.001*EKin);
	} else {
		dE = (RAlpha_Eq * Isotopes(Time,5) + Act_Isotopes - RH_Norm*dAtmos - ETotal)/(1 - RH_Norm + 0.001*EKin);
	}
} else if (Isotype == "D"){
	dAtmos = Met(Time,9) * RAlpha_Eq - EEq * AtmosphericShift;
	if (func_SVC (Time, Met, Lakevolumes, Lake, Interpolationtype) == 0){
		dE = (RAlpha_Eq * Isotopes(Time,6) + Act_Isotopes - RH_Norm*dAtmos - ETotal)/(1 - RH_Norm + 0.001*EKin);
	} else {
		dE = (RAlpha_Eq * Isotopes(Time,4) + Act_Isotopes - RH_Norm*dAtmos - ETotal)/(1 - RH_Norm + 0.001*EKin);
	}
}
//Limit = (RH_Norm * dAtmos + ETotal)/(RH_Norm - 0.001*ETotal);
//rerun the isotope calculation if feedback exists.
if (Feedback != 0){
	//Rcpp::Rcout << "Running Feedback" << std::endl;
	if (Isotype == "18O"){
		dAtmos = (Met(Time,8) * RAlpha_Eq - EEq * AtmosphericShift) * (1-Feedback) + Feedback * dE;
		if (func_SVC (Time, Met, Lakevolumes, Lake, Interpolationtype) == 0){
			dE = (RAlpha_Eq * Isotopes(Time,7) + Act_Isotopes - RH_Norm*dAtmos - ETotal)/(1 - RH_Norm + 0.001*EKin);
		} else {
			dE = (RAlpha_Eq * Isotopes(Time,5) + Act_Isotopes - RH_Norm*dAtmos - ETotal)/(1 - RH_Norm + 0.001*EKin);
		}
	} else if (Isotype == "D"){
		dAtmos = (Met(Time,9) * RAlpha_Eq - EEq * AtmosphericShift) * (1-Feedback) + Feedback * dE;
		if (func_SVC (Time, Met, Lakevolumes, Lake, Interpolationtype) == 0){
			dE = (RAlpha_Eq * Isotopes(Time,6) + Act_Isotopes  - RH_Norm*dAtmos - ETotal)/(1 - RH_Norm + 0.001*EKin);
		} else {
			dE = (RAlpha_Eq * Isotopes(Time,4) + Act_Isotopes - RH_Norm*dAtmos - ETotal)/(1 - RH_Norm + 0.001*EKin);
		}
	}
}

return dE;
}

// [[Rcpp::export]]
double RH_Normalised(NumericMatrix Met, long Time, double Feedback, double Act_Water, double tLake){
double Esa = 6.108 * exp(17.27 * Met(Time,4)/(Met(Time,4)+237.7));
double Esw = 6.108 * exp(17.27 * tLake/(tLake+237.7));
double RH_Norm;
	RH_Norm = func_FeedbackRH(Met, Time, Feedback) * Esa/Esw; //change 0 to Feedback to add feedback to evap.
	RH_Norm = RH_Norm / Act_Water; // Apply activity of water.
	if (RH_Norm > 0.95){ //Keep RH_norm below 0.95 to stop the fractionation equation breaking.
	RH_Norm = 0.95;
	}
	return RH_Norm;
}