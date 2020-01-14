#include <Rcpp.h>
#include "Hydrology_functions.hpp"
#include "Hypsographic_Calcs.hpp"
#include "ChimbleGW.hpp"
#include "Stratification.hpp"
using namespace Rcpp;

long double const pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899863;

// [[Rcpp::export]]
double func_CAe (double CA, NumericMatrix Lake, long Time) {
double CAe = CA - Lake(Time,6);
return CAe;
}
//daily checked

// [[Rcpp::export]]
NumericVector func_FpLake(long Time, NumericMatrix Met, NumericMatrix Lake, double Timestep){
double Fp;
double Fp18O;
double FpD;
if (Met(Time,3) > 0){
	Fp = Met(Time,3)*0.001*Lake(Time,6)*Timestep;
	Fp18O = Fp*Met(Time,8);
	FpD = Fp*Met(Time,9);
} else {
	Fp = 0;
	Fp18O = 0;
	FpD = 0;
}
NumericVector FpVals = NumericVector::create(Named("Fp") = Fp, Named("Fp18O") = Fp18O, Named("FpD") = FpD);
	return FpVals;
}
//daily checked

// [[Rcpp::export]]
NumericVector func_Fsm (double CA, long Time, NumericMatrix Lake,NumericMatrix Met, NumericMatrix Isotopes, std::string Datarate, double Timestep){
double snowmelt;
double Fsm;
double Fsm18O;
double FsmD;
if (Datarate == "Monthly"){
	snowmelt = 0.021; 
	} else {
	snowmelt = 0.00069; 
	}

if (Lake(Time, 11) > snowmelt * (Met(Time, 4)+2) * func_CAe(CA, Lake, Time)*Timestep && Met(Time, 4) > -2){
	Fsm = snowmelt * (Met(Time, 4)+2) * func_CAe(CA, Lake, Time)*Timestep;
	Fsm18O = Fsm * Isotopes(Time, 13); 
	FsmD = Fsm * Isotopes(Time, 12);
}else if (Lake(Time, 11) <= snowmelt * (Met(Time, 4)+2) * func_CAe(CA, Lake, Time)*Timestep && Met(Time, 4) > -2) {
	Fsm = Lake(Time, 11);
	Fsm18O = Fsm * Isotopes(Time, 13); 
	FsmD = Fsm * Isotopes(Time, 12);
} else {
	Fsm = 0;
	Fsm18O = 0; 
	FsmD = 0;
}
NumericVector FsmVals = NumericVector::create(Named("Fsm") = Fsm, Named("Fsm18O") = Fsm18O, Named("FsmD") = FsmD);
return FsmVals;
}
//daily checked

// [[Rcpp::export]]
NumericVector func_Fext (long Time, NumericMatrix Met, NumericMatrix Isotopes, double Timestep){
double Fext;
double Fext18O;
double FextD;
	Fext = Met(Time,11) * Timestep;
	Fext18O = Fext * Isotopes(Time, 17);
	FextD = Fext * Isotopes(Time, 16);
	NumericVector FextVals = NumericVector::create(Named("Fext") = Fext, Named("Fext18O") = Fext18O, Named("FextD") = FextD);
	return FextVals;
}
//daily checked

// [[Rcpp::export]]
NumericVector func_Fgw (NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakesedimentgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double Timestep, double ETexdepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, long Time, NumericMatrix Lake, NumericMatrix GWcatchmentgrid, NumericMatrix Lakevolumes, NumericMatrix GWlakekgrid, NumericMatrix Flux, NumericMatrix LakeDepths, double Csr, NumericMatrix Met, NumericMatrix Isotopes, std::string GroundwaterModule, List GWStats, NumericMatrix GWNumbers, int GWLakeID){
double Fgw = 0;
double Fgw18O;
double FgwD;
double FgwSos;
double FgwDos;
int StartDay=0;
int EndDay=0;
NumericMatrix FgwMatrix;
if (GroundwaterModule == "No"){
		Fgw = Met(Time,16)* Timestep;
		Fgw18O = Fgw*Met(Time,17);
		FgwD = Fgw*Met(Time,18);
		FgwSos = Csr*Lake(Time,7)*Timestep;
		FgwDos = Csr*Lake(Time,8)*Timestep;
	} else {
		if (Time == Met.nrow()-1){
			StartDay = func_DayofYear(Met(Time-1,1), Met(Time-1,2));
			EndDay = func_DayofYear(Met(Time,1), Met(Time,2));
			} else {
			StartDay = func_DayofYear(Met(Time,1), Met(Time,2));
			EndDay = func_DayofYear(Met(Time+1,1), Met(Time+1,2));
		}
		//updated to avoid error on last flux calc step. 10/09/19
		//int StartDay = func_DayofYear(Met(Time,1), Met(Time,2));
		//int EndDay = func_DayofYear(Met(Time+1,1), Met(Time+1,2));
		if (StartDay > EndDay){ //end of year
		EndDay += 365;
		}
	//Rcpp::Rcout << "Start and End " << StartDay << ", " <<EndDay << std::endl;
		double GWRuntime = EndDay - StartDay;
		double GWTimestep = GWRuntime;
		//Call the groundwater function.
		//Note. Check that it is using 2 matrixes, rather than just updating the heads matrix.
		//Prep the GW model
		func_LakeDepths(Lake, Time, LakeDepths);
		func_ClassGrid(Time, Topogrid, LakeDepths, GWcatchmentgrid, Classgrid, Lakevolumes);
		func_RechargeGrid (Time, Cellxdim, Cellydim, GWRechargegrid, Classgrid, Flux, Timestep);
		func_ETGrid (Time, GWETgrid, Flux); 
		//Lake conductance is not required, but it's efficient. We'll leave it in as we may want to update the sediment thickness as the lake evolves.
		func_LakeConductance(Lakecondgrid, Cellxdim, Cellydim, GWlakesedimentgrid, GWlakekgrid);
		
		//Rcpp::Rcout << "Topogrid: " << Topogrid(0,1) << std::endl;
		//Rcpp::Rcout << "LakeDepths: " << LakeDepths(0,3) << std::endl;
		//Rcpp::Rcout << "GWRechargegrid: " << GWRechargegrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWPumpinggrid: " << GWPumpinggrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWKgrid: " << GWKgrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWSgrid: " << GWSgrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWbasegrid: " << GWbasegrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWheadgrid: " << GWheadgrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWboundarygrid: " << GWboundarygrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWlakesedimentgrid: " << GWlakesedimentgrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWETgrid: " << GWETgrid(0,1) << std::endl;
		//Rcpp::Rcout << "Classgrid: " << Classgrid(0,1) << std::endl;
		//Rcpp::Rcout << "GWTimestep: " << GWTimestep << std::endl;
		//Rcpp::Rcout << "GWRuntime: " << GWRuntime << std::endl;
		//Run the groundwater Module.
		List FgwList = Headcalc(Topogrid, GWRechargegrid, GWPumpinggrid, GWKgrid, GWSgrid, GWbasegrid, GWheadgrid, GWboundarygrid, GWlakesedimentgrid, Lakecondgrid, GWETgrid, Classgrid, GWTimestep, GWRuntime, ETexdepth, LakeDepths, Cellxdim, Cellydim, Accuracy, Maxiterations, GWNumbers);
		NumericMatrix FgwMatrix = FgwList(0);
		NumericVector Stats = FgwList(1);
		//Rcpp::Rcout << "GW Flux" << - FgwMatrix(0,1) << std::endl;
		//Rcpp::Rcout << "Turning all back on" << std::endl;
		int ActiveLake = GWLakeID - 1;
		Fgw = -FgwMatrix(ActiveLake,1) + FgwMatrix(ActiveLake,3) ; //make positive (Surface and in-seepage)
		Fgw18O = Fgw * Isotopes(Time, 17);
		FgwD = Fgw * Isotopes(Time, 16);
		double SurfaceBias = pow((1 - Met(Time,15)/Lake(Time,4)),2);
		FgwSos = FgwMatrix(ActiveLake,2) * (1- SurfaceBias); 
		FgwDos = FgwMatrix(ActiveLake,2) * SurfaceBias;
		//dump the stats
		for (int i = 0; i < 14; i++){
		GWStats(i) = Stats(i);
		}
	}
	NumericVector FgwVals = NumericVector::create(Named("Fgw") = Fgw, Named("Fgw18O") = Fgw18O, Named("FgwD") = FgwD, Named("FgwSos") = FgwSos, Named("FgwDos") = FgwDos);
	return FgwVals;
}


// [[Rcpp::export]]
NumericVector func_FpCatchment(long Time, NumericMatrix Met, double CA, NumericMatrix Lake, double Timestep){
double Fr;
double Fr18O;
double FrD;
if (Met(Time,3) > 0){
	Fr = Met(Time,3)*0.001*func_CAe(CA, Lake, Time)*Timestep;
	Fr18O = Fr*Met(Time,8);
	FrD = Fr*Met(Time,9);
} else {
	Fr = 0;
	Fr18O = 0;
	FrD = 0;
}
NumericVector FpCatchmentVals = NumericVector::create(Named("Fr") = Fr, Named("Fr18O") = Fr18O, Named("FrD") = FrD);
	return FpCatchmentVals;
}
//daily checked

// [[Rcpp::export]]
int func_DayofYear(int Month, int Day){
// ignores leap years. Shouldn't matter for most applications.
int Days;
switch (Month){
	case 1: Days = 0; break;
	case 2 : Days = 31;break;
	case 3 : Days = 59;break;
	case 4 : Days = 90;break;
    case 5 : Days = 120;break;
    case 6 : Days = 151;break;
    case 7 : Days = 181;break;
    case 8 : Days = 212;break;
    case 9 : Days = 243;break;
    case 10 : Days = 273;break;
    case 11 : Days = 304;break;
    case 12 : Days = 334;break;
    }
Days = Days + Day;
return Days;
}


// [[Rcpp::export]]
int CDaysInMonth(int Month, int Year){
int leap = 0;
int Days;
if (Year % 4 == 0 && (Year % 100 != 0 || Year % 400 == 0)){
leap = 1;
}
switch (Month){
	case 1: Days = 31; break;
	case 2 : Days = 28 + leap;break;
	case 3 : Days = 31;break;
	case 4 : Days = 30;break;
    case 5 : Days = 31;break;
    case 6 : Days = 30;break;
    case 7 : Days = 31;break;
    case 8 : Days = 31;break;
    case 9 : Days = 30;break;
    case 10 : Days = 31;break;
    case 11 : Days = 30;break;
    case 12 : Days = 31;break;
    }
return Days;
}

// [[Rcpp::export]]
double func_MonthFraction(long Time, NumericMatrix Met, std::string Datarate, double Timestep){
double Monthfrac;
if (Datarate == "Monthly"){
	CharacterVector Rows = rownames(Met);
	float CurrTimestep = atof(Rows[Time]);
	Monthfrac = Met(Time,1) + fmod(CurrTimestep, 1) * Timestep *10; 
} else {
Monthfrac = Met(Time,1) + Met(Time,2)/30; 
}
return Monthfrac;
}
//Daily checked

// [[Rcpp::export]]
double func_Ra(long Time, double LatitudeRadians, NumericMatrix Met){
	//There is also an error in equation 2 (118 * N * etc). This yields unrealistic values for tropical regions. Therefore func.Ra has been rewritten (below) to use the FAO calculations - Allen 1998 FAO Chapter 3. eq 28. Previous version of function included below.
/*double func_Ra(long Time, double LatitudeRadians, NumericMatrix Met, std::string Datarate, double Timestep){
double Ra;
double N;
N = 4 * LatitudeRadians * sin(0.53 * func_MonthFraction(Time, Met, Datarate, Timestep)
 - 1.65) + 12;
if (fabs(LatitudeRadians) > 23.5*pi/180){ 
Ra = 3 * N * sin(0.131*N - 0.95 * fabs(LatitudeRadians));
} else if (fabs(LatitudeRadians) <= 23.5 * pi/180){
Ra = 118*(pow(N,0.2))*sin(0.131*N - (0.2 * fabs(LatitudeRadians)));
} else {
Ra = 999;
}
return Ra;
}
//Daily checked
*/

//Values for Ra for different latitudes are given in Annex 2 (Table 2.6). These values represent Ra on the 15th day of each month. These values deviate from values that are averaged over each day of the month by less than 1% for all latitudes during non-frozen periods and are included for simplicity of calculation. These values deviate slightly from the values in the Smithsonian Tables. For the winter months in latitudes greater than 55° (N or S), the equations for Ra have limited validity. Reference should be made to the Smithsonian Tables to assess possible deviations. - FAO - Allen 1998.
//double checked against https://inis.iaea.org/search/search.aspx?orig_q=RN:15063848
int DayofYear;
double Dr;
double SolarDec;
double Sunset;
double Ra;
	DayofYear = func_DayofYear(Met(Time,1), Met(Time,2));
	Dr = 1+0.033*cos(2*pi*DayofYear/365);
	SolarDec = 0.409*sin(2*pi*DayofYear/365-1.39);
	Sunset = acos(-tan(LatitudeRadians)*tan(SolarDec));
	Ra = 24*60/pi * 0.082 * Dr * (Sunset * sin(LatitudeRadians) * sin(SolarDec)+cos(LatitudeRadians)*cos(SolarDec)*sin(Sunset));
return Ra;
}


//double func_ClearSky(long Time, double LatitudeRadians, NumericMatrix Met){
//Update: This function may not be needed. Instead we'll use the clearness index to account for cloudiness.
//We use the Haurwitz model, which is reported to yield good values with only the need for zenith angle.
//Haurwitz, B., 1945: INSOLATION IN RELATION TO CLOUDINESS AND CLOUD DENSITY. J. Meteor., 2, 154–166, https://doi.org/10.1175/1520-0469(1945)002<0154:IIRTCA>2.0.CO;2 
// see https://github.com/sandialabs/MATLAB_PV_LIB/blob/master/pvl_clearsky_haurwitz.m for details on unit conversion.
//for model validation
//Badescu, V. (1998). Verification of some very simple clear and cloudy sky models to evaluate global solar irradiance. Solar Energy, 61(4), 251-264.
//Reno, M.J., Hansen, C.W. and Stein, J.S., 2012. Global horizontal irradiance clear sky models: Implementation and analysis. SANDIA report SAND2012-2389.
//}



// [[Rcpp::export]]
double func_E(long Time, NumericMatrix Met, double ALBlake, double PWF, double LatitudeRadians, std::string Datarate, double Timestep, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, double Feedback){
double E = 0;
//calculates the evaporation over the timestep as a monthly or daily rate. M/day, or M/month.
// Get activity of water to account for salinity.
double Act_Water;
if (ChemConcentrations(1,12) != 0){
	Act_Water = ChemConcentrations(1,13);
	//Rcpp::Rcout << "Using SL_Act_Water value of " << Act_Water << std::endl;
} else {
	Act_Water = ChemConcentrations(3,13);
	//Rcpp::Rcout << "Using DL_Act_Water value of " << Act_Water << std::endl;
}

if (StratificationModule == "Yes"){
	NumericVector Heat = func_LatSensHeat (Met(Time,4) + Met(Time,10), Met(Time,4), Met(Time,7), 0, 0, Met(Time,14), func_FeedbackRH(Met, Time, Feedback), NeutralDragCoeff);
	E = Heat(2); //mm per day
	E = E*Act_Water;
} else {
	if (Met(Time, 4) <= 0) {
		E = 0;
	} else {
		E = 0.051 * (1 - ALBlake) * Met(Time,6) *pow((Met(Time, 4) + 9.5),0.5) - (2.4 * pow(Met(Time,6)/func_Ra(Time, LatitudeRadians, Met),2)) + 0.052 * (Met(Time, 4) + 20) * (1 - func_FeedbackRH(Met, Time, 0)/Act_Water) * (PWF - 0.38 + 0.54*Met(Time,7)); //change 0 to Feedback to add feedback to RH
	}
}
E = E * 0.001; //M per day
if (Datarate == "Monthly"){
	E = E * CDaysInMonth(Met(Time,1), Met(Time,0));
}
	//Rcpp::Rcout << "CDaysInMonth(Met(Time,1), Met(Time,0))" << CDaysInMonth(Met(Time,1), Met(Time,0)) << std::endl;
return E;
}
//Daily checked

// [[Rcpp::export]]
double func_PET(long Time, NumericMatrix Met, double ALBearth, double PWF, double LatitudeRadians, std::string Datarate, double Timestep){
double PET;
//Activity of water not used in this calculation
if (Met(Time, 4) <= 0) {
PET = 0;
} else {
PET = 0.051 * (1 - ALBearth)*  Met(Time,6) *pow((Met(Time, 4) + 9.5),0.5) - (2.4 * pow(Met(Time,6)/func_Ra(Time, LatitudeRadians, Met),2)) + 0.048 * (Met(Time, 4) + 20) * (1 - func_FeedbackRH(Met, Time, 0)) * (0.5 + 0.536*Met(Time,7));//change 0 to Feedback to add feedback to RH
}
PET = PET * 0.001; 
if (Datarate == "Monthly"){
	PET = PET * CDaysInMonth(Met(Time,1), Met(Time,0));
}
return PET;
}
//Daily checked

// [[Rcpp::export]]
NumericVector func_Evap(double CA, long Time, NumericMatrix Met, NumericMatrix Lake,NumericMatrix Isotopes, float AWCss, float AWCds, double KcDS, double LatitudeRadians, double KcSS, double ALBearth, double ALBlake, std::string Datarate, double PWF, double Timestep){
double SSETC;
double DSETC;
double SSETa;
double DSETa;
double Fsse;
double Fdse;
double Fsse18O;
double Fdse18O;
double FsseD;
double FdseD;

SSETC = Lake(Time,9)/(func_CAe(CA, Lake, Time) * AWCss);
DSETC = Lake(Time,10)/(func_CAe(CA, Lake, Time) * AWCds);
SSETa = func_CAe(CA, Lake, Time) * func_PET(Time, Met, ALBearth, PWF, LatitudeRadians, Datarate, Timestep) * SSETC * KcSS;
DSETa = func_CAe(CA, Lake, Time) * func_PET(Time, Met, ALBearth, PWF, LatitudeRadians, Datarate, Timestep) * DSETC * KcDS;

if (Lake(Time,9) > SSETa * Timestep){
	Fsse = SSETa * Timestep;
	} else {
	Fsse = Lake(Time,9);
}
Fsse18O = Fsse * Isotopes(Time,9);
FsseD = Fsse * Isotopes(Time,8);

if (Lake(Time,10) > DSETa * Timestep){ 
	Fdse = DSETa * Timestep;
	} else {
	Fdse = Lake(Time,10);
}
Fdse18O = Fdse * Isotopes(Time,11);
FdseD = Fdse * Isotopes(Time,10);

NumericVector EvapVals = NumericVector::create(Named("Fsse") = Fsse, Named("Fdse") = Fdse, Named("Fsse18O") = Fsse18O, Named("FsseD") = FsseD, Named("Fdse18O") = Fdse18O, Named("FdseD") = FdseD);

return EvapVals;
}
//Daily checked

// [[Rcpp::export]]
NumericVector func_Soil(long Time, float AWCss, float AWCds, double CA, NumericMatrix Lake, NumericMatrix Met, NumericMatrix Isotopes, double KcDS, double  KcSS, double LatitudeRadians,float RunoffRatio, double ALBearth, double ALBlake, double PWF, std::string Datarate, double Timestep){
	double Fssi;
	double Fssi18O;
	double FssiD;
	double Fssd;
	double Fssd18O;
	double FssdD;
	double Fdsd;
	double Fdsd18O;
	double FdsdD;
	double Fro;
	double Fro18O;
	double FroD;
	double SSmax = func_CAe(CA, Lake, Time)* AWCss;
	double DSmax = func_CAe(CA, Lake, Time)* AWCds;
	
	double Influx = func_FpCatchment(Time, Met, CA, Lake, Timestep)[0] + func_Fsm(CA,Time,Lake,Met,Isotopes,Datarate,Timestep)[0];
	double SSflux = Lake(Time, 9) + Influx - func_Evap(CA, Time, Met, Lake, Isotopes, AWCss, AWCds, KcDS, LatitudeRadians, KcSS, ALBearth, ALBlake, Datarate, PWF,Timestep)[0];
	double SSexcess = SSflux - SSmax;
	double DSflux = Lake(Time,10) - func_Evap(CA, Time, Met, Lake, Isotopes, AWCss, AWCds, KcDS, LatitudeRadians, KcSS, ALBearth, ALBlake, Datarate, PWF,Timestep)[1];
	
double Influx18O = func_FpCatchment(Time, Met, CA, Lake, Timestep)[1] + func_Fsm(CA,Time,Lake,Met,Isotopes,Datarate,Timestep)[1];
double InfluxD = func_FpCatchment(Time, Met, CA, Lake, Timestep)[2] + func_Fsm(CA,Time,Lake,Met,Isotopes,Datarate,Timestep)[2];

if (Influx == 0) { 
	Fssi = 0;
	Fssi18O = 0;
	FssiD = 0;
	Fssd = 0;
	Fssd18O = 0;
	FssdD = 0;
	Fdsd = 0;
	Fdsd18O = 0;
	FdsdD = 0;
	Fro = 0;
	Fro18O = 0;
	FroD = 0;
	} else if (SSmax >= SSflux){ 
	Fssi = Influx;
	Fssi18O = Influx18O;
	FssiD = InfluxD;
	Fssd = 0;
	Fssd18O = 0;
	FssdD = 0;
	Fdsd = 0;
	Fdsd18O = 0;
	FdsdD = 0;
	Fro = 0;
	Fro18O = 0;
	FroD = 0;
	} else if (SSmax < SSflux && DSmax >= DSflux + SSexcess*RunoffRatio){
	Fssi = SSmax - Lake(Time,9) + func_Evap(CA, Time, Met, Lake, Isotopes, AWCss, AWCds, KcDS, LatitudeRadians, KcSS, ALBearth, ALBlake, Datarate, PWF,Timestep)[0];
	Fssd = SSexcess*RunoffRatio;
	Fdsd = 0;
	Fro = SSexcess*(1-RunoffRatio);
	Fssi18O = Influx18O * Fssi/(Fssi+Fssd+Fdsd+Fro);
	Fssd18O = Influx18O * Fssd/(Fssi+Fssd+Fdsd+Fro);
	Fdsd18O = Influx18O * Fdsd/(Fssi+Fssd+Fdsd+Fro);
	Fro18O = Influx18O * Fro/(Fssi+Fssd+Fdsd+Fro);
	FssiD = InfluxD * Fssi/(Fssi+Fssd+Fdsd+Fro);
	FssdD = InfluxD * Fssd/(Fssi+Fssd+Fdsd+Fro);
	FdsdD = InfluxD * Fdsd/(Fssi+Fssd+Fdsd+Fro);
	FroD = InfluxD * Fro/(Fssi+Fssd+Fdsd+Fro);
	} else if (SSmax < SSflux && DSmax < (DSflux + SSexcess*RunoffRatio)){
	Fssi = SSmax - Lake(Time,9) + func_Evap(CA, Time, Met, Lake, Isotopes, AWCss, AWCds, KcDS, LatitudeRadians, KcSS, ALBearth, ALBlake, Datarate, PWF,Timestep)[0];
	Fssd = DSmax - Lake(Time,10) + func_Evap(CA, Time, Met, Lake, Isotopes, AWCss, AWCds, KcDS, LatitudeRadians, KcSS, ALBearth, ALBlake, Datarate, PWF,Timestep)[1];
	Fdsd = SSexcess*RunoffRatio - Fssd;
	Fro = SSexcess*(1-RunoffRatio);
	Fssi18O = Influx18O * Fssi/(Fssi+Fssd+Fdsd+Fro);
	Fssd18O = Influx18O * Fssd/(Fssi+Fssd+Fdsd+Fro);
	Fdsd18O = Influx18O * Fdsd/(Fssi+Fssd+Fdsd+Fro);
	Fro18O = Influx18O * Fro/(Fssi+Fssd+Fdsd+Fro);
	FssiD = InfluxD * Fssi/(Fssi+Fssd+Fdsd+Fro);
	FssdD = InfluxD * Fssd/(Fssi+Fssd+Fdsd+Fro);
	FdsdD = InfluxD * Fdsd/(Fssi+Fssd+Fdsd+Fro);
	FroD = InfluxD * Fro/(Fssi+Fssd+Fdsd+Fro);
	} else {
	Rcpp::stop("Error in func_Soil conditions");
	}
	
	NumericVector Soil = NumericVector::create(Named("Fssi") = Fssi,Named("Fssd") = Fssd, Named("Fdsd") = Fdsd, Named("Fro")= Fro, Named("Fssi18O") = Fssi18O, Named("FssiD") = FssiD, Named("Fssd18O") = Fssd18O, Named("FssdD") = FssdD, Named("Fdsd18O") = Fdsd18O, Named("FdsdD") = FdsdD, Named("Fro18O") = Fro18O, Named("FroD") = FroD);
	return (Soil);
}

// [[Rcpp::export]]
NumericVector func_Fsf(double CA, NumericMatrix Lake, long Time, NumericMatrix Met, double Timestep){
double Fsf;
double Fsf18O;
double FsfD;
if (Met(Time,4) <= 0){
	Fsf = Met(Time,3)*0.001*func_CAe(CA, Lake, Time)*Timestep;
	Fsf18O = Fsf * Met(Time,8);
	FsfD = Fsf * Met(Time,9);
} else {
	Fsf = 0;
	Fsf18O = 0;
	FsfD = 0;
}
NumericVector Snow = NumericVector::create(Named("Fsf") = Fsf, Named("Fsf18O") = Fsf18O, Named("FsfD") = FsfD);
return Snow;
}

// [[Rcpp::export]]
NumericVector func_Fin(double CA, NumericMatrix Lake, long Time, NumericMatrix Isotopes, double Cin, double Timestep){
double Fin;
double Fin18O;
double FinD;
if (Lake (Time,12) > 0.001*func_CAe(CA, Lake, Time)) { 
Fin = Lake (Time,12) * Cin * Timestep;
Fin18O = Fin * Isotopes(Time,15);
FinD = Fin * Isotopes(Time,14);
}else {
Fin = 0;
Fin18O = 0;
FinD = 0;
}
NumericVector Inflow = NumericVector::create(Named("Fin") = Fin, Named("Fin18O") = Fin18O, Named("FinD") = FinD);
return Inflow;
}

// [[Rcpp::export]]
double func_SVC(long Time, NumericMatrix Met, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype){
double Basemix;
double SVC;
double SLval;
//little jiggery pokery here. At the last timestep, we're outside the Met file, so we'll use the value from the previous timestep.
if (Time == Met.nrow()){
SLval = Met(Time-1,15);
} else {
SLval = Met(Time,15);
}

if (Interpolationtype == "Loess"){
	Basemix = 0.6;
	} else {
	Basemix = 0;
}
if (SLval == 0) {
	SVC = 0;
//	Rcpp::Rcout << "SL = 0" << std::endl;
	} else if (func_Lakedepth(Lake(Time,5), Lakevolumes) - min(Lakevolumes(_,0)) - SLval <= Basemix){
	SVC = Lake(Time,5);
//	Rcpp::Rcout << "SL > 0" << std::endl;
	}else {
//	Rcpp::Rcout << "Third option" << std::endl;
	SVC = Lake(Time,5) - func_Lakevolume(Lake(Time,4) - SLval, Lakevolumes);
} 
return SVC;
	//Rcpp::Rcout << "SVC output" << std::endl;
}

// [[Rcpp::export]]
NumericVector func_Fdlm(long Time, NumericMatrix Isotopes, NumericMatrix Met, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype){
double Fdlm;
double Fdlm18O;
double FdlmD;
Fdlm = func_SVC(Time, Met, Lakevolumes, Lake, Interpolationtype) - Lake(Time,7);
if (Fdlm < 0){
	Fdlm18O = Fdlm*Isotopes(Time, 5);
	FdlmD = Fdlm*Isotopes(Time, 4);
} else { //use values from deep lake as flux flowing opposite way.
	Fdlm18O = Fdlm*Isotopes(Time, 7);
	FdlmD = Fdlm*Isotopes(Time, 6);
}
NumericVector FdlmVals = NumericVector::create(Named("Fdlm") = Fdlm, Named("Fdlm18O") = Fdlm18O, Named("FdlmD") = FdlmD);
return FdlmVals;
}


// [[Rcpp::export]]
double func_KcLake(long Time, double KcLwin, double KcLsum, NumericMatrix Met,std::string Datarate, double Timestep){
double KcAmplitude;
double KSeasonOffset = 0.8;
double KcShift;
double KcLake;
KcAmplitude = (KcLwin - KcLsum)/2;
KcShift = (KcLwin + KcLsum)/2;
KcLake = KcShift + sin(((func_MonthFraction(Time, Met, Datarate, Timestep)/12) + KSeasonOffset) *pi *2)*KcAmplitude;
return KcLake;
}
//daily checked

// [[Rcpp::export]]
NumericVector func_Fe(long Time, double Timestep, NumericMatrix Lake, NumericMatrix Met, double ALBlake, NumericMatrix Isotopes, double PWF, double LatitudeRadians, std::string Datarate, double KcLwin, double KcLsum, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, double Feedback){
	double Fe;
	double Fe18O;
	double FeD;
	if (StratificationModule == "Yes"){
		Fe = func_E(Time, Met, ALBlake, PWF, LatitudeRadians, Datarate, Timestep, ChemConcentrations, StratificationModule, NeutralDragCoeff, Feedback) * Lake(Time,6) * Timestep;
	} else {
		Fe = func_E(Time, Met, ALBlake, PWF, LatitudeRadians, Datarate, Timestep, ChemConcentrations, StratificationModule, NeutralDragCoeff, Feedback) * Lake(Time,6) * func_KcLake(Time, KcLwin, KcLsum, Met,Datarate, Timestep) * Timestep;
	}
	Fe18O = Fe * Isotopes(Time,18);
	FeD = Fe*Isotopes(Time,19);
	NumericVector FeVals = NumericVector::create(Named("Fe") = Fe, Named("Fe18O") = Fe18O, Named("FeD") = FeD);
	return FeVals;
}
//daily should be sorted in func_E.

// [[Rcpp::export]]
double func_FeedbackRH(NumericMatrix Met, long Time, double Feedback){
double RHModified = Met(Time,5)/100 * (1-Feedback) + Feedback;
if (RHModified > 1){
RHModified = 1;
}
return RHModified;
}

// [[Rcpp::export]]
double func_Fof(NumericMatrix Lake, NumericMatrix Lakevolumes, long Time){
double Overflow;
if (Lake(Time,7) + Lake(Time,8) > Lakevolumes(0,1)){
Overflow = (Lake(Time,7) + Lake(Time,8)) - Lakevolumes(0,1);
} else {
Overflow = 0;
}
return Overflow;
}