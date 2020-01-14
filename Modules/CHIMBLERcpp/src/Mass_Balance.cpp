#include <Rcpp.h>
#include "Hydrology_functions.hpp"
#include "Hypsographic_Calcs.hpp"
#include "Isotope_Calcs.hpp"
#include "Mass_Balance.hpp"
#include "ChimbleGW.hpp"
#include "Salinity_calcs.hpp"
#include "Stratification.hpp"

using namespace Rcpp;
using namespace std;

//#define FpCFlux FpC(0)
#define FpC(x) Flux(Time+(x),4)
#define FpL(x) Flux(Time+(x),5)
#define Fe(x) Flux(Time+(x),7)
#define Fsm(x) Flux(Time+(x),6)
#define Fext(x) Flux(Time+(x),20) 
#define Fdlm(x) Flux(Time+(x),10)
#define Fssi(x) Flux(Time+(x),11)
#define Fssd(x) Flux(Time+(x),13)
#define Fdsd(x) Flux(Time+(x),15)
#define Fro(x) Flux(Time+(x),16)
#define Fsse(x) Flux(Time+(x),12)
#define Fdse(x) Flux(Time+(x),14)
#define Fin(x) Flux(Time+(x),17)
#define Fsf(x) Flux(Time+(x),18)
#define Fgw(x) Flux(Time+(x),19)
#define Fdos(x) Flux(Time+(x),9)
#define Fsos(x) Flux(Time+(x),8)
#define Fgw_Fhead(x) Flux(Time+(x),22)
#define Fgw_Ext(x) Flux(Time+(x),23)
#define Fgw_Recharge(x) Flux(Time+(x),24)
#define Fgw_ET(x) Flux(Time+(x),25)
#define Fof(x) Flux(Time+(x),21)
#define EvapRate(x) Flux(Time+(x),26)
#define PETRate(x) Flux(Time+(x),27)
#define GWErr(x) Flux(Time+(x),28)

#define RESsl(x) Lake(Time+(x),7)
#define RESdl(x) Lake(Time+(x),8)
#define RESss(x) Lake(Time+(x),9)
#define RESds(x) Lake(Time+(x),10)
#define RESsp(x) Lake(Time+(x),11)
#define RESin(x) Lake(Time+(x),12)
#define RESgw(x) Lake(Time+(x),13)

#define LakeDepth(x) Lake(Time+(x),4)
#define LakeVolume(x) Lake(Time+(x),5)
#define LakeArea(x) Lake(Time+(x),6)

#define RESsl_D(x) Isotopes(Time+(x),4)
#define RESsl_18O(x) Isotopes(Time+(x),5)
#define RESdl_D(x) Isotopes(Time+(x),6)
#define RESdl_18O(x) Isotopes(Time+(x),7)
#define RESss_D(x) Isotopes(Time+(x),8)
#define RESss_18O(x) Isotopes(Time+(x),9)
#define RESds_D(x) Isotopes(Time+(x),10)
#define RESds_18O(x) Isotopes(Time+(x),11)
#define RESsp_D(x) Isotopes(Time+(x),12)
#define RESsp_18O(x) Isotopes(Time+(x),13)
#define RESin_D(x) Isotopes(Time+(x),14)
#define RESin_18O(x) Isotopes(Time+(x),15)
#define GW_D(x) Isotopes(Time+(x),16)
#define GW_18O(x) Isotopes(Time+(x),17)
#define E_D(x) Isotopes(Time+(x),18)
#define E_18O(x) Isotopes(Time+(x),19)
#define PET_D(x) Isotopes(Time+(x),20)
#define PET_18O(x) Isotopes(Time+(x),21)


// [[Rcpp::export]]
std::string StrSplitFront (std::string StringToSplit, std::string Delim){
//std::string StdStringToSplit = as<std::string>(StringToSplit);
int DelimLoc = 0;
for (int i= 0; i<StringToSplit.size(); i++){
	if (StringToSplit.substr(i,1) == Delim){
	DelimLoc = i;
	}
}
std::string FrontPart = StringToSplit.substr(0,DelimLoc);
return FrontPart;
}

// [[Rcpp::export]]
std::string StrSplitBack (std::string StringToSplit, std::string Delim){
//std::string StdStringToSplit = as<std::string>(StringToSplit);
int DelimLoc = 0;
for (int i= 0; i<StringToSplit.size(); i++){
	if (StringToSplit.substr(i,1) == Delim){
	DelimLoc = i;
	}
}
std::string BackPart = StringToSplit.substr(DelimLoc+1,(StringToSplit.size()-(DelimLoc+1)));
return BackPart;
}

// [[Rcpp::export]]
void func_FluxFpCatchment(NumericMatrix Flux, NumericMatrix Met, long Time, double CA, NumericMatrix Lake, double Timestep){
FpC(0) = func_FpCatchment(Time, Met, CA, Lake, Timestep)[0];
}

// [[Rcpp::export]]
void func_FluxFpLake(NumericMatrix Flux, long Time, NumericMatrix Met, NumericMatrix Lake, double Timestep){
FpL(0) = func_FpLake(Time, Met, Lake, Timestep)[0];
}

// [[Rcpp::export]]
void func_FluxFe (NumericMatrix Flux, long Time, double Timestep, NumericMatrix Lake, NumericMatrix Met, double ALBlake, NumericMatrix Isotopes, double PWF, double LatitudeRadians, std::string Datarate, double KcLwin, double KcLsum, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, double Feedback){
Fe(0) = func_Fe(Time, Timestep, Lake, Met, ALBlake, Isotopes, PWF, LatitudeRadians, Datarate, KcLwin, KcLsum, ChemConcentrations, StratificationModule, NeutralDragCoeff, Feedback)[0];
}

// [[Rcpp::export]]
void func_FluxFsm (NumericMatrix Flux, double CA, long Time, NumericMatrix Lake,NumericMatrix Met, NumericMatrix Isotopes, std::string Datarate, double Timestep){
Fsm(0) = func_Fsm(CA,Time,Lake,Met,Isotopes,Datarate,Timestep)[0];
}

// [[Rcpp::export]]
void func_FluxFext (NumericMatrix Flux, long Time, NumericMatrix Met, NumericMatrix Isotopes, double Timestep){
	Fext(0) = func_Fext(Time, Met, Isotopes, Timestep)[0];
}

// [[Rcpp::export]]
void func_FluxFdlm (NumericMatrix Flux, long Time, NumericMatrix Isotopes, NumericMatrix Met, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype){
Fdlm(0) = func_Fdlm(Time, Isotopes, Met, Lakevolumes, Lake, Interpolationtype)[0];
}

// [[Rcpp::export]]
void func_FluxSoil (NumericMatrix Flux, long Time, float AWCss, float AWCds, double CA, NumericMatrix Lake, NumericMatrix Met, NumericMatrix Isotopes, double KcDS, double  KcSS, double LatitudeRadians,float RunoffRatio, double ALBearth, double ALBlake, double PWF, std::string Datarate, double Timestep){
Fssi(0) = func_Soil(Time, AWCss, AWCds, CA, Lake, Met, Isotopes, KcDS, KcSS, LatitudeRadians,RunoffRatio,ALBearth,ALBlake,PWF, Datarate, Timestep)[0];
Fssd(0) = func_Soil(Time, AWCss, AWCds, CA, Lake, Met, Isotopes, KcDS, KcSS, LatitudeRadians,RunoffRatio,ALBearth,ALBlake,PWF, Datarate, Timestep)[1];
Fdsd(0) = func_Soil(Time, AWCss, AWCds, CA, Lake, Met, Isotopes, KcDS, KcSS, LatitudeRadians,RunoffRatio,ALBearth,ALBlake,PWF, Datarate, Timestep)[2];
Fro(0) = func_Soil(Time, AWCss, AWCds, CA, Lake, Met, Isotopes, KcDS, KcSS, LatitudeRadians,RunoffRatio,ALBearth,ALBlake,PWF, Datarate, Timestep)[3];
}

// [[Rcpp::export]]
void func_FluxSoilEvap (NumericMatrix Flux, double CA, long Time, NumericMatrix Met, NumericMatrix Lake,NumericMatrix Isotopes, float AWCss, float AWCds, double KcDS, double LatitudeRadians, double KcSS, double ALBearth, double ALBlake, std::string Datarate, double PWF, double Timestep){
Fsse(0) = func_Evap(CA, Time, Met, Lake, Isotopes, AWCss, AWCds, KcDS, LatitudeRadians, KcSS, ALBearth, ALBlake, Datarate, PWF,Timestep)[0];
Fdse(0) = func_Evap(CA, Time, Met, Lake, Isotopes, AWCss, AWCds, KcDS, LatitudeRadians, KcSS, ALBearth, ALBlake, Datarate, PWF,Timestep)[1];
}

// [[Rcpp::export]]
void func_FluxFin (NumericMatrix Flux, double CA, NumericMatrix Lake, long Time, NumericMatrix Isotopes, double Cin, double Timestep){
Fin(0) = func_Fin(CA, Lake, Time, Isotopes, Cin, Timestep)[0];
}

// [[Rcpp::export]]
void func_FluxFsf (NumericMatrix Flux, double CA, NumericMatrix Lake, long Time, NumericMatrix Met, double Timestep){
Fsf(0) = func_Fsf(CA, Lake, Time, Met, Timestep)[0];
}

// [[Rcpp::export]]
void func_FluxFgw (NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakesedimentgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double Timestep, double ETexdepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, long Time, NumericMatrix Lake, NumericMatrix GWcatchmentgrid, NumericMatrix Lakevolumes, NumericMatrix GWlakekgrid, NumericMatrix Flux, NumericMatrix LakeDepths, double Csr, NumericMatrix Met, NumericMatrix Isotopes, std::string GroundwaterModule, List GWStats, NumericMatrix GWNumbers, int GWLakeID){
NumericVector GWVals = func_Fgw (Topogrid, GWRechargegrid, GWPumpinggrid, GWKgrid, GWSgrid, GWbasegrid, GWheadgrid, GWboundarygrid, GWlakesedimentgrid, Lakecondgrid, GWETgrid, Classgrid, Timestep, ETexdepth, Cellxdim, Cellydim, Accuracy, Maxiterations, Time, Lake, GWcatchmentgrid, Lakevolumes, GWlakekgrid, Flux, LakeDepths, Csr, Met, Isotopes, GroundwaterModule, GWStats, GWNumbers, GWLakeID);
Fgw(0) = GWVals[0];
Fsos(0) = GWVals[3];
Fdos(0) = GWVals[4];
	if (GroundwaterModule == "Yes"){
		double FHeadFlows_3 = GWStats(3);
		double FHeadFlows_7 = GWStats(7);
		double FHeadFlows_6 = GWStats(6);
		double LakeFlows = -Fgw(0) + Fsos(0) + Fdos(0);
		Fgw_Fhead(0) = -FHeadFlows_3 + FHeadFlows_7 + FHeadFlows_6 - LakeFlows;
		Fgw_Ext(0) = GWStats(5);
		Fgw_Recharge(0) = GWStats(4);
		Fgw_ET(0) = GWStats(10);
		GWErr(0) = GWStats(2);
	}
}

// [[Rcpp::export]]
void func_FluxFof (NumericMatrix Flux, NumericMatrix Lake, NumericMatrix Lakevolumes, long Time){
Fof(0) = func_Fof(Lake, Lakevolumes, Time);
}

// [[Rcpp::export]]
void func_FluxE (NumericMatrix Flux, long Time, NumericMatrix Met, double ALBlake, double PWF, double LatitudeRadians, std::string Datarate, double Timestep, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, double Feedback){
EvapRate(0) = func_E(Time, Met, ALBlake, PWF, LatitudeRadians, Datarate, Timestep, ChemConcentrations, StratificationModule, NeutralDragCoeff, Feedback);
}

// [[Rcpp::export]]
void func_FluxPET (NumericMatrix Flux, long Time, NumericMatrix Met, double ALBearth, double PWF, double LatitudeRadians, std::string Datarate, double Timestep){
PETRate(0) = func_PET(Time, Met, ALBearth, PWF, LatitudeRadians, Datarate, Timestep);
}
// [[Rcpp::export]]
void func_Fluxes(double KcSS, double ALBearth, double ALBlake, double CA, double Cin, double Csr, double KcDS, double KcLsum, double KcLwin, double LatitudeRadians, float RunoffRatio, double PWF, double Timestep, float AWCds, float AWCss, long Time, NumericMatrix Flux, NumericMatrix Isotopes, NumericMatrix Lake, NumericMatrix Met, NumericMatrix Lakevolumes, std::string Datarate, std::string Interpolationtype, std::string GroundwaterModule, NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakesedimentgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double ETexdepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, NumericMatrix GWcatchmentgrid, NumericMatrix GWlakekgrid, NumericMatrix LakeDepths, List GWStats, NumericMatrix GWNumbers, NumericMatrix ChemConcentrations, std::string StratificationModule, double NeutralDragCoeff, int GWLakeID, double Feedback){
func_FluxFpCatchment(Flux, Met, Time, CA, Lake, Timestep);
func_FluxFpLake(Flux, Time, Met, Lake, Timestep);
func_FluxFe (Flux, Time, Timestep, Lake, Met, ALBlake, Isotopes, PWF, LatitudeRadians, Datarate, KcLwin, KcLsum, ChemConcentrations, StratificationModule, NeutralDragCoeff, Feedback);
func_FluxFsm (Flux, CA, Time, Lake,Met, Isotopes, Datarate, Timestep);
func_FluxSoil (Flux, Time, AWCss, AWCds, CA, Lake, Met, Isotopes, KcDS,  KcSS, LatitudeRadians, RunoffRatio, ALBearth, ALBlake, PWF, Datarate, Timestep);
func_FluxSoilEvap (Flux, CA, Time, Met, Lake,Isotopes, AWCss, AWCds, KcDS, LatitudeRadians, KcSS, ALBearth, ALBlake, Datarate, PWF, Timestep);
func_FluxFin (Flux, CA, Lake, Time, Isotopes, Cin, Timestep);
func_FluxFsf (Flux, CA, Lake, Time, Met, Timestep);
func_FluxFext (Flux, Time, Met, Isotopes, Timestep);
func_FluxFof (Flux, Lake, Lakevolumes, Time);
func_FluxE (Flux, Time, Met, ALBlake, PWF, LatitudeRadians, Datarate, Timestep, ChemConcentrations, StratificationModule, NeutralDragCoeff, Feedback);
func_FluxPET (Flux, Time, Met, ALBearth, PWF, LatitudeRadians, Datarate, Timestep);
//run groundwater last. It requires PET and FluxSoil to be defined.
func_FluxFgw (Topogrid, GWRechargegrid, GWPumpinggrid, GWKgrid, GWSgrid, GWbasegrid, GWheadgrid, GWboundarygrid, GWlakesedimentgrid, Lakecondgrid, GWETgrid, Classgrid, Timestep, ETexdepth, Cellxdim, Cellydim, Accuracy, Maxiterations, Time, Lake, GWcatchmentgrid, Lakevolumes, GWlakekgrid, Flux, LakeDepths, Csr, Met, Isotopes, GroundwaterModule, GWStats, GWNumbers, GWLakeID);
}

//******************************************************************************************************
//Hydrology mass balance
// [[Rcpp::export]]
void func_Hydrology(NumericMatrix Flux, NumericMatrix Lake, NumericMatrix Isotopes, long Time, NumericMatrix Met, NumericMatrix Lakevolumes, std::string Interpolationtype, std::string DirectExternalFlows, std::string GroundwaterModule, List GWStats){


	double RunoffPlusExt = Fro(-1);
	double ExtFlows = Fext(-1);
if (DirectExternalFlows == "No"){ //these should be temporary variables so they don't mess up isotopes.
	RunoffPlusExt = Fext(-1) + Fro(-1);
	ExtFlows = 0;
}


	//ResSL and ResDL - these, and following fluxes follow normal mass balance rules. Fluxes applied to reservoirs at time-1, to give a new reservoir value at time. Fdlm is then applied.
	//This leaves just runoff and external flows passing through the RESin reservoir.
	
	if (Met(Time-1,15) == 0){
		RESsl(0) = RESsl(-1); //Fdlm will sort out any excess volume.
		RESdl(0) = RESdl(-1)+Fgw(-1)+FpL(-1)+Fin(-1) + ExtFlows - Fe(-1) - Fdos(-1) - Fsos(-1) - Fof(-1);
	} else if (RESsl(-1) < Fsos(-1) + Fe(-1) + Fof(-1) - FpL(-1) - Fin(-1)-ExtFlows){
		double ExcessFlux =  (Fsos(-1) + Fe(-1) + Fof(-1) - FpL(-1) - Fin(-1) - ExtFlows - RESsl(-1));
		double ExcessFluxE = ExcessFlux * (Fe(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
		double ExcessFluxSos = ExcessFlux * (Fsos(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
		double ExcessFluxOf = ExcessFlux * (Fof(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
		RESsl(0) = RESsl(-1) + FpL(-1) + Fin(-1) + ExtFlows - (Fsos(-1) - ExcessFluxSos) - (Fe(-1) - ExcessFluxE) - (Fof(-1) - ExcessFluxOf); //should = 0.
		RESdl(0) = RESdl(-1) + Fgw(-1) - Fdos(-1) - ExcessFluxE - ExcessFluxSos - ExcessFluxOf;
	} else {
		RESsl(0) = RESsl(-1) + FpL(-1) + Fin(-1) + ExtFlows - Fsos(-1) - Fe(-1) - Fof(-1);
		RESdl(0) = RESdl(-1) + Fgw(-1) - Fdos(-1);
	}

//ResSS & ResDS
RESss(0) = RESss(-1) + Fssi(-1) - Fsse(-1);
RESds(0) = RESds(-1) + Fssd(-1) - Fdse(-1);


//Groundwater is always run into the lake. (Flux = 19)
//If groundwater module is active, Fdsd (Flux = 15) flows to the groundwater reservoir (in hydrology functions) and should be disabled here.
//ResIn
if (GroundwaterModule == "No"){
	RESin(0) = RESin(-1) + RunoffPlusExt + Fdsd(-1) - Fin(-1);
} else {
	RESin(0) = RESin(-1) + RunoffPlusExt - Fin(-1); //Fdsd is routed to groundwater through the groundwater function in hydrology.
}

//ResSP
RESsp(0) = RESsp(-1) + Fsf(-1) - Fsm(-1);

//ResGw is defined by the GW calculations. (volume applied in Hydrology functions)
if (GroundwaterModule == "Yes"){
RESgw(0) = GWStats(1);
}

//Depth, Area and Volume
Lake(Time, 5) = RESsl(0) + RESdl(0);
if (LakeVolume(0) > 0){
	LakeDepth(0) = func_Lakedepth(LakeVolume(0), Lakevolumes);
	LakeArea(0) = func_Lakearea(LakeVolume(0), Lakevolumes);
} else {
	LakeVolume(0) = 0;
	LakeDepth(0) = 0;
	LakeArea(0) = 0;
}

}

//******************************************************************************************************
//Isotopes mass balance
// [[Rcpp::export]]
void func_IsotopeFluxes (NumericMatrix Flux, long Time, NumericMatrix Isotopes, NumericMatrix Met, double AtmosphericShift, double Turbulence, double Feedback, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype, std::string DirectExternalFlows, NumericMatrix FluxD, NumericMatrix Flux18O, std::string GroundwaterModule, List GWStats, NumericMatrix GWNumbers, NumericMatrix ChemConcentrations, double Theta){
    E_D(-1) = func_dE(Time-1, "D", Isotopes, Met, AtmosphericShift, Turbulence, Feedback, Lakevolumes, Lake, Interpolationtype, ChemConcentrations, Theta);	
	E_18O(-1) = func_dE(Time-1, "18O", Isotopes, Met, AtmosphericShift, Turbulence, Feedback, Lakevolumes, Lake, Interpolationtype, ChemConcentrations, Theta);

double FpCatchmentD = FpC(-1) * Met(Time-1,9);
double FpLakeD = FpL(-1) * Met(Time-1,9);
double FsmD = Fsm(-1) * RESsp_D(-1);
double FeD = Fe(-1) * E_D(-1);
double FsosD = Fsos(-1) * RESsl_D(-1);
double FdosD = Fdos(-1) * RESdl_D(-1);
//double FdlmD;
double FssiD;
double FsseD = Fsse(-1) * RESss_D(-1); // no fractionation for evapotranspiration
double FssdD;
double FdseD = Fdse(-1) * RESds_D(-1); // no fractionation for evapotranspiration
double FdsdD;
double FroD;
double FinD = Fin(-1) * RESin_D(-1);
double FsfD = Fsf(-1) * Met(Time-1,9);
double FgwD = Fgw(-1) * GW_D(-1);//may need a caveat depending on groundwater module setting. Should be fixed now. Ext and GW flows separated.
double FextD = Fext(-1) * Met(Time-1,13);
double FofRESslD = RESsl_D(-1);
double FofRESdlD = RESdl_D(-1);

double FpCatchment18O = FpC(-1) * Met(Time-1,8);
double FpLake18O = FpL(-1) * Met(Time-1,8);
double Fsm18O = Fsm(-1) * RESsp_18O(-1);
double Fe18O = Fe(-1) * E_18O(-1);
double Fsos18O = Fsos(-1) * RESsl_18O(-1);
double Fdos18O = Fdos(-1) * RESdl_18O(-1);
//double Fdlm18O;
double Fssi18O;
double Fsse18O = Fsse(-1) * RESss_18O(-1);
double Fssd18O;
double Fdse18O = Fdse(-1) * RESds_18O(-1);
double Fdsd18O;
double Fro18O;
double Fin18O = Fin(-1) * RESin_18O(-1);
double Fsf18O = Fsf(-1) * Met(Time-1,8);
double Fgw18O = Fgw(-1) * GW_18O(-1); //may need a caveat depending on groundwater module setting. Should be fixed now. Ext and GW flows separated.
double Fext18O = Fext(-1) * Met(Time-1,12);
double FofRESsl18O = RESsl_18O(-1);
double FofRESdl18O = RESdl_18O(-1);

//calc groundwater fluxes
double GWF_HeadD = Fgw_Fhead(-1)*GW_D(-1);
double GWExtD = Fgw_Ext(-1)*GW_D(-1); //these may need to be updated for incoming fluxes.
double GW_ETD = Fgw_ET(-1) * GW_D(-1);
double GWF_Head18O = Fgw_Fhead(-1)*GW_18O(-1);
double GWExt18O = Fgw_Ext(-1)*GW_18O(-1); //these may need to be updated for incoming fluxes.
double GW_ET18O = Fgw_ET(-1) * GW_18O(-1);
double GW_RechargeD; // requires Fdsd. Covered in next section
double GW_Recharge18O;

if (Fsm(-1) + FpL(-1) > 0){ //Separate flows based on volume
	FssiD = (FsmD + FpCatchmentD) * Fssi(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
	FssdD = (FsmD + FpCatchmentD) * Fssd(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
	FdsdD = (FsmD + FpCatchmentD) * Fdsd(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
	FroD = (FsmD + FpCatchmentD) * Fro(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
	Fssi18O = (Fsm18O + FpCatchment18O) * Fssi(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
	Fssd18O = (Fsm18O + FpCatchment18O) * Fssd(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
	Fdsd18O = (Fsm18O + FpCatchment18O) * Fdsd(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
	Fro18O = (Fsm18O + FpCatchment18O) * Fro(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
	GW_RechargeD = Fgw_Recharge(-1)*((FsmD + FpCatchmentD)/(Fsm(-1)+FpC(-1))); //mix incoming sources 
	GW_Recharge18O = Fgw_Recharge(-1)*((Fsm18O + FpCatchment18O)/(Fsm(-1)+FpC(-1))); //mix incoming sources 
} else {
	FssiD = 0;
	FssdD = 0;
	FdsdD = 0;
	FroD = 0;
	Fssi18O = 0;
	Fssd18O = 0;
	Fdsd18O = 0;
	Fro18O = 0;
	GW_RechargeD = 0;
	GW_Recharge18O = 0;

}

//Sort out direct vs indirect external flows.
double RunoffPlusExt = Fro(-1);
double ExtFlows = Fext(-1);
double Fro_extD = FroD;
double Fro_ext18O = Fro18O;
if (DirectExternalFlows == "No"){
	RunoffPlusExt = Fext(-1) + Fro(-1);
	FroD = FextD + FroD;
	Fro18O = Fext18O + Fro18O;
	ExtFlows = 0;
	FextD = 0;
	Fext18O = 0;
}

//Fill out flux matrix (except Fdlm and Fof)
FluxD(Time-1,4) = FpCatchmentD;
FluxD(Time-1,5) = FpLakeD;
FluxD(Time-1,6) = FsmD;
FluxD(Time-1,7) = FeD;
FluxD(Time-1,8) = FsosD;
FluxD(Time-1,9) = FdosD;
FluxD(Time-1,11) = FssiD;
FluxD(Time-1,12) = FsseD;
FluxD(Time-1,13) = FssdD;
FluxD(Time-1,14) = FdseD;
FluxD(Time-1,15) = FdsdD;
FluxD(Time-1,16) = FroD;
FluxD(Time-1,17) = FinD;
FluxD(Time-1,18) = FsfD;
FluxD(Time-1,19) = FgwD;
FluxD(Time-1,20) = FextD;
FluxD(Time-1,21) = 0;//added later
FluxD(Time-1,22) = GWF_HeadD;
FluxD(Time-1,23) = GWExtD;
FluxD(Time-1,24) = GW_RechargeD;
FluxD(Time-1,25) = GW_ETD;
FluxD(Time-1,26) = 0;
FluxD(Time-1,27) = 0;

Flux18O(Time-1,4) = FpCatchment18O;
Flux18O(Time-1,5) = FpLake18O;
Flux18O(Time-1,6) = Fsm18O;
Flux18O(Time-1,7) = Fe18O;
Flux18O(Time-1,8) = Fsos18O;
Flux18O(Time-1,9) = Fdos18O;
Flux18O(Time-1,11) = Fssi18O;
Flux18O(Time-1,12) = Fsse18O;
Flux18O(Time-1,13) = Fssd18O;
Flux18O(Time-1,14) = Fdse18O;
Flux18O(Time-1,15) = Fdsd18O;
Flux18O(Time-1,16) = Fro18O;
Flux18O(Time-1,17) = Fin18O;
Flux18O(Time-1,18) = Fsf18O;
Flux18O(Time-1,19) = Fgw18O;
Flux18O(Time-1,20) = Fext18O;
Flux18O(Time-1,21) = 0; //added later
Flux18O(Time-1,22) = GWF_Head18O;
Flux18O(Time-1,23) = GWExt18O;
Flux18O(Time-1,24) = GW_Recharge18O;
Flux18O(Time-1,25) = GW_ET18O;
Flux18O(Time-1,26) = 0;
Flux18O(Time-1,27) = 0;

//Isotope mass balance
if (Met(Time-1,15) == 0){
	RESsl_D(0) = RESsl_D(-1);
	RESsl_18O(0) = RESsl_18O(-1);
	RESdl_D(0) = (RESdl(-1)*RESdl_D(-1) + FgwD + FpLakeD + FinD + FextD - FeD - FdosD - FsosD - FofRESdlD*Fof(-1))/RESdl(0);
	RESdl_18O(0) = (RESdl(-1)*RESdl_18O(-1)+Fgw18O+FpLake18O+Fin18O+Fext18O - Fe18O - Fdos18O - Fsos18O - FofRESdl18O*Fof(-1))/RESdl(0);
	FluxD(Time-1,21) = FofRESdlD*Fof(-1);
	Flux18O(Time-1,21) = FofRESdl18O*Fof(-1);

} else if (RESsl(-1) < Fsos(-1) + Fe(-1) + Fof(-1) - FpL(-1) - Fin(-1) - ExtFlows){
	double ExcessFlux =  (Fsos(-1) + Fe(-1) + Fof(-1) - FpL(-1) - Fin(-1) - ExtFlows - RESsl(-1));
	double ExcessFluxE = ExcessFlux * (Fe(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
	double ExcessFluxSos = ExcessFlux * (Fsos(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
	double ExcessFluxOf = ExcessFlux * (Fof(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
	RESsl_D(0) = 0;//(RESsl(-1)*RESsl_D(-1) + FpLakeD + FinD - (Fsos(-1) - ExcessFluxSos)*RESsl_D(-1) - (Fe(-1) - ExcessFluxE)*E_18O(-1) - (Fof(-1) - ExcessFluxOf)*FofRESslD)/RESsl(0);
	RESsl_18O(0) = 0;//(RESsl(-1)*RESsl_D(-1) + FpLakeD + FinD - (Fsos(-1) - ExcessFluxSos)*RESsl_D(-1) - (Fe(-1) - ExcessFluxE)*E_18O(-1) - (Fof(-1) - ExcessFluxOf)*FofRESslD)/RESsl(0);
	RESdl_D(0) = (RESdl(-1)*RESdl_D(-1) + FgwD - FdosD - ExcessFluxE*E_D(-1) - ExcessFluxSos*RESdl_D(-1) - ExcessFluxOf*FofRESdlD)/RESdl(0);
	RESdl_18O(0) = (RESdl(-1)*RESdl_18O(-1) + Fgw18O - Fdos18O - ExcessFluxE*E_18O(-1) - ExcessFluxSos*RESdl_18O(-1) - ExcessFluxOf*FofRESdl18O)/RESdl(0);
	//Fill out flux matrix
	FluxD(Time-1,8) = (Fsos(-1) - ExcessFluxSos)*FofRESslD + ExcessFluxSos*FofRESdlD;
	FluxD(Time-1,21) = (Fof(-1) - ExcessFluxOf)*FofRESslD + ExcessFluxOf*FofRESdlD;
	Flux18O(Time-1,8) = (Fsos(-1) - ExcessFluxSos)*FofRESsl18O + ExcessFluxSos*FofRESdl18O;
	Flux18O(Time-1,21) = (Fof(-1) - ExcessFluxOf)*FofRESsl18O + ExcessFluxOf*FofRESdl18O;
} else {
	RESsl_D(0) = (RESsl(-1)*RESsl_D(-1) + FpLakeD + FextD + FinD - FsosD - FeD - Fof(-1)*FofRESslD)/RESsl(0);
	RESsl_18O(0) = (RESsl(-1)*RESsl_18O(-1) + FpLake18O + Fext18O + Fin18O - Fsos18O - Fe18O - Fof(-1)*FofRESsl18O)/RESsl(0);
	RESdl_D(0) = (RESdl(-1)*RESdl_D(-1) + FgwD - FdosD)/RESdl(0);
	RESdl_18O(0) = (RESdl(-1)*RESdl_18O(-1) + Fgw18O - Fdos18O)/RESdl(0);
	FluxD(Time-1,21) = Fof(-1)*FofRESslD;
	Flux18O(Time-1,21) = Fof(-1)*FofRESsl18O;
}

//ResSS & ResDS
if (RESss(0) == 0){
	RESss_D(0) = 0;
	RESss_18O(0) = 0;
} else {
	RESss_D(0) = (RESss(-1)*RESss_D(-1) + FssiD - FsseD)/RESss(0);
	RESss_18O(0) = (RESss(-1)*RESss_18O(-1) + Fssi18O - Fsse18O)/RESss(0);
}
if (RESds(0) == 0){
	RESds_D(0) = 0;
	RESds_18O(0) = 0;
} else {
	RESds_D(0) = (RESds(-1)*RESds_D(-1) + FssdD - FdseD)/RESds(0);
	RESds_18O(0) = (RESds(-1)*RESds_18O(-1) + Fssd18O - Fdse18O)/RESds(0);
}

	//ResIn
	if (RESin(0) == 0){
	RESin_D(0) = 0;
	RESin_18O(0) = 0;
	} else {
		if (GroundwaterModule == "No"){
		RESin_D(0) = (RESin(-1)*RESin_D(-1) + Fro_extD + FdsdD - FinD)/RESin(0);
		RESin_18O(0) = (RESin(-1)*RESin_18O(-1) + Fro_ext18O + Fdsd18O - Fin18O)/RESin(0);
		} else {
		RESin_D(0) = (RESin(-1)*RESin_D(-1) + Fro_extD - FinD)/RESin(0);//FdsdD removed as per hydrology mass balance
		RESin_18O(0) = (RESin(-1)*RESin_18O(-1) + Fro_ext18O - Fin18O)/RESin(0); //Fdsd18O removed as per hydrology mass balance
		}
	}

//ResSP
if (RESsp(0) == 0){
	RESsp_D(0) = 0;
	RESsp_18O(0) = 0;
	} else {
	RESsp_D(0) = (RESsp(-1)*RESsp_D(-1) +FsfD - FsmD)/RESsp(0);
	RESsp_18O(0) = (RESsp(-1)*RESsp_18O(-1) +Fsf18O - Fsm18O)/RESsp(0);
}


//ResGW
	if (GroundwaterModule == "Yes"){
	GW_D(0) = (RESgw(-1)*GW_D(-1) + GW_RechargeD - FgwD + FsosD + FdosD + GWExtD - GWF_HeadD - GW_ETD)/RESgw(0); 
	GW_18O(0) = (RESgw(-1)*GW_18O(-1) +  GW_Recharge18O - Fgw18O + Fsos18O + Fdos18O + GWExt18O - GWF_Head18O - GW_ET18O)/RESgw(0); 
	}
}
//******************************************************************************************************
	//Chemistry mass balance
// [[Rcpp::export]]
void func_ChemistryFluxes (NumericMatrix Chemistry, long Time, NumericMatrix Flux, NumericMatrix Lake, NumericMatrix Met, std::string DirectExternalFlows, std::string GroundwaterModule, List GWStats, NumericMatrix GWNumbers){
std::string Delim = "_";
std::string	ResType;
std::string	ChemType;

//work on chemistry ion by ion. Each ion is a block of 9 columns (1 for each reservoir, plus P + Ext) in the chemistry file.
for (int i=4; i < Chemistry.ncol(); i = i+9){
    //StringVector ChemColumns = colnames(Chemistry);
    //std::string ColNameString = as<std::string>(ChemColumns(i));
    //Rcpp::Rcout << "ColNameString" << ColNameString << std::endl;
    //ResType = StrSplitFront(ColNameString,Delim);
	//ChemType = StrSplitBack(ColNameString,Delim);
	
	//increment the precip and external values. May need to look at a way to use a timeseries for these at some stage.
	Chemistry(Time,i+7) = Chemistry(Time-1,i+7);
	Chemistry(Time,i+8) = Chemistry(Time-1,i+8);	
	
	double FpCatchmentChem = FpC(-1) * Chemistry(Time-1,i+7);
	double FpLakeChem = FpL(-1) * Chemistry(Time-1,i+7);
	double FsmChem = Fsm(-1) * Chemistry(Time-1,i+4);
	double FeChem = 0;//Fe(-1) * Chemistry(Time-1,i+7);
	double FsosChem = Fsos(-1) * Chemistry(Time-1,i);
	double FdosChem = Fdos(-1) * Chemistry(Time-1,i+1);
	//double FdlmChem ;
	double FssiChem;
	double FsseChem = 0;//Fsse(-1) * Chemistry(Time-1,i+7);
	double FssdChem ;
	double FdseChem = 0;//Fdse(-1) * Chemistry(Time-1,i+7);
	double FdsdChem ;
	double FroChem ;
	double FinChem = Fin(-1) * Chemistry(Time-1,i+5);
	double FsfChem = Fsf(-1) * Chemistry(Time-1,i+7);
	double FgwChem = Fgw(-1) * Chemistry(Time-1,i+6); //may need a caveat depending on groundwater module setting. Should be fixed now. Ext and GW flows separated.
	double FextChem = Fext(-1) * Chemistry(Time-1,i+8);
	double FofRESslChem = Chemistry(Time-1,i);
	double FofRESdlChem = Chemistry(Time-1,i+1);

		if (Fsm(-1) + FpL(-1) > 0){
			FssiChem = (FsmChem + FpCatchmentChem) * Fssi(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
			FssdChem = (FsmChem + FpCatchmentChem) * Fssd(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
			FdsdChem = (FsmChem + FpCatchmentChem) * Fdsd(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));
			FroChem = (FsmChem + FpCatchmentChem) * Fro(-1)/(Fssi(-1) + Fssd(-1) + Fdsd(-1) + Fro(-1));

			} else {
				FssiChem = 0;
				FssdChem = 0;
				FdsdChem = 0;
				FroChem = 0;
		}

	double FextFlow = Fext(-1);
	double Fro_ext = Fro(-1);
	double Fro_extChem = FroChem;

		if (DirectExternalFlows == "No"){
			Fro_ext = FextFlow + Fro_ext;
			Fro_extChem = FextChem + FroChem;
			FextFlow = 0;
			FextChem = 0;
		}
		
//Sort out groundwaters
double GWF_HeadChem = Fgw_Fhead(-1)* Chemistry(Time-1,i+6);
double GWExtChem = Fgw_Ext(-1)*Chemistry(Time-1,i+6); //these may need to be updated for incoming fluxes.
double GW_RechargeChem = Fgw_Recharge(-1)*Chemistry(Time-1,i+7); //mix incoming sources 
double GW_ETChem = Fgw_ET(-1) * Chemistry(Time-1,i+6);

	//Chemistry mass balance
		if (Met(Time-1,15) == 0){
			Chemistry(Time,i) = Chemistry(Time-1,i);
			Chemistry(Time,i+1) = (RESdl(-1)*Chemistry(Time-1,i+1)+ FgwChem + FpLakeChem + FinChem + FextChem - FeChem - FdosChem - FsosChem - FofRESdlChem*Fof(-1))/RESdl(0);
		} else if (RESsl(-1) < Fsos(-1) + Fe(-1) + Fof(-1) - FpL(-1) - Fin(-1) - FextFlow){
			double ExcessFlux =  (Fsos(-1) + Fe(-1) + Fof(-1) - FpL(-1) - Fin(-1) - FextFlow - RESsl(-1));
			double ExcessFluxE = ExcessFlux * (Fe(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
			double ExcessFluxSos = ExcessFlux * (Fsos(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
			double ExcessFluxOf = ExcessFlux * (Fof(-1)/(Fe(-1) + Fsos(-1) + Fof(-1)));
			Chemistry(Time,i) = 0;//(RESsl(-1)*RESsl_D(-1) + FpLakeChem + FinChem - (Fsos(-1) - ExcessFluxSos)*RESsl_D(-1) - (Fe(-1) - ExcessFluxE)*E_18O(-1) - (Fof(-1) - ExcessFluxOf)*FofRESslD)/RESsl(0);
			Chemistry(Time,i+1) = (RESdl(-1)*Chemistry(Time-1,i+1) + FgwChem - FdosChem - ExcessFluxE * Chemistry(Time-1,i+7) - ExcessFluxSos*Chemistry(Time-1,i+1) - ExcessFluxOf*FofRESdlChem)/RESdl(0);
		} else {
			Chemistry(Time,i) = (RESsl(-1)*Chemistry(Time-1,i) + FpLakeChem + FextChem + FinChem - FsosChem - FeChem - Fof(-1)*FofRESslChem)/RESsl(0);
			Chemistry(Time,i+1) = (RESdl(-1)*Chemistry(Time-1,i+1) + FgwChem - FdosChem)/RESdl(0);
		}

	//ResSS & ResDS
		if (RESss(0) == 0){
			Chemistry(Time,i+2) = 0;
			} else {
			Chemistry(Time,i+2) = (RESss(-1)*Chemistry(Time-1,i+2) + FssiChem - FsseChem)/RESss(0);
		}
		if (RESds(0) == 0){
			Chemistry(Time,i+3) = 0;
			} else {
			Chemistry(Time,i+3) = (RESds(-1)*Chemistry(Time-1,i+3) + FssdChem - FdseChem)/RESds(0);
		}

	//ResIn
		if (RESin(0) == 0){
			Chemistry(Time,i+5) = 0;
			} else {
				if (GroundwaterModule == "No"){
				Chemistry(Time,i+5) = (RESin(-1)*Chemistry(Time-1,i+5) + Fro_extChem + FdsdChem - FinChem)/RESin(0);
				} else {
				Chemistry(Time,i+5) = (RESin(-1)*Chemistry(Time-1,i+5) + Fro_extChem - FinChem)/RESin(0);
				}
		}

	//ResSP
		if (RESsp(0) == 0){
			Chemistry(Time,i+4) = 0;
			} else {
			Chemistry(Time,i+4)= (RESsp(-1)*Chemistry(Time-1,i+4) +FsfChem - FsmChem)/RESsp(0);
		}
		
	//ResGW
		if (GroundwaterModule == "Yes"){
		Chemistry(Time,i+6) = (RESgw(-1)*Chemistry(Time-1,i+6) + GW_RechargeChem - FgwChem + FsosChem + FdosChem + GWExtChem -GWF_HeadChem -GW_ETChem)/RESgw(0); 
		} else {
		Chemistry(Time,i+6) = Chemistry(Time-1,i+6);
		}	
	}	
}



//Fdlm - operates at current timestep.
//******************************************************************************************************
// [[Rcpp::export]]
void func_LakeMixer(NumericMatrix Flux, NumericMatrix Lake, NumericMatrix Chemistry, NumericMatrix Isotopes, long Time, NumericMatrix Met, NumericMatrix Lakevolumes, std::string Interpolationtype, NumericMatrix FluxD, NumericMatrix Flux18O){
//Fdlm - operates at current timestep.
func_FluxFdlm (Flux, Time, Isotopes, Met, Lakevolumes, Lake, Interpolationtype); //This is calculated after the fluxes are all applied, as it's a condition to be met for each timestep.
double Fdlm;
double FdlmD;
double Fdlm18O;
int StratType;

Fdlm = Fdlm(0);
if (Fdlm < 0){
FdlmD = Fdlm*RESsl_D(0);
Fdlm18O = Fdlm*RESsl_18O(0);
} else {
FdlmD = Fdlm*RESdl_D(0);
Fdlm18O = Fdlm*RESdl_18O(0);
}

FluxD(Time,10) = FdlmD;
Flux18O(Time,10) = Fdlm18O;

//Isotopes
if (LakeVolume(0) == 0){
	StratType = 1;
	RESsl_D(0) = 0;
	RESsl_18O(0) = 0;
	RESdl_D(0) = 0;
	RESdl_18O(0) = 0;
} else if (RESsl(0) + Fdlm == 0){
	StratType = 2;
	RESsl_D(0) = 0;
	RESsl_18O(0) = 0;
	RESdl_D(0) = (RESdl(0)*RESdl_D(0) - FdlmD)/(RESdl(0) - Fdlm);
	RESdl_18O(0) = (RESdl(0)*RESdl_18O(0) - Fdlm18O)/(RESdl(0) - Fdlm);
} else if (RESdl(0) - Fdlm(0) <= 0.01){
	StratType = 3;
	RESsl_D(0) = (RESsl(0)*RESsl_D(0) + FdlmD)/(RESsl(0) + Fdlm);
	RESsl_18O(0) = (RESsl(0)*RESsl_18O(0) + Fdlm18O)/(RESsl(0) + Fdlm);
	RESdl_D(0) = 0;
	RESdl_18O(0) = 0;
} else {
	StratType = 4;
	RESsl_D(0) = (RESsl(0)*RESsl_D(0) + FdlmD)/(RESsl(0) + Fdlm);
	RESsl_18O(0) = (RESsl(0)*RESsl_18O(0) + Fdlm18O)/(RESsl(0) + Fdlm);
	RESdl_D(0) = (RESdl(0)*RESdl_D(0) - FdlmD)/(RESdl(0) - Fdlm);
	RESdl_18O(0) = (RESdl(0)*RESdl_18O(0) - Fdlm18O)/(RESdl(0) - Fdlm);
	}
	
//Chemistry
MixChems (Lake, Time, Flux, Chemistry, StratType);
	
//Hydrology
if (RESdl(0) - Fdlm(0) <= 0.01){
	RESdl(0) = 0;
	} else {
	RESsl(0) = RESsl(0) + Fdlm;
	RESdl(0) = RESdl(0) - Fdlm;
}
}

void MixChems(NumericMatrix Lake, long Time, NumericMatrix Flux, NumericMatrix Chemistry, int StratType){
	double Fdlm = Fdlm(0);
	double FdlmChem;

for (int i=4; i < Chemistry.ncol(); i = i+9){
	if (Fdlm < 0){
		FdlmChem = Fdlm*Chemistry(Time,i);
		} else {
		FdlmChem = Fdlm*Chemistry(Time,i+1);
	}

	switch(StratType){
		case 1: {
			Chemistry(Time,i) = Chemistry(Time,i);
			Chemistry(Time,i+1) = Chemistry(Time,i+1);
		break;
		}
	
		case 2: {
			Chemistry (Time,i) = 0;
			Chemistry (Time, i+1) = (RESdl(0)*Chemistry(Time,i+1) - FdlmChem)/(RESdl(0) - Fdlm);
		break;
		}
	
		case 3: {
			Chemistry (Time,i) = (RESsl(0)*Chemistry(Time,i) + FdlmChem)/(RESsl(0) + Fdlm);
			Chemistry (Time, i+1) = 0;
		break;
		}

		case 4: {
			Chemistry (Time,i) = (RESsl(0)*Chemistry(Time,i) + FdlmChem)/(RESsl(0) + Fdlm);
			Chemistry (Time, i+1) = (RESdl(0)*Chemistry(Time,i+1) - FdlmChem)/(RESdl(0) - Fdlm);
		break;
		}
	}
}
}


// && (Met(Time,1) == ScenarioEvents(i,1)) && (Met(Time,2) == ScenarioEvents(i,2))
// [[Rcpp::export]]
void MassBalance (double KcSS, double ALBearth, double ALBlake, double CA, double Cin, double Csr, double KcDS, double KcLsum, double KcLwin, double LatitudeRadians, float RunoffRatio, double PWF, double Timestep, float AWCds, float AWCss, long Time, NumericMatrix Flux, NumericMatrix Isotopes, NumericMatrix Lake, NumericMatrix Met, NumericMatrix Lakevolumes, std::string Datarate, std::string Interpolationtype, std::string GroundwaterModule, NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakesedimentgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double ETexdepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, NumericMatrix GWcatchmentgrid, NumericMatrix GWlakekgrid, NumericMatrix LakeDepths, NumericMatrix Chemistry, std::string DirectExternalFlows, double AtmosphericShift, double Turbulence, double Feedback, NumericMatrix FluxD, NumericMatrix Flux18O, List GWStats, NumericMatrix GWNumbers, List ScenarioEventsList, DataFrame ScenarioEvents, NumericMatrix ChemConcentrations, NumericMatrix ChemParams, double NeutralDragCoeff, double SWExtCoef, double ThermoclineTargetThickness, double ThermoclineMaximum, std::string StratificationModule, NumericMatrix Stratification, int GWLakeID, double Theta){
	//Prep scenario data for use.
	IntegerVector Year = ScenarioEvents["Year"];
	IntegerVector Month = ScenarioEvents["Month"];
	IntegerVector Day = ScenarioEvents["Day"];
	CharacterVector Names = ScenarioEvents["Name"];
	
//Start mass balance
for (long i = 1; i <= Met.nrow(); i++){
	Time = i;
	
//Check for events		
	if (Met(Time, 19) == 1){
	//Rcpp::Rcout << "Date checked" << Met(Time,0) << "/" << Met(Time,1) << "/" << Met (Time,2) << std::endl;
		for (int j = 0; j < ScenarioEvents.nrow() ;j++){
		//Rcpp::Rcout << "Parsing" << Names(j) << std::endl;
		if ((Met(Time,0) == Year(j)) && (Met(Time,1) == Month(j)) && (Met(Time,2) == Day(j))){
				std::string Label = Rcpp::as<std::string>(Names(j));
				if (Names(j) == "KcSS") {KcSS = ScenarioEventsList(j);}
				if (Names(j) == "ALBearth") {ALBearth  = ScenarioEventsList(j);}
				if (Names(j) == "ALBlake") {ALBlake  = ScenarioEventsList(j);}
				if (Names(j) == "CA") {CA  = ScenarioEventsList(j);}
				if (Names(j) == "Cin") {Cin  = ScenarioEventsList(j);}
				if (Names(j) == "Csr") {Csr  = ScenarioEventsList(j);}
				if (Names(j) == "KcDS") {KcDS  = ScenarioEventsList(j);}
				if (Names(j) == "KcLsum") {KcLsum  = ScenarioEventsList(j);}
				if (Names(j) == "KcLwin") {KcLwin  = ScenarioEventsList(j);}
				if (Names(j) == "LatitudeRadians") {LatitudeRadians  = ScenarioEventsList(j);}
				if (Names(j) == "RunoffRatio") {RunoffRatio  = ScenarioEventsList(j);}
				if (Names(j) == "PWF") {PWF  = ScenarioEventsList(j);}
				if (Names(j) == "PWF") {PWF  = ScenarioEventsList(j);}
				if (Names(j) == "NeutralDragCoeff") {NeutralDragCoeff  = ScenarioEventsList(j);}
				if (Names(j) == "SWExtCoef") {SWExtCoef  = ScenarioEventsList(j);}
				if (Names(j) == "AWCss") {AWCss  = ScenarioEventsList(j);}
				if (Names(j) == "MetData") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "HypsoData") {Lakevolumes  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "Interpolationtype") {Interpolationtype  = as<std::string>(ScenarioEventsList(j));}
				if (Names(j) == "Topogrid") {Topogrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWRechargegrid") {GWRechargegrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWPumpinggrid") {GWPumpinggrid = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWKgrid") {GWKgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWSgrid") {GWSgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWbasegrid") {GWbasegrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWheadgrid") {GWheadgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWboundarygrid") {GWboundarygrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWlakesedimentgrid") {GWlakesedimentgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "Lakecondgrid") {Lakecondgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWETgrid") {GWETgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "Classgrid") {Classgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "ETexdepth") {ETexdepth  = ScenarioEventsList(j);}
				if (Names(j) == "Cellxdim") {Cellxdim  = ScenarioEventsList(j);}
				if (Names(j) == "Cellydim") {Cellydim  = ScenarioEventsList(j);}
				if (Names(j) == "Accuracy") {Accuracy  = ScenarioEventsList(j);}
				if (Names(j) == "Maxiterations") {Maxiterations  = ScenarioEventsList(j);}
				if (Names(j) == "GWcatchmentgrid") {GWcatchmentgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWlakekgrid") {GWlakekgrid  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "LakeDepths") {LakeDepths  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "DirectExternalFlows") {DirectExternalFlows  = as<std::string>(ScenarioEventsList(j));}
				if (Names(j) == "AtmosphericShift") {AtmosphericShift  = ScenarioEventsList(j);}
				if (Names(j) == "Turbulence") {Turbulence  = ScenarioEventsList(j);}
				if (Names(j) == "Feedback") {Feedback  = ScenarioEventsList(j);}
				if (Names(j) == "StratData") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "GWData") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "PMod") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "TaMod") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "TwMod") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "WSMod") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "RHMod") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "RsMod") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}				
				if (Names(j) == "GWMod") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
				if (Names(j) == "Stratmod") {Met  = as<NumericMatrix>(ScenarioEventsList(j));}
			//Rcpp::Rcout << "AWCss " << AWCss << std::endl;
			//Rcpp::Rcout << "AWCds " << AWCds << std::endl;
			//Rcpp::Rcout << "MetP " << Met(Time,3) << std::endl;
			}
		}
	}
		
//Run Massbal
		SalinityCalcs (Time-1, ChemConcentrations, Chemistry, ChemParams);
		if (StratificationModule == "Yes" && Time < Met.nrow()){ //We have to stop te Stratification module early, as it looks forward one timestep.
		func_RunStratification (Time-1,  Met, Lake, ChemConcentrations, Stratification, LatitudeRadians, NeutralDragCoeff, SWExtCoef, ThermoclineTargetThickness, ThermoclineMaximum, Feedback);
		}
		func_Fluxes(KcSS, ALBearth, ALBlake, CA, Cin, Csr, KcDS, KcLsum, KcLwin, LatitudeRadians, RunoffRatio, PWF, Timestep, AWCds, AWCss, Time-1, Flux, Isotopes, Lake, Met, Lakevolumes, Datarate, Interpolationtype, GroundwaterModule, Topogrid, GWRechargegrid, GWPumpinggrid, GWKgrid, GWSgrid, GWbasegrid, GWheadgrid, GWboundarygrid, GWlakesedimentgrid, Lakecondgrid, GWETgrid, Classgrid, ETexdepth, Cellxdim, Cellydim, Accuracy, Maxiterations, GWcatchmentgrid, GWlakekgrid, LakeDepths, GWStats, GWNumbers, ChemConcentrations, StratificationModule, NeutralDragCoeff, GWLakeID, Feedback);
		func_Hydrology(Flux, Lake, Isotopes, Time, Met, Lakevolumes, Interpolationtype, DirectExternalFlows, GroundwaterModule, GWStats);
		func_IsotopeFluxes (Flux, Time, Isotopes, Met, AtmosphericShift, Turbulence, Feedback, Lakevolumes, Lake, Interpolationtype, DirectExternalFlows, FluxD, Flux18O, GroundwaterModule, GWStats, GWNumbers, ChemConcentrations, Theta);
		func_ChemistryFluxes (Chemistry, Time, Flux, Lake, Met, DirectExternalFlows, GroundwaterModule, GWStats, GWNumbers);
		func_LakeMixer (Flux, Lake, Chemistry, Isotopes, Time, Met, Lakevolumes, Interpolationtype, FluxD, Flux18O);
	}
	//and run Lakemixer once more for the final timestep. At this stage we're outside the scope of the Met file but Fdlm has been calculated.
	
}