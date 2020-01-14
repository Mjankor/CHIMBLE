#include <Rcpp.h>
#include "Hydrology_functions.hpp"
#include "Stratification.hpp"


//notes
//While some ice calcs are included, not all are completed yet
//When ice is included in the model, add an ice and snow column to the Lake matrix to store ice and snow thicknesses for each timestep.


using namespace Rcpp;

long double const pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899863;

//Lake specific parameters ********************************************
int StratTimestep = 86400;	//Timestep for the stratification module (in seconds)
double Obliquity = 23.4;      // obliquity
double Latitude = -6.30;      // latitude (negative for South)
double Longitude = 29.5;      // longitude (negative for West)
//double NeutralDragCoeff = .001;      //*** neutral drag coefficient 1.8 HAD 1.7GISS 1.2CCSM - Moved to preferences
//double SWExtCoef = 0.25;     //*** shortwave extinction coefficient (1/m) was 0.065 //moved to preferences
double Frac_Advected = 0.3;      // fraction of advected air over lake
double AlbedoSlush = 0.4;      // albedo of melting snow
double AlbedoSnow = 0.7;      // albedo of non-melting snow

//Simulation specific parameters **************************************
int ModelSpinupYears = 30;      // number of years for spinup
bool BoundaryFlag = 0;      // true for explict boundry layer computations; presently only for sigma coord climate models
double Sigma = 0.9925561;      // sigma level for boundary flag
bool MassBal_flag = 0;      // true for variable lake depth
bool VarIceflag = 0;      // true for variable ice cover
bool VarSalinity = 0;      // true for variable salinity
double MetHeight = 2;    // height of met inputs

//Other parameters DO NOT CHANGE without good reason for doing so******
double LayerThickness = 1.0;      // vertical layer thickness in m
double WaterSurfRoughness=0.0004;      // water surface roughness length
double StefanBoltzmannDel = 0.0000000567;      // s-b constant
double PureWaterDensity = 1000;      // density of water (g) (rhowat)
double SnowDensity = 330;      // density of snow (g)
double IceDensity = 917;      // density of ice
double SnowCutOff = 2.2;      // temp at which precip is snow
double SpecHeatCapAir = 1004.64;   // specific heat capacity dry air
double SGDryAir = 287.05;      // specific gas constant dry air
double SGWaterVapour = 461.495;  //specific gas constant water vapor
double SpecHeatCapVapor = 1810;     // specific heat capacity water vapor
double LatentHeatWater = 2450000;    //latent heat of vaporization water
double LatentHeatIce = 2500000;      //latent heat of vaporization ice
double LWEmmisivity = 0.97;      // longwave emmisivity
double LatentHeatFusion = 334000;      // latent heat of fusion
double SurfaceThickness = 1;      // surface thickness over which the flux exchange occurs in
double MinIceThick= 0.01;      // min ice thick in meters
double  MinIceFraction = 0.02;      // min ice fraction
double PollardSubIce = 86400.;      // D. Pollard sub-ice time constant
double HeatCapIce = 4200.;      // heat capacity of ice
double ConductivityWater = 0.58;      // conduction coefficinet for water
double SWSurfaceAbsorbed = 0.4;      // frac of solar rad abs in water surface layer
double AlbWeightSW = 0.5; // sw weighting for surface albedo
double AlbWeightLW = 0.5; // lw weighting for surface albedo
float SeaIceAlbSW = 0.6; // sw albedo for sea ice
float SeaIceAlbLW = 0.4; // lw albedo for sea ice
double VonKarman = 0.4; //vonkarman constant
double DiffusivityWater = 0.000000138889; // thermal molecular diffusivity of water m2/s
double CritSnowDepth = 0.05;      // Snow height for albedo
double IsotopeExchangeDepth  = 1.0;    // fixed depth (m) over which isotope exchange occurs
double LamdSWIce = 3.0;      // extinction coef for sw thru ice (lamisw)
double LamdLWIce = 20.0;     // extinction coef for lw thru ice (lamilw)
double LamdSWSnow = 3.0;      // extinction coef for sw thru snow (lamssw)
double LamdLWSnow = 20.0;      // extinction coef for lw thru snow (lamslw)
double VisibleFraction = 0.7;      // fraction of light in visible band for ice model (afrac1)
double IRFraction = 0.3;      // fraction of light in infrared band for ice model (afrac2)
double CondIce = 2.3;      // thermal conductivity of ice (condi)
double CondSnow = 0.31;     // thermal conductivity of snow (conds)
double Gravity = 9.80616;      // gravity (grav)
double LakeAlbedoLW = 0.03;     // Was 0.3 (Hostetler 1990, eq 24) - Revised to 0.03 (Henderson Sellers 1986)


// [[Rcpp::export]]
double func_Density (double LakeTemperature, double Salinity){
double WaterDensity; // (rhow)
double A;
double B;
double C;
//calculate density based on Dickson 2007. Guide to Best Practices for Ocean CO2 Measurements, Section 4.2. Assumes seawater. Returns the delta of density from 1000. eg density = 1024, returns 24.
//checked against density.f90
//Temp = C, Salinity = g/KG
WaterDensity = 999.842594 + 6.793952 * pow(10,-2) * LakeTemperature - 9.095290* pow(10,-3) * pow(LakeTemperature,2) + 1.001685*pow(10,-4) * pow(LakeTemperature,3) - 1.120083*pow(10,-6) * pow(LakeTemperature,4) + 6.536332*pow(10,-9) * pow(LakeTemperature,5);

A = 0.824493 - 4.0899* pow(10,-3) * LakeTemperature + 7.6438 * pow(10,-5) * pow(LakeTemperature,2) - 8.2467 * pow(10, -7) * pow(LakeTemperature, 3) + 5.3875 * pow(10,-9) * pow(LakeTemperature,4);

B = -5.72466 * pow(10,-3) + 1.0227 * pow(10,-4) * LakeTemperature - 1.6546 * pow(10,-6) * pow(LakeTemperature,2); 

C = 4.8314 * pow (10,-4);

WaterDensity = WaterDensity + A*Salinity + B * pow(Salinity, 1.5) + C * pow(Salinity,2);
WaterDensity = WaterDensity - 1000;
return WaterDensity;
}

// [[Rcpp::export]]
double func_LakeDrag (double LakeTemperature, double AirTemperature, double WindSpeed, double NeutralDragCoeff){
//checked against LakeDrag.f90
double RichardsonNeutralBulk;
double Drag;
double MinimumDrag = 0.0006;//was 0.0006

if (BoundaryFlag != 0){
MetHeight = 2;
}

//temps should be in K for this calculation
LakeTemperature = LakeTemperature + 273.15;
AirTemperature = AirTemperature + 273.15;

if (LakeTemperature/AirTemperature <= 1){ //unstable conditions
	RichardsonNeutralBulk = (MetHeight * Gravity * (1 - LakeTemperature/AirTemperature))/
	(pow(WindSpeed,2) + .01);
	} else { //stable conditions
	RichardsonNeutralBulk = (MetHeight * Gravity * (1 - LakeTemperature/AirTemperature))/
	(pow(WindSpeed,2) + 1);
}

if (RichardsonNeutralBulk < 0){ //unstable conditions
	Drag = NeutralDragCoeff * (1 + 24.5 * pow(-NeutralDragCoeff * RichardsonNeutralBulk,0.5));
	} else {
	Drag = NeutralDragCoeff / (1 + 11.5 * RichardsonNeutralBulk);
}

if (0.25 * NeutralDragCoeff > MinimumDrag){
	MinimumDrag = 0.25 * NeutralDragCoeff;
}

if (Drag < MinimumDrag){
	Drag = MinimumDrag;
}

return Drag;
}

// [[Rcpp::export]]
double func_LakeAlbedo (int Julian, double LatitudeRadians){
double Albedo_Lake;
//checked against albedo.f90
	if (LatitudeRadians < 0){ //adjust for southern hemisphere
	Albedo_Lake = 0.08 + 0.02 * sin (2 * pi * (Julian+182.5)/365 + pi/2);
	} else {
	Albedo_Lake = 0.08 + 0.02 * sin (2 * pi * Julian/365 + pi/2);
	}
return Albedo_Lake;
}

// [[Rcpp::export]]
double func_IceAlbedo (double AirTemperature, double FreezingTemp){
//checked against albedo.f90
//temps are agnostic. C or K.
double TempDelta;
double Albedo_Ice;
 //From BATS, calculates fragmented albedos (direct and diffuse) in wavelength regions split at 0.7um.
TempDelta = AirTemperature - FreezingTemp;
if (TempDelta < 0){
	TempDelta = 0;
}
if (TempDelta > 20){
	TempDelta = 20;
}

Albedo_Ice = AlbWeightSW * (SeaIceAlbSW - 0.011 * TempDelta) + AlbWeightLW * (SeaIceAlbLW - 0.0245 * TempDelta);
return Albedo_Ice;
}

// [[Rcpp::export]]
double func_SnowAlbedo (bool MeltFlag){
//checked against albedo.f90
double Albedo_Snow;
	if (MeltFlag == 1){
	Albedo_Snow = 0.4; //slush
	} else {
	Albedo_Snow = 0.7; //snow
	}
return Albedo_Snow;
}


// [[Rcpp::export]]
double func_Teten(double Temperature){
double Teten_A;
double Teten_B;
double TetenC = 610.78; //gives result in Pa.
double VapourPressure;
if (Temperature >= 0){
	Teten_A = 17.269; 
	Teten_B = 237.3; // or if using K for temp -35.86;
	} else {
	Teten_A = 21.874;
	Teten_B = 265.5; // or if using K for temp -7.66;
}
VapourPressure = TetenC * exp(Teten_A * (Temperature)/(Temperature + Teten_B));
return VapourPressure;
}



// [[Rcpp::export]]
NumericVector func_LatSensHeat (double LakeTemperature, double AirTemperature, double WindSpeed, double LakeIceThickness, double FreezingTemp, double SurfacePressure, double RH, double NeutralDragCoeff){

//RH and SurfacePressure can be acquired from the Met file.
//RH as fraction
//Checked against latsens.f90.
//Watch units. Temp in C used, with a +273.15 as required. Millibars for input pressure, with adjustment as needed via *100. Humidity input as a percent. These should be the same as the MET file. 
//Unlike env.f90, the calculation for the latent heat flux is calculated within the function. Evaporation rate is still calculated using Penman's in the hydrological functions file.

double LakeSatVapourPressure;
double LakeSpecHumidity;
double AirSatVapourPressure;
double AirVapourPressure;
double AirSpecHumidity;
double RelativeHumidity;
double HumidityGradient; 
double PartialVapourPressure;
double DryAirPressure;
double AirDensity;
double Teten_A;
double Teten_B;
double TetenC = 610.78; //gives result in Pa.
double Evap;
double Sensible;
double Latent;

double Drag = func_LakeDrag(LakeTemperature, AirTemperature, WindSpeed, NeutralDragCoeff);

if (LakeIceThickness == 0 && LakeTemperature > FreezingTemp){
	Teten_A = 17.269; 
	Teten_B = 237.3; // or if using K for temp -35.86;
	} else {
	Teten_A = 21.874;
	Teten_B = 265.5; // or if using K for temp -7.66;
}
//moved from datain to here. All in pascals
AirSatVapourPressure = TetenC * exp(Teten_A * (AirTemperature)/(AirTemperature + Teten_B)); 
LakeSatVapourPressure = TetenC * exp(Teten_A * (LakeTemperature)/(LakeTemperature + Teten_B));
 //return in Mb. Temp in C. 
// calculate air density - CRC Standard Mathematical Tables and Formulas Page 694.
//Also https://www.vaisala.com/sites/default/files/documents/Humidity_Conversion_Formulas_B210973EN-F.pdf
AirVapourPressure = RH * AirSatVapourPressure;
AirSpecHumidity = 0.622 * (AirVapourPressure/(SurfacePressure*100 - 0.378*AirVapourPressure)); 
LakeSpecHumidity = 0.622 * (LakeSatVapourPressure/(SurfacePressure*100 - 0.378*LakeSatVapourPressure));

RelativeHumidity = AirSpecHumidity/LakeSpecHumidity; //as a fraction
PartialVapourPressure = RelativeHumidity * LakeSatVapourPressure;
DryAirPressure = SurfacePressure*100 - PartialVapourPressure; //All in Pascals.

AirDensity = DryAirPressure/(SGDryAir * (AirTemperature + 273.15)) + PartialVapourPressure/(SGWaterVapour*(AirTemperature + 273.15)); //eq 4a - https://wahiduddin.net/calc/density_altitude.htm.
HumidityGradient = AirSpecHumidity - LakeSpecHumidity;
Evap = -(Drag * WindSpeed * AirDensity) * HumidityGradient; // From Small et al 1999, eq, 2.
Latent = -Evap*LatentHeatWater;
Sensible = Drag * WindSpeed * AirDensity * SpecHeatCapAir * (AirTemperature - LakeTemperature);
NumericVector Heat = NumericVector::create(Named("Latent") = Latent, Named("Sensible") = Sensible, Named("DailyEvapRate") = Evap*StratTimestep); //Latent, Sensible and Evap (in mm/timestep)
return Heat;
}

// [[Rcpp::export]]
double func_SpecificHeat (double LakeTemperature, double Salinity){
double PureWaterSpecHeat;
double LakeSpecificHeat;
//Checked against Specheat.f90
//from https://www.coaps.fsu.edu/RVSMDC/marine_workshop3/docs/Appendixes_A,B,C.pdf 
//Salinity is in PSU from 0 to 40; This could be a problem. PSU = PPT, but range is limited. Temp is in C. Salinity is in g.kg. Output is in J/Kg.K
PureWaterSpecHeat = 4217.4 - 3.720283 * LakeTemperature + 0.1412855 * pow(LakeTemperature,2) - 2.654387*pow(10,-3) * pow(LakeTemperature,3) + 2.093236 * pow(10,-5) * pow(LakeTemperature,4);
LakeSpecificHeat = PureWaterSpecHeat + Salinity * (-7.64444 + 0.10728*LakeTemperature - 1.384 * pow(10,-3) * pow(LakeTemperature,2)) + pow(Salinity,1.5) * (0.177 - 4.08 * pow(10,-3)*LakeTemperature + 5.35 * pow(10,-5)*pow(LakeTemperature,2));

return LakeSpecificHeat;
}

// [[Rcpp::export]]
NumericVector TridiagSolver (int NumSystems, int Columns, int NumUnknowns, NumericMatrix SubdiagMatrix, NumericMatrix MainDiagMatrix, NumericMatrix SuperDiagMatrix, NumericMatrix RHS){

//check code in R
/*
a <- matrix (nrow = 3, ncol = 1, data = -1)
b <- matrix (nrow = 3, ncol = 1, data = 3)
c <- matrix (nrow = 3, ncol = 1, data = -1)
RHS <- matrix(nrow = 3, ncol = 1)
RHS[1,1] = -1
RHS[2,1] = 7
RHS[3,1] = 7
TridiagSolver(1,1,nrow(RHS),a,b,c,RHS)
Answers = 0.952, 3.857, 3.619;
https://www.slideshare.net/ParidhiMadurar/thomas-algorithm
*/

int i=0;
int j=0;
int Shortdim = NumUnknowns -1;
int Shortdim2 = 0;
//unlike the Fortran code this is based on this expects a row for each value, and a column for the "Ydim".
NumericMatrix alpha(NumUnknowns, Columns);
NumericMatrix gamma(NumUnknowns, Columns);
NumericMatrix Solution(NumUnknowns, Columns);

//start with lu decompositions;
//initialise boundaries
for (j = 0; j < NumSystems; j++){
	alpha(0,j) = 1/MainDiagMatrix(0,j);
	gamma(0,j) = SuperDiagMatrix(0,j) * alpha(0,j);
}

//fill in 
for (i=1 ; i < Shortdim ; i++){
	for (j = 0; j < NumSystems; j++){
		alpha(i,j) = 1/(MainDiagMatrix(i,j) - SubdiagMatrix(i,j)*gamma(i-1,j));
		gamma (i,j) = SuperDiagMatrix(i,j)*alpha(i,j);
	}
}   

for (j = 0; j < NumSystems; j++){
	Solution(0,j) = RHS(0,j)*alpha(0,j);
}

for (i=1 ; i < Shortdim ; i++){
	for (j = 0; j < NumSystems; j++){
	Solution(i,j) = (RHS(i,j) - SubdiagMatrix(i,j)*Solution(i-1,j))*alpha(i,j);
	}
}


for (j = 0; j < NumSystems; j++){
	Solution(Shortdim,j) = (RHS(Shortdim,j) - SubdiagMatrix(Shortdim,j)* Solution(Shortdim-1,j))/(MainDiagMatrix(Shortdim,j) - SubdiagMatrix(Shortdim,j)*gamma(Shortdim-1, j));
}

for (i = 0;i < Shortdim; i ++) {
	Shortdim2 = Shortdim - (i+1);
	for (j=0;j < NumSystems; j ++){
	Solution(Shortdim2,j) = Solution(Shortdim2,j) - gamma(Shortdim2,j)*Solution(Shortdim2+1,j);
	}
}	
	return Solution;
}


//NumericVector func_EddyDiffusion (double LakeIceThickness, double Windspeed, NumericMatrix Stratification, NumericVector LayerTemp, NumericVector LayerDensity, double LakeDepth, NumericVector Salinity, double LatitudeRadians){
//Update - moved into the temp_profile equation to cut down on the doubling up of loops

// [[Rcpp::export]]
NumericVector func_IceSWRad (double SW, double LakeIceThickness, double SnowThickness){
//Calculates terms needed for surface energy balance in presence of lake ice and snow equations based on Patterson and Hamblin (1988) Limnology and Oceanography v. 33 p. 323, eq7.
double Conductance = (SnowThickness*CondIce + LakeIceThickness*CondSnow)/(CondIce*CondSnow); //this looks wrong. I'd expect it Ice and Snow to be coupled, rather than split.
double a = (1-exp(-LamdSWSnow*SnowThickness)) / (CondSnow*LamdLWSnow);
double b = exp(-LamdSWSnow*SnowThickness)*(1-exp(-LamdSWIce*LakeIceThickness))/(CondIce*LamdSWIce);
double c = (1-exp(LamdLWSnow*SnowThickness))/(CondSnow*LamdLWSnow);
double d = exp(-LamdSWSnow*SnowThickness)*(1-exp(-LamdLWIce*LakeIceThickness))/(CondIce*LamdLWIce);
double val = SW * VisibleFraction*(a+b) + SW * IRFraction*(c+d);
double val2 = SW * VisibleFraction*(1-exp(-(LamdSWSnow * SnowThickness+LamdSWIce * LakeIceThickness))) - SW * IRFraction*(1-exp(-(LamdLWSnow * SnowThickness + LamdLWIce * LakeIceThickness)));

NumericVector IceRad = NumericVector::create(Named("Conductance") = Conductance, Named("val") = val, Named("val2") = val2);
return IceRad;
}


// [[Rcpp::export]]
double func_FreezingPoint(double SurfacePressure, double Salinity){
//calculates freezing point, in Gill 1982, from Millero 1978
double PressureDecibars = SurfacePressure/10000; //convert Millibars to Decibars(as per env.F90)
double FreezingPoint = -0.0575*Salinity + .001710523*pow(Salinity,1.5) - 0.0002154996*pow(Salinity,2) - 0.000753*PressureDecibars;
return FreezingPoint;
}


// [[Rcpp::export]]
double func_OldDownwardLongwave(double Ra, double Rs, double AirTemperature, double WaterVaporPressure){
//This function has been replaced with an updated version using Prata and Crawfords technique. Sugita and Brutsaert tend to underestimate longwave. 
//Based on eq 10 from SUGITA M. & BRUTSAERT W. 1993. Cloud effect in the estimation of instantaneous downward longwave radiation. Water Resources Research. 29, 599-605.
//validated SANTOS C. A. C. D., SILVA B. B. D., RAO T. V. R., SATYAMURTY P. & MANZI A. O. 2011. Downward longwave radiation estimates for clear-sky conditions over northeast Brazil. Revista Brasileira de Meteorologia. 26, 443-450. 
// Description of revised coefficients (and validation) DUARTE H. F., DIAS N. L. & MAGGIOTTO S. R. 2006. Assessing daytime downward longwave radiation estimates for clear and cloudy skies in Southern Brazil. Agricultural and Forest Meteorology. 139, 171-181.
//Temp needs to be K. Water vapor pressure in Pa.
//Using the method described in Prata, A. J. (1996). A new long-wave formula for estimating downward clear-sky radiation at the surface. Q. J . R. Meteorol. Soc., 122, 1127-1151.
double ClearnessIndex = Rs/Ra;
double ClearskyEmissivity = 0.714*pow((WaterVaporPressure/(AirTemperature+273.15)),0.0687);
double ClearskyLongWave = ClearskyEmissivity*StefanBoltzmannDel * pow(AirTemperature+273.15,4);
double DownwardLongWave = 1.02*ClearskyLongWave * pow(ClearnessIndex, -0.0227);
return DownwardLongWave;
}


// [[Rcpp::export]]
double func_DownwardLongwave(long Time, double LatitudeRadians, NumericMatrix Met){
//The Clearsky longwave is calculated uing the method described in Prata, A. J. (1996). A new long-wave formula for estimating downward clear-sky radiation at the surface. Q. J . R. Meteorol. Soc., 122, 1127-1151.
//The effect of clouds is calculated using the method from Crawford, T. M., and C. E. Duchon (1999). An improved parameterization for estimating effective atmospheric emissivity for use in calculating daytime downwelhng longwave radiation, .1. .4pp/. Meteorol., 48, 474 480.
//This is one of the most accurate longwave cloudy sky models as determined by: Flerchinger. G. N., W. Xaio, D. Marks. T. J. Sauer, and Q. Yu (2009), Comparison of algorithms for incoming atmospheric long-wave radiation, Water Resour Res., 45, W03423, doi: 10.1 029/2008WR0073 94.

//Clearsky Longwave.
double WaterMass;
double OpticalAirMass;
double ClearskyEmissivity;
double CloudyEmissivity;
double AirVapourPressure;
double DownwardLongWave;
AirVapourPressure = Met(Time,5)/100 * func_Teten(Met(Time,4));
WaterMass = 0.465*(AirVapourPressure/(Met(Time,4)+273.15));
ClearskyEmissivity = 1 - (1+WaterMass)*exp(-(pow(1.2 + 3 * WaterMass,0.5)));

//Next we calculate the ClearSky Shortwave radiation to determine a cloudiness index (referred to as solar index in literature).
double SinAlpha;
double SolarDec;
int DayofYear;
double SolarConstant = 1360;
double Noon = 12; //as we're calculating whole day Rs, the "noon" time can be any value.
double TrTpg;
double Tw;
double Ta;
double HourlyClrSWRad = 0;
double DailyClrSWRad = 0;
double SolarIndex;

	DayofYear = func_DayofYear(Met(Time,1), Met(Time,2));
	SolarDec = 0.409*sin(2*pi*DayofYear/365-1.39);
	
	for (int i=0;i<24;i++){ //average over 1 day
		SinAlpha = sin(LatitudeRadians) * sin(SolarDec)+cos(LatitudeRadians)*cos(SolarDec)*cos(pi*(i - Noon)/12);
		OpticalAirMass = 35*SinAlpha*pow((1224*pow(SinAlpha,2)+1),-0.5);
		Tw = 1-0.077*pow(WaterMass*OpticalAirMass,0.3);
		TrTpg = 1.021 - 0.084*pow(OpticalAirMass*(0.000949*Met(Time,14)+0.051),0.5);
		Ta = pow(0.932,OpticalAirMass);
		HourlyClrSWRad = SolarConstant * SinAlpha * Tw * TrTpg * Ta;
		if (HourlyClrSWRad > 0){
		DailyClrSWRad += HourlyClrSWRad;
		}
	}
	
	DailyClrSWRad = (DailyClrSWRad / 24)/1000000 * 24*60*60; //averaged and converted to MJ/M2/Day
	
//Solar index calc
SolarIndex = Met(Time,6)/DailyClrSWRad;
CloudyEmissivity = (1-SolarIndex) + SolarIndex * ClearskyEmissivity;
DownwardLongWave = CloudyEmissivity * StefanBoltzmannDel * pow(Met(Time,4)+273.15,4);
return DownwardLongWave;
}

// [[Rcpp::export]]
double func_LongWave (double LakeIceThickness, double WaterTemperature, double TempIce, double Ra, double Rs, double AirTemperature, double RH, NumericMatrix Met, long Time, double LatitudeRadians){
	double UpwardLWRad;
	double DownwardLWRad;
	double LakeNetLWRad;
	double AirVapourPressure;
	AirVapourPressure = RH * func_Teten(AirTemperature);
	DownwardLWRad = func_DownwardLongwave(Time, LatitudeRadians, Met);
if (LakeIceThickness == 0){
	UpwardLWRad = LWEmmisivity * StefanBoltzmannDel * pow(WaterTemperature+273.15,4); //as W/m2
} else {
	UpwardLWRad = LWEmmisivity * StefanBoltzmannDel * pow(TempIce+273.15,4);	//as W/m2

}
DownwardLWRad = DownwardLWRad*(1 - LakeAlbedoLW); //Already in W/m2. Albedo from Hostetler 1990 eq, 23)
LakeNetLWRad = (DownwardLWRad - UpwardLWRad);
return LakeNetLWRad; //as W/m2
}



// [[Rcpp::export]]
void func_TempProfile (double LakeIceThickness, double LakeDepth, double LakeSurfaceArea, double LakeTemperature, double AirTemperature, double WindSpeed, double SurfacePressure, double RH, double Salinity, double TempIce, double Ra, double Rs, double LatitudeRadians, int Julian, NumericMatrix Stratification, double NeutralDragCoeff, double SWExtCoef, NumericMatrix Met, long Time){
//test code func_TempProfile(0,161.1,5520860,23,25,2,1020,75,0,0,40.65991,22.419,-0.6673218, 2,Stratification)

int LayerIndexStart = 0;
double t1=0;
double ShortwaveRad =0; 	
double UpperArea;
double LowerArea;
double LayerThickness;
double SurfaceThickness;
double AverageLayerArea;
double LayerTopDepth;
double LayerBottomDepth;
double CNExtra;
//eddy variables
double ks; 
double ws;
double z;
double Po = 1;
double radmax = 40000;
double BruntVaisala; //N2
double DensityDelta;
double rad;
double RichardsonNumber;

double FreezingTemp = func_FreezingPoint(SurfacePressure, Salinity);
NumericVector LayerDeltas (Stratification.nrow());
NumericVector Heat = func_LatSensHeat (LakeTemperature, AirTemperature, WindSpeed, LakeIceThickness, FreezingTemp, SurfacePressure, RH, NeutralDragCoeff);

double LakeNetLWRad = func_LongWave (LakeIceThickness, LakeTemperature, TempIce, Ra, Rs, AirTemperature, RH, Met, Time, LatitudeRadians);  // converted to Mj/m2/day to W/m2 in the function. ;
	//The stratification matrix has rows for each layer, and columns for LowerDepth, UpperDepth, Thickness, PreviousTemp, PreviousD, Previous18O, PreviousSalinity, Temp, D, 18O, Salinity.
	//Area in the matrix is the area of the lower level. 
		//There is going to be some quantisation of the model as the top layer is defined by a surface thickness parameter, rather than the thickness between the lake surface and the first complete layer. //In the future we'll need to assess the effect of having a variable thickness top layer.
		//Update. I'm changing this to use the surface layer thickness.
	
 	//area notes
 	//	 ktop = 999 - 100 = 899
 	//	 area_1 =(area(k+ktop) + area(k+ktop)) /2. = (899+1 + 899+1)/2 Surface
    //	area_2 =(area(k+ktop) + area(k+ktop+1)) / 2. = (899+1 + 899+2)/2 Surf-L1
    //  top = (surf+(k-2)*dz) = 1, 2dz, 3dz
	//   bot = (surf+(k-1)*dz) = 2, 3dz
	//   area_1 =(area(k+ktop-1) + area(k+ktop)) / 2 =  (899 + 1)+(899 + 2)/2 Upper area
	//   area_2 =(area(k+ktop) + area(k+ktop+1)) / 2. = (899 + 2)+(899+3)/2 Lower area

	//The stratification matrix is a matrix going from the lake minimum to maximum height. Therefore, the values used will typically ignore the first few rows of the matrix, and run from the first layer below the lake surface height and the bottom of the lake. Here we define that starting index and the distances between the midpoint of each layer (zhalf) and the height between the centre of each layer (LayerDeltas). 
	for (int i = 0; i < Stratification.nrow(); i++){
		LayerDeltas (i) = (Stratification(i,2) + Stratification(i+1,2))/2; //reset the layer thicknesses. //Layerdeltas are from i to i+1;
		if (LakeDepth < Stratification(i,1) + Stratification(i,2)/2){; // we subtract a half layer so the surface layer is a near approximation to the layer thickness. eg: Lake height of 161.1 would then have a surface layer of 1.1m, rather than 0.1m.
		LayerIndexStart = i; //define uppermost complete layer.
		}
	}
	LakeDepth = Stratification(LayerIndexStart,1);//need this to enforce a standard thickness surface layer. Remove this once variable thickness layers is possible.

	//Define surface layer thickness
	SurfaceThickness = LakeDepth - Stratification(LayerIndexStart,0);
	LayerDeltas(LayerIndexStart) = (SurfaceThickness + Stratification(LayerIndexStart+1,2))/2;
	LayerDeltas(LayerDeltas.size()-1) = 1; //remove for variable layers
	
	//Backup a bunch of data to the previous timestep and calculate updated densities and specific heat for each layer.
	for (int i = LayerIndexStart; i < Stratification.nrow(); i++){
	Stratification(i,3) = Stratification(i,7); //transfer current values to previous timestep. May use #defines to make these easier to identify later.
	Stratification(i,4) = Stratification(i,8); // D
	Stratification(i,5) = Stratification(i,9); // 18O
	Stratification(i,6) = Stratification(i,10); //
	//and call density and specific heat for each layer
	Stratification (i,11) = func_Density (Stratification(i,7), Stratification(i,10));
	Stratification (i,12) = func_SpecificHeat (Stratification(i,7), Stratification(i,10));
	}

//Prep data for Eddy Calcs
//Eddy uses difference in density over depth, before the eddy diffusion values are used in the temp_profile calculation. However, they both use the same for loops. 
	if (LakeIceThickness != 0){ //no ice, so set diffusivity to water diffusivity as a starting point.
		for (int k = 0; k < Stratification.nrow(); k++){
			Stratification(k,13) = 	DiffusivityWater;
		}
	} else {
	//do open water calculations
	//from Henderson Sellers 1985 eq16
	if (WindSpeed < 0.5){ //apparently here to stop ks = NAN
		WindSpeed = 0.5;
	}
	ks = 6.6*sqrt(fabs(sin(LatitudeRadians))) * pow(WindSpeed,-1.84);
	ws = 0.0012 * WindSpeed;

//Do Eddy calculations
//Eddy results are checked against Eddy.f90. Results are the same when z = LakeDepth - Stratification(k,0); , However, there are minor differences due to layer placement. In our model, layers run from (eg) 100 to 101 height giving an average layer depth of (assuming lake surface of 101) 0.5. In contrast, Eddy.f90 uses a layer centred on the integer, so 101 surface, with a layer from 100.5 to 99.5, and a z to the layer of 1m.
	for (int k = LayerIndexStart;k < Stratification.nrow(); k++){ //may need to be depth -1
		//Maybe eq 30 in Henderson Seller 1985.
		DensityDelta = (Stratification(k+1,11)-Stratification(k,11))/LayerDeltas(k); //check layers go right way
		BruntVaisala = (DensityDelta/(1000 + Stratification(k,11)))*Gravity;
		z = LakeDepth - (Stratification(k,1)+Stratification(k,0))/2; //average depth of layer
		//z = LakeDepth - Stratification(k,0); alternative (and slightly wrong) to match fortran model setup. 

		if (z * ks/ws > 40) {
			rad = radmax; // apparently to avoid a NAN.
		} else {
			rad = 1 + 40 *BruntVaisala * pow(VonKarman * z,2)/(pow(ws,2) * exp(-2*ks*z));
			if (rad > radmax){
				rad = radmax;
			}
		}
		if (rad < 1){ // so RichardsonNumber lower limit = 0;
			rad = 1; 
		}
		RichardsonNumber = (-1+sqrt(rad))/20; 
		if ((ks*z/ws) > 40){
			Stratification(k,13) = DiffusivityWater; //for upper layers
		} else {
			Stratification(k,13) = DiffusivityWater + VonKarman*ws*z*Po*exp(-ks*z)/(1+37*pow(RichardsonNumber,2));

		}
			Stratification(Stratification.nrow()-1,13)  = DiffusivityWater; //bottom layers
			Stratification(Stratification.nrow()-2,13) = Stratification(Stratification.nrow()-1,13); //apparently necessary for CN solution. 
	}
	}
// End Eddy calculations.


 //************************************************************************
//prep vectors for the Crank Nicholson and the Tridiag solver.
NumericMatrix RHS(Stratification.nrow() - LayerIndexStart,1);
NumericMatrix DiagMatrix(Stratification.nrow() - LayerIndexStart,1);
NumericMatrix SuperDiagMatrix(Stratification.nrow() - LayerIndexStart,1);
NumericMatrix SubDiagMatrix(Stratification.nrow() - LayerIndexStart,1);

 //************************************************************************
//Do first layer (surface layer)

//Modify the surface to initial layer to account for surface layer thickness,
UpperArea = LakeSurfaceArea; //current surface area 
LowerArea = Stratification(LayerIndexStart+1,14); //surface at bottom of layer.
AverageLayerArea = (UpperArea+LowerArea)/2; //mid layer area

if (LakeIceThickness == 0){ //open water
 ShortwaveRad = (Rs*1000000/(24*60*60))*(1-func_LakeAlbedo (Julian, LatitudeRadians));//convert Mj/m2/day to W/m2
 t1 = ShortwaveRad*SWSurfaceAbsorbed + (1 - SWSurfaceAbsorbed)*ShortwaveRad * (1-exp(-SWExtCoef*SurfaceThickness)) * LowerArea/UpperArea + (LakeNetLWRad + Heat(0) + Heat(1));// * UpperArea/UpperArea (which goes to 1);
 } else {
 //Do ice calc here
 }

CNExtra = 0.5 * LowerArea/UpperArea * (Stratification(LayerIndexStart,13)/LayerDeltas(LayerIndexStart)) * StratTimestep/SurfaceThickness * (Stratification(LayerIndexStart+1,7) - Stratification(LayerIndexStart,7));


RHS(0,0) = Stratification(LayerIndexStart,7) + t1*StratTimestep/((1000+Stratification(LayerIndexStart,11))*Stratification(LayerIndexStart,12)*SurfaceThickness) + CNExtra;

//Calc the A and B vectors
SuperDiagMatrix(0,0) = -0.5 * (Stratification(LayerIndexStart,13)/LayerDeltas(LayerIndexStart) * StratTimestep/SurfaceThickness * LowerArea/UpperArea);
DiagMatrix(0,0) = 1-SuperDiagMatrix(0,0);
 
 //************************************************************************
//Do rest of the water column;
for (int i = LayerIndexStart + 1;i < Stratification.nrow()-1; i++){

	//calculate critical dimensions for each layer
	LayerTopDepth = LakeDepth - Stratification(i,1);
	LayerBottomDepth = LakeDepth - Stratification(i,0);
	UpperArea = Stratification(i,14);
	LowerArea = Stratification(i+1,14);
	AverageLayerArea = (UpperArea+LowerArea)/2;
	LayerThickness = Stratification(i,2); //z
	
	if (LakeIceThickness == 0){
		t1 = (1-SWSurfaceAbsorbed)*ShortwaveRad * ((UpperArea*exp(-SWExtCoef*LayerTopDepth) - LowerArea * exp(-SWExtCoef*LayerBottomDepth))/AverageLayerArea);

		CNExtra = 0.5 * 1/AverageLayerArea * (((Stratification(i,13)/LayerDeltas(i)) * (StratTimestep/LayerThickness) * (Stratification(i+1,7) - Stratification(i,7))) * LowerArea - ((Stratification(i-1,13)/LayerDeltas(i-1))*(StratTimestep/LayerThickness) * (Stratification(i,7) - Stratification(i-1,7)))*UpperArea);

		RHS(i-LayerIndexStart,0) = Stratification(i,7) + t1*StratTimestep/((1000+Stratification(i,11)) * Stratification(i,12) * LayerThickness) + CNExtra;
	} else {
		// do ice calcs
	}

//Calc the A, B and C vectors
SuperDiagMatrix(i-LayerIndexStart,0) = -0.5 * (Stratification(i,13)/LayerDeltas(i)) * StratTimestep/LayerThickness * LowerArea/AverageLayerArea;
SubDiagMatrix(i-LayerIndexStart,0) = -0.5 * (Stratification(i-1,13)/LayerDeltas(i-1)) * StratTimestep/LayerThickness * UpperArea/AverageLayerArea;
DiagMatrix(i-LayerIndexStart,0) = 1 - SuperDiagMatrix(i-LayerIndexStart,0) - SubDiagMatrix(i-LayerIndexStart,0);
}

 //************************************************************************
//Do the lower layer;
int i = Stratification.nrow()-1;
LayerTopDepth = LakeDepth - Stratification(i,1);
LayerBottomDepth = LakeDepth - Stratification(i,0);
UpperArea = Stratification(i,14);
LowerArea = 0; //Check this. Divide by zero hurts.
AverageLayerArea = (UpperArea+LowerArea)/2;
LayerThickness = Stratification(i,2); //z
if (LakeIceThickness == 0){
t1 = (1-SWSurfaceAbsorbed)*ShortwaveRad * UpperArea*exp(-SWExtCoef*LayerTopDepth)/AverageLayerArea;
} else {
//do ice calcs.
}


CNExtra = 0.5* UpperArea/AverageLayerArea * (Stratification(i-1,13)/LayerDeltas(i-1)) * (StratTimestep/LayerThickness) * (Stratification(i,7) - Stratification(i-1,7));
RHS(i-LayerIndexStart,0) = Stratification(i,7) + t1*StratTimestep/((1000+Stratification(i,11)) * Stratification(i,12) * LayerThickness) + CNExtra;

//missing c and b calcs here.
SubDiagMatrix(i-LayerIndexStart,0) = -0.5 * (Stratification(i,13)/LayerDeltas(i)) * StratTimestep/LayerThickness * UpperArea/AverageLayerArea;
DiagMatrix(i-LayerIndexStart,0) = 1 - SubDiagMatrix(i-LayerIndexStart,0);



//Calculate Arrays for the tridiagonal matrix.
	//SubDiagMatrix(SubDiagMatrix.size(),0) = -0.5 * (Stratification(i,13)/LayerDeltas(i))* StratTimestep/LayerThickness * UpperArea/AverageLayerArea;
	//DiagMatrix(DiagMatrix.size(),0) = 1-SubDiagMatrix(SubDiagMatrix.size(),0); //these may need a -1 for c++ indexing.
NumericVector Solution = TridiagSolver (1, 1, (Stratification.nrow() - LayerIndexStart),  SubDiagMatrix, DiagMatrix, SuperDiagMatrix, RHS); //Errrrg. Watch the argument order.
	for (i=LayerIndexStart; i < Stratification.nrow(); i ++){
	Stratification(i,7) = Solution(i - LayerIndexStart); //update temps
	Stratification(i,11) = func_Density (Stratification(i,7), Stratification(i,10));
	Stratification (i,12) = func_SpecificHeat (Stratification(i,7), Stratification(i,10)); //update temp, density and specheat.
	}

}

// [[Rcpp::export]]
void func_ConvectionMixing (NumericMatrix Stratification, double LakeDepth){
int LayerIndexStart;
int MixStart;
double CombinedVol = 0;
double CombinedHeat = 0;
double CombinedHeatVol = 0;
double LayerHeatVol = 0;
double AveHeat = 0;
double AveSalinity = 0;
double AveDensity = 0;
int InstabilityDepth;
double CombinedSalinity = 0;
bool RestartLoopFlag = 1;
double UpperArea;
double LowerArea;
//Specheat and density done at the end of temp_profile.

	for (int i = 0; i < Stratification.nrow(); i++){
		if (LakeDepth < Stratification(i,1) + Stratification(i,2)/2){; 
			LayerIndexStart = i; //define uppermost complete layer.
		} 
	}

while (RestartLoopFlag == 1){
MixStart = LayerIndexStart;
CombinedHeat = 0;
CombinedHeatVol = 0;
CombinedSalinity = 0;
CombinedVol = 0;
RestartLoopFlag = 0;

	for (int j = LayerIndexStart; j < Stratification.nrow()-1; j++){ //define top of instability
		if (Stratification(j,11) > Stratification(j+1,11)){
		//define highest unstable layer;
		MixStart = j; 
		AveDensity = Stratification(j,11);	
		//toggle the restart loop flag. We assume that if any instability is found, that the loop will have to be rerun.
		RestartLoopFlag = 1;
		//define inital values for the layer at the top of the instability
		//areas of initial layer
		UpperArea = Stratification(j,14);
		LowerArea = Stratification(j,14);
		//define starting values
		LayerHeatVol = Stratification(j,2) * (1000+Stratification(j,11)) * Stratification(j,12) * (UpperArea + LowerArea)/2; //upper layer
		CombinedHeat = LayerHeatVol * Stratification(j,7);
		CombinedHeatVol = LayerHeatVol;
		CombinedSalinity = Stratification(j,2) * Stratification(j,10) * (UpperArea + LowerArea)/2;
		CombinedVol = Stratification(j,2) * (UpperArea + LowerArea)/2;

		//Run the averaging loop over the zone of instability, working from the highest layer downwards.
		for (int i = MixStart; i < Stratification.nrow()-1; i++){		
			if (Stratification(i+1,11) < AveDensity){
				//define areas for active layer
				UpperArea = Stratification(i+1,14);
				if (i+1 == Stratification.nrow()){
					LowerArea = 0;
				} else {
					LowerArea = Stratification(i+1,14);
				}
	
				LayerHeatVol = Stratification(i+1,2) * (1000+Stratification(i+1,11)) * Stratification(i+1,12) * (UpperArea + LowerArea)/2; //upper layer
				CombinedHeat += LayerHeatVol * Stratification(i+1,7);
				CombinedHeatVol += LayerHeatVol;
				CombinedSalinity += Stratification(i+1,2) * Stratification(i+1,10) * (UpperArea + LowerArea)/2;
				CombinedVol += Stratification(i+1,2) * (UpperArea + LowerArea)/2;

				//calculate averages
				AveHeat = CombinedHeat/CombinedHeatVol;
				AveSalinity = CombinedSalinity/CombinedVol;
				AveDensity = func_Density (AveHeat, AveSalinity);
				InstabilityDepth = i+1;
			} else {
			//Do something?
			}
		}
		
		//add new mixed layer averages back into the main matrix.
		for (int i = MixStart; i <= InstabilityDepth; i++){
			Stratification (i,7) = AveHeat;
			Stratification (i,11) = AveDensity;
			Stratification (i,10) = AveSalinity;
			Stratification (i,12) = func_SpecificHeat (Stratification(i,7), Stratification(i,10));
		}

		} else { //if no instabilities are found
		//RestartLoopFlag = 0;
		}
	}
	
	}
}

// [[Rcpp::export]]
double func_Interpolate(double StartValue, double EndValue, int Days, int CurrentDay){
//Days should be the from t+1 to t+n, where n = number of days (so if going from 1st to 6th of Jan, Days = 5 (2,3,4,5,6). Current day should be t+Current day.
double InterpVal = (EndValue - StartValue)/Days * (CurrentDay) + StartValue;
return InterpVal;
}

// [[Rcpp::export]]
double func_CleanAirAct (NumericMatrix Stratification, double LakeDepth){
//This function adds or removes layers from the stratification matrix to represent current lake level. If the lake level is higher than the current stratification layer, then the top stratification layer is duplicated upwards until the new lake level is reached.
//If the lake is lower, then layers above the lake level are reset to 0.
int LayerIndexStart = 0;
int CurrentStratificationTop = 0;
double TopLayerTemp;

	//define current layer index for current lake depth.
	for (int i = 0; i < Stratification.nrow(); i++){ 
		if (LakeDepth < Stratification(i,1) + Stratification(i,2)/2){; 
			LayerIndexStart = i; //define uppermost complete layer.
		}
	}
	//define layer index for top of current stratification matrix.
	for (int i = 0; i < Stratification.nrow(); i++){ 
		if (Stratification(i,7) != 0){; 
			CurrentStratificationTop = i;
			break; //define uppermost complete layer.
		}
	}

	if (CurrentStratificationTop > LayerIndexStart){
		for (int i = LayerIndexStart; i < CurrentStratificationTop; i++){
		Stratification(i,3) = Stratification(CurrentStratificationTop,3);
		Stratification(i,4) = Stratification(CurrentStratificationTop,4);
		Stratification(i,5) = Stratification(CurrentStratificationTop,5);
		Stratification(i,6) = Stratification(CurrentStratificationTop,6);
		Stratification(i,7) = Stratification(CurrentStratificationTop,7);
		Stratification(i,8) = Stratification(CurrentStratificationTop,8);
		Stratification(i,9) = Stratification(CurrentStratificationTop,9);
		Stratification(i,10) = Stratification(CurrentStratificationTop,10);
		Stratification(i,11) = Stratification(CurrentStratificationTop,11);
		Stratification(i,12) = Stratification(CurrentStratificationTop,12);
		Stratification(i,13) = Stratification(CurrentStratificationTop,13);
		}
		TopLayerTemp = Stratification(CurrentStratificationTop,7);
	} else {	
		for (int i = 0; i < LayerIndexStart; i++){
		Stratification(i,3) = 0;
		Stratification(i,4) = 0;
		Stratification(i,5) = 0;
		Stratification(i,6) = 0;
		Stratification(i,7) = 0;
		Stratification(i,8) = 0;
		Stratification(i,9) = 0;
		Stratification(i,10) = 0;
		Stratification(i,11) = 0;
		Stratification(i,12) = 0;
		Stratification(i,13) = 0;
		}
		TopLayerTemp = Stratification(LayerIndexStart,7);
	}
return TopLayerTemp;
}


// [[Rcpp::export]]
void func_RunStratification (long Time, NumericMatrix Met, NumericMatrix Lake, NumericMatrix ChemConcentrations, NumericMatrix Stratification, double LatitudeRadians, double NeutralDragCoeff, double SWExtCoef, double ThermoclineTargetThickness, double ThermoclineMaximum) {
//Get current values
double LakeDepth = Lake(Time, 4);
double LakeSurfaceArea = Lake(Time,6);
double LakeTemperature;
double AirTemperature;
double WindSpeed;
double SurfacePressure;
double RH;
double Rs;
double Ra;
int Julian;
int StartDay;
int EndDay;
int Days;
double RaStart = func_Ra(Time, LatitudeRadians, Met);
double RaEnd = func_Ra(Time+1, LatitudeRadians, Met);
double LakeIceThickness = 0;
double TempIce = 0;
double Salinity = ChemConcentrations(1,12);//assuming that RESsl is always present.
int LayerIndexStart = 0;
double ThermoclineDepth = 0;

//Number of days
StartDay = func_DayofYear(Met(Time,1), Met(Time,2));
EndDay = func_DayofYear(Met(Time+1,1), Met(Time+1,2));
	if (StartDay > EndDay){ //end of year
		EndDay += 365;
	}
	//Rcpp::Rcout << "Start and End " << StartDay << ", " <<EndDay << std::endl;

Days = EndDay - StartDay;
	for (int i = 0; i < Days; i++){
		AirTemperature = func_Interpolate(Met(Time,4), Met(Time+1,4), Days, i);
		WindSpeed = func_Interpolate(Met(Time,7), Met(Time+1,7), Days, i);
		SurfacePressure = func_Interpolate(Met(Time,14), Met(Time+1,14), Days, i);
		RH = func_Interpolate(Met(Time,5), Met(Time+1,5), Days, i)/100;
		Rs = func_Interpolate(Met(Time,6), Met(Time+1,6), Days, i);
		Ra = func_Interpolate(RaStart, RaEnd, Days, i);
		Julian = StartDay + i;
		if (Julian > 365){
			Julian = Julian - 365;
		}
	
		//Check for change in lake level vs the stratification matrix. If the lake level has dropped, then just reset any values above lake level to out of bounds. If lake level has risen to a new layer, duplicate the next valid lowest layer.
		//Also, grab current lake temp value.
		LakeTemperature = func_CleanAirAct (Stratification, LakeDepth);

		//Rcpp::Rcout << "TempProf Vars " << LakeIceThickness << ", " << LakeDepth << ", " << LakeSurfaceArea << ", " << LakeTemperature << ", " << AirTemperature << ", " << WindSpeed << ", " << SurfacePressure << ", " << RH << ", " << Salinity << ", " << TempIce << ", " << Ra << ", " << Rs << ", " << LatitudeRadians << ", " << Julian << std::endl;

		//Update Stratification Matrix
		func_TempProfile (LakeIceThickness, LakeDepth, LakeSurfaceArea, LakeTemperature, AirTemperature, WindSpeed, SurfacePressure, RH, Salinity, TempIce, Ra, Rs, LatitudeRadians, Julian, Stratification, NeutralDragCoeff, SWExtCoef, Met, Time);
		
		//And run the convection mixing.
		func_ConvectionMixing (Stratification, LakeDepth);
	}
	
	//Add final value for lake temp and stratification depth into the Met File.
	//still using temp offsets for the lake. Consider changing to direct temp values in the future. Currently seems to be only used in the isotope fractionation equation.
	Met(Time+1,10) = func_CleanAirAct (Stratification, LakeDepth) - Met(Time+1,4);
	
	//find the Thermocline
	for (int i = 0; i < Stratification.nrow(); i++){ 
		if (LakeDepth < Stratification(i,1) + Stratification(i,2)/2){; 
			LayerIndexStart = i; //define uppermost complete layer.
		}
	}
	
	for (int i = LayerIndexStart; i < Stratification.nrow() - ThermoclineTargetThickness; i++){ 
		if (Stratification(i, 7) - Stratification(i+ThermoclineTargetThickness, 7) >= ThermoclineMaximum){
			// This will pluck out the highest series of layers that match the characteristics of the thermocline we're anticipating.
			ThermoclineDepth = Stratification(LayerIndexStart, 1) - Stratification(ceil(i+(ThermoclineTargetThickness)/2),0);
			//Rcpp::Rcout << "LayerIndex and Thermo Depth " << Stratification(LayerIndexStart, 1) << ", " << Stratification(ceil((i+ThermoclineTargetThickness)/2),0) << std::endl;

			break; 
		} else {
			ThermoclineDepth = 0;
		}
	}
	//Rcpp::Rcout << "LayerIndex and Thermo Depth " << LayerIndexStart << ", " <<ThermoclineDepth << std::endl;
	Met(Time+1,15) = ThermoclineDepth;
}