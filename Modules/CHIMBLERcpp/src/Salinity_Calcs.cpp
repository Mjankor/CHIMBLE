#include <Rcpp.h>

using namespace Rcpp;

#define RESsl_NaMol	ChemConcentrations(1,0)
#define RESsl_MgMol	ChemConcentrations(1,1)
#define RESsl_CaMol	ChemConcentrations(1,2)
#define RESsl_KMol	ChemConcentrations(1,3)
#define RESsl_ClMol	ChemConcentrations(1,4)
#define RESsl_BrMol	ChemConcentrations(1,5)
#define RESsl_NaCl	ChemConcentrations(1,6)
#define RESsl_MgCl2	ChemConcentrations(1,7)
#define RESsl_CaCl2	ChemConcentrations(1,8)
#define RESsl_KCl	ChemConcentrations(1,9)
#define RESsl_NaBr	ChemConcentrations(1,10)

#define RESdl_NaMol	ChemConcentrations(3,0)
#define RESdl_MgMol	ChemConcentrations(3,1)
#define RESdl_CaMol	ChemConcentrations(3,2)
#define RESdl_KMol	ChemConcentrations(3,3)
#define RESdl_ClMol	ChemConcentrations(3,4)
#define RESdl_BrMol	ChemConcentrations(3,5)
#define RESdl_NaCl	ChemConcentrations(3,6)
#define RESdl_MgCl2	ChemConcentrations(3,7)
#define RESdl_CaCl2	ChemConcentrations(3,8)
#define RESdl_KCl	ChemConcentrations(3,9)
#define RESdl_NaBr	ChemConcentrations(3,10)

#define RESsl_Density	ChemConcentrations(1,11)
#define RESdl_Density	ChemConcentrations(3,11)
#define RESsl_Salinity	ChemConcentrations(1,12)
#define RESdl_Salinity	ChemConcentrations(3,12)
#define RESsl_ActWater	ChemConcentrations(1,13)
#define RESdl_ActWater	ChemConcentrations(3,13)

#define NaCl_MolMass	ChemConcentrations(4,6)
#define MgCl2_MolMass	ChemConcentrations(4,7)
#define CaCl2_MolMass	ChemConcentrations(4,8)
#define KCl_MolMass	ChemConcentrations(4,9)
#define NaBr_MolMass	ChemConcentrations(4,10)

// [[Rcpp::export]]
void CalcMolarities(double Time, NumericMatrix Chemistry, NumericMatrix ChemConcentrations){
for (int i = 0; i < 6; i++){
//RESslChemIndex & RESdlChemIndex need to decrease by 1 as the index is created in R.
	int RESslChemIndex = ChemConcentrations (0,i)-1;
	int RESdlChemIndex = ChemConcentrations (2,i)-1;
	ChemConcentrations (1,i) = (Chemistry (Time, RESslChemIndex) / 1000)/ChemConcentrations (4,i);
	ChemConcentrations (3,i) = (Chemistry (Time, RESdlChemIndex) / 1000)/ChemConcentrations (4,i);
}
} 

// [[Rcpp::export]]
void EstimateSaltMolarities(NumericMatrix ChemConcentrations){
//This routine just balance out any differences in ion concentrations. Sum of positive - sum of negative. The cations are adjusted to match the anions.
//Need a conditional in case of 0 volume reservoir.
	if (RESsl_NaMol + 2*RESsl_MgMol + 2*RESsl_CaMol + RESsl_KMol != 0){
		double RESsl_IonDiff = 1-((RESsl_NaMol + 2*RESsl_MgMol + 2*RESsl_CaMol + RESsl_KMol - RESsl_ClMol - RESsl_BrMol)/(RESsl_ClMol + RESsl_BrMol));
		RESsl_NaMol = RESsl_NaMol *RESsl_IonDiff;
		RESsl_MgMol = RESsl_MgMol *RESsl_IonDiff;
		RESsl_CaMol = RESsl_CaMol *RESsl_IonDiff;
		RESsl_KMol = RESsl_KMol *RESsl_IonDiff;
	}
	if (RESdl_NaMol + 2*RESdl_MgMol + 2*RESdl_CaMol + RESdl_KMol != 0){
		double RESdl_IonDiff = 1-((RESdl_NaMol + 2*RESdl_MgMol + 2*RESdl_CaMol + RESdl_KMol - RESdl_ClMol - RESdl_BrMol)/(RESdl_ClMol + RESdl_BrMol));
		RESdl_NaMol = RESdl_NaMol *RESdl_IonDiff;
		RESdl_MgMol = RESdl_MgMol *RESdl_IonDiff;
		RESdl_CaMol = RESdl_CaMol *RESdl_IonDiff;
		RESdl_KMol = RESdl_KMol *RESdl_IonDiff;
	}
}

// [[Rcpp::export]]
void CalcSalts(NumericMatrix ChemConcentrations){
//This routine calcs the salt molarity for the 5 keys salts considered for activity of water.
RESsl_NaCl = RESsl_NaMol - RESsl_BrMol;
RESsl_MgCl2 = RESsl_MgMol;
RESsl_CaCl2 = RESsl_CaMol;
RESsl_KCl = RESsl_KMol;
RESsl_NaBr = RESsl_BrMol;

RESdl_NaCl = RESdl_NaMol - RESdl_BrMol;
RESdl_MgCl2 = RESdl_MgMol;
RESdl_CaCl2 = RESdl_CaMol;
RESdl_KCl = RESdl_KMol;
RESdl_NaBr = RESdl_BrMol;
}

// [[Rcpp::export]]
void CalcSaltDensity (NumericMatrix ChemConcentrations, NumericMatrix ChemParams){
//This routine estimates the density of the salt solution. Full disassociation is assumed.
//first estimate 1L of liquid density. (Sum of salts + water).
//Need a conditional in case of 0 volume reservoir.
double InterpolationMainStep = 0;
double InterpolationPartStep = 0;

	if (RESsl_NaMol + 2*RESsl_MgMol + 2*RESsl_CaMol + RESsl_KMol != 0){
		double RESsl_DensitySalts = RESsl_NaCl * NaCl_MolMass + RESsl_MgCl2*MgCl2_MolMass + RESsl_CaCl2*CaCl2_MolMass + RESsl_KCl*KCl_MolMass + RESsl_NaBr*NaBr_MolMass;
		RESsl_Density = (1000 + RESsl_DensitySalts)/1000;
		// and estimate initial salinity estimates
		RESsl_Salinity = RESsl_DensitySalts/RESsl_Density;
		//Compare density to density in parameters.
		InterpolationMainStep = floor(RESsl_Salinity/25);
		InterpolationPartStep = RESsl_Salinity/25 - InterpolationMainStep;
		RESsl_Density = ChemParams(InterpolationMainStep,1) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,1) - ChemParams(InterpolationMainStep,1));
		//Now revise salinity based on updated density.
		RESsl_Salinity = RESsl_DensitySalts*1000/RESsl_Density;
	} else {
	RESsl_Density = 0;
	RESsl_Salinity = 0;
	}

	if (RESdl_NaMol + 2*RESdl_MgMol + 2*RESdl_CaMol + RESdl_KMol != 0){
		double RESdl_DensitySalts = RESdl_NaCl * NaCl_MolMass + RESdl_MgCl2*MgCl2_MolMass + RESdl_CaCl2*CaCl2_MolMass + RESdl_KCl*KCl_MolMass + RESdl_NaBr*NaBr_MolMass;
		RESdl_Density = (1000 + RESdl_DensitySalts)/1000;
		RESdl_Salinity = RESdl_DensitySalts/RESdl_Density;
		InterpolationMainStep = floor(RESdl_Salinity/25);
		InterpolationPartStep = RESdl_Salinity/25 - InterpolationMainStep;
		RESdl_Density = ChemParams(InterpolationMainStep,1) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,1) - ChemParams(InterpolationMainStep,1));
		RESdl_Salinity = RESdl_DensitySalts*1000/RESdl_Density;
	} else {
	RESdl_Density = 0;
	RESdl_Salinity = 0;
	}
}

// [[Rcpp::export]]
void CalcWaterActivity (NumericMatrix ChemConcentrations, NumericMatrix ChemParams){
double InterpolationMainStep = 0;
double InterpolationPartStep = 0;
double NaClAW = 0;
double MgCl2AW = 0;
double CaCl2AW = 0;
double KClAW = 0;
double NaBrAW = 0;
	if (RESsl_NaMol + 2*RESsl_MgMol + 2*RESsl_CaMol + RESsl_KMol != 0){
		InterpolationMainStep = floor(RESsl_Salinity/25);
		InterpolationPartStep = RESsl_Salinity/25 - InterpolationMainStep;
		double ResSL_TotalMolality = RESsl_NaCl+RESsl_MgCl2+RESsl_CaCl2+RESsl_KCl+RESsl_NaBr;
		NaClAW = ChemParams(InterpolationMainStep, 2) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,2) - ChemParams(InterpolationMainStep,2));
		MgCl2AW = ChemParams(InterpolationMainStep, 3) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,3) - ChemParams(InterpolationMainStep,3));
		CaCl2AW = ChemParams(InterpolationMainStep, 4) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,4) - ChemParams(InterpolationMainStep,4));
		KClAW = ChemParams(InterpolationMainStep, 5) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,5) - ChemParams(InterpolationMainStep,5));
		NaBrAW = ChemParams(InterpolationMainStep, 6) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,6) - ChemParams(InterpolationMainStep,6));
		RESsl_ActWater = (RESsl_NaCl * NaClAW + RESsl_MgCl2 * MgCl2AW + RESsl_CaCl2 * CaCl2AW + RESsl_KCl * KClAW + RESsl_NaBr * NaBrAW)/ResSL_TotalMolality;
	} else {
	RESsl_ActWater = 1;
	}
	//and for the RESdl
	if (RESdl_NaMol + 2*RESdl_MgMol + 2*RESdl_CaMol + RESdl_KMol != 0){
		InterpolationMainStep = floor(RESdl_Salinity/25);
		InterpolationPartStep = RESdl_Salinity/25 - InterpolationMainStep;
		double ResDL_TotalMolality = RESdl_NaCl+RESdl_MgCl2+RESdl_CaCl2+RESdl_KCl+RESdl_NaBr;
		NaClAW = ChemParams(InterpolationMainStep, 2) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,2) - ChemParams(InterpolationMainStep,2));
		MgCl2AW = ChemParams(InterpolationMainStep, 3) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,3) - ChemParams(InterpolationMainStep,3));
		CaCl2AW = ChemParams(InterpolationMainStep, 4) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,4) - ChemParams(InterpolationMainStep,4));
		KClAW = ChemParams(InterpolationMainStep, 5) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,5) - ChemParams(InterpolationMainStep,5));
		NaBrAW = ChemParams(InterpolationMainStep, 6) + InterpolationPartStep*(ChemParams(InterpolationMainStep+1,6) - ChemParams(InterpolationMainStep,6));
		RESdl_ActWater = (RESdl_NaCl * NaClAW + RESdl_MgCl2 * MgCl2AW + RESdl_CaCl2 * CaCl2AW + RESdl_KCl * KClAW + RESdl_NaBr * NaBrAW)/ResDL_TotalMolality;
	} else {
	RESdl_ActWater = 1;
	}
}

// [[Rcpp::export]]
void SalinityCalcs (double Time, NumericMatrix ChemConcentrations, NumericMatrix Chemistry, NumericMatrix ChemParams){
	CalcMolarities(Time, Chemistry, ChemConcentrations);
	EstimateSaltMolarities(ChemConcentrations);
	CalcSalts(ChemConcentrations);
	CalcSaltDensity (ChemConcentrations, ChemParams);
	CalcWaterActivity (ChemConcentrations, ChemParams);
}
