#include <Rcpp.h>
using namespace Rcpp;
#include "Hypsographic_Calcs.hpp"
#include "Hydrology_functions.hpp"

// [[Rcpp::export]]
List Headcalc(NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakebedgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double Timestep, double Runtime, double ETexdepth, NumericMatrix Lakedepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations, NumericMatrix GWNumbers) {


int nrow = Topogrid.nrow();
int ncol = Topogrid.ncol();
int time = 0;
double Hold = 0;
double maxdiff = Accuracy+1;
int itercount = 0;
double outseepage = 0;
double inflows = 0;
double outflows = 0; 
double surfaceflows = 0; // Flows into the catchments
double massinit = 0;
double massfinal = 0;
double Extflows = 0;
double RechargeFlows = 0;
double ETflux = 0;
double Fixedheadflows = 0;
double ETfluxall = 0;
//NumericMatrix GWNumbers (9,5); 
NumericMatrix Qgrid (nrow,ncol); //remove in future
NumericMatrix Sgrid (nrow,ncol);
NumericMatrix RHSbase (nrow,ncol);
NumericMatrix HCOFbase (nrow,ncol);
NumericMatrix RHS (nrow,ncol);
NumericMatrix HCOF (nrow,ncol);
NumericMatrix Aqthick (nrow,ncol);
NumericMatrix GWnewheadgrid = clone(GWheadgrid);
NumericMatrix Flowsgrid (nrow,ncol);
NumericMatrix CW (nrow,ncol);
NumericMatrix CE (nrow,ncol);
NumericMatrix CS (nrow,ncol);
NumericMatrix CN (nrow,ncol);
NumericMatrix Laplace (nrow,ncol);
NumericMatrix Lakecoords (nrow*ncol,4); //x,y,depth,flow
NumericMatrix Landcoords (nrow*ncol,2); 
NumericMatrix Allcoords (nrow*ncol,2); 
// debug stuff below
NumericMatrix diff (nrow,ncol);
bool Overlandflows = TRUE;

// ************************ New method based on spreadsheet model **********************
// *****************************************************************************
// ================= Before iteration loop (Q, S, HCOF, RHS) ===================
// *****************************************************************************

// ******************** Add lake numbers to the GWNumbers matrix ************************
// And some column names.
//colnames (GWNumbers) = CharacterVector::create("Lake ID", "Seepage_in", "Seepage_out", "Surface Flows", "Net_flows");
for (int i = 0; i < GWNumbers.nrow(); i++){
GWNumbers(i,0) = i+1;
GWNumbers(i,1) = 0;
GWNumbers(i,2) = 0;
GWNumbers(i,3) = 0;
GWNumbers(i,4) = 0;

}

// ************************ Define all lake cells **********************
//This is no longer needed as the lake now uses Cauchy boundaries. Therefore only a few cells will use fixed head boundaries, and they can be caught with an if statement. 
//However, this is still useful to cut down on HCOF and RHS calcs for the lake, and to define the areas of "confined" aquifers.
//Remaining cells can be added to a 2 column matrix, which will then be used for the for loops in the solution calculation. This will cut down the number of cells requiring calculation and enforce the boundary conditions.
// First two columns (0 and 1) are grid coordinates, Third column (2) is lake depth. Fourth column (3) is flux to and from the lake cell via the Cauchy boundary for each timestep. 
std::fill(Lakecoords.begin(), Lakecoords.end(), NumericVector::get_na());
//Then start a counter
int Lakecounter = 0;
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
	if (GWboundarygrid(i,j) == 0 && Classgrid(i,j) < 10){
	Lakecoords (Lakecounter,0) = i;
	Lakecoords (Lakecounter,1) = j;
	Lakecoords (Lakecounter,2) = Lakedepth(Classgrid(i,j)-1,3);
	//Needs the -1 to account for C++ starting the Lakedepth matrix at 0.
	Lakecounter++;
	
	//and we use counter further down (in the solution) to stop running into the NA cells.
	}
	}
}
//Cells on the land (not fixed head). 
std::fill(Landcoords.begin(), Landcoords.end(), NumericVector::get_na());
//Then start a counter
int Landcounter = 0;
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
	if (GWboundarygrid(i,j) == 0 && Classgrid(i,j) >= 10){
	Landcoords (Landcounter,0) = i;
	Landcoords (Landcounter,1) = j;
	Landcounter++;
	//and we use counter further down (in the solution) to stop running into the NA cells.
	}
	}
}
//Cells needed for solution
int Allcounter = 0;
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
	if (GWboundarygrid(i,j) == 0){
	Allcoords (Allcounter,0) = i;
	Allcoords (Allcounter,1) = j;
	Allcounter++;
	//and we use counter further down (in the solution) to stop running into the NA cells.
	}
	}
}

// The above divides up the cells into lake cells and land cells, and a remained of fixed head boundary cells. This makes it easy to apply specific calculations to zones. Lake cells calculate for lake seepage, whereas land cells calculate for ET, recharge and wells.


// ************************ Define starting heads **********************

//Define starting head heights, based on existing
//Strange behaviour here. This assigns the wrong way - making GWheadgrid the same as GWnewheadgrid. Maybe something to do with proxy models.
// see http://stackoverflow.com/questions/30344257/rcpp-inconsistent-behavior-with-proxy-model
//http://stackoverflow.com/questions/11300048/rcpp-pass-by-reference-vs-by-value/11300707#11300707
//GWnewheadgrid = clone(GWheadgrid);// this version preserves headgrid
// This version updates headgrid
// Moved up to declaration area.

// ************************ Define starting aquifer thickness **********************

	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
		Aqthick(i,j) = (GWheadgrid(i,j) - GWbasegrid(i,j));
		if (Aqthick(i,j) < 0){ //just in case head height drops below the base of the aquifer.
		Aqthick(i,j) = 0.1;
	}
	}
	}


// ************************ Define combined fluxes **********************
//Only applied to land cells. Recharge and Evap of lake done in CHIMBLE.


for (int row = 0; row < Landcounter; row++){
			int i = Landcoords(row,0);
			int j = Landcoords(row,1);
			Qgrid (i,j) = GWRechargegrid (i,j)*Cellxdim*Cellydim + GWPumpinggrid(i,j);
			RechargeFlows += GWRechargegrid (i,j)*Cellxdim*Cellydim * Timestep;
			Extflows += GWPumpinggrid(i,j) * Timestep; //double check volume per timestep vs daily rate.
	}
	
// And apply pumping to lake cells too. This is to allow for seepage through the vents below the lake to deeper aquifers, eg: perched lakes.
for (int row = 0; row < Lakecounter; row++){
	int i = Lakecoords(row,0);
	int j = Lakecoords(row,1);
		Qgrid (i,j) = GWPumpinggrid(i,j);
		Extflows += GWPumpinggrid(i,j) * Timestep;
	}
//return Qgrid;

// ************************ Define cell storage **********************
// Single matrix
// Requires GWSgrid, Cellxdim, Cellydim
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
			Sgrid (i,j) = GWSgrid (i,j) * Cellxdim * Cellydim;
	}
	}
//return Sgrid;


// ************************ Define initial mass balance *****************************
//Calculate Mass Balance
//make sure there is no mass above the land surface.
// While it may appear that we should treat this as two separate aquifers, instead we are treating it as one single unconfined aquifer with a big low permeability blob in the middle that can remove or add water
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
			massinit += (GWheadgrid(i,j) - GWbasegrid(i,j))*Sgrid(i,j);
	}
	}
	
	// ************************ Define HCOF **********************
// Requires cell storage matrix, Timestep
// Should be updated with each iteration with head dependent nodes.
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
		HCOFbase(i,j) = -GWSgrid (i,j) * Cellxdim * Cellydim/Timestep ;
	}
	}
//return HCOF;

// **************************************************************************************
//================== Inside main loop (Conductance, Laplace, Solution) ==================
// **************************************************************************************
	while (time < Runtime){
	//reset variables
	Hold = 0;
	maxdiff = Accuracy+1;
	itercount = 0;
	

	
	// ********************************** Define RHS **********************************
	// Requires Combined flux matrix, Cell Storage matrix, GWheadgrid, Timestep
	// Should be updated with each iteration with head dependent nodes.
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
		RHSbase (i,j) = -Qgrid(i,j) - (GWSgrid (i,j) * Cellxdim * Cellydim * GWheadgrid(i,j))/Timestep ;
	//below option uses matrixes from above code.
	}
	}
	//return RHSbase;
	

	// *********************************************************************************
	//============= Inside iteration loop (Conductance, Laplace, Solution) =============
	// *********************************************************************************
		while (maxdiff > Accuracy && itercount < Maxiterations){
		//Reset variables
			maxdiff = 0;


	// ************************ Define aquifer thickness **********************
	// needed for each iteration
	//This treats land areas as unconfined, but the cells beneath the lake as confined, as per A2016 V5
	for (int row = 0; row < Landcounter; row++){
	int i = Landcoords(row,0);
	int j = Landcoords(row,1);
	Aqthick(i,j) = GWnewheadgrid(i,j) - GWbasegrid(i,j);
	if (Aqthick(i,j) < 0){ //just in case head height drops below the base of the aquifer.
		Aqthick(i,j) = 0.1;
	}
	}

	// Split this into two runs, as it saves us using an if/or statement for all cells
	for (int row = 0; row < Lakecounter; row++){
	int i = Lakecoords(row,0);
	int j = Lakecoords(row,1);
	Aqthick(i,j) = (Topogrid(i,j) - GWlakebedgrid(i,j)) - GWbasegrid(i,j);
	}

	// ************************ Define Lake cells HCOF and RHS **********************
	// here we can use the counter routine as only lake cells require the calculation
	//need to modify this to account for multiple lakes.
	for (int row = 0; row < Lakecounter; row++){
	int i = Lakecoords(row,0);
	int j = Lakecoords(row,1);
		if (GWnewheadgrid(i,j) > Topogrid(i,j) - GWlakebedgrid(i,j)){
		HCOF(i,j) = HCOFbase(i,j) - Lakecondgrid(i,j);
		RHS(i,j) = RHSbase(i,j) - Lakecondgrid(i,j) * Lakecoords(row,2);
		Lakecoords(row,3) = Lakecondgrid(i,j)*(Lakecoords(row,2) - GWnewheadgrid(i,j));
		} else {
		RHS(i,j) = RHSbase(i,j) - Lakecondgrid(i,j) * (Lakecoords(row,2) - (Topogrid(i,j) - GWlakebedgrid(i,j)));
		HCOF(i,j) = HCOFbase(i,j);
		//And update the flux 
		Lakecoords(row,3) = Lakecondgrid(i,j)*(Lakecoords(row,2) - (Topogrid(i,j) - GWlakebedgrid(i,j)));
		}
	}

	// ************************** Define ET HCOF and RHS ****************************	
	ETflux = 0;
	//Rcpp::Rcout << " ETGrid Value: " << GWETgrid(20,20) << std::endl;
	//Rcpp::Rcout << "Confirmed ETGrid Value: " << GWETgrid(20,20) << std::endl;
	for (int row = 0; row < Landcounter; row++){
	int i = Landcoords(row,0);
	int j = Landcoords(row,1);
		if (GWnewheadgrid(i,j) > Topogrid (i,j)){
			RHS(i,j) = RHSbase(i,j) + GWETgrid(i,j)*Cellxdim*Cellydim;
			HCOF(i,j) = HCOFbase(i,j);
			//And add it to the scalar value of ETflux.
			ETflux += GWETgrid(i,j)*Cellxdim*Cellydim;
		} else if (GWnewheadgrid(i,j) < Topogrid (i,j) - ETexdepth) {
			HCOF(i,j) = HCOFbase(i,j);
			RHS(i,j) = RHSbase(i,j);
		} else {
			HCOF(i,j) = HCOFbase(i,j) - (GWETgrid(i,j)*Cellxdim*Cellydim/ETexdepth);
			RHS(i,j) = RHSbase(i,j) - (-GWETgrid(i,j)*Cellxdim*Cellydim + (GWETgrid(i,j)*Cellxdim*Cellydim *(Topogrid(i,j)/ETexdepth)));
			ETflux += GWETgrid(i,j)*Cellxdim*Cellydim* (GWnewheadgrid(i,j)-(Topogrid (i,j) - ETexdepth))/ETexdepth;
		}
	}

		// ************************ Define conductance **********************
		// 4 matrixes for conductance in each direction
		// Requires GWKgrid, Cellxdim, Cellydim, GWnewheadgrid
		// Use edge of model to define 0 conductance and a no-flow.
		// Should be updated with each iteration.
		//These have been updated to reflect full harmonic mean as per A2016 V5
		//Conductance West
		//Future thing to test - Combine these, and the Laplace core, into a bunch of values calculated for each single point during the solution. That will cut down on the for loops, and the additional matrixes. May give a small speed boost.
		

			for (int j = 1; j < ncol; j++){
			for (int i = 0; i < nrow; i++){
			CW(i,j) = (2*GWKgrid(i,j-1) * Aqthick(i,j-1) * GWKgrid(i,j) * Aqthick(i,j))/(GWKgrid(i,j-1) * Aqthick(i,j-1) + GWKgrid(i,j) * Aqthick(i,j))*Cellydim/Cellxdim;
			}
			}

			for (int j = 0; j < ncol-1; j++){
			for (int i = 0; i < nrow; i++){
			CE(i,j) = (2*GWKgrid(i,j+1) * Aqthick(i,j+1) * GWKgrid(i,j) * Aqthick(i,j))/(GWKgrid(i,j+1) * Aqthick(i,j+1) + GWKgrid(i,j) * Aqthick(i,j))*Cellydim/Cellxdim;
			}
			}

			for (int j = 0; j < ncol; j++){
			for (int i = 0; i < nrow-1; i++){
			CS(i,j) = (2*GWKgrid(i+1,j) * Aqthick(i+1,j) * GWKgrid(i,j) * Aqthick(i,j))/(GWKgrid(i+1,j) * Aqthick(i+1,j) + GWKgrid(i,j) * Aqthick(i,j))*Cellxdim/Cellydim;
			}
			}

			for (int j = 0; j < ncol; j++){
			for (int i = 1; i < nrow; i++){
			CN(i,j) = (2*GWKgrid(i-1,j) * Aqthick(i-1,j) * GWKgrid(i,j) * Aqthick(i,j))/(GWKgrid(i-1,j) * Aqthick(i-1,j) + GWKgrid(i,j) * Aqthick(i,j))*Cellxdim/Cellydim;
			}
			}

		// ************************ Define Laplace Core **************************
		//Requires Conductances and GWnewheadgrid
		// Should be updated with each iteration
			for (int j = 1; j < ncol-1; j++){
			for (int i = 1; i < nrow-1; i++){
			Laplace(i,j) = CW(i,j)*GWnewheadgrid(i,j-1) + CE(i,j)*GWnewheadgrid(i,j+1) + CS(i,j)*GWnewheadgrid(i+1,j)+CN(i,j)*GWnewheadgrid(i-1,j);
			}
			}

		//The use of conductance results in 0's along the edges of the model. This means that even referencing outside the matrix is possible as it multiplies by 0. Therefore this next 4 loops aren't required and we can just parse the entire matrix at once. :D
		//Cancel the above. Using data from just outside the matrix occasionally results in incorrect data being used, if the matrix is put into memory next to another object.
		
			//edges
			for (int j = 1;j < ncol-1; j++){
			Laplace(0,j) = CW(0,j)*GWnewheadgrid(0,j-1) + CE(0,j)*GWnewheadgrid(0,j+1) + CS(0,j)*GWnewheadgrid(1,j);
			}

			for (int j = 1;j < ncol-1; j++){
			Laplace(nrow-1,j) = CW(nrow-1,j)*GWnewheadgrid(nrow-1,j-1) + CE(nrow-1,j)*GWnewheadgrid(nrow-1,j+1) + CN(nrow-1,j)*GWnewheadgrid(nrow-2,j);
			}

			for (int i = 1;i < nrow-1; i++){
			Laplace(i,0) = CE(i,0)*GWnewheadgrid(i,1) + CS(i,0)*GWnewheadgrid(i+1,0)+CN(i,0)*GWnewheadgrid(i-1,0);
			}

			for (int i = 1;i < nrow-1; i++){
			Laplace(i,ncol-1) = CW(i,ncol-1)*GWnewheadgrid(i,ncol-2) + CS(i,ncol-1)*GWnewheadgrid(i+1,ncol-1)+CN(i,ncol-1)*GWnewheadgrid(i-1,ncol-1);
			}
		
			//corners
			Laplace(0,0) = CE(0,0)*GWnewheadgrid(0,1) + CS(0,0)*GWnewheadgrid(1,0);
			Laplace(nrow-1,0) = CE(nrow-1,0)*GWnewheadgrid(nrow-1,1) + CN(nrow-1,0)*GWnewheadgrid(nrow-2,0);
			Laplace(0,ncol-1) = CS(0,ncol-1)*GWnewheadgrid(1,ncol-1) + CW(0,ncol-1)*GWnewheadgrid(0,ncol-2);
			Laplace(nrow-1,ncol-1) = CW(nrow-1,ncol-1)*GWnewheadgrid(nrow-1,ncol-2) + CN(nrow-1,ncol-1)*GWnewheadgrid(nrow-2,ncol-1);

			// ************************ Define Solution *****************************
			// Requires Laplace Core, RHS, HCOF, Conductances	
			for (int row = 0; row < Allcounter; row++){
			int i = Allcoords(row,0);
			int j = Allcoords(row,1);
			Hold = GWnewheadgrid(i,j);
			//grab previous value for comparison.
			GWnewheadgrid(i,j) = (RHS(i,j) - Laplace(i,j)) / (-CW(i,j)-CE(i,j)-CS(i,j)-CN(i,j)+HCOF(i,j));
				if (fabs(GWnewheadgrid(i,j) - Hold) > maxdiff){
				maxdiff = fabs(GWnewheadgrid(i,j) - Hold);
				}
			}

		itercount++;
		}
		
		for (int j = 1; j < ncol-1; j++){//added to prevent runaway condition when pumping draws a cell down to <0 thickness. 25/09/19
		for (int i = 1; i < nrow-1; i++){
			if (GWnewheadgrid(i,j) <= GWbasegrid(i,j)){ 
			GWnewheadgrid(i,j) = GWbasegrid(i,j) + 0.1;
			}
		}
		}
	// ******************************************************************
	//======================= End iteration loop ========================
	// ******************************************************************
	// ************************ Update Flows *****************************
	// May not be needed to crunch the entire matrix. We only need the flows in the lake area. Boundary = 0
	// Having said that, full thing may be useful for mass balance.
	// As flows are calculated based on the head heights at the end of the timestep, they represent the average flow per time unit over that timestep. Therefore flows have to be multiplied by the timestep.

			for (int j = 1; j < ncol-1; j++){
			for (int i = 1; i < nrow-1; i++){
			Flowsgrid(i,j) = Timestep*(CW(i,j)*(GWnewheadgrid(i,j-1)-GWnewheadgrid(i,j)) + CE(i,j)*(GWnewheadgrid(i,j+1)-GWnewheadgrid(i,j)) + CN(i,j)*(GWnewheadgrid(i-1,j)-GWnewheadgrid(i,j)) + CS(i,j)*(GWnewheadgrid(i+1,j)-GWnewheadgrid(i,j))) + Flowsgrid(i,j);

			}
			}
		
			//edges
			for (int j = 1;j < ncol-1; j++){
			Flowsgrid(0,j) = Timestep*(CW(0,j)*(GWnewheadgrid(0,j-1)-GWnewheadgrid(0,j)) + CE(0,j)*(GWnewheadgrid(0,j+1)-GWnewheadgrid(0,j)) + CS(0,j)*(GWnewheadgrid(1,j)-GWnewheadgrid(0,j))) + Flowsgrid(0,j);
			}

			for (int j = 1;j < ncol-1; j++){
			Flowsgrid(nrow-1,j) = Timestep*(CW(nrow-1,j) *(GWnewheadgrid(nrow-1,j-1)-GWnewheadgrid(nrow-1,j)) + CE(nrow-1,j)*(GWnewheadgrid(nrow-1,j+1)-GWnewheadgrid(nrow-1,j)) + CN(nrow-1,j)*(GWnewheadgrid(nrow-2,j)-GWnewheadgrid(nrow-1,j))) + Flowsgrid(nrow-1,j);
			}

			for (int i = 1;i < nrow-1; i++){
			Flowsgrid(i,0) = Timestep*(CE(i,0)*(GWnewheadgrid(i,1)-GWnewheadgrid(i,0)) + CN(i,0)*(GWnewheadgrid(i-1,0)-GWnewheadgrid(i,0)) + CS(i,0)*(GWnewheadgrid(i+1,0)-GWnewheadgrid(i,0))) + Flowsgrid(i,0);
			}

			for (int i = 1;i < nrow-1; i++){
			Flowsgrid(i,ncol-1) = Timestep*(CW(i,ncol-1)*(GWnewheadgrid(i,ncol-2)-GWnewheadgrid(i,ncol-1)) + CN(i,ncol-1)*(GWnewheadgrid(i-1,ncol-1)-GWnewheadgrid(i,ncol-1)) + CS(i,ncol-1)*(GWnewheadgrid(i+1,ncol-1)-GWnewheadgrid(i,ncol-1))) + Flowsgrid(i,ncol-1);
			}
		
			//corners
			Flowsgrid(0,0) = Timestep*(CE(0,0)*(GWnewheadgrid(0,1)-GWnewheadgrid(0,0)) + CS(0,0)*(GWnewheadgrid(1,0)-GWnewheadgrid(0,0))) + Flowsgrid(0,0);
			
			Flowsgrid(nrow-1,0) = Timestep*(CE(nrow-1,0)*(GWnewheadgrid(nrow-1,1)-GWnewheadgrid(nrow-1,0)) + CN(nrow-1,0)*(GWnewheadgrid(nrow-2,0)-GWnewheadgrid(nrow-1,0))) + Flowsgrid(nrow-1,0);
			
			Flowsgrid(0,ncol-1) = Timestep*(CW(0,ncol-1)*(GWnewheadgrid(0,ncol-2)-GWnewheadgrid(0,ncol-1)) + CS(0,ncol-1)*(GWnewheadgrid(1,ncol-1)-GWnewheadgrid(0,ncol-1))) + Flowsgrid(0,ncol-1);
			
			Flowsgrid(nrow-1,ncol-1) = Timestep*(CW(nrow-1,ncol-1)*(GWnewheadgrid(nrow-1,ncol-2)-GWnewheadgrid(nrow-1,ncol-1)) + CN(nrow-1,ncol-1)*(GWnewheadgrid(nrow-2,ncol-1)-GWnewheadgrid(nrow-1,ncol-1))) + Flowsgrid(nrow-1,ncol-1);

				//Rcpp::Rcout << "Lake number: " << GWNumbers(0,0) << std::endl;				
				//Rcpp::Rcout << "Check on inflow: " << GWNumbers(0,1) << std::endl;		
				//Rcpp::Rcout << "Check on outseepage: " << GWNumbers(0,2) << std::endl;

// ************************ Quantify external flows balance *****************************
//		for (int row = 0; row < Landcounter; row++){
//		int i = Landcoords(row,0);
//		int j = Landcoords(row,1);
//		extflows += Qgrid(i,j)*Timestep;
//		}
//Now done during the formation of Qgrid so as to separate recharge from Ext Flows.		
		
// ************************ Quantify lake in and outseepage *****************************
// This is cumulative over the runtime of the model. \
// ***Check this. May be tallying up all lake flows as Bullen Merri flows. Cancel that. This works as intended, with a tally for each lake.
// With updated model it needs to separate out inflows and outflows. 
for (int row = 0; row < Lakecounter; row++){
	int i = Lakecoords(row,0);
	int j = Lakecoords(row,1);
	int lakenum = (Classgrid(i,j));
	if (Lakecoords(row,3) < 0){
	GWNumbers(lakenum-1,1) += Lakecoords(row,3)*Timestep;
	} else {
	GWNumbers(lakenum-1,2) += Lakecoords(row,3)*Timestep;
	}
	}

// ******************************* Quantify ET loss**************************************
// This is cumulative over the runtime of the model. 
//Not calculated for lake or fixed head cells
		ETfluxall += ETflux*Timestep;
		//Rcpp::Rcout << "ETflux " << ETflux << std::endl;
		//Rcpp::Rcout << "ETfluxall " << ETfluxall << std::endl;



// **************************** Define Overland Flows *********************************
if (Overlandflows){
		//Rcpp::Rcout << "Overland flows active." << std::endl;
// If head heights are above the topography, then they are treated as overland flows and are removed from the catchment either into the lakes or lost through ET/drainage. In future we could look at better solutions (seepage faces, flows to downhill cells), but for now we're using a simple approach.
//This section could be merged with the above section (lake in and outseepage), however, as the above section only applies to the lake cells it may be faster to keep it separate.
for (int row = 0; row < Landcounter; row++){
			int i = Landcoords(row,0);
			int j = Landcoords(row,1);
	
			if (Classgrid (i,j) < 20 && GWnewheadgrid (i,j) > Topogrid(i,j)){
			GWNumbers((Classgrid (i,j) - 11),3) += (GWnewheadgrid(i,j) - Topogrid(i,j))*Sgrid(i,j);
			GWnewheadgrid (i,j) = Topogrid (i,j);
			} 
			if (Classgrid(i,j) == 20 && GWnewheadgrid (i,j) > Topogrid(i,j)){
			outflows += (GWnewheadgrid(i,j) - Topogrid(i,j))*Sgrid(i,j);
			GWnewheadgrid (i,j) = Topogrid (i,j);
			}	
		}
} else {
	//Rcpp::Rcout << "Overland flows disabled" << std::endl;
}



// ***************************** Finalise Main Loop *************************************
		//Rcpp::Rcout << "Classgrid@100,100: " << Classgrid(100,100) << std::endl;
		//Rcpp::Rcout << "Qgrid@100,100: " << Qgrid(100,100) << std::endl;
		//Rcpp::Rcout << "Sgrid@100,100: " << Sgrid(100,100) << std::endl;
		//Rcpp::Rcout << "Laplace@100,100: " << Laplace(100,100) << std::endl;
		//Rcpp::Rcout << "Topogrid@100,100: " << Topogrid(100,100) << std::endl;
		//Rcpp::Rcout << "Head@100,100: " << GWheadgrid(100,100) << std::endl;
		//Rcpp::Rcout << "NewHead@100,100: " << GWnewheadgrid(100,100) << std::endl;
		//Update time
		//Rcpp::Rcout << "Runtime " << time << std::endl;
		time = time+Timestep;

		//Update headheights for new timestep

		for (int j = 0; j < ncol; j++){
		for (int i = 0; i < nrow; i++){
		GWheadgrid(i,j) = GWnewheadgrid(i,j);
		}
		}


}
// **************************************************************************************
//================== End main loop (Conductance, Laplace, Solution) ==================
// **************************************************************************************

//	return GWheadgrid;

// ************************ Flows to/from fixed head cells *****************************

	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
			if (GWboundarygrid(i,j) == 1){
			Fixedheadflows += Flowsgrid(i,j);
			}
	}
	}
		
// ************************ Define final mass balance *****************************
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
			massfinal += (GWheadgrid(i,j) - GWbasegrid(i,j))*Sgrid(i,j);
	}
	}
// ************************ Calculate Mass Balance Error *****************************
// Calculate all flows into and out of all lakes for mass balance
for (int i = 0; i < GWNumbers.nrow(); i++){
	inflows += GWNumbers(i,1);
	outseepage += GWNumbers (i,2);
	surfaceflows += GWNumbers (i,3);
	GWNumbers (i,4) = -GWNumbers(i,1) - GWNumbers(i,2) + GWNumbers (i,3);
}

//Rcpp::Rcout << "Fixed Head, outseepage, inflows " << Fixedheadflows << ", " << outseepage << ", " << inflows << std::endl;


// ***Modify this to account for all lake and catchment flows
// See page 7 - http://www.sspa.com/sites/default/files/images/stories/documents/MT3D-notes_mass-budget.pdf
	//double masserror = 100*(massinit - massfinal + Fixedheadflows + Extflows + RechargeFlows + outseepage + inflows + surfaceflows - ETfluxall - outflows)/((Fixedheadflows + Extflows + RechargeFlows + outseepage + inflows + surfaceflows - ETfluxall - outflows));
	
	double masserror = 100*((massinit + Extflows +  RechargeFlows + outseepage ) - (massfinal + Fixedheadflows + inflows + surfaceflows - ETfluxall - outflows))/(0.5*((massinit + Extflows +  RechargeFlows + outseepage) + (massfinal + Fixedheadflows + inflows + surfaceflows - ETfluxall - outflows)));
	
	//double masserror = 100*(massinit + Extflows + outseepage + inflows) - (massfinal +Fixedheadflows + surfaceflows + ETfluxall + outflows)/(0.5*(massinit + Extflows + outseepage + inflows)+(massfinal +Fixedheadflows + surfaceflows + ETfluxall + outflows));
	//double masserror = 100*( Extflows - Fixedheadflows + outseepage + inflows - surfaceflows - ETfluxall - outflows)/(Fixedheadflows + Extflows + outseepage + inflows - surfaceflows - ETfluxall - outflows); Modflow variant - Chapter A1 - seems to assume steady state?
	

double rec=GWRechargegrid(1,1);
double K=GWKgrid(1,1);

NumericVector Stats = NumericVector::create(Named("Massinit")=massinit, Named("Massfinal")=massfinal, Named("Mass error")=masserror, Named("F_Head flux")=Fixedheadflows, Named("Recharge") = RechargeFlows, Named("Extflows")=Extflows, Named("Lake seepage(in)")=inflows, Named("Lake seepage(out)") = outseepage, Named("Surface flows")=surfaceflows, Named("Non catchment outflows")=outflows, Named("ET")=ETfluxall, Named("Timestep")=Timestep, Named("Runtime")=Runtime, Named("Recharge")=rec, Named("Conductivity")=K);

// ***Modify this to account for all and catchment flows
List GW = List::create(
Named("Flows") = GWNumbers,
Named("Stats") = Stats);
return GW;     
} // End c function


// [[Rcpp::export]]
void func_ClassGrid(long Time, NumericMatrix Topogrid, NumericMatrix LakeDepths, NumericMatrix GWcatchmentgrid, NumericMatrix Classgrid, NumericMatrix Lakevolumes){
int nrow = Topogrid.nrow();
int ncol = Topogrid.ncol();
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
//process through main grid
		if (GWcatchmentgrid(i,j) != 0){ //if a catchment
			for (int k = 0; k < LakeDepths.nrow(); k++){
				if (GWcatchmentgrid(i,j) == k+1){
					if (Topogrid(i,j) - LakeDepths(k,3) <= 0){
						Classgrid(i,j) = k+1;
					} else {
						Classgrid(i,j) = k+11;
				}
			}
			}
		} else {
			Classgrid(i,j) = 20;
		}
	}
	}
}


/* //old method - updated to use LakeDepths matrix Feb2019.
double Lake1Depth = func_Lakedepth(Lake(Time,5), Lakevolumes);
double Lake2Depth = func_Lakedepth(Lake(Time,5), Lakevolumes) - 40; //hack for gnotuk
//Rcpp::Rcout << "Lake1depth " << Lake1Depth << std::endl;
//Rcpp::Rcout << "Lake2depth " << Lake2Depth << std::endl;
for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
		if (GWcatchmentgrid(i,j) == 1){
			if (Topogrid(i,j) - Lake1Depth <= 0){
				Classgrid(i,j) = 1 ;
				//Rcpp::Rcout << "Classgrid set to 1 at" << i << "," << j << std::endl;
			} else {
				Classgrid(i,j) = 11;
				//Rcpp::Rcout << "Classgrid set to 11 at" << i << "," << j << std::endl;
			}	
		} else if (GWcatchmentgrid(i,j) == 2){
			if (Topogrid(i,j) - Lake2Depth <= 0){
				Classgrid(i,j) = 2;
				//Rcpp::Rcout << "Classgrid set to 2 at" << i << "," << j << std::endl;
			} else {
				Classgrid(i,j) = 12;
				//Rcpp::Rcout << "Classgrid set to 12 at" << i << "," << j << std::endl;
			}	
		} else {
		Classgrid(i,j) = 20;
		//Rcpp::Rcout << "Classgrid set to 20 at" << i << "," << j << std::endl;
		}
	}
}
}
*/

// [[Rcpp::export]]
void func_LakeConductance(NumericMatrix Lakecondgrid, int Cellxdim, int Cellydim, NumericMatrix GWlakesedimentgrid, NumericMatrix GWlakekgrid){
int nrow = GWlakesedimentgrid.nrow();
int ncol = GWlakesedimentgrid.ncol();
int area = Cellxdim * Cellydim;
	for (int j = 0; j < ncol; j++){
		for (int i = 0; i < nrow; i++){
		if (GWlakesedimentgrid(i,j) != 0){
		Lakecondgrid(i,j) = GWlakekgrid(i,j) * area / GWlakesedimentgrid(i,j);
		} else {
		Lakecondgrid(i,j) = 1000000; //set very high conductivity so that sediments are essentially ignored. 1000000 is equivalent to 1m/d over .01cm sediment.
		}
	}
}
}

// [[Rcpp::export]]
void func_RechargeGrid (long Time, int Cellxdim, int Cellydim, NumericMatrix GWRechargegrid, NumericMatrix Classgrid, NumericMatrix Flux, double Timestep){
//This has to run after influx is calculated as it uses calculated flux values from the mass balance. Either that or it needs modification to run the func_soil routine again.
// Also requires that func_Classgrid runs first so it uses the latest lake/catchment areas.
// Currently uses just deep soil drainage. 
// Note. We convert this to a rate and it's then converted back to a volume in the GW module. This is because most people are more familiar with mm/day of recharge.
int nrow = GWRechargegrid.nrow();
int ncol = GWRechargegrid.ncol();
int Area = Cellxdim * Cellydim;
int CatchmentCells = 0;
	for (int j = 0; j < ncol; j++){ // count cells in catchment
		for (int i = 0; i < nrow; i++){
			if (Classgrid(i,j) == 11){
			CatchmentCells++;
			}
		}
	}
	//Rcpp::Rcout << "Catchment Cells = " << CatchmentCells << std::endl;
double RechargeRate = Flux(Time,15)/(Timestep*CatchmentCells*Area*CDaysInMonth(Flux(Time,1), Flux(Time,0))); //flux is a volumetric amount for the catchment. Divided by catchment area and days in month to give a rate and assumed rate for entire groundwater area.
	for (int j = 0; j < ncol; j++){
		for (int i = 0; i < nrow; i++){
		GWRechargegrid(i,j) = RechargeRate;
		}
	}
}

// [[Rcpp::export]]
void func_ETGrid (long Time, NumericMatrix GWETgrid, NumericMatrix Flux){
int nrow = GWETgrid.nrow();
int ncol = GWETgrid.ncol();
double ETRate = Flux(Time,27)/CDaysInMonth(Flux(Time,1), Flux(Time,0)); //Convert to daily rate
	//Rcpp::Rcout << "Time: " << Time <<" ETGrid Value: " << ETRate << std::endl;
	for (int j = 0; j < ncol; j++){
		for (int i = 0; i < nrow; i++){
		GWETgrid(i,j) = ETRate;
		}
	}
	//Rcpp::Rcout << "Confirmed: " << Time <<" ETGrid Value: " << GWETgrid(20,20) << std::endl;
}

// [[Rcpp::export]]
void func_LakeDepths(NumericMatrix Lake, long Time, NumericMatrix LakeDepths){
	double CurrentLakedepth = Lake(Time,4);
	//Rcpp::Rcout << "Current Lake Depth " << Lake(Time,4) << std::endl;
	for (int i = 0; i < LakeDepths.nrow(); i++){
		if (CurrentLakedepth + LakeDepths(i,1) > LakeDepths(i,2)){
		LakeDepths(i,3) = CurrentLakedepth + LakeDepths(i,1);
		//Rcpp::Rcout << "Setting lake " << LakeDepths(i,0) << " to " << CurrentLakedepth + LakeDepths(i,1) << std::endl;
		} else {
		LakeDepths(i,3) = LakeDepths(i,2);
		}
	}	//return LakeDepths;
}
