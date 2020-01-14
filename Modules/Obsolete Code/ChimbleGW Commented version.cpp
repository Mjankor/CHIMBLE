#include <Rcpp.h>
using namespace Rcpp;
#include "Hypsographic_Calcs.hpp"
#include "Hydrology_functions.hpp"

// [[Rcpp::export]]
List Headcalc(NumericMatrix Topogrid, NumericMatrix GWRechargegrid, NumericMatrix GWPumpinggrid, NumericMatrix GWKgrid, NumericMatrix GWSgrid, NumericMatrix GWbasegrid, NumericMatrix GWheadgrid, NumericMatrix GWboundarygrid, NumericMatrix GWlakebedgrid, NumericMatrix Lakecondgrid, NumericMatrix GWETgrid, NumericMatrix Classgrid, double Timestep, double Runtime, double ETexdepth, NumericMatrix Lakedepth, int Cellxdim, int Cellydim, double Accuracy, int Maxiterations) {

//Experiments suggest that lumping all the functions into one horrible big code block will be faster. (0.112 for functioned C file, 0.012 for this original, linear C++ file). Ideally we would break this up into functions, but it appears each function call copies the matrix, which can be a huge file. 
// Grid names as per R module
// Version notes 
// 19/6/17 Gnotuk hack removed. Output will no longer be a NumericVector. Instead it will be a list structure, consisting of a matrix for lake and catchment surface flows and a numeric vector for mass balance and model parameters. The matrix is going to consists of lake name, inflow, outflow, catchment surface flow and combined flow cumulative for each model simulation. The vector will contain cumulative values for Massinit, Massfinal, F_Head flux, Extflows, Masserror, ET, Timestep, Runtime, Recharge, Conductivity.

// Require Pumping grid , Topogrid, GWKgrid, GWSgrid, GWbasegrid, Recharge grid, GWheadgrid, GWboundarygrid, Classgrid, Timestep, Cellxdim, Cellydim
Rcpp::Rcout << "Passed to here " << std::endl;

int nrow = Topogrid.nrow();
int ncol = Topogrid.ncol();
int time = 0;
double Hold = 0;
double maxdiff = Accuracy+1;
int itercount = 0;
//The below are the five variables we want to calculate
double outseepage = 0;
double inflows = 0;
double outflows = 0; // total flows lost outside the catchment (head > topo)
double surfaceflows = 0; // Flows into the catchments
//double Alllakeflows[10] = 0;
double massinit = 0;
double massfinal = 0;
double extflows = 0;
double ETflux = 0;
double Fixedheadflows = 0;
double ETfluxall = 0;
NumericMatrix Allflows (9,5); //columns are lake name, lake inflows, lake outseepage, catchment surface flow and combined (net) flows.
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
NumericMatrix Lakecoords (nrow*ncol,4); //matrix is cell coords (0,1), depth (2) and flux per timestep (3) 
NumericMatrix Landcoords (nrow*ncol,2); 
NumericMatrix Allcoords (nrow*ncol,2); 
// debug stuff below
NumericMatrix diff (nrow,ncol);
bool Overlandflows = FALSE;

// ************************ New method based on spreadsheet model **********************
// *****************************************************************************
// ================= Before iteration loop (Q, S, HCOF, RHS) ===================
// *****************************************************************************

// ******************** Add lake numbers to the Allflows matrix ************************
// And some column names.
colnames (Allflows) = CharacterVector::create("Lake ID", "Seepage_in", "Seepage_out", "Surface Flows", "Net_flows");
for (int i = 0; i < Allflows.nrow(); i++){
Allflows(i,0) = i+1;
Allflows(i,1) = 0;
Allflows(i,2) = 0;
Allflows(i,3) = 0;
Allflows(i,4) = 0;
//This just puts lake numbers 1 - 9 in the first column of the Allflows matrix and fills the rest of the matrix with 0.
}

// ************************ Define all lake cells **********************
//This is no longer needed as the lake now uses Cauchy boundaries. Therefore only a few cells will use fixed head boundaries, and they can be caught with an if statement. 
//However, this is still useful to cut down on HCOF and RHS calcs for the lake, and to define the areas of "confined" aquifers.
//Remaining cells can be added to a 2 column matrix, which will then be used for the for loops in the solution calculation. This will cut down the number of cells requiring calculation and enforce the boundary conditions.
// First two columns (0 and 1) are grid coordinates, Third column (2) is lake depth. Fourth column (3) is flux to and from the lake cell via the Cauchy boundary for each timestep. 
//First fill our coord matrix with NA
std::fill(Lakecoords.begin(), Lakecoords.end(), NumericVector::get_na());
//Then start a counter
int Lakecounter = 0;
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
	if (GWboundarygrid(i,j) == 0 && Classgrid(i,j) < 10){
	Lakecoords (Lakecounter,0) = i;
	Lakecoords (Lakecounter,1) = j;
	Lakecoords (Lakecounter,2) = Lakedepth(Classgrid(i,j)-1,1);
	//Needs the -1 to account for C++ starting the Lakedepth matrix at 0.
	Lakecounter++;
	
	//and we use counter further down (in the solution) to stop running into the NA cells.
	}
	}
}
//Cells on the land (not fixed head). Not this limits the use of fixed head and fixed flow combined cells.
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
	}
	}


// ************************ Define combined fluxes **********************
//Only applied to land cells. Recharge and Evap of lake done in CHIMBLE.

/* This method is slightly slower than the twin for loops
int ncells = nrow*ncol;
for (i = 0; i < ncells; i++){
Qgrid [i] = GWRechargegrid [i]*Cellxdim*Cellydim+ GWPumpinggrid[i];
}
*/
for (int row = 0; row < Landcounter; row++){
			int i = Landcoords(row,0);
			int j = Landcoords(row,1);
			Qgrid (i,j) = GWRechargegrid (i,j)*Cellxdim*Cellydim + GWPumpinggrid(i,j);
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
	//below option uses matrixes from above code.
	//	HCOF (i,j) = -Sgrid(i,j)/Timestep ;
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
	
//	if (time == 3599){
//	return GWheadgrid;
//	}
	
	// ********************************** Define RHS **********************************
	// Requires Combined flux matrix, Cell Storage matrix, GWheadgrid, Timestep
	// Should be updated with each iteration with head dependent nodes.
	for (int j = 0; j < ncol; j++){
	for (int i = 0; i < nrow; i++){
		RHSbase (i,j) = -Qgrid(i,j) - (GWSgrid (i,j) * Cellxdim * Cellydim * GWheadgrid(i,j))/Timestep ;
	//below option uses matrixes from above code.
	//		RHS (i,j) = -Qgrid(i,j) - Sgrid(i,j) * GWheadgrid(i,j)/Timestep ;
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
		//And update the flux (flux to/from lake per unit time)
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

		//update iteration count
		//Rcpp::Rcout << "Node Height for iter " << (GWheadgrid(2,2)) << std::endl;
		//Rcpp::Rcout << "Iteration counter " << itercount << std::endl;
		//Rcpp::Rcout << "Maxdiff " << maxdiff << std::endl;
		itercount++;

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
			// + Qgrid(i,j); // Qgrid removed as rainfall on the lake will be included in the CHIMBLE mass balance.
			}
			}
		
			//edges
			for (int j = 1;j < ncol-1; j++){
			Flowsgrid(0,j) = Timestep*(CW(0,j)*(GWnewheadgrid(0,j-1)-GWnewheadgrid(0,j)) + CE(0,j)*(GWnewheadgrid(0,j+1)-GWnewheadgrid(0,j)) + CS(0,j)*(GWnewheadgrid(1,j)-GWnewheadgrid(0,j))) + Flowsgrid(0,j);// + Qgrid(0,j);
			}

			for (int j = 1;j < ncol-1; j++){
			Flowsgrid(nrow-1,j) = Timestep*(CW(nrow-1,j) *(GWnewheadgrid(nrow-1,j-1)-GWnewheadgrid(nrow-1,j)) + CE(nrow-1,j)*(GWnewheadgrid(nrow-1,j+1)-GWnewheadgrid(nrow-1,j)) + CN(nrow-1,j)*(GWnewheadgrid(nrow-2,j)-GWnewheadgrid(nrow-1,j))) + Flowsgrid(nrow-1,j);// + Qgrid(nrow-1,j);
			}

			for (int i = 1;i < nrow-1; i++){
			Flowsgrid(i,0) = Timestep*(CE(i,0)*(GWnewheadgrid(i,1)-GWnewheadgrid(i,0)) + CN(i,0)*(GWnewheadgrid(i-1,0)-GWnewheadgrid(i,0)) + CS(i,0)*(GWnewheadgrid(i+1,0)-GWnewheadgrid(i,0))) + Flowsgrid(i,0);// + Qgrid(i,0);
			}

			for (int i = 1;i < nrow-1; i++){
			Flowsgrid(i,ncol-1) = Timestep*(CW(i,ncol-1)*(GWnewheadgrid(i,ncol-2)-GWnewheadgrid(i,ncol-1)) + CN(i,ncol-1)*(GWnewheadgrid(i-1,ncol-1)-GWnewheadgrid(i,ncol-1)) + CS(i,ncol-1)*(GWnewheadgrid(i+1,ncol-1)-GWnewheadgrid(i,ncol-1))) + Flowsgrid(i,ncol-1);// + Qgrid(i,ncol-2);
			}
		
			//corners
			Flowsgrid(0,0) = Timestep*(CE(0,0)*(GWnewheadgrid(0,1)-GWnewheadgrid(0,0)) + CS(0,0)*(GWnewheadgrid(1,0)-GWnewheadgrid(0,0))) + Flowsgrid(0,0);// + Qgrid(0,0);
			
			Flowsgrid(nrow-1,0) = Timestep*(CE(nrow-1,0)*(GWnewheadgrid(nrow-1,1)-GWnewheadgrid(nrow-1,0)) + CN(nrow-1,0)*(GWnewheadgrid(nrow-2,0)-GWnewheadgrid(nrow-1,0))) + Flowsgrid(nrow-1,0);// + Qgrid(nrow-1,0);
			
			Flowsgrid(0,ncol-1) = Timestep*(CW(0,ncol-1)*(GWnewheadgrid(0,ncol-2)-GWnewheadgrid(0,ncol-1)) + CS(0,ncol-1)*(GWnewheadgrid(1,ncol-1)-GWnewheadgrid(0,ncol-1))) + Flowsgrid(0,ncol-1);// + Qgrid(0,ncol-1);
			
			Flowsgrid(nrow-1,ncol-1) = Timestep*(CW(nrow-1,ncol-1)*(GWnewheadgrid(nrow-1,ncol-2)-GWnewheadgrid(nrow-1,ncol-1)) + CN(nrow-1,ncol-1)*(GWnewheadgrid(nrow-2,ncol-1)-GWnewheadgrid(nrow-1,ncol-1))) + Flowsgrid(nrow-1,ncol-1);// + Qgrid(ncol-1,ncol-1);

				//Rcpp::Rcout << "Lake number: " << Allflows(0,0) << std::endl;				
				//Rcpp::Rcout << "Check on inflow: " << Allflows(0,1) << std::endl;		
				//Rcpp::Rcout << "Check on outseepage: " << Allflows(0,2) << std::endl;

// ************************ Quantify external flows balance *****************************
		for (int row = 0; row < Landcounter; row++){
		int i = Landcoords(row,0);
		int j = Landcoords(row,1);
		extflows += Qgrid(i,j)*Timestep;
		}
		
		
// ************************ Quantify lake in and outseepage *****************************
// This is cumulative over the runtime of the model. \
// ***Check this. May be tallying up all lake flows as Bullen Merri flows. Cancel that. This works as intended, with a tally for each lake.
// With updated model it needs to separate out inflows and outflows. 
for (int row = 0; row < Lakecounter; row++){
	int i = Lakecoords(row,0);
	int j = Lakecoords(row,1);
	int lakenum = (Classgrid(i,j));
	if (Lakecoords(row,3) < 0){
	Allflows(lakenum-1,1) += Lakecoords(row,3)*Timestep;
	} else {
	Allflows(lakenum-1,2) += Lakecoords(row,3)*Timestep;
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
		//if class grid suggests catchment, then it goes is partioned based on the catchment to Allflows(Classgrid(i,j) - 11, 3), else it goes into outflows (double)
			if (Classgrid (i,j) < 20 && GWnewheadgrid (i,j) > Topogrid(i,j)){
			Allflows((Classgrid (i,j) - 11),3) += (GWnewheadgrid(i,j) - Topogrid(i,j))*Sgrid(i,j);
			GWnewheadgrid (i,j) = Topogrid (i,j);
			} 
			if (Classgrid(i,j) == 20 && GWnewheadgrid (i,j) > Topogrid(i,j)){
			outflows += (GWnewheadgrid(i,j) - Topogrid(i,j))*Sgrid(i,j);
			GWnewheadgrid (i,j) = Topogrid (i,j);
			}
			
/*			if (Classgrid (i,j) == 11 && GWnewheadgrid (i,j) > Topogrid(i,j)){
			inflows += (GWnewheadgrid(i,j) - Topogrid(i,j))*Sgrid(i,j);
			GWnewheadgrid (i,j) = Topogrid (i,j);
			}
			if (Classgrid (i,j) == 20 && GWnewheadgrid (i,j) > Topogrid(i,j)){
			outflows += (GWnewheadgrid(i,j) - Topogrid(i,j))*Sgrid(i,j);
			GWnewheadgrid (i,j) = Topogrid (i,j);
			}
			//hack for gnotuk
			if (Classgrid (i,j) == 12 && GWnewheadgrid (i,j) > Topogrid(i,j)){
			gnotukflows += (GWnewheadgrid(i,j) - Topogrid(i,j))*Sgrid(i,j);
			GWnewheadgrid (i,j) = Topogrid (i,j);
			}				
*/		}
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
for (int i = 0; i < Allflows.nrow(); i++){
	inflows += Allflows(i,1);
	outseepage += Allflows (i,2);
	surfaceflows += Allflows (i,3);
	Allflows (i,4) = -Allflows(i,1) - Allflows(i,2) + Allflows (i,3);
}

// ***Modify this to account for all lake and catchment flows
	double masserror = 100*(massinit - massfinal - Fixedheadflows + extflows + outseepage + inflows - surfaceflows - ETfluxall - outflows)/( Fixedheadflows + extflows + outseepage + inflows - surfaceflows - ETfluxall - outflows);
	
	//double masserror = 100*(massinit - massfinal - Fixedheadflows + extflows + Alllakeflows[1] + Alllakeflows[2] - gnotukflows  - outflows - inflows - ETfluxall)/(-Fixedheadflows + extflows - Alllakeflows[1] - Alllakeflows[2]  + outflows + inflows + ETfluxall);

//Rcpp::Rcout << "Massdelta " << (massinit - massfinal - Fixedheadflows + extflows + Alllakeflows[1] + Alllakeflows[2]  - outflows - inflows - ETfluxall) << std::endl;
// Data to return as vector
//Overland flows lost, overland to lake, lake flows, mass balance before, mass balance after. 
//If we decide to modify this to include P flows and conductance through the lake floor, the Lake flows can be modified to show this new result.
// GWheadgrid is modified during the C++ function run.

double rec=GWRechargegrid(1,1);
double K=GWKgrid(1,1);

NumericVector Stats = NumericVector::create(Named("Massinit")=massinit, Named("Massfinal")=massfinal, Named("Mass error")=masserror, Named("F_Head flux")=Fixedheadflows, Named("Extflows")=extflows, Named("Lake seepage(in)")=inflows, Named("Lake seepage(out)") = outseepage, Named("Surface flows")=surfaceflows, Named("Non catchment outflows")=outflows, Named("ET")=ETfluxall, Named("Timestep")=Timestep, Named("Runtime")=Runtime, Named("Recharge")=rec, Named("Conductivity")=K);

// ***Modify this to account for all and catchment flows
List GW = List::create(
Named("Flows") = Allflows,
Named("Stats") = Stats);//, Named("inflows") = inflows, Named("outflows") = outflows;
// removed - Named("Masserror")=masserror, Named("Gnotukflows")=Alllakeflows[2], Named("Gnotuksurfflows")=gnotukflows,
return GW;
                                 
} // End c function


// [[Rcpp::export]]
void func_ClassGrid(long Time, NumericMatrix Topogrid, NumericMatrix Lake, NumericMatrix GWcatchmentgrid, NumericMatrix Classgrid, NumericMatrix Lakevolumes){
int nrow = Topogrid.nrow();
int ncol = Topogrid.ncol();
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
		Lakecondgrid(i,j) = 0;
		}
	}
}
}

// [[Rcpp::export]]
void func_RechargeGrid (long Time, int Cellxdim, int Cellydim, NumericMatrix GWRechargegrid, NumericMatrix Classgrid, NumericMatrix Flux){
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
double RechargeRate = Flux(Time,15)/(CatchmentCells*Area); //flux is a volumetric amount for the catchment. Divided by catchment area to give a rate and assumed rate for entire groundwater area.
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
double ETRate = Flux(Time,23)/CDaysInMonth(Flux(Time,1), Flux(Time,0)); //Convert to daily rate
	for (int j = 0; j < ncol; j++){
		for (int i = 0; i < nrow; i++){
		GWETgrid(i,j) = ETRate;
		}
	}
}
