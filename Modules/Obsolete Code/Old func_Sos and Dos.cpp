// [[Rcpp::export]]
NumericVector func_Fsos(long Time, double Csr, NumericMatrix Lake, double Timestep, NumericMatrix Isotopes){
	double Fsos;
	double Fsos18O;
	double FsosD;
	Fsos = Csr*Lake(Time,7)*Timestep;
	Fsos18O = Fsos * Isotopes(Time,5);
	FsosD =  Fsos * Isotopes(Time,4);
	NumericVector FsosVals = NumericVector::create(Named("Fsos") = Fsos, Named("Fsos18O") = Fsos18O, Named("FsosD") = FsosD);
	return FsosVals;
}

// [[Rcpp::export]]
NumericVector func_Fdos(long Time, double Csr, NumericMatrix Lake, double Timestep, NumericMatrix Isotopes){
	double Fdos;
	double Fdos18O;
	double FdosD;
	Fdos = Csr*Lake(Time,8)*Timestep;
	Fdos18O = Fdos * Isotopes(Time,7);
	FdosD =  Fdos * Isotopes(Time,6);
	NumericVector FdosVals = NumericVector::create(Named("Fdos") = Fdos, Named("Fdos18O") = Fdos18O, Named("FdosD") = FdosD);
	return FdosVals;
}