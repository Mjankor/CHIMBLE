// [[Rcpp::export]]
void func_IsotopeFluxes (NumericMatrix Flux, long Time, std::string Isotype, NumericMatrix Isotopes, NumericMatrix Met, double AtmosphericShift, double Turbulence, double Feedback, NumericMatrix Lakevolumes, NumericMatrix Lake, std::string Interpolationtype, std::string Delim){
for (int i=4; i < Isotopes.ncol(); i++){
    StringVector Isocolumns = colnames(Isotopes);
    std::string ColNameString = as<std::string>(Isocolumns(i));
    //Rcpp::Rcout << "Front = " << StrSplitFront(ColNameString,Delim) << std::endl;
    //Rcpp::Rcout << "Back = " << StrSplitBack(ColNameString,Delim) << std::endl;
 
    Isotopes(Time-1,18) = func_dE(Time, "18O", Isotopes, Met, AtmosphericShift, Turbulence, Feedback, Lakevolumes, Lake, Interpolationtype)
    Isotopes(Time-1,19) = func_dE(Time, "D", Isotopes, Met, AtmosphericShift, Turbulence, Feedback, Lakevolumes, Lake, Interpolationtype)
    //Isotopes(Time-1,20) = Met(Time-1,8);
    //Isotopes(Time-1,21) = Met(Time-1,9);
 
    if (StrSplitBack(ColNameString,Delim) == "D"){
        switch(StrSplitFront(ColNameString,Delim)){
            case "RESsl" :
            if (Lake(Time,7) == 0){
                Isotopes(Time,i) = 0;           
            } else {
                if (Met(Time,11) == 0){
                    Isotopes(Time,i) = 0; //Lake(Time-1,7);
                //  Lake(Time,8) = Lake(Time-1,8)+Flux(Time-1,19)+Flux(Time-1,5)+Flux(Time-1,17) - Flux(Time-1,7) - Flux(Time-1,9) - Flux(Time-1,8) - Flux(Time-1,20);
                } else if (Lake(Time-1,7) < Flux(Time-1,8) + Flux(Time-1,7) + Flux(Time-1,20) - Flux(Time-1,5) - Flux(Time-1,17)){
                    double ExcessFluxE = (Flux(Time-1,8) + Flux(Time-1,7) + Flux(Time-1,20) - Flux(Time-1,5) - Flux(Time-1,17) - Lake(Time-1,7)) * (Flux(Time-1,7)/(Flux(Time-1,7) + Flux(Time-1,8)));
                    double ExcessFluxSos = (Flux(Time-1,8) + Flux(Time-1,7) + Flux(Time-1,20) - Flux(Time-1,5) - Flux(Time-1,17) - Lake(Time-1,7)) * (Flux(Time-1,8)/(Flux(Time-1,7) + Flux(Time-1,8)));
                    Isotopes(Time,i) = 0;//(Lake(Time-1,7)*Isotopes(Time-1,4) + Flux(Time-1,5)*Met(Time-1,9) + Flux(Time-1,17)*Isotopes(Time-1,14) - (Flux(Time-1,8) - ExcessFluxSos)*Isotopes(Time-1,4) - (Flux(Time-1,7) - ExcessFluxE)*Isotopes(Time-1,19))/Lake(Time,7);
                    //Lake(Time,8) = Lake(Time-1,8) + Flux(Time-1,19) - Flux(Time-1,9) - ExcessFluxE - ExcessFluxSos;
                } else {
                    Isotopes(Time,i) = (Lake(Time-1,7)*Isotopes(Time-1,4) + Flux(Time-1,5)*Met(Time-1,9) + Flux(Time-1,17)*Isotopes(Time-1,14) - Flux(Time-1,8)*Isotopes(Time-1,4) - Flux(Time-1,7)*Isotopes(Time-1,19) - Flux(Time-1,20)*Isotopes(Time-1,4))/Lake(Time,7);
                    //Lake(Time,8) = Lake(Time-1,8) + Flux(Time-1,19) - Flux(Time-1,9);
                }
            }
            break;
            case "RESdl" :
            if (Lake(Time,8) == 0){
                Isotopes(Time,i) = 0;           
            } else {
                if (Met(Time,11) == 0){
                    //Isotopes(Time,i) = 0; //Lake(Time-1,7);
                    Isotopes(Time,i) = (Lake(Time-1,8)*Isotopes(Time-1,6)+Flux(Time-1,19)*Isotopes(Time-1,16)+Flux(Time-1,5)*Met(Time-1,9)+Flux(Time-1,17)*Isotopes(Time-1,14) - Flux(Time-1,7)*Isotopes(Time-1,19) - Flux(Time-1,9)*Isotopes(Time-1,6) - Flux(Time-1,8)*Isotopes(Time-1,4) - Flux(Time-1,20)*Isotopes(Time-1,6))/Lake(Time,8);
                } else if (Lake(Time-1,7) < Flux(Time-1,8) + Flux(Time-1,7) + Flux(Time-1,20) - Flux(Time-1,5) - Flux(Time-1,17)){
                    double ExcessFluxE = (Flux(Time-1,8) + Flux(Time-1,7) + Flux(Time-1,20) - Flux(Time-1,5) - Flux(Time-1,17) - Lake(Time-1,7)) * (Flux(Time-1,7)/(Flux(Time-1,7) + Flux(Time-1,8)));
                    double ExcessFluxSos = (Flux(Time-1,8) + Flux(Time-1,7) + Flux(Time-1,20) - Flux(Time-1,5) - Flux(Time-1,17) - Lake(Time-1,7)) * (Flux(Time-1,8)/(Flux(Time-1,7) + Flux(Time-1,8)));
                    //Isotopes(Time,i) = 0;//(Lake(Time-1,7)*Isotopes(Time-1,4) + Flux(Time-1,5)*Met(Time-1,9) + Flux(Time-1,17)*Isotopes(Time-1,14) - (Flux(Time-1,8) - ExcessFluxSos)*Isotopes(Time-1,4) - (Flux(Time-1,7) - ExcessFluxE)*Isotopes(Time-1,19))/Lake(Time,7);
                    Isotopes(Time,i) = (Lake(Time-1,8)*Isotopes(Time-1,6) + Flux(Time-1,19)*Isotopes(Time-1,16) - Flux(Time-1,9)*Isotopes(Time-1,6) - ExcessFluxE*Isotopes(Time-1,19) - ExcessFluxSos*Isotopes(Time-1,6)/Lake(Time,8);
                } else {
                    //Isotopes(Time,i) = (Lake(Time-1,7)*Isotopes(Time-1,4) + Flux(Time-1,5)*Met(Time-1,9) + Flux(Time-1,17)*Isotopes(Time-1,14) - Flux(Time-1,8)*Isotopes(Time-1,4) - Flux(Time-1,7)*Isotopes(Time-1,19) - Flux(Time-1,20)*Isotopes(Time-1,4))/Lake(Time,7);
                    Isotopes(Time,i) = (Lake(Time-1,8)*Isotopes(Time-1,6) + Flux(Time-1,19)*Isotopes(Time-1,16) - Flux(Time-1,9)*Isotopes(Time-1,6))/Lake(Time,8);
                }
            }
            break;
            case "RESss" :
            if (Lake(Time,9) == 0){
                Isotopes(Time,i) = 0;
            } else {
                double FssiD = 0;
                if (Flux(Time-1,4) + Flux(Time-1,6) > 0){
                    double FssiD = (Flux(Time-1,6)*Isotopes(Time-1,12) + Flux(Time-1,4)*Met(Time-1,9))*(Flux(Time-1,11)/(Flux(Time-1,11) + Flux(Time-1,13) + Flux(Time-1,15)+ Flux(Time-1,16)));
                } else {
                    FssiD = 0;
                }
                Isotopes(Time,i) = (Lake(Time-1,9)*Isotopes(Time-1,8) + FssiD - Flux(Time-1,12)*Isotopes(Time-1,8))/Lake(Time,9);
            }
            break;
            case "RESds" :
                if (Lake(Time,10) == 0){
                Isotopes(Time,i) = 0;
            } else {
                double FssdD = 0;
                if (Flux(Time-1,4) + Flux(Time-1,6) > 0){
                    double FssdD = (Flux(Time-1,6)*Isotopes(Time-1,12) + Flux(Time-1,4)*Met(Time-1,9))*(Flux(Time-1,13)/(Flux(Time-1,11) + Flux(Time-1,13) + Flux(Time-1,15)+ Flux(Time-1,16)));
                } else {
                    FssiD = 0;
                }
                Isotopes(Time,i) = (Lake(Time-1,10)*Isotopes(Time-1,10) + FssdD - Flux(Time-1,14)*Isotopes(Time-1,10)/Lake(Time,10);
            }
            break;
            case "RESsp" :
            if (Lake(Time,11) == 0{
                Isotopes(Time,i) = 0;
            } else {
                Isotopes(Time,i) = (Lake(Time-1,11)*Isotopes(Time-1,12) + Flux(Time-1,18)*Met(Time-1,9) - Flux(Time-1,6)*Isotopes(Time-1,12))/Lake(Time,11);
            }
            break;
            case "RESin" :
                if (Lake(Time,12) == 0{
                Isotopes(Time,i) = 0;
            } else {
                double FroD = 0;
                if (Flux(Time-1,4) + Flux(Time-1,6) > 0){
                    double FroD = (Flux(Time-1,6)*Isotopes(Time-1,12) + Flux(Time-1,4)*Met(Time-1,9))*(Flux(Time-1,13)/(Flux(Time-1,11) + Flux(Time-1,13) + Flux(Time-1,15)+ Flux(Time-1,16)));
                } else {
                    FssiD = 0;
                }
                 
                 
                Isotopes(Time,i) = (Lake(Time-1,12)*Isotopes(Time-1,14) + )/Lake(Time,12);
            }
            break;
            case "GW" :
             
     
            break;
        }
    }
}
}
//Fdlm - This should be done after isotopes and chemistry for the lake is done. Needs to be a two step process.