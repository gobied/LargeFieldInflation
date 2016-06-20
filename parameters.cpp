// parameters.cpp
// Defines the constructors for parameters.h

#include "parameters.h"

//gsr_params
gsr_params::gsr_params(std::vector<double>* Gp_in, std::vector<double>* lneta_in): Gp(Gp_in), lneta(lneta_in){
                SET_FAIL(fail);
		SET_FAIL(fail_lin);
                sp.lamda = 0;
                sp.c = 0;
		sp_lin.lamda = 0;
		sp_lin.c = 0;	

                lnetar = NAG_ALLOC(lneta->size() - 4, double);
		etar = NAG_ALLOC(lneta->size()-4, double);
                Gpr = NAG_ALLOC(Gp->size() - 4, double);
		scale = 1.;//(*(std::max_element(Gp->begin() + 2, Gp->begin() + Gp->size() - 3)))/3.06115;

                for(int ir=0; ir<(lneta->size() - 4); ir++){
                        Gpr[ir] = (Gp->rbegin()[2+ir])/scale;
                        lnetar[ir] = lneta->rbegin()[2+ir];
			etar[ir] = exp(lneta->rbegin()[2+ir]);
                } // end of for ir loop

		std::cout << "Min/Max lneta = " << lnetar[0] << ", " << lnetar[lneta->size() -5] << std::endl;

                nag_1d_spline_interpolant((Integer)lneta->size()-4, lnetar, Gpr, &sp, NAGERR_DEFAULT); 
		nag_1d_spline_interpolant((Integer)lneta->size()-4, etar, Gpr, &sp_lin, NAGERR_DEFAULT);
}

double gsr_params::interp_Gp(double lneta){
                nag_1d_spline_evaluate(lneta, &fit, &sp, NAGERR_DEFAULT); //interpolating in logspace but OK for now
                return fit;
}

double gsr_params::interp_Gp_lin(double lneta){
		nag_1d_spline_evaluate(exp(lneta), &fit_lin, &sp_lin, NAGERR_DEFAULT);
		return fit_lin;
}

//eta_params
eta_params::eta_params(std::vector<double>* a_in, std::vector<double>* H_in): a(a_in), H(H_in){
	//SET_FAIL(fail);
	sp.lamda = 0;
	sp.c = 0;

	aHinv = NAG_ALLOC(a->size(), double);
	pa = NAG_ALLOC(a->size(), double);

	for(int i=0; i<a->size(); i++){
		aHinv[i] = 1./((*a)[i] * (*H)[i]);
		pa[i] = (*a)[i];
	}

	nag_1d_spline_interpolant((Integer) a->size(), pa, aHinv, &sp, NAGERR_DEFAULT);

} // end of eta_params constructor

double eta_params::interp_aHinv(double lna_in){
	nag_1d_spline_evaluate(exp(lna_in), &fit, &sp, NAGERR_DEFAULT);
	return fit;
} // end of itnerp_aHinv
