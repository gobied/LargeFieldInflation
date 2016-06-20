// gsr.h
// Calculates Gmin, I0 and I1 and fills Delta2GSR

#ifndef approx_h
#define approx_h

#include "utils.h"
#include "constants.h"
#include "bkg.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

//NAG includes
#include <nag.h>
#include <nag_stdlib.h>
#include <nagd01.h>
#include <nagx01.h>

class approx {
    protected:
	bkg* mybkg;
	std::vector<double>* a;
	std::vector<double>* H;
	std::vector<double>* eta;
	std::vector<double>* lneta;
	std::vector<double>* G;
	std::vector<double>* Gp;
	std::vector<double>* ks;
	int kindex;
	std::vector<double> I0;
	std::vector<double> I1;
	std::vector<double> Delta2GSR;
	std::vector<double> Delta2NLA;

	//for G Interpolation
	double* lnetar = 0;
	double* Gr = 0;
	//static NagError Gpr_fail;
	Nag_Spline Gr_sp;
	double Gr_fit = 0.;

	//NLA stuff
	std::vector<double> B;
	std::vector<double> cosphi;
	std::vector<double> deltaI0;
	std::vector<double> deltaI1;
	double ns = 0.96;
	double nsm1 = -0.04;
	double phi_bar = -1.*M_PI*(1.-ns)/2.;
	double cosphibar = cos(phi_bar);
	double sinphibar = sin(phi_bar);
	double kp = 0.05;
	double I1_bar = M_PI*(1. - ns)/(2.*sqrt(2));
	double fstar = -21000.; //get this from the bkg but OK for now
	


	void nag_qp_free(Nag_QuadProgress& qp);

    public:
	//approx(std::vector<double>* a_in, std::vector<double>* H_in, std::vector<double>* eta_in, std::vector<double>* G_in, std::vector<double>* Gp_in, std::vector<double>* ks_in): a(a_in), H(H_in), eta(eta_in), G(G_in), Gp(Gp_in), ks(ks_in) {}	
	approx(bkg* mybkg_in, std::vector<double>* ks_in): mybkg(mybkg_in), a(mybkg_in->geta()), H(mybkg_in->getH()), eta(mybkg_in->geteta()), lneta(mybkg_in->getlneta()), G(mybkg_in->getG()), Gp(mybkg_in->getGp()), ks(ks_in) {
		Gr_sp.lamda = 0;
		Gr_sp.c = 0;

		lnetar = NAG_ALLOC(lneta->size()-4, double);
		Gr = NAG_ALLOC(lneta->size()-4, double);

		for(int i=0; i<lneta->size()-4; i++){
			Gr[i] = (G->rbegin()[i+2]);
			lnetar[i] = lneta->rbegin()[i+2];
		}

		std::cout << "MIN/MAX lneta in G " << lnetar[0] << "   " << lnetar[lneta->size() - 5] << std::endl;

		nag_1d_spline_interpolant((Integer)lneta->size()-4, lnetar, Gr, &Gr_sp, NAGERR_DEFAULT);


	} // end of constructor	

	void writeGSR(std::string name, bool suf);
	void writeNLA(std::string name, bool suf);

	double interp_G(double lneta){
		nag_1d_spline_evaluate(lneta, &Gr_fit, &Gr_sp, NAGERR_DEFAULT);
		return Gr_fit;
	}

	//Consider putting approx_init() here but OK for now
	static double NAG_CALL I0integrand(double lneta_in, Nag_Comm *comm){
		gsr_params* p = (gsr_params *)comm->p;
		//nag_1d_spline_evaluate(lneta_in, &(p->fit), p->sp, NAGERR_DEFAULT);
		//std::cout << "THIS IS IN I0 " << ((p->get_kc()) * exp(lneta_in)) << std::endl;
		return (win::W((p->get_kc()) * exp(lneta_in))*(p->interp_Gp_lin(lneta_in)));
	} // end I0 integrand

	static double NAG_CALL I1integrand(double lneta_in, Nag_Comm *comm){
		gsr_params* p = (gsr_params *)comm->p;
		//nag_1d_spline_evaluate(lneta_in, &(p->fit), p->sp, NAGERR_DEFAULT);
		return (win::X((p->get_kc()) * exp(lneta_in))*(p->interp_Gp_lin(lneta_in)));
	} // end I1 integrand

	static double NAG_CALL I0integrand_osc(double lneta_in, Nag_User *comm){
		gsr_params* p = (gsr_params *)comm->p;
		//nag_1d_spline_evaluate(lneta_in, &(p->fit), p->sp, NAGERR_DEFAULT);
		return (win::W((p->get_kc()) * exp(lneta_in))*(p->interp_Gp_lin(lneta_in)));
	} // end I0 integrand

	static double NAG_CALL I1integrand_osc(double lneta_in, Nag_User *comm){
		gsr_params* p = (gsr_params *)comm->p;
		//nag_1d_spline_evaluate(lneta_in, &(p->fit), p->sp, NAGERR_DEFAULT);
		return (win::X((p->get_kc()) * exp(lneta_in))*(p->interp_Gp_lin(lneta_in)));
	} // end I1 integrand

	static double NAG_CALL I0integrand_osc2(double eta_in, Nag_User *comm){
		gsr_params* p = (gsr_params *)comm->p;
		//nag_1d_spline_evaluate(lneta_in, &(p->fit), p->sp, NAGERR_DEFAULT);
		return (win::W((p->get_kc()) * eta_in)*(p->interp_Gp_lin(log(eta_in)))/eta_in);
	} // end I0 integrand

	static double NAG_CALL I1integrand_osc2(double eta_in, Nag_User *comm){
		gsr_params* p = (gsr_params *)comm->p;
		//nag_1d_spline_evaluate(lneta_in, &(p->fit), p->sp, NAGERR_DEFAULT);
		return (win::X((p->get_kc()) * eta_in)*(p->interp_Gp_lin(log(eta_in)))/eta_in);
	} // end I1 integrand


	static double NAG_CALL I0integrand_gen(double eta_in, Nag_User *comm){
		gsr_params* p = (gsr_params *)comm->p;
		//nag_1d_spline_evaluate(lneta_in, &(p->fit), p->sp, NAGERR_DEFAULT);
		return (win::W((p->get_kc()) * eta_in)*(p->interp_Gp_lin(log(eta_in)))/eta_in);
	} // end I0 integrand

	static double NAG_CALL I1integrand_gen(double eta_in, Nag_User *comm){
		gsr_params* p = (gsr_params *)comm->p;
		//nag_1d_spline_evaluate(lneta_in, &(p->fit), p->sp, NAGERR_DEFAULT);
		return (win::X((p->get_kc()) * eta_in)*(p->interp_Gp_lin(log(eta_in)))/eta_in);
	} // end I1 integrand

	//GSL Functions without the matching in eta
	static double gsl_I0integrand(double lneta_in, void* params){
		gsr_params* p = (gsr_params*)params;
		return (win::W((p->get_kc())*exp(lneta_in)) * (p->interp_Gp_lin(lneta_in))); // the 1.e-75 takes care of the matching
	}
	static double gsl_I1integrand(double lneta_in, void* params){
		gsr_params* p = (gsr_params*)params;
		return (win::X((p->get_kc())*exp(lneta_in)) * (p->interp_Gp_lin(lneta_in)));
	}


	void approx_calc();

	std::vector<double>* getDelta2GSR() { return &Delta2GSR; }
	std::vector<double>* getDelta2NLA() { return &Delta2NLA; }

};


#endif
