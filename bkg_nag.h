// bkg.h
// Defines functions for calculating the background phi, phiprime, H
// GO
#ifndef bkg_nag_h
#define bkg_nag_h

#include "utils.h"
#include "constants.h"
#include <cmath>
//#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
//#include <stdio.h>

//GSL Imports
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

//NAG Imports
#include <nag.h>
#include <nag_stdlib.h>
#include <nage01.h>
#include <nage02.h>
#include <nagd02.h>
#include <nagx01.h>

class bkg_nag {
    protected:
	double phi_init;
	double phip_init;
	double H_init;
	double lna0;
	double lna1;
	double matching;
	int nsteps;
	step_params p;
	std::vector<double> lna;
	std::vector<double> a;
	std::vector<double> phi;
	std::vector<double> phip;
	std::vector<double> H;
	std::vector<double> epsH;

	std::vector<double> eta;
	std::vector<double> lneta;
	std::vector<double> f;
	std::vector<double> dfdlneta;
	std::vector<double> ddfdlneta;
	std::vector<double> G;
	std::vector<double> Gp; //derivative wrt lneta

	//For Gp interpolation
	/*
	double* lnetar = 0;
	double* Gpr = 0;
	NagError Gpr_fail;
	Nag_Spline Gpr_sp;
	double Gpr_fit = 0.;
	*/

	static int psiprime(double lna, const double psi[], double psip[], void* params);

	//NAG functions
	static void NAG_CALL bkg_out(Integer neq, double *xsol, const double y[], Nag_User *comm);

	static void NAG_CALL psiprime_nag(Integer neq, double x, const double y[], double f[], Nag_User *comm);

	static double NAG_CALL epsHm1(Integer neq, double x, const double y[], Nag_User *comm);

    public:
	int set_inits(double phi_i, double lna_i, double lna_f, int n, step_params p_in){
		//SET_FAIL(Gpr_fail);
		//Gpr_sp.lamda = 0;
		//Gpr_sp.c = 0;
		//Add some form of error checking


		phi_init = phi_i;
		H_init = sqrt( pot::V(phi_init, p)/(3.*Mpl*Mpl) ); //SLOW-ROLL H
		phip_init = -1.*pot::Vp(phi_init, p)/(3.*H_init*H_init); //SLOW-ROLL phip

		std::cout<< "Inits = " << phi_init << "  " << phip_init <<"  " << H_init << "  " << std::endl;
		p = p_in;
		lna0 = lna_i;
		lna1 = lna_f;
		nsteps = n;
		return 0;
	}

	void bkg_nag_calc();

	void approx_init();

	std::vector<double>* getlna(){ return &lna; }
	std::vector<double>* geta(){ return &a; }
	std::vector<double>* getphi(){ return &phi; }
	std::vector<double>* getphip(){ return &phip; }
	std::vector<double>* getH(){ return &H; }
	std::vector<double>* getepsH() { return &epsH; }
	std::vector<double>* geteta() { return &eta; }
	std::vector<double>* getlneta() { return &lneta; }
	std::vector<double>* getG() { return &G; }
	std::vector<double>* getGp() { return &Gp; }

	/*double interp_Gp(const double lneta) { 
		nag_1d_spline_evaluate(lneta, &Gpr_fit, &Gpr_sp, &Gpr_fail); //Add some error checking but OK for now
		return Gpr_fit;
	} // end of interp_Gp*/
}; // end of class bkg_nag








#endif
