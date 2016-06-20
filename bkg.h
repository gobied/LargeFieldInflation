// bkg.h
// Defines functions for calculating the background phi, phiprime, H

#ifndef bkg_h
#define bkg_h

#include "utils.h"
#include "constants.h"
#include <cmath>
//#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
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

class bkg {
    protected:
	double phi_init;
	double phip_init;
	double H_init;
	double lna0;
	double lna1;
	double matching;
	int nsteps;
	std::vector<double> lna;
	std::vector<double> a;
	std::vector<double> lna_matched;
	std::vector<double> a_matched;
	std::vector<double> phi;
	std::vector<double> phip;
	std::vector<double> H;
	std::vector<double> epsH;
	std::vector<double> etaH;
	std::vector<double> delta2;

	std::vector<double> eta;
	std::vector<double> lneta;
	std::vector<double> f;
	std::vector<double> dfdlneta;
	std::vector<double> ddfdlneta;
	std::vector<double> G;
	std::vector<double> Gp; //derivative wrt lneta
	//std::vector<double> df_num;
	//std::vector<double> ddf_num;
	//std::vector<double> Gp_num; //numerical derivative

	//For Gp interpolation
	/*
	double* lnetar = 0;
	double* Gpr = 0;
	NagError Gpr_fail;
	Nag_Spline Gpr_sp;
	double Gpr_fit = 0.;
	*/

	static int psiprime(double lna, const double psi[], double psip[], void* params);

    public:
	step_params p;
	int indexetalo;
	int indexetahi;

	int set_inits(double phi_i, double lna_i, double lna_f, int n, step_params p_in){
		//SET_FAIL(Gpr_fail);
		//Gpr_sp.lamda = 0;
		//Gpr_sp.c = 0;
		//Add some form of error checking

		p = p_in;

		phi_init = phi_i;
		H_init = sqrt( pot::V(phi_init, p)/(3.*Mpl*Mpl) ); //SLOW-ROLL H
		phip_init = 0.;//-1.*pot::Vp(phi_init, p)/(3.*H_init*H_init); //SLOW-ROLL phip

		std::cout<< "Inits = " << phi_init << "  " << phip_init <<"  " << H_init << "  " << std::endl;
		std::cout << "IN BKG " << p.c << std::endl;
		lna0 = lna_i;
		lna1 = lna_f;
		nsteps = n;
		return 0;
	}

	void bkg_calc();

	void approx_init();

	void nag_approx_init();

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

	static double NAG_CALL etaintegrand(double lna_in, Nag_Comm *comm){
		eta_params* p = (eta_params*)comm->p;
		return (p->interp_aHinv(lna_in));
	} // end of etaintegrand

	/*double interp_Gp(const double lneta) { 
		nag_1d_spline_evaluate(lneta, &Gpr_fit, &Gpr_sp, &Gpr_fail); //Add some error checking but OK for now
		return Gpr_fit;
	} // end of interp_Gp*/
}; // end of class bkg








#endif
