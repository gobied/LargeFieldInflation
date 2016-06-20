// Parameters.h
// GO
//Defines a structure containing the parameters and reads the values from a file

#ifndef parameters_h
#define parameters_h

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <nag.h>
#include <nag_stdlib.h>
#include <nagd01.h>
#include <nagx01.h>
#include <nage01.h>
#include <nage02.h>

//class bkg;

struct step_params {
	double m, b, c, d;
	step_params(): m(7.126e-6), b(14.668), c(0.001505), d(0.02705) {}
}; //end of step params

struct bkg_nag_params {
	step_params* sparams;
	std::vector<double>* plna;
	std::vector<double>* pa;
	std::vector<double>* pphi;
	std::vector<double>* pphip;
	std::vector<double>* pH;
	std::vector<double>* pepsH;

	double lna1;
	double lna0;
	Integer istep;
	double lna_step = (lna1 - lna0)/(double(istep)+1.);

	bkg_nag_params(double lna0_in, double lna1_in, Integer istep_in, std::vector<double>* plna_in, std::vector<double>* pa_in, std::vector<double>* pphi_in, std::vector<double>* pphip_in, std::vector<double>* pH_in, std::vector<double>* pepsH_in, step_params* sparams_in): lna0(lna0_in), lna1(lna1_in), istep(istep_in), plna(plna_in), pa(pa_in), pphi(pphi_in), pphip(pphip_in), pH(pH_in), pepsH(pepsH_in), sparams(sparams_in) {}
}; //end of bkg_nag_params

struct muk_params {
        step_params* sparams;
        std::vector<double>* a;
        std::vector<double>* phi;
        std::vector<double>* phip;
        std::vector<double>* H;
	std::vector<double>* eta;
        std::vector<double>* ks; 

        int kindex;

        gsl_interp_accel* acc_phi;
        gsl_interp_accel* acc_phip;
        gsl_interp_accel* acc_H;

        gsl_spline* sp_phi;
        gsl_spline* sp_phip;
        gsl_spline* sp_H;

        muk_params(std::vector<double>* a_in, std::vector<double>* phi_in, std::vector<double>* phip_in, std::vector<double>* H_in, std::vector<double>* eta_in, std::vector<double>* ks_in, step_params* sp): a(a_in), phi(phi_in), phip(phip_in), H(H_in),eta(eta_in), ks(ks_in), sparams(sp)
        {   
                acc_phi = gsl_interp_accel_alloc();
                acc_phip = gsl_interp_accel_alloc();
                acc_H = gsl_interp_accel_alloc();

                sp_phi = gsl_spline_alloc(gsl_interp_cspline, phi->size());
                sp_phip = gsl_spline_alloc(gsl_interp_cspline, phip->size());
                sp_H = gsl_spline_alloc(gsl_interp_cspline, H->size());


                gsl_spline_init(sp_phi, &(*a)[0], &(*phi)[0], phi->size());
                gsl_spline_init(sp_phip, &(*a)[0], &(*phip)[0], phip->size());
                gsl_spline_init(sp_H, &(*a)[0], &(*H)[0], H->size());
        }   
        //fix interpolation so it happens before a is scaled (these value are too small!!!)
	double interp_phi(double aval) { return gsl_spline_eval(sp_phi, aval, acc_phi); }
	double interp_phip(double aval) { return gsl_spline_eval(sp_phip, aval, acc_phip); }
	double interp_H(double aval) { return gsl_spline_eval(sp_H, aval, acc_H); }
	int set_kindex(int newk) { kindex = newk; return 0;} 
};//end of muk_params

struct gsr_params{
	protected:
	std::vector<double>* Gp;
	std::vector<double>* lneta;
	
	double kc = 0.;

        double* lnetar = 0;
	double* etar = 0;
        double* Gpr = 0;
        NagError fail;
	NagError fail_lin;
        Nag_Spline sp;
	Nag_Spline sp_lin;
        double fit = 0.;
	double fit_lin = 0.;
	long double scale;
	public:
	gsr_params(std::vector<double>* Gp_in, std::vector<double>* lneta_in);

	double interp_Gp(double lneta);

	double interp_Gp_lin(double lneta);

	void set_kc(double kc_in) { kc = kc_in; }

	double get_kc() { return kc; }

	double get_scale(){ return scale; }

}; //end of gsr_params

struct eta_params{
	protected:
	std::vector<double>* a;
	std::vector<double>* H;

	double* aHinv; 
	double* pa;

	//NagError fail;
	Nag_Spline sp;
	double fit = 0.;

	public:
	eta_params(std::vector<double>* a_in, std::vector<double>* H_in);

	double interp_aHinv(double lna_in);

}; //end of eta_params

#endif
