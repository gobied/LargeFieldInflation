// muk.h
// Defines functions for calculating the mukhanov power spectrum
// GO
#ifndef muk_h
#define muk_h

#include "parameters.h"
#include "utils.h"
#include "constants.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

//GSL Imports 
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

class muk {
    protected:
	int nsteps;
	muk_params* mparams;

	const double u1init = 1.;
	const double u2init = 0.;
	const double u1pinit= 0.;
	const double u2pinit= 1.;

	std::vector<double> Delta2;

	static int modfprime(double lna, const double modf[], double modfp[], void* params);
    public:
	int set_inits(int n, muk_params* p_in){
		nsteps = n;
		mparams = p_in;
		std::cout << "MUK Inits = " << u1init << "  " << u2init << "  " << u1pinit << "  " << u2pinit << std::endl;

		return 0;
	}
	
	void muk_calc();

	std::vector<double>* getDelta2(){ return &Delta2; }

	void writetofile(std::string name, bool suf);

}; // end of class muk

#endif
