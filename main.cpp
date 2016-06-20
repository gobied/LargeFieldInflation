//main.cpp

#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <getopt.h>

#include "parameters.h"
#include "bkg.h"
#include "bkg_nag.h"
#include "muk.h"
#include "approx.h"

//GO
typedef std::chrono::high_resolution_clock Clock;

int main(int argc, char** argv){
	int bkgnsteps = 250000;
	std::cout << "Hello world!" << std::endl;
	/*	
	//NAG
	bkg_nag mybkg_nag;	
	step_params stepparams; // perhaps edit this to read from a file later but OK for now
	mybkg_nag.set_inits(19.9507458217414, -43.0, 150.0, bkgnsteps, stepparams);
	mybkg_nag.bkg_nag_calc();
	*/
	
	//GSL
	bkg mybkg;
	step_params stepparams; 
	if(argc > 1){
		int opt;
		double c_in, d_in;
		while( (opt = getopt(argc, argv, "c:d:")) != -1){
			switch(opt){
				case 'c':
					c_in = std::atof(optarg);
					stepparams.c = c_in;
					break;
				case 'd':
					d_in = std::atof(optarg);
					stepparams.d = d_in;
					break;
			}
		}
		// perhaps edit this to read from a file later but OK for now
	}

	mybkg.set_inits(19.9507458217414, -43.0, 150.0, bkgnsteps, stepparams);
	//mybkg.set_inits(19.9507458217414, -83.0, 150.0, bkgnsteps, stepparams);

	//mybkg.set_inits(15.9507458217414, -43.0, 0.0, bkgnsteps, stepparams);

	mybkg.bkg_calc();

	//Solving Mukhanov
	//Reading k values
	std::vector<double> ks;
	std::fstream ksfile("ks.txt", std::ios_base::in);
	double temp;
	while (ksfile >> temp){
		ks.push_back(temp);
	}
	ksfile.close();

	
	//Setting initial Conditions
	int muknsteps = 250000;

	muk_params mparams(mybkg.geta(), mybkg.getphi(), mybkg.getphip(), mybkg.getH(), mybkg.geteta(), &ks, &stepparams);

	
	//muk mymuk;
	//mymuk.set_inits(muknsteps, &mparams);
	//mymuk.muk_calc();
	//mymuk.writetofile("/n/home13/gobied/LargeFieldInflation_noconv/pows/muk", true); // this assumes you start calculating from the very first k!!

	mybkg.nag_approx_init();
	approx myapprox(&mybkg, &ks);	
	myapprox.approx_calc();
	myapprox.writeGSR("/n/home13/gobied/LargeFieldInflation_noconv/pows/gsr", true); // this assumes you start calculating from the very first k!!
	//myapprox.writeNLA("/n/home13/gobied/LargeFieldInflation_noconv/pows/nla2", true); // this assumes you start calculating from the very first k!!

	return 0;

}
