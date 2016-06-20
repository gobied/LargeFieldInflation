// bkg.cpp
// Defines functions for calculating the background phi, phiprime, H


#include "bkg.h"

int bkg::psiprime(double lna, const double psi[], double psip[], void* params){
	//psi is a vector whose elements are phi, phiprime, H
	(void)(lna);

        psip[0] = psi[1];
        psip[1] = ( (((psi[1]*psi[1]/(2.*Mpl*Mpl)) - 3.)*psi[1]) - (pot::Vp(psi[0], *(step_params*)params)/(psi[2]*psi[2])) );
        psip[2] = -1.*psi[2]*psi[1]*psi[1]/(2.*Mpl*Mpl); 

        return GSL_SUCCESS;
}

void bkg::bkg_calc(){
	printf("Calling bkg_calc\n");
	//step_params p;

	std::cout << "in BKG_CALC " << p.c << std::endl;

        gsl_odeiv2_system sys = {psiprime, nullptr, 3, &p};

        gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 0.001, 0., 1.e-13);

	double lna_step = (lna1 - lna0)/nsteps;
	double psi[3] = {phi_init, phip_init, H_init};
	double lna_i;
	int status;
	bool fillin1 = true;
	bool fillin2 = false;
	double factor;
	if(p.d>0.005){
		factor = 3.;
	} else if (p.d > 0.0005){
		factor = 4.;
	} else {
		factor = 8.;
	} // this is a bad solution but OK for now

	//Add initial conditions to bkg arrays
	lna.push_back(lna0);
	a.push_back(exp(lna0));
	lna_matched.push_back(lna0);
	a_matched.push_back(exp(lna0));
	phi.push_back(phi_init);
	phip.push_back(phip_init);
	H.push_back(H_init);
	epsH.push_back(phip_init*phip_init/2.);
	/*
	etaH.push_back( -1.*( ( (((phip_init*phip_init/2.) - 3.)*phip_init) - (pot::Vp(phi_init, p)/(H_init*H_init))) - (phip_init*phip_init*phip_init/2.) )/phip_init);
	delta2.push_back(( (pow(phip_init, 5.)/4.) 
		+ (( (pot::Vp(phi_init, p)/(H_init*H_init)) - (phip_init*((phip_init*phip_init/4.) - 3.)) )*phip_init*phip_init) 
		- ((3.*phip_init*phip_init/2.) * ( (((phip_init*phip_init/2.)-3.)*phip_init) - (pot::Vp(phi_init, p)/(H_init*H_init)) )) 
		+ ((3.*((phip_init*phip_init/2.)-1.)*( (((phip_init*phip_init/2.)-3.)*phip_init) - (pot::Vp(phi_init, p)/(H_init*H_init)) )) - (phip_init*phip_init*pot::Vp(phi_init, p)/(H_init*H_init)) - (phip_init*pot::Vpp(phi_init, p)/(H_init*H_init))) )
		/phip_init);
	*/
	//printf("%.9e  %.9e  %.9e  %.9e  %.9e  %.9e  %.9e  %.9e\n", lna[0], exp(lna[0]), phi[0], phip[0], H_init, phip[0]*phip[0]/2., etaH.back(), delta2.back());

	//Evolve initial conditions and fill bkg arrays until inflation ends
	for(int i=1; i<=nsteps; i++){	
		//std::cout <<(14.668 + (3.*p.d)) << "  " << (14.668 - (3.*p.d)) << "  " << phi.back() << std::endl;
		if ( phi.back() > (14.668 + (factor*p.d)) || phi.back() < (14.668 - (5.*p.d)) ){
			if(fillin2){
				indexetalo = i;
				fillin2 = false;
			}
			lna_i = lna0+lna_step;
		} else {
			if(fillin1){
				indexetahi = i;
				fillin1 = false;
				fillin2 = true;
			}
			lna_i = lna0 + (lna_step/30.); //make this factor a function of c and d (likely c/d^2) but OK for now
		}
		status = gsl_odeiv2_driver_apply(d, &lna0, lna_i, psi);

		if(status!=GSL_SUCCESS){
			printf("Error, bkg integrator return value=%d\n", status);
			break;
		}

		lna.push_back(lna0);
		a.push_back(exp(lna0));
		lna_matched.push_back(lna0);
		a_matched.push_back(exp(lna0));
		phi.push_back(psi[0]);
		phip.push_back(psi[1]);	
		H.push_back(psi[2]);
		epsH.push_back(psi[1]*psi[1]/(2.*Mpl*Mpl));
		/*
		etaH.push_back( -1.*( ( (((psi[1]*psi[1]/2.) - 3.)*psi[1]) - (pot::Vp(psi[0], p)/(psi[2]*psi[2]))) 
				- (psi[1]*psi[1]*psi[1]/2.) )
				/psi[1]);
		delta2.push_back(( (pow(psi[1], 5.)/4.) 
			+ (( (pot::Vp(psi[0], p)/(psi[2]*psi[2])) - (psi[1]*((psi[1]*psi[1]/4.) - 3.)) )*psi[1]*psi[1]) 
			- ((3.*psi[1]*psi[1]/2.) * ( (((psi[1]*psi[1]/2.)-3.)*psi[1]) - (pot::Vp(psi[0], p)/(psi[2]*psi[2])) )) 
			+ ((3.*((psi[1]*psi[1]/2.)-1.)*( (((psi[1]*psi[1]/2.)-3.)*psi[1]) - (pot::Vp(psi[0], p)/(psi[2]*psi[2])) )) - (psi[1]*psi[1]*pot::Vp(psi[0], p)/(psi[2]*psi[2])) - (psi[1]*pot::Vpp(psi[0], p)/(psi[2]*psi[2]))) )
			/psi[1]);
		*/

		//printf("%.9e  %.9e  %.9e  %.9e  %.9e  %.9e  %.9e  %.9e\n", lna0, exp(lna0), psi[0], psi[1], psi[2], epsH.back(), etaH.back(), delta2.back());


		if((phip.back()*phip.back()/(2.*Mpl*Mpl)) >= 1.){
			double Nefolds = lna.back() - lna[0]; //improve this by some interpolation but OK for now
			int index50ef;
			if(Nefolds < 50.){
				printf("Error, not enough efolds\n");
			} else {
				for(int i=0; i<lna.size(); i++){
					if((lna.back()-lna[i]) <= 50.){
						index50ef = i;
						break;
					}
				}
			}
			matching = log(0.05) - log(H[index50ef]) + 50. - lna.back();// + log(mpcconv);
			printf("aH 50 ef = %.9e\n", exp(lna[index50ef]) * H[index50ef]);
			printf("Nefolds before end = %.9e\n", (lna.back() - lna[index50ef]));
			printf("logH50ef, lnaend = %.9e  %.9e\n", log(H[index50ef]), lna.back());
			printf("matching value = %.9e   %.9e\n", exp(matching), matching);	
			printf("H_pivot = %.5e\n", H[index50ef]);	
			printf("alpha_e = %.5e %.5e\n", lna.back(), lna.rbegin()[1]);	
			printf("conversion factor = %.5e\n", mpcconv);	




			printf("lna0 and last = %.5e  %.5e\n", lna[0], lna.back());
			//Consider saving the matching condition to prevent dealing with large numbers
			std::transform(lna_matched.begin(), lna_matched.end(), lna_matched.begin(), std::bind2nd(std::plus<double>(), matching));
			std::transform(a_matched.begin(), a_matched.end(), a_matched.begin(), std::bind2nd(std::multiplies<double>(), exp(matching)));

			std::transform(lna.begin(), lna.end(), lna.begin(), std::bind2nd(std::plus<double>(), matching));
			std::transform(a.begin(), a.end(), a.begin(), std::bind2nd(std::multiplies<double>(), exp(matching)));

			printf("aH 50 ef = %.9e\n", exp(lna[index50ef]) * H[index50ef]);
			printf("Nefolds before end = %.9e\n", (lna.back() - lna[index50ef]));
			printf("logH50ef, lnaend = %.9e  %.9e\n", log(H[index50ef]), lna.back());

			printf("lna0 and last = %.5e  %.5e\n", lna[0], lna.back());

			printf("Inflation ended with %.5e efolds\n", Nefolds);
			
			break;
		} // end of inflation and matching		
		//printf("%.5e %.5e %.5e %.5e\n", lna0, psi[0], psi[1], psi[2]);
	} // end of integration loop

	//Calculating eta
	Nag_Comm Icomm;
	
	double etatemp;
	double etaerr;

	eta_params etap(&a, &H);
	double ruser[1] = {-1.0};
	Icomm.user = ruser;
	Icomm.p = (Pointer) &etap;

	//eta.push_back(0.); lneta.push_back(-1.*std::numeric_limits<double>::infinity());
	// will figure something out with these but OK for now

	eta.push_back(0.);
	lneta.push_back(-1.*std::numeric_limits<double>::infinity());
	for(int i=1; i<a.size();i++){
		//nag_1d_quad_vals((Integer) i, &lna[a.size()-i], &etaintegrand[a.size()-i], &etatemp, &etaerr, &etafail);
		nag_quad_1d_fin_smooth(etaintegrand, lna[lna.size() - 1 - i], lna[lna.size()-1] , 1.e-13, 1.e-13, &etatemp, &etaerr, &Icomm);
		eta.push_back(etatemp);
		lneta.push_back(log(etatemp));

		//printf("%.9e  %.9e  %.9e  %.9e  %.9e  %.9e  %.9e  %.9e  %.9e\n", lna[i], exp(lna[i]), lneta.back(), phi[i], phip[i], H[i], epsH[i], etaH[i], delta2[i]);

		//printf("%.15e  %.15e  %.15e  %.15e\n",lna[i], eta.back(), lneta.back(), etaerr);
	}
	std::reverse(eta.begin(), eta.end());
	std::reverse(lneta.begin(), lneta.end());

	//maybe recalculate refined bkg here

} // end of bkg_calc

void bkg::approx_init(){
	//Calculate eta, f, fp, fpp, G, and Gp for GSR and NL

	//Calculate eta, fix this to know what to do based on the size of lna but OK for now
	std::cout << "SIZES. a = " << a.size() << " H = " << H.size() << std::endl;
	std::vector<double>::reverse_iterator ait = a.rbegin();
	std::vector<double>::reverse_iterator Hit = H.rbegin();
	double h = lna[1] - lna[0];


	eta.push_back(0.);
	lneta.push_back(-1.*std::numeric_limits<double>::infinity());	
	//std::cout << lneta.rbegin()[0] << std::endl;
	eta.push_back(h*((1./(ait[0]*Hit[0])) + (1./(ait[1]*Hit[1])))/2.);
	lneta.push_back(log(eta.rbegin()[0]));
	//std::cout << lneta.rbegin()[0] << std::endl;

	//eta.push_back(h*((1./(ait[0]*Hit[0])) + 4.*(1./(ait[1]*Hit[1])) + (1./(ait[2]*Hit[2])))/3.); //Simpson's rule
	ait += 1;
	Hit += 1;

	for(ait; ait!=a.rend(); ait+=2, Hit+=2){
		//std::cout << ait[0] << "   " << ait[1] << "   " << ait[2] << std::endl;
		//std::cout << Hit[0] << "   " << Hit[1] << "   " << Hit[2] << std:: endl;
		eta.push_back(eta.back() + (h*((1./(ait[0]*Hit[0])) + (1./(ait[1]*Hit[1])))/2.));
		//std::cout <<"COMPLETED ONCE" << std::endl;

		lneta.push_back(log(eta.rbegin()[0]));
		//std::cout << lneta.back() << std::endl;
		eta.push_back(eta.rbegin()[1] + (h*((1./(ait[0]*Hit[0])) + 4.*(1./(ait[1]*Hit[1])) + (1./(ait[2]*Hit[2])))/3.));
		lneta.push_back(log(eta.rbegin()[0]));
		//std::cout << lneta.rbegin()[0] << std::endl;

	}
	std::reverse(eta.begin(), eta.end());
	std::reverse(lneta.begin(), lneta.end());


	//std::cout << "THIS IS ETA " << eta[1000] << "   " << eta[eta.size() - 25] << "   " << eta.size() << "  " << std::endl;
	//Calculate f, dfdlneta
	for(int i=0; i!=eta.size(); i++){
		double Hp = -1.*H[i]*phip[i]*phip[i]/(2.*Mpl*Mpl);
		double phipp = ( (((phip[i]*phip[i]/(2.*Mpl*Mpl)) - 3.)*phip[i]) - (pot::Vp(phi[i], p)/(H[i]*H[i])) );
		double phippp = (3.*phipp*((phip[i]*phip[i]/2.) - 1.)) - (phip[i]*phip[i]*pot::Vp(phi[i], p)/(H[i]*H[i])) - (phip[i]*pot::Vpp(phi[i],p)/(H[i]*H[i]));

		f.push_back(2.*M_PI*phip[i]*a[i]*eta[i]);
		dfdlneta.push_back(a[i]*H[i]*eta[i]*((2.*M_PI*phip[i]/H[i]) - (f[i]) - (2.*M_PI*phipp*a[i]*eta[i])));
		ddfdlneta.push_back(2.*M_PI*a[i]*a[i]*eta[i]*eta[i]*H[i]* 
			((phip[i]/(a[i]*H[i]*eta[i])) +
			((phippp + (3.*phipp) + (2.*phip[i]))*a[i]*H[i]*eta[i]) +
			((phipp + phip[i])*((a[i]*Hp*eta[i]) - 3.))));


		G.push_back(log(1./(f[i]*f[i])) + (2.*dfdlneta[i]/(3.*f[i])));
		Gp.push_back( 2.*((f[i]*ddfdlneta[i]) - (3.*f[i]*dfdlneta[i]) - (dfdlneta[i]*dfdlneta[i]))/(3.*f[i]*f[i]) );

		//Numerical fp, fpp, Gp
		//df_num.push_back(-1.*a[i]*H[i]*eta[i]*);
		//ddf_num.push_back();
		
		printf("%.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e\n", log(a[i]), a[i], phi[i], phip[i], H[i], eta[i], f[i], dfdlneta[i], ddfdlneta[i], G[i], Gp[i]);
		fflush(stdout);
		// df should be O(100) and ddf should be O(1)
	} // end of f, fp, fpp, G, Gp loop
	

	//lna   a   phi   phip       H                         eta   f   dfdlneta   ddfdlneta           G   Gp
	/*
	lnetar = NAG_ALLOC(lneta.size() - 4, double);
	Gpr = NAG_ALLOC(lneta.size() - 4, double);

	
	for(int i=0; i<(lneta.size() - 4); i++){
		Gpr[i] = Gp.rbegin()[2+i];
		lnetar[i] = lneta.rbegin()[2+i];
	} //end of for loop to fill Gpr and lnetar

	nag_1d_spline_interpolant((Integer)lneta.size() - 4, lnetar, Gpr, &Gpr_sp, &Gpr_fail);
	//Add some form of error checking but OK for now

	std::cout << "TESTIN INTERP " <<  interp_Gp((lneta[1000] + lneta[1001])/2.) << std::endl;
	*/
} // end of approx_init

void bkg::nag_approx_init(){
	std::cout << "SIZES. a = " << a.size() << " H = " << H.size() << std::endl;

	/*
	Nag_Comm Icomm;
	
	double etatemp;
	double etaerr;

	eta_params etap(&a, &H);
	static double ruser[1] = {-1.0};
	Icomm.user = ruser;
	Icomm.p = (Pointer) &etap;

	//eta.push_back(0.); lneta.push_back(-1.*std::numeric_limits<double>::infinity());
	// will figure something out with these but OK for now

	for(int i=0; i<a.size();i++){
		//nag_1d_quad_vals((Integer) i, &lna[a.size()-i], &etaintegrand[a.size()-i], &etatemp, &etaerr, &etafail);
		nag_quad_1d_fin_smooth(etaintegrand, lna[lna.size() - 1 - i], lna[lna.size()-1] , 1.e-13, 1.e-13, &etatemp, &etaerr, &Icomm);
		eta.push_back(etatemp);
		lneta.push_back(log(etatemp));

		//printf("%.15e  %.15e  %.15e  %.15e\n",lna[i], eta.back(), lneta.back(), etaerr);
	}
	std::reverse(eta.begin(), eta.end());
	std::reverse(lneta.begin(), lneta.end());
	*/

	double Hp, phipp, phippp;
	double iphi, iphip, ia, ieta, iH;
	double ft, dft, ddft;


	//for reconstruction test
	/*
	std::fstream Gpfile("pows/Gp_c0.001505_d0.02705_rate10.dat", std::ios_base::in);
	std::string Gptemp;
	std::string::size_type sz;
	while (Gpfile >> Gptemp){
		if (Gptemp == "-inf") {
			Gp.push_back(-1.*std::numeric_limits<double>::infinity());
		} else if (Gptemp == "nan"){
			Gp.push_back(std::nan(""));
		} else {
			Gp.push_back(std::stod(Gptemp, &sz));
		}
	}
	Gpfile.close();
	std::cout << "size of Gprime " << Gp.size() << std::endl;
	*/
	for(int i=0; i!=eta.size(); i++){
		iphi = phi[i];
		iphip = phip[i];
		ia = a[i];
		ieta = eta[i];
		iH = H[i];
		Hp = -1.*iH*iphip*iphip/(2.*Mpl*Mpl);
		phipp = ( (((iphip*iphip/(2.*Mpl*Mpl)) - 3.)*iphip) - (pot::Vp(iphi, p)/(iH*iH)) );
		phippp = (3.*phipp*((iphip*iphip/2.) - 1.)) - (iphip*iphip*pot::Vp(iphi, p)/(iH*iH)) - (iphip*pot::Vpp(iphi,p)/(iH*iH));


		ft = 2.*M_PI*iphip*ia*ieta;
		dft = ia*iH*ieta*((2.*M_PI*iphip/iH) - (ft) - (2.*M_PI*phipp*ia*ieta));
		ddft = 2.*M_PI*ia*ia*ieta*ieta*iH* 
                        ((iphip/(ia*iH*ieta)) +
                        ((phippp + (3.*phipp) + (2.*iphip))*ia*iH*ieta) +
                        ((phipp + iphip)*((ia*Hp*ieta) - 3.))); 



		f.push_back(ft);
		dfdlneta.push_back(dft);
		ddfdlneta.push_back(ddft);


		G.push_back(log(1./(ft*ft)) + (2.*dft/(3.*ft)));
		Gp.push_back( 2.*((ft*ddft) - (3.*ft*dft) - (dft*dft))/(3.*ft*ft) );

		//For reconstruction test comment the above line
		

		//Numerical fp, fpp, Gp
		//df_num.push_back(-1.*a[i]*H[i]*eta[i]*);
		//ddf_num.push_back();

		//printf("%.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e\n", log(a[i]), a[i], phi[i], phip[i], H[i], eta[i], f[i], dfdlneta[i], ddfdlneta[i], G[i], Gp[i]);
		//fflush(stdout);
		// df should be O(100) and ddf should be O(1)
	} // end of f, fp, fpp, G, Gp loop

	std::cout << "size of eta " << eta.size() << std::endl;

} //end of nag_approx_init()














