// bkg.cpp
// Defines functions for calculating the background phi, phiprime, H
//GO

#include "bkg_nag.h"

int bkg_nag::psiprime(double lna, const double psi[], double psip[], void* params){
	//psi is a vector whose elements are phi, phiprime, H
	(void)(lna);

        psip[0] = psi[1];
        psip[1] = ( (((psi[1]*psi[1]/(2.*Mpl*Mpl)) - 3.)*psi[1]) - (pot::Vp(psi[0], *(step_params*)params)/(psi[2]*psi[2])) );
        psip[2] = -1.*psi[2]*psi[1]*psi[1]/(2.*Mpl*Mpl); 


        return GSL_SUCCESS;
}

void NAG_CALL bkg_nag::bkg_out(Integer neq, double *lnasol, const double psi[], Nag_User *comm){
	bkg_nag_params* s = (bkg_nag_params*)comm->p;

	//Printing results
	//printf("%.9e  %.9e  %.9e  %.9e  %.9e  %.9e\n", *lnasol, exp(*lnasol), psi[0], psi[1], psi[2], psi[1]*psi[1]/2.);

	//Fill arrays
	s->plna->push_back(*lnasol);
	s->pa->push_back(exp(*lnasol));
	s->pphi->push_back(psi[0]);
	s->pphip->push_back(psi[1]);
	s->pH->push_back(psi[2]);
	s->pepsH->push_back(psi[1]*psi[1]/2.);

	//Updat lnasol
	*lnasol = s->lna1 - (double)((s->istep)*(s->lna_step));
	(s->istep)--;
}// end of bkg_out

void NAG_CALL bkg_nag::psiprime_nag(Integer neq, double lna, const double psi[], double psip[], Nag_User *comm){

        psip[0] = psi[1];
        psip[1] = ( (((psi[1]*psi[1]/(2.*Mpl*Mpl)) - 3.)*psi[1]) - (pot::Vp(psi[0], ((bkg_nag_params*)(comm->p))->sparams)/(psi[2]*psi[2])) );
        psip[2] = -1.*psi[2]*psi[1]*psi[1]/(2.*Mpl*Mpl); 

} // end of psiprime-nag

double NAG_CALL bkg_nag::epsHm1(Integer neq, double lna, const double phi[], Nag_User *comm){
	return (1.-((phi[1]*phi[1])/(2.*Mpl*Mpl)));

} //end of epsHm1

void bkg_nag::bkg_nag_calc(){
	printf("Calling bkg_nag_calc\n");
	Integer neq = 3;
	double psi[3] = {phi_init, phip_init, H_init};
	double lnatemp = lna0;
	Nag_User comm;
	NagError fail;

	bkg_nag_params nagparams(lna0, lna1, nsteps, &lna, &a, &phi, &phip, &H, &epsH);

	SET_FAIL(fail);
	comm.p = (Pointer) &nagparams;

	nag_ode_ivp_adams_gen(neq, psiprime_nag, &lnatemp, psi, nagparams.lna1, 1.e-13, Nag_Mixed, bkg_out, epsHm1, &comm, &fail);

	printf("Number of efolds = %.9e\n", lnatemp);

	//matching
	

	int index50ef;
	if(lnatemp < 50.){
		printf("Error, not enough efolds\n");
	} else {
		for(int i=0; i<lna.size(); i++){
			if((lna.back()-lna[i]) <= 50.){
				index50ef = i;
				break;
			}   
		}   
	}  
	double matching = log(0.05) - log(H[index50ef]) + 50. - lna.back();// + log(mpcconv);
	std::transform(lna.begin(), lna.end(), lna.begin(), std::bind2nd(std::plus<double>(), matching));
	std::transform(a.begin(), a.end(), a.begin(), std::bind2nd(std::multiplies<double>(), exp(matching)));

} // end of bkg_calc

void bkg_nag::approx_init(){
	//Calculate eta, f, fp, fpp, G, and Gp for GSR and NL
	
	//Calculate eta, fix this to know what to do based on the size of lna but OK for now
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
	//ait += 1;
	//Hit += 1;

	for(ait; ait!=a.rend(); ait+=2, Hit+=2){
		//std::cout << ait[0] << "   " << ait[1] << "   " << ait[2] << std::endl;
		//std::cout << Hit[0] << "   " << Hit[1] << "   " << Hit[2] << std:: endl;
		eta.push_back(eta.back() + (h*((1./(ait[0]*Hit[0])) + (1./(ait[1]*Hit[1])))/2.));
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

                printf("%.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e\n", log(a[i]), a[i], phi[i], phip[i], H[i], eta[i], f[i], dfdlneta[i], ddfdlneta[i], G[i], Gp[i]); 
                fflush(stdout);

	
		// df should be O(100) and ddf should be O(1)
	} // end of f, fp, fpp, G, Gp loop
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















