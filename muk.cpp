// muk.cpp
// Defines functions for calculating the power spectrum
// GO
#include "muk.h"

int muk::modfprime(double lna, const double modf[], double modfp[], void* params){
	//returns the derivative of modf which is
	// u1, u2, u1p, u2p

	muk_params* p = (muk_params*) params;

        double m = p->sparams->m;
        double b = p->sparams->b;
        double c = p->sparams->c;
        double d = p->sparams->d;

        double k = (*(p->ks))[p->kindex];
        double a = exp(lna);

        double p1 = p->interp_phi(a);
        double p2 = p->interp_phip(a);
        double p3 = p->interp_H(a);

        double Hp = -1.*p3*p2*p2/(2.*Mpl*Mpl);
        double phipp = ( (((p2*p2/(2.*Mpl*Mpl)) - 3.)*p2) - (pot::Vp(p1, *(p->sparams))/(p3*p3)) );

        double X = ((k*k)/(a*a*p3*p3)) - 2. + (4.*Hp*phipp/(p3*p2)) + (2*Hp*Hp/(p3*p3)) + (5.*Hp/p3) + (pot::Vpp(p1, *(p->sparams))/(p3*p3)); // NO MPCCONV!!!

        modfp[0] = modf[2];     //u1p = u1p_in
        modfp[1] = modf[3];     //u2p = u2p_in
        modfp[2] = -1.*(((Hp/p3)+1.)*modf[2]) - (X*modf[0]);
        modfp[3] = -1.*(((Hp/p3)+1.)*modf[3]) - (X*modf[1]);

	return GSL_SUCCESS;
}

void muk::muk_calc(){

	gsl_odeiv2_system muksys = {modfprime, nullptr, 4, mparams};
        gsl_odeiv2_driver* mukd = gsl_odeiv2_driver_alloc_y_new(&muksys, gsl_odeiv2_step_rkf45, 0.01, 0., 1.e-9); // perhaps make the tolerance adjustable later but OK for now

	double kmpc;
	double kc;
	double lna0muk;
	double lna1muk;

	double step_muk;
	double H_i;
	double a_i;
	double lna_i;
	int status;
	double z;
	double conv;

	double modf[4];

	std::vector<double> u1; 
	std::vector<double> u2; 
	std::vector<double> convtest;

	int indexlna0 = 0;
	int indexlna1 = 0;
	double etaminval;
	double etamaxval;

	int asize = mparams->a->size();

	//This is for chekcing only
	//FILE * u1modf;
	//FILE * u2modf;

	//end of checking
	int lenks = mparams->ks->size();

        for(int kin=0; kin <300; kin++){
                mparams->set_kindex(kin);
		// Make these depend on k but OK for now
		kmpc = (*(mparams->ks))[mparams->kindex];
		kc = kmpc;// * mpcconv; NO MPCCONV!!
                // lna0muk = log(kmpc*1.e20);
                // lna1muk = log(kmpc*1.e26);
                etaminval = 0.0001/kc;
		etamaxval = 1000./kc;
		for(int i=indexlna0; i!=asize; ++i){
			if( (*(mparams->eta))[i] < etamaxval ){ // use inheritance here, no need to derference
				indexlna0 = i;
				break;
			}
		}

		for(int i=indexlna1; i!=asize; ++i){
			if( (*(mparams->eta))[i] < etaminval ){
				indexlna1 = i;
				break;
			}
		}
		//Here I did not reset the values of indexlna becuse kc is always increasing but be careful otherwise; OK for now


		lna0muk = log((*(mparams->a))[indexlna0]);
		lna1muk = log((*(mparams->a))[indexlna1]);

		step_muk = (lna1muk-lna0muk)/nsteps; // change this since we changed the a steps near the feature;  OK for now
                H_i = mparams->interp_H(exp(lna0muk));
                a_i = exp(lna0muk); 

                modf[0] = u1init; modf[1] = u2init; modf[2] = u1pinit; modf[3] = u2pinit;

                u1.push_back(u1init);
                u2.push_back(u2init);


	/*	
		double factor;
		if( (mparams->sparams->d) >0.005){
			factor = 3.; 
		} else if ( (mparams->sparams->d) > 0.0005){
			factor = 4.; 
		} else {
			factor = 8.; 
		} // this is a bad solution but OK for now
	*/	
		//checking
		//if((kin > 130) && (kin < 140)){
		//	u1modf = fopen( ("modf/modf1_k" + std::to_string(kc) + ".dat").c_str() , "w");
		//	u2modf = fopen( ("modf/modf2_k" + std::to_string(kc) + ".dat").c_str() , "w");
		//}//end of checking


                for(int i=1; i<=nsteps; ++i){
			lna_i = lna0muk + step_muk;
		
	/* 
 			//Variable step
			if ( ( (mparams->interp_phi(exp(lna_i))) > (14.668 + (factor*(mparams->sparams->d))) ) || ( (mparams->interp_phi(exp(lna_i))) < (14.668 - (5.*(mparams->sparams->d)))) ) {
                        	lna_i = lna0muk + step_muk;
			} else {
				lna_i = lna0muk + (step_muk/30.);
			}
	*/
			
                        status = gsl_odeiv2_driver_apply(mukd, &lna0muk, lna_i, modf);
			

                        u1.push_back(modf[0]);
                        u2.push_back(modf[1]);

                        z = exp(lna0muk) * mparams->interp_phip(exp(lna0muk));
                        //double kc = (*(mparams->ks))[mparams->kindex] * mpcconv;
                        conv = ( (u1.back()*u1.back()/(2.*kc)) + (kc*u2.back()*u2.back()/(2.*a_i*a_i*H_i*H_i)) )/(z*z);

			//checking only
			//if((kin > 130) && (kin < 140)){
			//	fprintf(u1modf, "%.9e  %.9e  %.9e\n", lna0muk, modf[0], z);
			//	fprintf(u2modf, "%.9e  %.9e  %.9e\n", lna0muk, modf[1], z);
			//}
			//end of checking


                        convtest.push_back(conv);
                        //if((std::abs(convtest.back() - convtest.end()[-2]) < (0.00000001 * convtest.back()))){
                        //        break;
			//}
		} // end of modf integration loop

		//checking 
		//if((kin > 130) && (kin < 140)){
		//	std::cout << "kval = " << kc << std::endl;
		//	fclose(u1modf);
		//	fclose(u2modf);
		//}
		//end of checking


                //double kc = (*(mparams->ks))[mparams->kindex] * mpcconv;
                Delta2.push_back( (kc*kc*kc/(2*M_PI*M_PI))*convtest.back() );
                //std::cout << ((*(mparams->ks))[mparams->kindex]) << "  " << Delta2.back() << std::endl;

                std::vector<double>().swap(u1);
                std::vector<double>().swap(u2);
                std::vector<double>().swap(convtest);


	} // end of kindex loop


} // end of muk_calc

void muk::writetofile(std::string name, bool suf){
	std::string suffix = "";
	std::string suffix2 = "";
	if (suf)
		suffix = suffix + "_b" + std::to_string((mparams->sparams->b)) + "_c" + std::to_string((mparams->sparams->c)) + "_d" + std::to_string((mparams->sparams->d)) + ".dat";

	if (!suf)
		suffix2 = suffix2 + ".dat";
	std::ofstream myfile;
	myfile.open(name + suffix + suffix2);

	for(int i=0; i!=Delta2.size(); i++){
		myfile << ((*(mparams->ks))[i]) << "  " << Delta2[i] << std::endl;
	}

	myfile.close();


} // end of writetofile















