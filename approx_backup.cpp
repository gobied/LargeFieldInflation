// gsr.cpp
// Calculates Gmin, I0, I1 and combines them to get the power spectrum
//GO
#include "approx.h"

void approx::nag_qp_free(Nag_QuadProgress& qp){	
	NAG_FREE(qp.sub_int_beg_pts);
	NAG_FREE(qp.sub_int_end_pts);
	NAG_FREE(qp.sub_int_result);
	NAG_FREE(qp.sub_int_error);
}

void approx::approx_calc(){
	// Calculate indices first
	// Perhaps make these k-dependent? but OK for now
	//double etaminval = 1.e52; // at most
	//double etamaxval = 1.e66; // at leas
	//double etalow = 9.e58; // where step in Gp starts
	//double etahi = 4.e59; // where step in Gp ends

	//no matching values
	double etaminval;// = 1.;//0.0001; //1.e-20; //at most
	double etamaxval;// = 1.e10;//1.e-9; //at least
	double Gmin;
	std::cout << "INDICES " << lneta->size() << "  " << mybkg->indexetalo << "   " << mybkg->indexetahi << std::endl;
	double etalo = exp((*lneta)[mybkg->indexetalo]); //180.; // make a function to calculate this but OK for now
	double etahi = exp((*lneta)[mybkg->indexetahi]);//1200.; 
	std::cout << "etalo and etahi " << mybkg->indexetalo << "   " << etalo << "  " << mybkg->indexetahi << "   " << etahi <<  "   " << lneta->size() << std::endl;

	//Set Gmin
	//Gmin = (*G)[etaminindex];
	//std::cout << "THIS IS Gmin " << (*G)[etaminindex -1] << "   " << (*G)[etaminindex] << "   " << (*G)[etaminindex+1] << std::endl;

	double kmpc;
	double kc;
	//double h = log((*a)[1]) - log((*a)[0]);
	double I0temp;
	double I1temp;
	double I0int;
	double I1int;
	double I0err;
	double I1err;
	Nag_QuadProgress I0qp;
	Nag_QuadProgress I1qp;
	Nag_User Icomm_osc;
	Nag_Comm Icomm;
	static double ruser[1] = {-1.0};
	Icomm.user = ruser;
	NagError I0fail;
	NagError I1fail;
	gsr_params p(Gp, lneta);
	double iscale = p.get_scale();

	printf("iscale = %g\n",iscale);

	//find the peak
	int pmaxindex = std::max_element(Gp->begin()+mybkg->indexetahi, Gp->begin()+mybkg->indexetalo) - (Gp->begin() + mybkg->indexetahi);
	int pminindex = std::min_element(Gp->begin()+mybkg->indexetahi, Gp->begin()+mybkg->indexetalo) - (Gp->begin() + mybkg->indexetahi);
	double Gpmax = *(Gp->begin() + mybkg->indexetahi + pmaxindex);
	double Gpmin = *(Gp->begin() + mybkg->indexetahi + pminindex);

	std::cout << "Peak indices " << pmaxindex << "   " << pminindex << "   " << (*(Gp->begin() + mybkg->indexetahi + pmaxindex)) << "   " << Gpmax << "   " << (*(Gp->begin() + mybkg->indexetahi + pminindex)) << "   " << Gpmin << std::endl;

	int fineindexhi = -1;
	int fineindexlo = -1;

	for(int i = pmaxindex; i > 0; i--){
		if( (*(Gp->begin() + mybkg->indexetahi + i)) < (Gpmax/100.) ){
			fineindexhi = mybkg->indexetahi + i;
			break;
		}
	}

	for(int i = pminindex; i + mybkg->indexetahi < mybkg->indexetalo; i++){
		if( (*(Gp->begin() + mybkg->indexetahi + i)) > (Gpmin/100.) ){
			fineindexlo = mybkg->indexetahi + i;
			break;
		}
	}

	if (fineindexhi == -1)
		fineindexhi = mybkg->indexetahi;
	if (fineindexlo == -1)
		fineindexlo = mybkg->indexetalo;

	std::cout << "Fine indices " << fineindexlo << "   " << fineindexhi << "   " << ((*(Gp))[fineindexlo]) << "   " << ((*(Gp))[fineindexhi]) << "   " << exp((*lneta)[fineindexlo]) << "   " << exp((*lneta)[fineindexhi]) <<  std::endl;

	//Divide the peak into this many intervals
	int nosc = 1.;//40; // make this depend on c and d but OK for now
	double fineetalo = exp((*lneta)[fineindexlo]);
	double fineetahi = exp((*lneta)[fineindexhi]);
	double oscint = (fineetahi - fineetalo)/nosc;

	SET_FAIL(I0fail);
	SET_FAIL(I1fail);

		
	//NAG library implementation	
	for(Integer kin=0; kin<300; kin++){ //200 and above is case I
		kmpc = (*ks)[kin];
		kc = kmpc;// * mpcconv; no MPCCONV!! change the names.. but OK for now
		//lnetainterval = log(4.*M_PI/(2.*kc)); //4. is a guess but OK for now
		I0temp = 0.;
		I1temp = 0.;
		p.set_kc(kc);

		Icomm.p = (Pointer)&p;
		Icomm_osc.p = (Pointer)&p;

		//Calculate etaminval, etamaxval and Gmin
		etaminval = 0.001/kc;
		etamaxval = 1000./kc;
		Gmin = interp_G(log(etaminval));

/*
	//this is just for checking
        std::ofstream myfile_G;
	std::ofstream myfile_gsr;
        myfile_G.open("limits_Gpeak.dat");
	myfile_gsr.open("limits.dat");
	FILE * gpeakf;
	FILE * gsrerrf; 
	gpeakf = fopen("Gpeak.dat", "w");
	gsrerrf = fopen("trialgsr_err.dat", "w");
	
	int etamaxindex, etaminindex;
	double lneta_in;
	double interval = (etahi-etalo)/5000.;
	for(int i=0; i!=eta->size(); i++){
		if((*eta)[i] < etamaxval){
			etamaxindex = i;
			break;
		}
	}
	for(int i=etamaxindex; i!=eta->size(); i++){
		if((*eta)[i] < etaminval){
			etaminindex = i;
			break;
		}
	}
	for(int i=etamaxindex; i<etaminindex; i=i+1){
		//break;
		p.set_kc((*ks)[0]); // NO MPCCONV!!!
		lneta_in = (*lneta)[i];
		//std::cout << "STUPID TRIAL " << p.get_kc() << "   " << exp(lneta_in) << "   " << std::endl;
		fprintf(gsrerrf, "%.9e  %.9e  %.9e  %.9e  %.9e\n", lneta_in, (p.get_kc()*exp(lneta_in)), I0integrand_osc(lneta_in, &Icomm_osc), win::W((p.get_kc()) *exp(lneta_in)), p.interp_Gp(lneta_in));
		fflush(stdout);
	}

	for(int i=0; i<5000; i++){
		//break;
		p.set_kc((*ks)[0]);
		fprintf(gpeakf, "%.9e  %.9e\n", (log(etalo+(i*interval))),p.interp_Gp(log(etalo+(i*interval))));
	}
	std::cout << "Gp at etalo and etahi " << p.interp_Gp(log(etalo)) << "  " << p.interp_Gp(log(etahi)) << std::endl;

	for(int osci=0; osci<nosc; osci++){
		myfile_G << ((etalo + (osci*oscint))) << std::endl;
	}

	for(int osci=0; osci<(mybkg->indexetalo - mybkg->indexetahi); osci++){
		myfile_gsr << exp((*lneta)[(mybkg->indexetalo)-osci]) << std::endl;
	}


	myfile_G.close();
	myfile_gsr.close();
	fclose(gpeakf);
	fclose(gsrerrf);

	//end of chekcing code
*/	
		if (fineetalo > (1./kc)){
			//Smooth integral
			nag_quad_1d_fin_smooth(I0integrand, log(etaminval), log(1./kc), 1.e-8, 1.e-8, &I0int, &I0err, &Icomm);
			nag_quad_1d_fin_smooth(I1integrand, log(etaminval), log(1./kc), 1.e-8, 1.e-8, &I1int, &I1err, &Icomm);

			I0temp = I0temp + I0int;
			I1temp = I1temp + I1int;
			//std::cout << "Case 11 " << I0int << "(+/-)" <<  I0err << "   " << I1int<< "(+/-)" << I1err << std::endl;

			//Do oscillatory integral
			nag_1d_quad_osc_1(I0integrand_osc, log(1./kc), log(fineetalo), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
			nag_1d_quad_osc_1(I1integrand_osc, log(1./kc), log(fineetalo), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

			I0temp = I0temp + I0int;
			I1temp = I1temp + I1int;

			//std::cout << "Case 12 " << I0int << "(+/-)" <<  I0err << "   " << I1int << "(+/-)" << I1err << std::endl;

			nag_qp_free(I0qp);
			nag_qp_free(I1qp);

			//Fine integral
			for(int osci=0; osci<nosc; osci++){
				nag_1d_quad_osc_1(I0integrand_osc2, fineetalo + (osci*oscint), fineetalo + ((osci+1)*oscint), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
				nag_1d_quad_osc_1(I1integrand_osc2, fineetalo + (osci*oscint), fineetalo + ((osci+1)*oscint), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

				I0temp = I0temp + I0int;
				I1temp = I1temp + I1int;

				//std::cout << "Case 13 " << I0int << "(+/-)" <<  I0err << "   " << I1int << "(+/-)" << I1err << std::endl;

				nag_qp_free(I0qp);
				nag_qp_free(I1qp);
			}

		       	nag_1d_quad_osc_1(I0integrand_osc, log(fineetahi), log(etamaxval), 1.e-8, 1.e-8, 1.e7, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
                        nag_1d_quad_osc_1(I1integrand_osc, log(fineetahi), log(etamaxval), 1.e-8, 1.e-8, 1.e7, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);
			//std::cout << "Case 14 " << I0int << "(+/-)" <<  I0err << "   " << I1int << "(+/-)" << I1err << std::endl;

                        I0temp = I0temp + I0int;
                        I1temp = I1temp + I1int;
                        nag_qp_free(I0qp);
                        nag_qp_free(I1qp);
		

		} else if ( (fineetalo < (1./kc)) && (fineetahi > (1./kc)) ){
			//Smooth integral
			nag_quad_1d_fin_smooth(I0integrand, log(etaminval), log(fineetalo), 1.e-8, 1.e-8, &I0int, &I0err, &Icomm);
			nag_quad_1d_fin_smooth(I1integrand, log(etaminval), log(fineetalo), 1.e-8, 1.e-8, &I1int, &I1err, &Icomm);

			I0temp = I0temp + I0int;
			I1temp = I1temp + I1int;
			//std::cout << "Case 21 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;

			//Do oscillatory integral
			for(int osci=0; osci<nosc; osci++){
				nag_1d_quad_osc_1(I0integrand_osc2, fineetalo + (osci*oscint), fineetalo + ((osci+1)*oscint), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
				nag_1d_quad_osc_1(I1integrand_osc2, fineetalo + (osci*oscint), fineetalo + ((osci+1)*oscint), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

				I0temp = I0temp + I0int;
				I1temp = I1temp + I1int;
				//std::cout << "Case 22 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;

				nag_qp_free(I0qp);
				nag_qp_free(I1qp);
			}

                        nag_1d_quad_osc_1(I0integrand_osc, log(fineetahi), log(etamaxval), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
                        nag_1d_quad_osc_1(I1integrand_osc, log(fineetahi), log(etamaxval), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

                        I0temp = I0temp + I0int;
                        I1temp = I1temp + I1int;
                        nag_qp_free(I0qp);
                        nag_qp_free(I1qp);

		} else if (fineetahi < (1./kc)) {
			//Smooth integral (until the fine division)
			nag_quad_1d_fin_smooth(I0integrand, log(etaminval), (*lneta)[fineindexlo], 1.e-8, 1.e-8, &I0int, &I0err, &Icomm);
			nag_quad_1d_fin_smooth(I1integrand, log(etaminval), (*lneta)[fineindexlo], 1.e-8, 1.e-8, &I1int, &I1err, &Icomm);

			I0temp = I0temp + I0int;
			I1temp = I1temp + I1int;
			//std::cout << "Case 31 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;
			//Do oscillatory integral
			//Old way
	
			for(int osci=0; osci<nosc; osci++){
				nag_1d_quad_gen_1(I0integrand_gen, (fineetalo + (osci*oscint)), (fineetalo + ((osci+1)*oscint)), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
				nag_1d_quad_gen_1(I1integrand_gen, (fineetalo + (osci*oscint)), (fineetalo + ((osci+1)*oscint)), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

				I0temp = I0temp + I0int;
				I1temp = I1temp + I1int;
			//	std::cout << "Case 31.5 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;
				//printf("%.9e  %.9e\n", (etalo + (osci*oscint)), (etalo + ((osci+1)*oscint)));
				nag_qp_free(I0qp);
				nag_qp_free(I1qp);
			}
			//std::cout << "checking limits " << (fineetalo + ((nosc)*oscint)) << "  " << fineetahi << std::endl;
	
			//new way
/*		
			
			//nag_1d_quad_gen_1(I0integrand_gen, etalo, exp((*lneta)[mybkg->indexetalo]), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
			//nag_1d_quad_gen_1(I1integrand_gen, etalo, exp((*lneta)[mybkg->indexetalo]), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

			//I0temp = I0temp + I0int;
			//I1temp = I1temp + I1int;
			//std::cout << "Case 31.5 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;
			//nag_qp_free(I0qp);
			//nag_qp_free(I1qp);

			//Fine-division integral
			for(int osci=0; osci<(mybkg->indexetalo - mybkg->indexetahi); osci++){
				nag_1d_quad_gen_1(I0integrand_gen, exp((*lneta)[(mybkg->indexetalo)-osci]), exp((*lneta)[(mybkg->indexetalo) - osci - 1]), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
				nag_1d_quad_gen_1(I1integrand_gen, exp((*lneta)[(mybkg->indexetalo)+osci]), exp((*lneta)[(mybkg->indexetalo) - osci - 1]), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

				I0temp = I0temp + I0int;
				I1temp = I1temp + I1int;
				std::cout << "Case 32 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;

				nag_qp_free(I0qp);
				nag_qp_free(I1qp);
			}

			//nag_1d_quad_gen_1(I0integrand_gen, exp((*lneta)[mybkg->indexetahi]), etahi, 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
			//nag_1d_quad_gen_1(I1integrand_gen, exp((*lneta)[mybkg->indexetahi]), etahi, 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

			//I0temp = I0temp + I0int;
			//I1temp = I1temp + I1int;
			//std::cout << "Case 32.5 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;
			//nag_qp_free(I0qp);
			//nag_qp_free(I1qp);
*/			
	
			//Smooth integral
			nag_quad_1d_fin_smooth(I0integrand, (*lneta)[fineindexhi], log(1./kc), 1.e-8, 1.e-8, &I0int, &I0err, &Icomm);
			nag_quad_1d_fin_smooth(I1integrand, (*lneta)[fineindexhi], log(1./kc), 1.e-8, 1.e-8, &I1int, &I1err, &Icomm);

			I0temp = I0temp + I0int;
			I1temp = I1temp + I1int;
			//std::cout << "Case 33 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;

			//Do oscillatory integral
			nag_1d_quad_osc_1(I0integrand_osc, log(1./kc), log(etamaxval), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
			nag_1d_quad_osc_1(I1integrand_osc, log(1./kc), log(etamaxval), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

			I0temp = I0temp + I0int;
			I1temp = I1temp + I1int;
			//std::cout << "Case 34 " << I0int << "(+/-)" << I0err << "   " << I1int << "(+/-)" << I1err << std::endl;

			nag_qp_free(I0qp);
			nag_qp_free(I1qp);

		} // end of integration
		//Fill Delta2GSR and Delta2NLA...
		Delta2GSR.push_back( exp(Gmin + (iscale*I0temp) + log(1. + (iscale*iscale*I1temp*I1temp)/2.)) );
		//std::cout << kmpc << "   " << Delta2GSR.back() << "   " << Gmin << "   " << I0temp << "   " << I1temp << std::endl;


		//for(int i=etamaxindex; i<etaminindex; i++){
		//	break;
		//	lneta_in = (*lneta)[i];
		//	printf("%.9e  %.9e  %.9e  %.9e  %.9e  %.9e  %.9e\n", lneta_in, I0integrand(lneta_in, &Icomm), I0integrand_osc(lneta_in, &Icomm_osc),  win::W((p.get_kc()) *exp(lneta_in)), p.interp_Gp(lneta_in), gsl_I1integrand(lneta_in, &p), win::X((p.get_kc())*exp(lneta_in)));
		//}

		/*	
		nag_quad_1d_fin_smooth(I0integrand, (*lneta)[etaminindex], log(1/kc), 1.e-8, 1.e-8, &I0int, &I0err, &Icomm);
		nag_quad_1d_fin_smooth(I1integrand, (*lneta)[etaminindex], log(1/kc), 1.e-8, 1.e-8, &I1int, &I1err, &Icomm);

		I0temp = I0temp + I0int;
		I1temp = I1temp + I1int;


		//Do oscillatory integral
		nag_1d_quad_osc_1(I0integrand_osc, log(1/kc), (*lneta)[etamaxindex], 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
		nag_1d_quad_osc_1(I1integrand_osc, log(1/kc), (*lneta)[etamaxindex], 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

		I0temp = I0temp + I0int;
		I1temp = I1temp + I1int;

		nag_qp_free(I0qp);
		nag_qp_free(I1qp);
		*/
		/*		
		nag_quad_1d_fin_smooth(I0integrand, (*lneta)[etaminindex], log(etalo), 1.e-8, 1.e-8, &I0int, &I0err, &Icomm);
		nag_quad_1d_fin_smooth(I1integrand, (*lneta)[etaminindex], log(etalo), 1.e-8, 1.e-8, &I1int, &I1err, &Icomm);

		I0temp = I0temp + I0int;
		I1temp = I1temp + I1int;

		//std::cout << I0err << "   " << I1err << std::endl;
		//Do oscillatory integral
		nag_1d_quad_osc_1(I0integrand_osc, log(etalo), log(etahi), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
		nag_1d_quad_osc_1(I1integrand_osc, log(etalo), log(etahi), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);

		I0temp = I0temp + I0int;
		I1temp = I1temp + I1int;
		//std::cout << I0err << "   " << I1err << std::endl;

		nag_qp_free(I0qp);
		nag_qp_free(I1qp);

		nag_quad_1d_fin_smooth(I0integrand, log(etahi), log(1./kc), 1.e-8, 1.e-8, &I0int, &I0err, &Icomm);
		nag_quad_1d_fin_smooth(I1integrand, log(etahi), log(1./kc), 1.e-8, 1.e-8, &I1int, &I1err, &Icomm);
		//std::cout << I0err << "   " << I1err << std::endl;

		I0temp = I0temp + I0int;
		I1temp = I1temp + I1int;

		nag_1d_quad_osc_1(I0integrand_osc, log(1./kc), log(etamaxval), 1.e-6, 1.e-6, 1.e6, &I0int, &I0err, &I0qp, &Icomm_osc, &I0fail);
		nag_1d_quad_osc_1(I1integrand_osc, log(1./kc), log(etamaxval), 1.e-6, 1.e-6, 1.e6, &I1int, &I1err, &I1qp, &Icomm_osc, &I1fail);
		//std::cout << I0err << "   " << I1err << std::endl;

		I0temp = I0temp + I0int;
		I1temp = I1temp + I1int;

		nag_qp_free(I0qp);
		nag_qp_free(I1qp);
		*/

	//for(int kin=200; kin<501; kin++){
	//	break;
	//	kmpc = (*ks)[kin];
	//	kc = kmpc;//*mpcconv;
	//	p.set_kc(kc);

	//	I0.params = &p;
	//	I1.params = &p;

	//	gsl_integration_qags (&I0, (*lneta)[etaminindex], log(etalo), 0, 1e-4, 1e6, w, &gslI0, &gslErr0);
	//	gsl_integration_qags (&I1, (*lneta)[etaminindex], log(etalo), 0, 1e-4, 1e6, w, &gslI1, &gslErr1);
 	//	Delta2GSR.push_back( exp(Gmin + gslI0 + log(1. + gslI1*gslI1/2.)) );
        //      std::cout << kmpc << "   " << Delta2GSR.back() << "   " << I0temp << "   " << I0err << std::endl;
	//}// end of GSL kin loop


		/*
		//This was commented out
 		// Integral after 1/k
		for(int i=0; i<10; i++){
			nag_1d_quad_gen_1(I0integrand, log(1/kc) + i, log(1/kc) + i+1, 1.e-6, 1.e-6, 1e7, &I0int, &I0err, &I0qp, &Icomm, &I0fail);
			nag_1d_quad_gen_1(I1integrand, log(1/kc) + i, log(1/kc) + i+1, 1.e-6, 1.e-6, 1e7, &I1int, &I1err, &I1qp, &Icomm, &I1fail);
			I0temp = I0temp + I0int;
			I1temp = I1temp + I1int;

			nag_qp_free(I0qp);
			nag_qp_free(I1qp);
		}
		// until here (end of inner comment)
		*/
		//Consider making lnetamin and lnetamax depend on k but OK for now
		//nag_1d_quad_osc_1(I0integrand, log(etaminval), log(etamaxval), 1.e-3, 1.e-3, 1e5, &I0temp, &I0err, &I0qp, &Icomm, &I0fail);
		//nag_1d_quad_osc_1(I1integrand, log(etaminval), log(etamaxval), 1.e-3, 1.e-3, 1e5, &I1temp, &I1err, &I1qp, &Icomm, &I1fail);

		//I0.push_back(I0temp);
		//I1.push_back(I1temp/sqrt(2.));

	//double lnetainterval; //4. is a guess but OK for now

	//double gslI0, gslI1, gslErr0, gslErr1;
	//gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e8);
	//gsl_function I0;
	//gsl_function I1;
	//I0.function = &gsl_I0integrand;
	//I1.function = &gsl_I1integrand;
	//GSL Library implementation
	//print integrand for a check
	//double lneta_in;
	//std::cout << "keta range = " << ((*ks)[0] * etaminval) << "   " << ((*ks).rbegin()[0] * etamaxval) << std::endl;

	
	//for(int i=etamaxindex; i<etaminindex; i=i+1){
	//	break;
	//	p.set_kc((*ks)[400]); // NO MPCCONV!!!
	//	lneta_in = (*lneta)[i];
	//	//std::cout << "STUPID TRIAL " << p.get_kc() << "   " << exp(lneta_in) << "   " << std::endl;
	//	printf("%.9e  %.9e  %.9e  %.9e  %.9e\n", lneta_in, (p.get_kc()*exp(lneta_in)), gsl_I0integrand(lneta_in, &p), win::W((p.get_kc()) *exp(lneta_in)), p.interp_Gp(lneta_in));
	//	fflush(stdout);

	//}




	} //end of NAG kin loop
	



	/*	
	for(int kin=100; kin<101; kin++){
		kmpc = (*ks)[kin];
		kc = kmpc * mpcconv;

		//Calculate I0 and I1
		//Write this in an integrator method in utils but OK for now
		I0temp=0.;
		I1temp=0.;

		if(((etaminindex-etamaxindex)%2) == 0){ //Use Simpson's Rule throughout
			for(int i=etamaxindex; i!=etaminindex; i+=2){
				I0temp += h*( (2.*win::W(kc*(*eta)[i])*(*Gp)[i]/((*a)[i]*(*H)[i]*(*eta)[i])) + (4.*win::W(kc*(*eta)[i+1])*(*Gp)[i+1]/((*a)[i+1]*(*H)[i+1]*(*eta)[i+1])) )/3.;	
				I1temp += h*( (2.*win::X(kc*(*eta)[i])*(*Gp)[i]/((*a)[i]*(*H)[i]*(*eta)[i])) + (4.*win::X(kc*(*eta)[i+1])*(*Gp)[i+1]/((*a)[i+1]*(*H)[i+1]*(*eta)[i+1])) )/3.;	
			}

			I0temp += h*((-1.*win::W(kc*(*eta)[etamaxindex])*(*Gp)[etamaxindex]/((*a)[etamaxindex]*(*H)[etamaxindex]*(*eta)[etamaxindex])) + (win::W(kc*(*eta)[etaminindex])*(*Gp)[etaminindex]/((*a)[etaminindex]*(*H)[etaminindex]*(*eta)[etaminindex])))/3.; // correction for under/over counting in the above loop
			I1temp += h*((-1.*win::X(kc*(*eta)[etamaxindex])*(*Gp)[etamaxindex]/((*a)[etamaxindex]*(*H)[etamaxindex]*(*eta)[etamaxindex])) + (win::X(kc*(*eta)[etaminindex])*(*Gp)[etaminindex]/((*a)[etaminindex]*(*H)[etaminindex]*(*eta)[etaminindex])))/3.; // correction for under/over counting in the above loop


		} else { // Use trapezoidal rule for first interval and then Simpson's throughout
			I0temp += h*(( win::W(kc*(*eta)[etamaxindex])*(*Gp)[etamaxindex]/((*a)[etamaxindex]*(*H)[etamaxindex]*(*eta)[etamaxindex]) ) + ( win::W(kc*(*eta)[etamaxindex])*(*Gp)[etamaxindex]/((*a)[etamaxindex]*(*H)[etamaxindex]*(*eta)[etamaxindex]) ))/2.;
			I1temp += h*(( win::X(kc*(*eta)[etamaxindex])*(*Gp)[etamaxindex]/((*a)[etamaxindex]*(*H)[etamaxindex]*(*eta)[etamaxindex]) ) + ( win::X(kc*(*eta)[etamaxindex])*(*Gp)[etamaxindex]/((*a)[etamaxindex]*(*H)[etamaxindex]*(*eta)[etamaxindex]) ))/2.;
			std::cout << "THIS IS KETA " << (kc*(*eta)[etaminindex]) << "   " << (kc*(*eta)[etamaxindex]) << std::endl;
			for(int i=etamaxindex+1; i!=etaminindex; i+=2){
				I0temp += h*( (2.*win::W(kc*(*eta)[i])*(*Gp)[i]/((*a)[i]*(*H)[i]*(*eta)[i])) + (4.*win::W(kc*(*eta)[i+1])*(*Gp)[i+1]/((*a)[i+1]*(*H)[i+1]*(*eta)[i+1])) )/3.;
				I1temp += h*( (2.*win::X(kc*(*eta)[i])*(*Gp)[i]/((*a)[i]*(*H)[i]*(*eta)[i])) + (4.*win::X(kc*(*eta)[i+1])*(*Gp)[i+1]/((*a)[i+1]*(*H)[i+1]*(*eta)[i+1])) )/3.;
				std::cout << (*eta)[i] << "   " << win::W(kc*(*eta)[i]) << "   " << (*Gp)[i] << std::endl;
			}

			I0temp += h*((-1.*win::W(kc*(*eta)[etamaxindex+1])*(*Gp)[etamaxindex+1]/((*a)[etamaxindex+1]*(*H)[etamaxindex+1]*(*eta)[etamaxindex+1])) + (win::W(kc*(*eta)[etaminindex])*(*Gp)[etaminindex]/((*a)[etaminindex]*(*H)[etaminindex]*(*eta)[etaminindex])))/3.; // correction for under/over counting in the above loop
			I1temp += h*((-1.*win::X(kc*(*eta)[etamaxindex+1])*(*Gp)[etamaxindex+1]/((*a)[etamaxindex+1]*(*H)[etamaxindex+1]*(*eta)[etamaxindex+1])) + (win::X(kc*(*eta)[etaminindex])*(*Gp)[etaminindex]/((*a)[etaminindex]*(*H)[etaminindex]*(*eta)[etaminindex])))/3.; // correction for under/over counting in the above loop

		}
		

		//Write results to I0 and I1
		I0.push_back(I0temp);
		I1.push_back(I1temp/sqrt(2.));

		//Fill Delta2GSR and Delta2NLA...
		Delta2GSR.push_back( exp(Gmin + I0temp + log(1. + I1temp*I1temp/2.)) );
		std::cout << kmpc << "   " << Delta2GSR.back() << "   " << I0temp << std::endl;

	} // end of kin loop
	*/

} // end of approx_calc()


void approx::writetofile(std::string name, bool suf){
	std::string suffix = "";
	std::string suffix2 = "";

	if (suf)
		suffix = suffix + "_b" + std::to_string((mybkg->p.b)) + "_c" + std::to_string((mybkg->p.c)) + "_d" + std::to_string((mybkg->p.d)) + ".dat";

	if (!suf)
		suffix2 = suffix2 + ".dat";


        std::ofstream myfile;
        myfile.open(name + suffix + suffix2);

        for(int i=0; i!=Delta2GSR.size(); i++){
                myfile << ((*ks)[i]) << "  " << Delta2GSR[i] << std::endl;
        }   

        myfile.close();

} // end of write to file












