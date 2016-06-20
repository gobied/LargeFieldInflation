// utils.h
// Calculates potentials and their derivatives as well as window functions
// GO
#ifndef utils_h
#define utils_h

#include <cmath>
#include "parameters.h"

namespace pot {

inline double V(const double phi, const step_params& p){
	printf("params = %.5e %.5e %.5e %.5e\n", p.m, p.b, p.c, p.d);
	return ( p.m*p.m*phi*phi*(1+(p.c*tanh((phi- p.b)/p.d)))/2. );
}

inline double Vp(const double phi, const step_params& p){
	double sec2 = 1. - (tanh((phi - p.b)/p.d)*tanh((phi - p.b)/p.d));
        return (p.m*p.m*phi*( 1+(p.c*tanh((phi-p.b)/p.d)) ) + (p.c*p.m*p.m*phi*phi*sec2/(2.*p.d)));
}

inline double Vpp(const double phi, const step_params& p){
	double sec2 = 1. - (tanh((phi - p.b)/p.d)*tanh((phi - p.b)/p.d));
	return (p.m*p.m*(1 + (p.c*tanh((phi-p.b)/p.d)) ) + ( (2.*p.c*p.m*p.m*phi*sec2)/p.d ) - ( (p.c*p.m*p.m*phi*phi*sec2*tanh((phi-p.b)/p.d))/(p.d*p.d) ));
}

}

namespace win {

inline double W(double u){
	return ((3.*sin(2.*u)/(2.*u*u*u)) - (3.*cos(2.*u)/(u*u)) - (3.*sin(2.*u)/(2.*u)));
}

inline double X(double u){
	return (3.*( sin(u) - (u*cos(u)) )*( sin(u) - (u*cos(u)) )/(u*u*u));
}

}

#endif
