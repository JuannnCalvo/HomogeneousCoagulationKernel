/********************************************************************
Copyright (C) 2022  Juan Calvo Yagüe

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**********************************************************************/


#include "Header.h"


////
//B-Spline evaluations:

//Forms a linear combination of B-splines
double linear_Bcombination (double x, void *params){

	double sum=0.0;
	int i;
	struct bspline_comb_params *p;

	p= (struct bspline_comb_params *) params;

	gsl_bspline_eval(x,p->B,p->bw);

	for (i=0;i<(p->nfuncs);i++){
		sum=sum+gsl_vector_get(p->coefs,i)*gsl_vector_get(p->B,i);
	};

	return sum;

};


//Plain evaluation of a B-spline
double betaspline0 (double x, void *params){

	struct bspline_params0 *p;

	p= (struct bspline_params0 *) params;

	gsl_bspline_eval(x,p->B,p->bw);

	return gsl_vector_get(p->B,p->index); 
	
};


//Form the integrand for the L2 scalar product
double betapairproduct0 (double x, void *params){

	struct bspline_params0 pars1, pars2;
	struct bspline_prod0 *p;

	p= (struct bspline_prod0 *) params;

	pars1.index=(p->index_1);
	pars1.B=(p->B);
	pars1.bw=(p->bw);

	pars2.index=(p->index_2);
	pars2.B=(p->B);
	pars2.bw=(p->bw);

	return betaspline0(x,&pars1)*betaspline0(x,&pars2);
};


//Evaluates first derivative
double betaspline1 (double x, void *params){

	struct bspline_params012 *p;

	p= (struct bspline_params012 *) params;

	gsl_bspline_deriv_eval(x,3,p->B_derivs,p->bw);

	return gsl_matrix_get(p->B_derivs,p->index,1); 
	
};


//Evaluates second derivative
double betaspline2 (double x, void *params){

	struct bspline_params012 *p;

	p= (struct bspline_params012 *) params;

	gsl_bspline_deriv_eval(x,3,p->B_derivs,p->bw);

	return gsl_matrix_get(p->B_derivs,p->index,2); 
	
};


//Integrand for scalar product of first derivatives
double betapairproduct1 (double x, void *params){

	struct bspline_params012 pars1, pars2;
	struct bspline_prod012 *p;

	p= (struct bspline_prod012 *) params;

	pars1.index=(p->index_1);
	pars1.B_derivs=(p->B_derivs);
	pars1.bw=(p->bw);

	pars2.index=(p->index_2);
	pars2.B_derivs=(p->B_derivs);
	pars2.bw=(p->bw);

	return betaspline1(x,&pars1)*betaspline1(x,&pars2);
};


//Integrand for scalar product of second derivatives
double betapairproduct2 (double x, void *params){

	struct bspline_params012 pars1, pars2;
	struct bspline_prod012 *p;

	p= (struct bspline_prod012 *) params;

	pars1.index=(p->index_1);
	pars1.B_derivs=(p->B_derivs);
	pars1.bw=(p->bw);

	pars2.index=(p->index_2);
	pars2.B_derivs=(p->B_derivs);
	pars2.bw=(p->bw);

	return betaspline2(x,&pars1)*betaspline2(x,&pars2);
};


////
//Construction of the coefficient matrix:

//Form integrands for the discretized inverse problem
double Xintegrand (double xi, void *params){ 

	struct bspline_params0 pars_bspline;
	struct Xintegrand_params *p;
	double zeta;
	double term1=0.0;
	double term2=0.0;

	
	p= (struct Xintegrand_params *) params;

	zeta=(p->collocation_point);

	pars_bspline.index=(p->index_1);
	pars_bspline.B=(p->B);
	pars_bspline.bw=(p->bw);

	
	if(zeta/xi>-2.0*log(0.5*ABSERR)){
		term1=-exp(-0.5*zeta)/(PI*(1.0+xi)*sqrt(xi));
	}
	else{term1=(exp(-(1.0+xi)*zeta/(2.0*xi)) -exp(-0.5*zeta))/(PI*(1.0+xi)*sqrt(xi));};

	if(xi<sqrt(0.5*ABSERR)){
		term2=-zeta*exp(-0.5*zeta)/(2.0*PI*(1.0+xi)*sqrt(xi));
	}
	else{term2=(exp(-(1.0+xi)*0.5*zeta)-exp(-0.5*zeta))/(PI*(1.0+xi)*pow(xi,1.5));};


	return betaspline0(xi,&pars_bspline)*(term1+term2);			

};
