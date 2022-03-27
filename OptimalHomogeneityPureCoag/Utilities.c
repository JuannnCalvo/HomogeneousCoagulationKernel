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

//formula for the profile to be inverted
double ge(double zeta){ 

  return 4.0*exp(-2*zeta);//test case with the constant kernel
  
};


////
//B-Spline evaluations:

double linear_Bcombination (double x, void *params){

	double sum=0.0;
  double top_range;
	int i;
	struct bspline_comb_params *p;
  size_t top_index;

	p= (struct bspline_comb_params *) params;
  top_index=(((p->bw)->knots)->size)-1;
  top_range=gsl_vector_get((p->bw)->knots,top_index);
  //(bw->knots) is a gsl vector containing the nodes
  //(gsl_vector-> size) is a size_t

  if(x<=top_range){
	   gsl_bspline_eval(x,p->B,p->bw);

	   for (i=0;i<(p->nfuncs);i++){
		    sum=sum+gsl_vector_get(p->coefs,i)*gsl_vector_get(p->B,i);
	   };
  };//end_if (else return zero)

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

//Auxiliary integrands
double integrand(double y, void * params) { 
  struct inner_integral_params *p;
  double landa, xi, f;

  p= (struct inner_integral_params *) params;
  landa = (p->homogeneity_degree); 
  xi = (p->xi);
  
  f =pow(y,2+landa)*ge(y)*ge(y*xi);

  return f;
};



//Form integrands for the discretized inverse problem
double Xintegrand (double xi, void *params){ 

  struct bspline_params0 pars_bspline_kernel;
  struct Xintegrand_params *p;
  double landa; 
  double zeta;
  double integral1=0.0;
  double integral2=0.0;
  double error; 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(WSIZE);
  

  p= (struct Xintegrand_params *) params;

  pars_bspline_kernel.index=(p->index_1);
  pars_bspline_kernel.B=(p->B_reg);
  pars_bspline_kernel.bw=(p->bw_reg);

  landa=(p->homogeneity_degree);
  zeta=(p->collocation_point);
  
  //Create the functions to integrate
  struct inner_integral_params pars_int ={landa,xi};

  gsl_function Int;
  Int.function=&integrand;
  Int.params = &pars_int;

  
   gsl_integration_qag(&Int, 
                            zeta/(1.0+xi),  
                            zeta/(xi+SOFTENING), 
                            ABSERR, 
                            0,
                            WSIZE,
                            RULE1,
                            w, 
                            &integral1, 
                            &error
                            );

   if(xi>sqrt(ABSERR)){ 
      gsl_integration_qag(&Int, 
                            zeta/(1.0+xi),  
                            zeta, 
                            ABSERR, 
                            0,
                            WSIZE,
                            RULE1,
                            w, 
                            &integral2, 
                            &error
                            );
    }
    else{ 
      integral2=xi*pow(zeta,3+landa)*ge(zeta)*ge(xi*zeta)/(1+xi);
    };
  
    gsl_integration_workspace_free (w);
  return -betaspline0(xi,&pars_bspline_kernel)*(xi*integral1+integral2);

};




