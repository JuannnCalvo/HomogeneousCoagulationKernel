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


#include <math.h>
#include <stdlib.h> 
#include <string.h> 
#include <stdio.h> 
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_bspline.h> 
#include <gsl/gsl_math.h> 
#include <gsl/gsl_multifit.h> 


//////////
//Symbolic constants
/////////

#define SOFTENING 1e-6 // DEfault 1e-6 
#define NCOLLOCATIONPOINTS 100 //DEfault is 100
#define NCOEFFS 12 //number of B-splines used (number of fit coefficients)
#define NPOINTS_CURVE 1000 //number of points on L-curve and GCV curve 
#define WSIZE 200 //size of integration workspace
#define ABSERR 1e-8 //absolute error limit for the numerical integrator
#define RELERR 1e-8 //relative error limit for the numerical integrator
#define RULE1 1 //15-point rule for the numerical integrator
#define RULE6 6 //most-exhaustive-point rule for the numerical integrator



/////
//Structs
/////


/////
//B-Spline evaluations:

struct bspline_params0{ //for plain evaluation
	int index;
	gsl_bspline_workspace *bw; 
	gsl_vector *B;	
};

struct bspline_comb_params{ //for a linear combination
	gsl_vector *coefs;
	gsl_bspline_workspace *bw; 
	gsl_vector *B;
	int nfuncs;
};

struct bspline_params012{ //for evaluation with derivatives of order 1 & 2
	int index;
	gsl_bspline_workspace *bw; 
	gsl_matrix *B_derivs;	
};

struct bspline_prod0{ //For the scalar product in L2
	int index_1;
	int index_2;
	gsl_bspline_workspace *bw; 
	gsl_vector *B;	
};

struct bspline_prod012{ //For the scalar product in H2
	int index_1;
	int index_2;
	gsl_bspline_workspace *bw; 
	gsl_matrix *B_derivs;	
};


///
//Construction of the coefficient matrix:

struct inner_integral_params{
  double homogeneity_degree;
  double xi;
  struct bspline_comb_params *params_fit;
};

struct Xintegrand_params{
	int index_1; //bspline index
	double collocation_point; //the z_i
	gsl_bspline_workspace *bw_reg;  //pertaining kernel basis
	gsl_vector *B_reg;	//pertaining kernel basis
	struct bspline_comb_params *pars_data; 
			//Pointer to struc containing (pertaining data basis):
			//gsl_vector *coefs;
			//gsl_bspline_workspace *bw_fit; 
			//gsl_vector *B_fit;
			//int nfuncs; 
	double homogeneity_degree; //homogeneity exponent
};


/////
//Functions and integrands
/////



double linear_Bcombination (double x, void *params);

double betaspline0 (double x, void *params);	

double betapairproduct0 (double x, void *params);

double betaspline1 (double x, void *params);	

double betapairproduct1 (double x, void *params);

double betaspline2 (double x, void *params);	

double betapairproduct2 (double x, void *params);

double integrand(double x, void * params); 

double Xintegrand (double xi, void *params);






