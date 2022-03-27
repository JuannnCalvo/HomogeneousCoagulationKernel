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


extern double nu_homo_uno;//=2.0;
//specific for the borderline case of homogeneity one


double residual(double homogeneity_degree, void *params){

////////////////
//Declarations
/////////////////
size_t i, j;
double nu;
double residual=-1.0;


//Control the range of homogeneity_degree
if(homogeneity_degree<1.0){
  	nu=-1.0/(1.0-homogeneity_degree); 
}
else if(homogeneity_degree==1.0){
	printf("Kernel homogeneity equals one, reading extra input parameter \n");
	nu=nu_homo_uno;
}
else{
	return(1e10); //Too large homogeneity degree for the implementation
};



////////////////
//Constructing log-spaced collocation points on the data range

double collocation_points[NCOLLOCATIONPOINTS]; //(the z_i's) 

for(i=0;i<NCOLLOCATIONPOINTS;i++){ //ranging from 2^-8 to 4 aprox. 
	collocation_points[i]=pow(2.0,-8.0+i/10.0);
};


////////////////////////
//Linear system: data structures
gsl_matrix *coeff_matrix =gsl_matrix_alloc(NCOLLOCATIONPOINTS,NCOEFFS);
gsl_vector *data_vector =gsl_vector_alloc(NCOLLOCATIONPOINTS);

gsl_matrix *regularizer_matrix =gsl_matrix_calloc(NCOEFFS,NCOEFFS);
gsl_matrix *regularizer_matrix_0 =gsl_matrix_alloc(NCOEFFS,NCOEFFS);
gsl_matrix *regularizer_matrix_1 =gsl_matrix_alloc(NCOEFFS,NCOEFFS);
gsl_matrix *regularizer_matrix_2 =gsl_matrix_alloc(NCOEFFS,NCOEFFS);


///////////////
//B-splines: regularizing matrix and kernel expansion

//nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 
const size_t nbreak = NCOEFFS - 2; 

gsl_bspline_workspace *bw_reg; 
gsl_vector *B_reg;  
gsl_matrix *B_derivs;

//allocate a cubic bspline workspace (k = 4) 
bw_reg = gsl_bspline_alloc(4, nbreak);
B_reg = gsl_vector_alloc(NCOEFFS);
B_derivs= gsl_matrix_alloc(NCOEFFS,4);
//use uniform breakpoints on [0, 1] 
gsl_bspline_knots_uniform(0.0, 1.0, bw_reg);

//Data structures to construct regularizer matrix:
struct bspline_prod0 pars_bscalarprod0={0,0,bw_reg,B_reg};
//L2 product: {index 1, index 2, b-workspace, B-vector}

struct bspline_prod012 pars_bscalarprod012={0,0,bw_reg,B_derivs};
//H2 product (sort of): {index 1, index 2, b-workspace, B-matrix}


///////////////////
//Parameters for the integrator 
double result, error; 
gsl_integration_workspace * w_int = gsl_integration_workspace_alloc(WSIZE);

//Integral structures for the coefficient matrix 
struct Xintegrand_params pars_Xintegrand={0, //bspline index
										0.0, //the z_i
										bw_reg,//pertaining kernel basis
										B_reg, //pertaining kernel basis
                                        homogeneity_degree //homogeneity exponent
										};

//Defining the function structure to integrate for the coefficient matrix
gsl_function X_integrands; 
X_integrands.function=&Xintegrand;
X_integrands.params = &pars_Xintegrand;

//Same thing for the regularizer matrix
gsl_function W_integrands0; //regularizer, zeroth order 
W_integrands0.function=&betapairproduct0;
W_integrands0.params = &pars_bscalarprod0;

gsl_function W_integrands1; //regularizer, first derivatives
W_integrands1.function=&betapairproduct1;
W_integrands1.params = &pars_bscalarprod012;

gsl_function W_integrands2; //regularizer, second derivatives
W_integrands2.function=&betapairproduct2;
W_integrands2.params = &pars_bscalarprod012;


////////////////
//L-curve and GCV-curve handlers
gsl_multifit_linear_workspace *w_lin = gsl_multifit_linear_alloc(NCOLLOCATIONPOINTS, NCOEFFS);

gsl_vector *kernel_expansion_coefs_lcurve = gsl_vector_alloc(NCOEFFS); //regularized solution (L-curve) 
gsl_vector *reg_param = gsl_vector_alloc(NPOINTS_CURVE);
gsl_vector *rho = gsl_vector_alloc(NPOINTS_CURVE); // residual norms 
gsl_vector *eta = gsl_vector_alloc(NPOINTS_CURVE); // solution norms 

double rnorm, snorm;
double lambda_l; // optimal regularization parameter L-curve
size_t reg_idx; // index of optimal lambda //The actual weight is lambda^2
double rcond; // reciprocal condition number of X

///////////////////
//end-decls
//////////////////////

/////////////
//Construct data vector (rhs of the linear system)
for(i=0;i<
	NCOLLOCATIONPOINTS
	;i++){ 
    gsl_vector_set(data_vector,i,
                    nu*collocation_points[i]*collocation_points[i]
                        *ge(collocation_points[i])
                        );
    
};

////////////////////
//Construct coefficient matrix
//nrows= #collocation points, ncols=#basis splines


for(i=0;i<NCOLLOCATIONPOINTS;i++){
    pars_Xintegrand.collocation_point=collocation_points[i];

    for(j=0;j<NCOEFFS;j++){

        pars_Xintegrand.index_1=j;
        X_integrands.params = &pars_Xintegrand;
        gsl_integration_qag (&X_integrands, //ref to the function to integrate
					0, //left end
					1, //right end
					ABSERR, //absolute error limit
					0, //relative error limit (we compute to the specified absolute error)
					WSIZE, //max workspace size
					RULE1, //(15)-point rule
					w_int, //workspace
					&result, 
					&error); //estimate for absolute error
       
        gsl_matrix_set(coeff_matrix,i,j,result);

    };
};



/////////////
//Construct regularizing matrix

for(i=0;i<NCOEFFS;i++){
	pars_bscalarprod0.index_1=i;

	for(j=0;j<NCOEFFS;j++){

		pars_bscalarprod0.index_2=j;
		W_integrands0.params = &pars_bscalarprod0;

		gsl_integration_qag (&W_integrands0, 
					0, 
					1, 
					ABSERR, 
					RELERR, 
					WSIZE, 
					RULE1,
					w_int, 
					&result, 
					&error);
	
		gsl_matrix_set(regularizer_matrix_0,i,j,result);

    };
};


//Include first and second order derivatives as well: 

for(i=0;i<NCOEFFS;i++){
	pars_bscalarprod012.index_1=i;

	for(j=0;j<NCOEFFS;j++){

		pars_bscalarprod012.index_2=j;
		W_integrands1.params = &pars_bscalarprod012;

		gsl_integration_qag (&W_integrands1, 
					0, 
					1, 
					ABSERR, 
					RELERR, 
					WSIZE, 
					RULE1,
					w_int, 
					&result, 
					&error);
	
		gsl_matrix_set(regularizer_matrix_1,i,j,result);

    };
};


for(i=0;i<NCOEFFS;i++){
	pars_bscalarprod012.index_1=i;

	for(j=0;j<NCOEFFS;j++){

		pars_bscalarprod012.index_2=j;
		W_integrands2.params = &pars_bscalarprod012;

		gsl_integration_qag (&W_integrands2, 
					0, 
					1, 
					ABSERR, 
					RELERR, 
					WSIZE, 
					RULE1,
					w_int, 
					&result, 
					&error);
	
		gsl_matrix_set(regularizer_matrix_2,i,j,result);

    };
};

gsl_matrix_add (regularizer_matrix,regularizer_matrix_0);
gsl_matrix_add (regularizer_matrix,regularizer_matrix_1); 
gsl_matrix_add (regularizer_matrix,regularizer_matrix_2);


/////////////
//Passing to standard form and factorizing

//Cholesky-decompose regularizing matrix
//(the lower part of regularizer_matrix is overwritten)
gsl_linalg_cholesky_decomp1 (regularizer_matrix); 

/////////////

//Compute coeff_matrix times the inverse of the transpose of the Cholesky factor 
//(as we regularize and we are not in standard form)
gsl_blas_dtrsm(CblasRight,
		CblasLower,
		CblasTrans,
		CblasNonUnit,
		1.0,
		regularizer_matrix,
		coeff_matrix);
//The result is overwritten in coef_matrix
//Now we are in standard form

////Perform SVD for the coefficient matrix 
//With the ad-hoc routine that goes fine for the GSL solver
gsl_multifit_linear_svd(coeff_matrix, w_lin); 
//coeff_matrix is handled as const
//The factors are stored in the workspace w_lin

rcond = gsl_multifit_linear_rcond(w_lin);


/////////////////////////////
// calculate L-curve and find its corner 
gsl_multifit_linear_lcurve(data_vector, 
			reg_param, //sampled values of lambda
			rho, //residual norms 
			eta, //solution norms
			w_lin); //workspace

gsl_multifit_linear_lcorner(rho, eta, &reg_idx);

// store optimal regularization parameter 
lambda_l = gsl_vector_get(reg_param, reg_idx);

///////////////
// Solve using the optimal lambda_l. First solve the problem in standard form:
gsl_multifit_linear_solve(lambda_l, //regularization weight is lambda^2
			coeff_matrix, 
			data_vector, 
			kernel_expansion_coefs_lcurve, 
			&rnorm, //residual norm 
			&snorm, //solution norm
			w_lin);

//Get back to original form
gsl_blas_dtrsv (CblasLower,
		CblasTrans,
		CblasNonUnit,
		regularizer_matrix,
		kernel_expansion_coefs_lcurve);



/////////////////////////////////

///////////////////////////
residual=gsl_vector_get(rho, reg_idx);
//////////////////////
//Free stuff


gsl_bspline_free(bw_reg); 
gsl_vector_free(B_reg); 
gsl_matrix_free(B_derivs);
gsl_multifit_linear_free(w_lin);

gsl_vector_free(kernel_expansion_coefs_lcurve);
gsl_vector_free(reg_param);
gsl_vector_free(rho);
gsl_vector_free(eta);
gsl_integration_workspace_free (w_int);
gsl_matrix_free(coeff_matrix);
gsl_vector_free (data_vector);

gsl_matrix_free(regularizer_matrix);
gsl_matrix_free(regularizer_matrix_0);
gsl_matrix_free(regularizer_matrix_1);
gsl_matrix_free(regularizer_matrix_2);

return residual;

} //end residual