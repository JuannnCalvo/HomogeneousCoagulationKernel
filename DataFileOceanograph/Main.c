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


int main (int argc, char **argv){

////////////////
//Declarations
/////////////////
size_t i, j;

//Handling model's constants 
double homogeneity_degree=0.0; 
double c_S=1.0; //prefactor for Stokes' law
double x_min=0.01; //cutoff for the inflow term
double mu=1.0; //prefactor for the inflow term



/////////
//Input data files
FILE *DATA_FILE, *PARS_FILE; 
gsl_vector *data_abcissae, *data_ordinates;

//Read input parameters
if(argc!=3){
        printf("Error, incorrect argument number\n");
        printf("Use as: \n exec_name \n input_data_file \n input_pars_file \n");
        return(1);
};

if((PARS_FILE=fopen(argv[2],"rt"))==NULL){
        printf("Error: input_pars_file could not be opened\n");
        return(1) ;
};

fscanf(PARS_FILE,"%lf %lf %lf %lf",&homogeneity_degree,
					   &x_min,
					   &mu,
					   &c_S);

fclose(PARS_FILE);

//Abort when the homogeneity degree is not suitable
if(homogeneity_degree>=1.0){
  printf("Too large homogeneity degree for the implementation, aborting execution \n");
  exit(1);
};


//////
if((DATA_FILE=fopen(argv[1],"rt"))==NULL){
        printf("Error: input_data_file could not be opened\n");
        return(1) ;
};

//Read entry file 
//Assume column entry format (no "\n" after the last datum)
//We use GSL vectors to be able to use 
//the GSL routine to perform lest squares fits to the input data
    int c;
    long count, npoints=0;
    
    //Compute number of spatial slots
    while((c=getc(DATA_FILE))!=EOF){if(c=='\n'){npoints++;};};
    npoints++; 
    rewind(DATA_FILE);

    data_abcissae = gsl_vector_alloc(npoints);
    data_ordinates = gsl_vector_alloc(npoints);

    //Read the info
    for(count=0;count<npoints;count++){
        fscanf(DATA_FILE,"%lf %lf",
          gsl_vector_ptr(data_abcissae,count),
          gsl_vector_ptr(data_ordinates,count)
          );
    };
    
    fclose(DATA_FILE);

//Abort when the number of data points is low
if(npoints<NCOEFFS){
  printf("Too few points, aborting execution \n");
  exit(1);
};

double top_range=gsl_vector_get(data_abcissae,npoints-1);
double bottom_range=gsl_vector_get(data_abcissae,0);

//////////////////////////
//B-splines I: data filtering

//nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 
const size_t nbreak = NCOEFFS - 2; 

gsl_bspline_workspace *bw_fit; 
//allocate a cubic bspline workspace (k = 4) 
bw_fit = gsl_bspline_alloc(4, nbreak);

//We shall use uniform breakpoints on [0, top_range] 
gsl_bspline_knots_uniform(0.0,top_range, bw_fit);


/////////////////////
//Interpolate input data with the constructed B-splines

gsl_vector *B_fit, *coefs_fit;
B_fit = gsl_vector_alloc(NCOEFFS);
coefs_fit=gsl_vector_alloc(NCOEFFS);

//Zero inicialization
for(i=0;i<NCOEFFS;i++){
    gsl_vector_set(coefs_fit,i,0.0);
};

gsl_matrix *X,*cov;
X=gsl_matrix_alloc(npoints,NCOEFFS);
cov=gsl_matrix_alloc(NCOEFFS,NCOEFFS);

gsl_multifit_linear_workspace *mv_fit;
mv_fit=gsl_multifit_linear_alloc(npoints,NCOEFFS);

double xi, chisq=-1.0;
for (i = 0; i < npoints; ++i){
    xi=gsl_vector_get(data_abcissae,i);
    gsl_bspline_eval(xi,B_fit,bw_fit);

    for (j= 0; j < NCOEFFS; ++j){
        gsl_matrix_set(X,i,j,gsl_vector_get(B_fit,j));
    };
};

gsl_multifit_linear(X,
                    data_ordinates,
                    coefs_fit,
                    cov,
                    &chisq,  
                    mv_fit);


///////
//Output reconstructed curve:
struct bspline_comb_params params_fit ={coefs_fit,bw_fit,B_fit,NCOEFFS};

printf("Reconstructed data curve:\n");
for (i = 0; i < npoints; ++i){
	printf("%f %f\n",gsl_vector_get(data_abcissae,i),
        linear_Bcombination (gsl_vector_get(data_abcissae,i),&params_fit)
        );
};

printf("\n \n");


////////////////
//Constructing log-spaced collocation points on the data range
double collocation_points[NCOLLOCATIONPOINTS]; 
long k; 
double stride;

if(bottom_range<=SOFTENING){
	stride=-log2((bottom_range+SOFTENING)/top_range)/(NCOLLOCATIONPOINTS-1.0);
}
else{
	stride=-log2(bottom_range/top_range)/(NCOLLOCATIONPOINTS-1.0);
};

collocation_points[0]=bottom_range;
for(k = 1; k < NCOLLOCATIONPOINTS-1; ++k){
    collocation_points[k]=top_range*pow(2.0,-(NCOLLOCATIONPOINTS-1-k)*stride); 
};
collocation_points[NCOLLOCATIONPOINTS-1]=top_range;


////////////////////////
//Linear system: data structures
gsl_matrix *coeff_matrix =gsl_matrix_alloc(NCOLLOCATIONPOINTS,NCOEFFS);
gsl_vector *data_vector =gsl_vector_alloc(NCOLLOCATIONPOINTS);

gsl_matrix *regularizer_matrix =gsl_matrix_calloc(NCOEFFS,NCOEFFS);
gsl_matrix *regularizer_matrix_0 =gsl_matrix_alloc(NCOEFFS,NCOEFFS);
gsl_matrix *regularizer_matrix_1 =gsl_matrix_alloc(NCOEFFS,NCOEFFS);
gsl_matrix *regularizer_matrix_2 =gsl_matrix_alloc(NCOEFFS,NCOEFFS);


///////////////
//B-splines II: regularizing matrix 
gsl_bspline_workspace *bw_reg; 
gsl_vector *B_reg;  
gsl_matrix *B_derivs;

//allocate a cubic bspline workspace (k = 4) 
bw_reg = gsl_bspline_alloc(4, nbreak);
B_reg = gsl_vector_alloc(NCOEFFS);
B_derivs= gsl_matrix_alloc(NCOEFFS,4);
//use uniform breakpoints on [0, 1] 
gsl_bspline_knots_uniform(0.0, 1.0, bw_reg);

//Strucutures to construct regularizer matrix:
struct bspline_prod0 pars_bscalarprod0={0,0,bw_reg,B_reg};
//L2 product: {index 1, index 2, b-workspace, B-vector}

struct bspline_prod012 pars_bscalarprod012={0,0,bw_reg,B_derivs};
//H2 product (sort of): {index 1, index 2, b-workspace, B-matrix}


///////////////////
//Parameters for the integrator 
double result, error; 
gsl_integration_workspace * w = gsl_integration_workspace_alloc(WSIZE);

//Integral structures for the coefficient matrix 
struct Xintegrand_params pars_Xintegrand={0, //bspline index
										1.0, //the x_i
										bw_reg,
										B_reg,
                                        &params_fit, //pointer to data fit struct (linear combination of B-splines)
										homogeneity_degree //homogeneity exponent
										};

//Defining the function structure to integrate for the coefficient matrix
gsl_function X_integrands; 
X_integrands.function=&Xintegrand;
X_integrands.params = &pars_Xintegrand;

//Defining the function structure for Stokes' law
gsl_function stokes_integrand; 
stokes_integrand.function=&StokesIntegrand;
stokes_integrand.params = &params_fit;

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
gsl_vector *kernel_expansion_coefs_gcv = gsl_vector_alloc(NCOEFFS); //regularized solution (GCV) 
gsl_vector *reg_param = gsl_vector_alloc(NPOINTS_CURVE);
gsl_vector *rho = gsl_vector_alloc(NPOINTS_CURVE); // residual norms 
gsl_vector *eta = gsl_vector_alloc(NPOINTS_CURVE); // solution norms 
gsl_vector *G = gsl_vector_alloc(NPOINTS_CURVE); // GCV function values 

double rnorm, snorm;
double lambda_l; // optimal regularization parameter L-curve
double lambda_gcv; // optimal regularization parameter GCV
double G_gcv; // G(lambda_gcv) 
size_t reg_idx; // index of optimal lambda //The actual weight is lambda^2
double rcond; // reciprocal condition number of X

gsl_vector *vector_h,*vector_gamma;
vector_h = gsl_vector_alloc(NCOEFFS);
vector_gamma = gsl_vector_alloc(NCOEFFS);

///////////////////
//end-decls
//////////////////////


/////////////
//Construct data vector 
double aux;
for(i=0;i<NCOLLOCATIONPOINTS;i++){ 
	//We compute Stokes' term separately
	gsl_integration_qag(&stokes_integrand, //ref to the function to integrate
                            0, //left end
                            collocation_points[i], //right end
                            1e-3, //absolute error limit 
                            0, //relative error limit 
                            WSIZE, //max workspace size
                            RULE6,
                            w, //workspace
                            &result, 
                            &error); 

	aux=-0.5*mu*pow(min(collocation_points[i],x_min),2)+c_S*result;
	gsl_vector_set(data_vector,i,aux); 
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
					1e-3,//ABSERR, //absolute error limit
					0, //relative error limit (we compute to the specified absolute error)
					WSIZE, //max workspace size
					RULE1, //(15)-point rule
					w, //workspace
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
					w, 
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
					w, 
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
					w, 
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
//coeff_matrix is handled as constt
//The factors are stored in the workspace w_lin

rcond = gsl_multifit_linear_rcond(w_lin);
fprintf(stdout, "matrix condition number = %e\n\n", 1.0 / rcond);


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

////////////////////
//Output various infos (L-curve strategy)

chisq = pow(rnorm, 2.0) + pow(lambda_l * snorm, 2.0);

fprintf(stdout, "\n=== Regularized fit (L-curve) ===\n");
fprintf(stdout, "optimal lambda: %g\n", lambda_l);
fprintf(stdout, "residual norm = %g\n", rnorm); 
fprintf(stdout, "solution norm = %g\n", snorm); 
fprintf(stdout, "chisq/dof = %g\n", chisq / (NCOLLOCATIONPOINTS - NCOEFFS));

struct bspline_comb_params params_linearcomb_lcurve ={kernel_expansion_coefs_lcurve,bw_reg,B_reg,NCOEFFS};
//{coeff list, b-workspace, B-vector, number of functions to combine}

fprintf(stdout, "reconstructed kappa:\n");
for (i = 0; i < NCOLLOCATIONPOINTS+1; ++i){
	printf("%f %f\n", (1.0/NCOLLOCATIONPOINTS)*i,
		linear_Bcombination ((1.0/NCOLLOCATIONPOINTS)*i,&params_linearcomb_lcurve));
};


///////////////////////
// calculate GCV curve and find its minimum 

gsl_multifit_linear_gcv(data_vector, 
			reg_param, //sampled values of lambda
			G, //values of G(lambda)
			&lambda_gcv, //optimal lambda
			&G_gcv, //minimum of the G-curve
			w_lin); //workspace

// regularize with lambda_gcv. Solve the problem in standard form:
gsl_multifit_linear_solve(lambda_gcv, //regularization weight is lambda^2
			coeff_matrix,
			data_vector,
			kernel_expansion_coefs_gcv,
			&rnorm,
			&snorm, 
			w_lin);

/////////////////
//Get back to original form
gsl_blas_dtrsv (CblasLower,
		CblasTrans,
		CblasNonUnit,
		regularizer_matrix,
		kernel_expansion_coefs_gcv);


///////////////////
//Output various infos (GCV-curve strategy)
chisq = pow(rnorm, 2.0) + pow(lambda_gcv * snorm, 2.0);

fprintf(stdout, "\n=== Regularized fit (GCV) ===\n");
fprintf(stdout, "optimal lambda: %g\n", lambda_gcv); //regularization weight is lambda^2
fprintf(stdout, "residual norm = %g\n", rnorm);
fprintf(stdout, "solution norm = %g\n", snorm);
fprintf(stdout, "chisq/dof = %g\n", chisq / (NCOLLOCATIONPOINTS - NCOEFFS));

struct bspline_comb_params params_linearcomb_gcv ={kernel_expansion_coefs_gcv,bw_reg,B_reg,NCOEFFS};
//{coeff list, b-workspace, B-vector, number of functions to combine}

fprintf(stdout, "reconstructed kappa:\n");
for (i = 0; i < NCOLLOCATIONPOINTS+1; ++i){
	printf("%f %f\n", (1.0/NCOLLOCATIONPOINTS)*i,
		linear_Bcombination ((1.0/NCOLLOCATIONPOINTS)*i,&params_linearcomb_gcv)
		);
};


////////////////
//Output both curves

printf("\n\nL-curve and GCV-curve (reg_param | rho | eta | G):\n");
printf("To plot the L-curve: (x=rho,y=eta)\n");
printf("To plot the GCV-curve: (x=reg_param , y= G)\n");
// output L-curve and GCV curve 
for (i = 0; i < NPOINTS_CURVE; ++i){
	printf("%e %e %e %e\n",
		gsl_vector_get(reg_param, i),
		gsl_vector_get(rho, i),
		gsl_vector_get(eta, i),
		gsl_vector_get(G, i)
		);
};

// output L-curve corner point 
printf("\n\nL-curve corner point: %.15lf %lf\n",
	gsl_vector_get(rho, reg_idx),gsl_vector_get(eta, reg_idx));

// output GCV curve corner minimum 
printf("\n\nGCV curve corner minimum: %e %e\n",lambda_gcv,G_gcv);

printf("\n"); 


/////////////////////////////////

///////////////////////////

//////////////////////
//Free stuff

free(data_abcissae);
free(data_ordinates);
gsl_bspline_free(bw_fit); 
gsl_vector_free(B_fit); 
gsl_vector_free(coefs_fit); 
gsl_bspline_free(bw_reg); 
gsl_vector_free(B_reg); 
gsl_matrix_free(B_derivs);
gsl_matrix_free(X);
gsl_matrix_free(cov);
gsl_multifit_linear_free(mv_fit);
gsl_multifit_linear_free(w_lin);

gsl_vector_free(kernel_expansion_coefs_lcurve);
gsl_vector_free(kernel_expansion_coefs_gcv);
gsl_vector_free(reg_param);
gsl_vector_free(rho);
gsl_vector_free(eta);
gsl_vector_free(G);
gsl_vector_free(vector_h);
gsl_vector_free(vector_gamma);
gsl_integration_workspace_free (w);
gsl_matrix_free(coeff_matrix);
gsl_vector_free (data_vector);

gsl_matrix_free(regularizer_matrix);
gsl_matrix_free(regularizer_matrix_0);
gsl_matrix_free(regularizer_matrix_1);
gsl_matrix_free(regularizer_matrix_2);

return 0;

}
