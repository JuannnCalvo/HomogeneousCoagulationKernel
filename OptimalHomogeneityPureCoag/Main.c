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



double nu_homo_uno=2.0;

int main (int argc, char **argv){

	double homogeneity_degree=0.25; //starting value
	//Some starting values may yield an error "endpoints do not enclose a minimum" 
	//Arguably this will happen when the starting value provides a worse score 
	//than the ends of the initial search interval 
	//(specific issue for Brent's minimization method implementation)

	int status;
	int iter = 0, max_iter = 100;

	const gsl_min_fminimizer_type *T; 
	gsl_min_fminimizer *s;

 	double homo_min = -3.0, homo_max = 1.0;
 	
  	

  	gsl_function F;

	F.function = &residual;
	F.params = 0;
	
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc (T);
	gsl_min_fminimizer_set (s, &F, homogeneity_degree, homo_min, homo_max);

	printf ("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "min", "err(est)");
	printf ("%5d [%.7f, %.7f] %.7f %.7f\n", 
		iter, 
		homo_min, 
		homo_max, 
		homogeneity_degree, 
		homo_max - homo_min);


	do {
		iter++;
		status = gsl_min_fminimizer_iterate (s);

		homogeneity_degree = gsl_min_fminimizer_x_minimum (s);
 		homo_min = gsl_min_fminimizer_x_lower (s);
		homo_max = gsl_min_fminimizer_x_upper (s);

		status= gsl_min_test_interval (homo_min, 
										homo_max, 
										0.001, //absolute tolerance
										0.0 //relative tolerance
										);

		if (status == GSL_SUCCESS){
  			printf ("Converged:\n");
  		};	
		
		printf ("%5d [%.7f, %.7f] " "%.7f  %.7f\n",
			iter, homo_min, homo_max, homogeneity_degree, homo_max - homo_min);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	//end do-while

  	gsl_min_fminimizer_free (s);

	return status;
}
