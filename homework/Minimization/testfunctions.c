#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

//Til opgave A//

double Rosenbrock(gsl_vector* vec){
	double x=gsl_vector_get(vec,0);
	double y=gsl_vector_get(vec,1);
	double value=(1-x)*(1-x)+100*(y-x*x)*(y-x*x);
return value;
}

double Himmelblau(gsl_vector* vec){
	double x=gsl_vector_get(vec,0);
	double y=gsl_vector_get(vec,1);
	double value=(x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
return value;
}
