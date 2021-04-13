#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>

double kvrod(double x){
	double y=sqrt(x);
return y;
}

double kvrod_2(double x){
	double y=4*sqrt(1-x*x);
return y;
}
