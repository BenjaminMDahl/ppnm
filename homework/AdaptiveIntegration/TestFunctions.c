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

double kvrod_2_cc(double x){
	double arg=(double) 0.5*cos(x)+0.5;
	double y=kvrod_2(arg)*sin(x)*0.5;
return y;
}

double inv_kvrod(double x){
	double y=1/sqrt(x);
return y;
}

double inv_kvrod_cc(double x){ //Den transformerede af inv_kvrod
	double arg=(double) 0.5*cos(x)+0.5;
	double y=inv_kvrod(arg)*sin(x)*0.5;
return y;
}

double ln_kvrod(double x){
	double y=log(x)/sqrt(x);
return y;}


double ln_kvrod_cc(double x){
	double arg=(double) 0.5*cos(x)+0.5;
	double y=ln_kvrod(arg)*sin(x)*0.5;
return y;}


