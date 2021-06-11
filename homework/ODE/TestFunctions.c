#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>

//ODE y''=-y
void odecos(double t,gsl_vector* y, gsl_vector* dydt){
	double y1=gsl_vector_get(y,1);
	double y2=-gsl_vector_get(y,0);
	gsl_vector_set(dydt,0,y1);
	gsl_vector_set(dydt,1,y2);
}

//ODE: y'=1
void line(double t,gsl_vector* y, gsl_vector* dydt){
	double y0=1;
	gsl_vector_set(dydt,0,y0);}

//ODE: SIR-model
void epidemic(double t,gsl_vector* y, gsl_vector* dydt){
	double N=5500000, Tr=7, Tc=2;
	double y0=-(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc);
	double y1=(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc)-gsl_vector_get(y,1)/Tr;
	double y2=gsl_vector_get(y,1)/Tr;
	gsl_vector_set(dydt,0,y0);				//dS/dt=-(I*S)/(N*T_c)
	gsl_vector_set(dydt,1,y1);		//dI/dt=(I*S)/(N*T_c)-I/T_r
	gsl_vector_set(dydt,2,y2);								//dR/dt=I/T_r
}


