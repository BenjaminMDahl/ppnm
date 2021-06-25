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


//ODE: Newtonian gravitational three-body problem(inspireret af eksemplet på blackboard)

void threebody(double t,gsl_vector* X, gsl_vector* dXdt){
	//Pakker startgættet ud og definere konstanter
	double M1=1,M2=1,M3=1,G=1;
	double x1,y1,x2,y2,x3,y3;
	x1=gsl_vector_get(X,0);y1=gsl_vector_get(X,1);
	x2=gsl_vector_get(X,2);y2=gsl_vector_get(X,3);
	x3=gsl_vector_get(X,4);y3=gsl_vector_get(X,5);
	double vx1,vy1,vx2,vy2,vx3,vy3;
	vx1=gsl_vector_get(X,6);vy1=gsl_vector_get(X,7);
	vx2=gsl_vector_get(X,8);vy2=gsl_vector_get(X,9);
	vx3=gsl_vector_get(X,10);vy3=gsl_vector_get(X,11);

	//Relative positioner/afstande
	double dx12=x2-x1, dy12=y2-y1;
	double dx13=x3-x1, dy13=y3-y1;
	double dx23=x3-x2, dy23=y3-y2;

	double r12=sqrt(dx12*dx12+dy12*dy12);
	double r13=sqrt(dx13*dx13+dy13*dy13);
	double r23=sqrt(dx23*dx23+dy23*dy23);

	//Fra newtons 2. lov
	double f12=M1*M2*G/r12/r12;
	double f13=M1*M3*G/r13/r13;
	double f23=M2*M3*G/r23/r23;

	//Fra opgave formuleringen
	double dvx1dt=1./M1*( f12*dx12/r12+f13*dx13/r13);
	double dvy1dt=1./M1*( f12*dy12/r12+f13*dy13/r13);
	double dvx2dt=1./M2*(-f12*dx12/r12+f23*dx23/r23);
	double dvy2dt=1./M2*(-f12*dy12/r12+f23*dy23/r23);
	double dvx3dt=1./M3*(-f13*dx13/r13-f23*dx23/r23);
	double dvy3dt=1./M3*(-f13*dy13/r13-f23*dy23/r23);

	//Laver andenordens diff.ligningen
	gsl_vector_set(dXdt,0,vx1); gsl_vector_set(dXdt,1,vy1);
	gsl_vector_set(dXdt,2,vx2); gsl_vector_set(dXdt,3,vy2);
	gsl_vector_set(dXdt,4,vx3); gsl_vector_set(dXdt,5,vy3);
	gsl_vector_set(dXdt,6,dvx1dt); gsl_vector_set(dXdt,7,dvy1dt);
	gsl_vector_set(dXdt,8,dvx2dt); gsl_vector_set(dXdt,9,dvy2dt);
	gsl_vector_set(dXdt,10,dvx3dt); gsl_vector_set(dXdt,11,dvy3dt);
}
