#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>


//Til opgave A//


void test_f1(gsl_vector* x,gsl_vector* f){
double x0=gsl_vector_get(x,0);
gsl_vector_set(f,0,x0*x0);

}

void test_f2(gsl_vector* x_vec,gsl_vector* f){
double x=gsl_vector_get(x_vec,0);
double y=gsl_vector_get(x_vec,1);
double z=gsl_vector_get(x_vec,2);

gsl_vector_set(f,0,x+z);
gsl_vector_set(f,1,y+z);
gsl_vector_set(f,2,z+x+y);
}

void test_f3(gsl_vector* x_vec,gsl_vector* f){
double x=gsl_vector_get(x_vec,0);
double y=gsl_vector_get(x_vec,1);


gsl_vector_set(f,0,(2-2*x+400*x*(y-x*x))*(-1));
gsl_vector_set(f,1,200*y-200*x*x);
}

//Til opgave B//

int OdeDriver(	void f(double,gsl_vector*,gsl_vector*), 		// right-hand-side of dy/dt=f(t,y)
	double a,			                     	// the start-point a
	double b,                     				// the end-point of the integration
	gsl_vector* ya,                    				// y(a)
	gsl_vector* yb,                     			// y(b) to be calculated
	double h,                     				// initial step-size
	double acc,                   				// absolute accuracy goal
	double eps,
	char *program);                    			// relative accuracy goal



static double E; //Energy
void diff(double r,gsl_vector* u,gsl_vector* dudr){
	gsl_vector_set(dudr,0,gsl_vector_get(u,1));
	gsl_vector_set(dudr,1,2*(-1/r-E)*gsl_vector_get(u,0));
}

void test_diff(gsl_vector* x,gsl_vector* f){
	E=gsl_vector_get(x,0);
	double a=0.00001,b=8,acc=0.00001,eps=0.00001,h=0.01;

	gsl_vector* ya=gsl_vector_alloc(2);
	gsl_vector* yb=gsl_vector_alloc(2);
	gsl_vector_set(ya,0,a-a*a);
	gsl_vector_set(ya,1,1-2*a);

	OdeDriver(diff,a,b,ya,yb,h,acc,eps,"kvant.txt");//problemet er at diff i ode ikke kan tage et 4. argument i form af energien, som skal kunne variere.
	gsl_vector_set(f,0,gsl_vector_get(yb,0));
	gsl_vector_free(ya);gsl_vector_free(yb);
}
