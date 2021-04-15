#include<math.h>
#include<assert.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>

//Vi tager udgangspunkt i et åbent interval og bruger 4 punkter, så vi kan bruge trapezium og rectangle rule fra bogen (vægtene fra sætning (49) og (50) i integrations kapitelet)

double quad(double f(double), double a, double b, double delta, double eps, double f2, double f3, int nrec, int *record)
{
(*record)++;
assert(nrec<1000000);// tjekker at vi ikke får for mange recursions
double f1=(double) f(a+(b-a)/6);
double f4=(double) f(a+5*(b-a)/6);

double Q = (2*f1+f2+f3+2*f4)*(b-a)/6; // Baseret på (49) fra bogen kapitel Numerical integration
double q = (f1+f2+f3+f4)*(b-a)/4; // Baseret på (50) fra bogen kapitel numerical integration
double err=fabs(Q-q);
if (err < delta+eps*fabs(Q)) return Q;
else return quad(f,a,(a+b)/2,delta/sqrt(2),eps,f1,f2,nrec+1,record)+quad(f,(a+b)/2,b,delta/sqrt(2),eps,f3,f4,nrec+1,record);
}


double integrate(double f(double), double a, double b, double delta, double eps, int *record) // Dette første trin er for at kunne genbruge punkterne løbende
{
	(*record)=0;
	int nrec=0;
	double f2=(double) f(a+2*(b-a)/6);
	double f3=(double) f(a+4*(b-a)/6);
	return quad(f,a,b,delta,eps,f2,f3,nrec,record);
}


	// her prøver vi med transformation //

/*


void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);


double quad_CC(double f(double), double a, double b, double delta, double eps, double f2, double f3, int nrec)
{
assert(nrec<1000000);

// Vægtene bestemmes via min egen linsolve funktion
gsl_matrix* A=gsl_matrix_alloc(4,4);
gsl_matrix* R=gsl_matrix_alloc(4,4);

double a1=(double) (a+b)/2+(b-a)/2*cos(M_PI/6), a2=(double) (a+b)/2+(b-a)/2*cos(M_PI/6*2), a3=(double) (a+b)/2+(b-a)/2*cos(M_PI/6*4), a4=(double) (a+b)/2+(b-a)/2*cos(M_PI/6*5);
gsl_matrix_set(A,0,0,sin(M_PI/6)), gsl_matrix_set(A,0,1,sin(M_PI/6*2)), gsl_matrix_set(A,0,2,sin(M_PI/6*4)), gsl_matrix_set(A,0,3,sin(M_PI/6*5));
gsl_matrix_set(A,1,0,sin(M_PI/6)*a1), gsl_matrix_set(A,1,1,sin(M_PI/6*2)*a2), gsl_matrix_set(A,1,2,sin(M_PI/6*4)*a3), gsl_matrix_set(A,1,3,sin(M_PI/6*5)*a4);
gsl_matrix_set(A,2,0,sin(M_PI/6)*a1*a1), gsl_matrix_set(A,2,1,sin(M_PI/6*2)*a2*a2), gsl_matrix_set(A,2,2,sin(M_PI/6*4)*a3*a3), gsl_matrix_set(A,2,3,sin(M_PI/6*5)*a4*a4);
gsl_matrix_set(A,3,0,sin(M_PI/6)*a1*a1*a1), gsl_matrix_set(A,3,1,sin(M_PI/6*2)*a2*a2*a2), gsl_matrix_set(A,3,2,sin(M_PI/6*4)*a3*a3*a3), gsl_matrix_set(A,3,3,sin(M_PI/6*5)*a4*a4*a4);

gsl_vector* B=gsl_vector_alloc(4);
gsl_vector* w=gsl_vector_alloc(4);

gsl_vector_set(B,0,2), gsl_vector_set(B,1,a+b), gsl_vector_set(B,2,2/3*(a*a+a*b+b*b)), gsl_vector_set(B,3,1/2*(a+b)*(a*a+b*b));

GS_decomp(A,R);
GS_solve(A,R,B,w);
double w1=gsl_vector_get(w,0), w2=gsl_vector_get(w,1), w3=gsl_vector_get(w,2), w4=gsl_vector_get(w,3);

gsl_matrix_free(A);gsl_matrix_free(R);gsl_vector_free(B);gsl_vector_free(w);

	double f1=(double) f((a+b)/2+(b-a)/2*cos(M_PI/6))*sin(M_PI/6)*(b-a)/2;
	double f4=(double) f((a+b)/2+(b-a)/2*cos(M_PI/6*5))*sin(M_PI/6*5)*(b-a)/2;

double Q = (w1*f1+w2*f2+w3*f3+w4*f4); double q = (f1+f2+f3+f4)*M_PI/4;
double err=fabs(Q-q);
if (err < delta+eps*fabs(Q)) return Q;
else return quad_CC(f,a,(a+b)/2,delta/sqrt(2),eps,f1,f2,nrec+1)+quad_CC(f,(a+b)/2,b,delta/sqrt(2),eps,f3,f4,nrec+1);
}


double integrate_CC(double f(double), double a, double b, double delta, double eps)
{
	assert(b>a);
	int nrec=0;
	double f2=(double) f((a+b)/2+(b-a)/2*cos(M_PI/6*2))*sin(M_PI/6*2)*(b-a)/2;
	printf("f2=%g\n",f2);
	double f3=(double) f((a+b)/2+(b-a)/2*cos(M_PI/6*4))*sin(M_PI/6*4)*(b-a)/2;
	return quad_CC(f,a,b,delta,eps,f2,f3,nrec);
}




*/
