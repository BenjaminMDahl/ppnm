#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_sort_vector.h>
#include<assert.h>

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);

void secular_equation_guess(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps){
assert(D->size>1); // Vi arbejder kun med matricer ikke skalar

	// Vi starter med at opdatere D, så hvis den p'te indgang i u ikke er 0, updateres D og u sættes lig 0.
double dstart=gsl_vector_get(D,p);
double ustart=gsl_vector_get(u,p);
gsl_vector_set(D,p,dstart+2*ustart);
gsl_vector_set(u,p,0);

	// Det karakteristiske polynomium dannes
	void Se(gsl_vector* v, gsl_vector* f){
		double x=gsl_vector_get(v,0);
		double S=x-gsl_vector_get(D,p);
		double P=1;
		for(int k=0;k<p;k++){
			double dk=gsl_vector_get(D,k);
			double uk=gsl_vector_get(u,k);
			S=S+uk*uk/(dk-x);
			P=P*(dk-x);}
		for(int k=p+1;k<u->size;k++){
			double dk=gsl_vector_get(D,k);
			double uk=gsl_vector_get(u,k);
			S=S+uk*uk/(dk-x);
			P=P*(dk-x);}
		gsl_vector_set(f,0,S*P);}


	// Da vi bruger vores newton metode lavet i kurset, skal alt leveres i vektorer, da metoden er lavet til at kunne klare problemer af højere dimension.
gsl_vector* g=gsl_vector_alloc(1);
for(int i=0;i<x->size;i++){
		double xi=gsl_vector_get(x,i);
		gsl_vector_set(g,0,xi);
		newton(Se,g,eps);
		xi=gsl_vector_get(g,0);
		gsl_vector_set(x,i,xi);}

gsl_vector_free(g);
}

