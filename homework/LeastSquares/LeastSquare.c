#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<assert.h>

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void backsub(gsl_matrix* R, gsl_vector* c);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* I);

void least_square(int n, int m, double* x,double* y, double* dy,double f(int,double),gsl_vector* c,gsl_matrix* dc){
	//makes A=f_k(x)//
	gsl_matrix* A=gsl_matrix_alloc(n,m);
	gsl_matrix* R=gsl_matrix_alloc(m,m);
	for(int i=0; i< A->size1; i++)
		for(int k=0; k<A->size2; k++)
		{
		double Aik=f(k,x[i])/dy[i];
			gsl_matrix_set(A,i,k,Aik);
		}
	//makes b=y/dy//
	gsl_vector* b=gsl_vector_alloc(n);
	for(int i=0; i< b->size; i++){
		double bi=y[i]/dy[i];
		gsl_vector_set(b,i,bi);}
	GS_decomp(A,R);
	GS_solve(A,R,b,c);

	// inverser R og gemmer i dc til variance//
	gsl_matrix* Rr=gsl_matrix_alloc(m,m);
	GS_decomp(R,Rr);
	GS_inverse(R,Rr,dc);

	gsl_matrix_free(A);gsl_matrix_free(R);gsl_vector_free(b); gsl_matrix_free(Rr);
}
