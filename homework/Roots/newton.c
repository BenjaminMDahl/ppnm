#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<assert.h>

void GS_solve(gsl_matrix* A, gsl_matrix* R, gsl_vector* b, gsl_vector* x);


void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){

	// gør diverse gsl matricer og vektorer klar
	int n=x->size;double dx=sqrt(eps);
	gsl_matrix* J=gsl_matrix_alloc(n,n);
	gsl_matrix* R=gsl_matrix_alloc(n,n);
	gsl_vector* fx=gsl_vector_alloc(n);
	gsl_vector* dfx=gsl_vector_alloc(n);
	gsl_vector* fy=gsl_vector_alloc(n);
	gsl_vector* xstep=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_calloc(n);

	// Herunder dannes jacobianen jævnført ligning (7) i kapitlet om nonlinear equations og min GS_solver kaldes på den
	while(1){
		f(x,fx);
		for(int i=0;i<n;i++){
			double xi=gsl_vector_get(x,i);
			gsl_vector_set(x,i,xi+dx);
			f(x,dfx);
			gsl_vector_sub(dfx,fx);
			gsl_vector_scale(dfx,1/dx);
			for(int j=0;j<n;j++)gsl_matrix_set(J,j,i,gsl_vector_get(dfx,j));
			double dxi=gsl_vector_get(x,i);
			gsl_vector_set(x,i,dxi-dx);}
		gsl_vector_scale(fx,-1);//-f(x)
		GS_solve(J,R,fx,xstep);

		// Her tjekkes om det funde x er tilfredsstillende
		double s=2;
		gsl_vector_memcpy(y,x);
		while(1){
			s/=2;
			gsl_vector_scale(xstep,s);
			gsl_vector_add(y,xstep);
			gsl_vector_scale(xstep,1/s);
			f(y,fy);
//		printf("fy=%g\n",gsl_blas_dnrm2(fy));
//		printf("faktorfx=%g\n",(1-s/2)*gsl_blas_dnrm2(fx));
			if(gsl_blas_dnrm2(fy)<(1-s/2)*gsl_blas_dnrm2(fx))break;
//		printf("s=%g\n",s);
			if(s<0.01)break;
		}
		gsl_vector_memcpy(x,y);gsl_vector_memcpy(fx,fy);
//	printf("step norm=%g\n",gsl_blas_dnrm2(xstep));
		if(gsl_blas_dnrm2(xstep)<dx) break;
//	printf("value norm=%g\n",gsl_blas_dnrm2(fx));
		if(gsl_blas_dnrm2(fx)<eps) break;
	}
	gsl_matrix_free(J);gsl_matrix_free(R);
	gsl_vector_free(fx);gsl_vector_free(dfx);gsl_vector_free(y);gsl_vector_free(fy);
}
