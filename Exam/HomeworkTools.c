#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<assert.h>

void GS_decomp(gsl_matrix* A, gsl_matrix* R){
	int i,j,k;
	//Vi laver A om til Q og Danner R //
	for(i=0; i< A->size2; i++){
		double lenAi=0;
		for(j=0;j<A->size1; j++){  //Finder indreprodukt af søjle ai med sig selv
			double ai=gsl_matrix_get(A,j,i);
			lenAi+=ai*ai;
		}
		for(j=0; j<A->size1; j++){
			gsl_matrix_set(A,j,i,gsl_matrix_get(A,j,i)/sqrt(lenAi));
		}
		gsl_matrix_set(R,i,i,sqrt(lenAi));
		for(j=i+1;j<A->size2; j++){
			double qiaj=0;
			for(k=0; k<A->size1; k++){  //Finder indreprodukt mellem ak og nye ai(qi)
				double qi=gsl_matrix_get(A,k,i);
				double aj=gsl_matrix_get(A,k,j);
				qiaj+=qi*aj;
			}
			for(k=0; k<A->size1; k++){
				gsl_matrix_set(A,k,j,gsl_matrix_get(A,k,j)-qiaj*gsl_matrix_get(A,k,i));
			}
			gsl_matrix_set(R,i,j,qiaj);
		}
	}
}


void backsub(gsl_matrix* R, gsl_vector* c){
	for(int i=c->size-1; i>=0; i--){
		double s=gsl_vector_get(c,i);
		for(int k=i+1; k< R->size1; k++)s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));
	}

}

void GS_solve(gsl_matrix* A, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	GS_decomp(A,R);
	gsl_blas_dgemv(CblasTrans,1,A,b,0,x); //x=Q^(T)*b
	backsub(R,x);
}


void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* I){
	int n=Q->size1; assert(Q->size1==Q->size2); //Tjekker at Q og dermed A er en square matrix
	gsl_vector* bi=gsl_vector_calloc(n);
	gsl_vector* xi=gsl_vector_alloc(n);

	for(int i=0;i<n;i++){
		gsl_vector_set(bi,i,1);
		if(i>0)gsl_vector_set(bi,i-1,0);
		GS_solve(Q,R,bi,xi);
		for(int j=0;j<n;j++)gsl_matrix_set(I,j,i,gsl_vector_get(xi,j));
	}
	gsl_vector_free(bi);gsl_vector_free(xi);
}

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
	int Nmax=0;
	while(1){
		Nmax++;
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

		// Her tjekkes om det fundne x er tilfredsstillende
		double s=2;
		gsl_vector_memcpy(y,x);
		while(1){
			s/=2;
			gsl_vector_scale(xstep,s);
			gsl_vector_add(y,xstep);
			gsl_vector_scale(xstep,1/s);
			f(y,fy);
			if(gsl_blas_dnrm2(fy)<(1-s/2)*gsl_blas_dnrm2(fx))break;
			if(s<0.01)break;
		}
		gsl_vector_memcpy(x,y);gsl_vector_memcpy(fx,fy);
		if(gsl_blas_dnrm2(xstep)<dx) break;
		if(gsl_blas_dnrm2(fx)<eps) break;
		if(Nmax>1000) break;
	}
	gsl_matrix_free(J);gsl_matrix_free(R);
	gsl_vector_free(fx);gsl_vector_free(dfx);gsl_vector_free(y);gsl_vector_free(fy);
}
