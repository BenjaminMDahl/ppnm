#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<assert.h>


void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}

void matrix_print(char s[], gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0;i< A->size1 ;i++){							// Note til selv size1=vertical, size2=horisontal
		for(int j=0;j< A->size2 ;j++)printf("%10g ",gsl_matrix_get(A,i,j));
		printf("\n");}
	printf("\n");
}

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
		for(int k=i+1; k< R->size2; k++)s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));
	}

}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x); //x=Q^(T)*b
	backsub(R,x);
}

double fun(int i, double x){
	switch(i){
		case  0: return 1  ; break;
		case  1: return x  ; break;
		default: return NAN;
		}
	}

void least_square(int n, int m, double* x,double* y, double* dy,double f(int,double),gsl_vector* c){
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
	gsl_matrix_free(A);gsl_matrix_free(R);gsl_vector_free(b);
}

int main(){
	double x[]={1,2,3,4,6,9,10,13,15};
	double y[]={117,100,88,72,53,29.5,25.2,15.2,11.1};
	int n=sizeof x / sizeof x[0];
	double dy[n],ln_y[n],dln_y[n];
	for(int i=0;i<n;i++){
	dy[i]=y[i]/20;
	ln_y[i]=log(y[i]);
	dln_y[i]=dy[i]/y[i];}
	printf("data:x y log(y) dy dlog(y)\n");
	for(int i=0;i<n;i++)printf("%6g %6g %6g %6g %6g\n",x[i],y[i],ln_y[i],dy[i],dln_y[i]);

	int m=2; //antal funktioner
	gsl_vector* c=gsl_vector_alloc(m);
	double (*f)(int,double);f=fun;		 //jeg kan ikke bare give least square fun direkte fra main, der skal være noget i main som peger på fun jeg kan give videre
	least_square(n,m,x,ln_y,dln_y,f,c);
	vector_print("her er mine c'er",c);

return 0;
}
