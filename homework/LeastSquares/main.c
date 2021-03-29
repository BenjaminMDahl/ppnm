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
		case  1: return -x  ; break;
		default: return NAN;
		}
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

int main(){
	double x[]={1,2,3,4,6,9,10,13,15};
	double y[]={117,100,88,72,53,29.5,25.2,15.2,11.1};
	int n=sizeof x / sizeof x[0];
	double dy[n],ln_y[n],dln_y[n];
	for(int i=0;i<n;i++){
	dy[i]=y[i]/20;
	ln_y[i]=log(y[i]);
	dln_y[i]=dy[i]/y[i];}

	int m=2; //antal funktioner
	gsl_vector* c=gsl_vector_alloc(m);
	gsl_matrix* dc=gsl_matrix_alloc(m,m);
	double (*f)(int,double);f=fun;		 //jeg kan ikke bare give least square fun direkte fra main, der skal være noget i main som peger på fun jeg kan give videre
	least_square(n,m,x,ln_y,dln_y,f,c,dc);
	printf("#index -1: Resultater som ikke skal plottes\n");
	vector_print("her er mine c'er",c);
	printf("Fra overstående c[1] får vi følgende halverignstid\n 224Ra_1/2=%7g \n\n",log(2)/gsl_vector_get(c,1));
	matrix_print("Her fås min covariance matrice for C'erne",dc);
	printf("Fra overstående kan man udlede at usikkerheden på min funde halveringstid til\n sigma_224Ra_1/2=%7g \n",sqrt(gsl_matrix_get(dc,1,1))/(gsl_vector_get(c,1)*gsl_vector_get(c,1)));
	printf("Det ses at usikkerheden er ret stor, og den man er sikker på i dag (3.627 af hvad jeg har kunne finde på nettet) er inden for de usikkerheder.\n");
	printf("\n\n");


	printf("#index 0: data:x y log(y) dy dlog(y)\n");
	for(int i=0;i<n;i++)printf("%6g %6g %6g %6g %6g\n",x[i],y[i],ln_y[i],dy[i],dln_y[i]);
	printf("\n\n");

	printf("#index 1: x-fit, y-fit, y-fit+sigma, y-fit-sigma, y-fit+2*sigma, y-fit-2*sigma\n");
	int N=200; double x_fit[N],y_fit[N],y_fitp1[N],y_fitm1[N],y_fitp2[N],y_fitm2[N];
	for(int i=0;i<N;i++){
	x_fit[i]=(double)(i*15.5)/N;
	y_fit[i]=gsl_vector_get(c,0)-x_fit[i]*gsl_vector_get(c,1);
	y_fitp1[i]=gsl_vector_get(c,0)+gsl_matrix_get(dc,0,0)-x_fit[i]*(gsl_vector_get(c,1)+gsl_matrix_get(dc,1,1));
	y_fitm1[i]=gsl_vector_get(c,0)-gsl_matrix_get(dc,0,0)-x_fit[i]*(gsl_vector_get(c,1)-gsl_matrix_get(dc,1,1));
	y_fitp2[i]=gsl_vector_get(c,0)+2*gsl_matrix_get(dc,0,0)-x_fit[i]*(gsl_vector_get(c,1)+2*gsl_matrix_get(dc,1,1));
	y_fitm2[i]=gsl_vector_get(c,0)-2*gsl_matrix_get(dc,0,0)-x_fit[i]*(gsl_vector_get(c,1)-2*gsl_matrix_get(dc,1,1));
	printf("%7g %7g %7g %7g %7g %7g\n",x_fit[i],y_fit[i],y_fitp1[i],y_fitm1[i],y_fitp2[i],y_fitm2[i]);}
	printf("\n\n");



return 0;
}
