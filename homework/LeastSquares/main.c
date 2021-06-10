#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<assert.h>


void vector_print(char s[], gsl_vector* v);
void matrix_print(char s[], gsl_matrix* A);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void backsub(gsl_matrix* R, gsl_vector* c);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* I);
void least_square(int n, int m, double* x,double* y, double* dy,double f(int,double),gsl_vector* c,gsl_matrix* dc);



double fun(int i, double x){
	switch(i){
		case  0: return 1  ; break;
		case  1: return -x  ; break;
		default: return NAN;
		}
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
