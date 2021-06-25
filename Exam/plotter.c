#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_eigen.h>
#include<assert.h>
#define RND (double)rand()/RAND_MAX*10


void make_rnd_vector(gsl_vector* v);
void secular_equation_default(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps,int m);

double plot_data(double x){
	double S=x-9+1/(1-x)+4/(2-x)+16/(4-x)+25/(5-x);
return S;
}

int main(){


	////plotning af egenværdier////
int n=5;
gsl_vector* D=gsl_vector_alloc(n);
gsl_vector* u=gsl_vector_alloc(n);
gsl_vector* x=gsl_vector_alloc(n);
int p=2;
double tal[]={1,2,3,4,5};
for(int i=0;i<n;i++){
	gsl_vector_set(D,i,tal[i]);
	gsl_vector_set(u,i,tal[i]);}
	secular_equation_default(D,u,p,x,0.00001,0);

// Først dataene til resultatet fra vores funktion
printf("#index 1:lambda_i     f(lambda_i)\n");
for(int i=0;i<x->size;i++){
	double xi=gsl_vector_get(x,i);
	double fi=plot_data(xi);
	printf("%10g %10g\n",xi,fi);}


printf("\n\n");
//Derefter data til at lave det karakteristiske polynomium
printf("#index 2:x_i    f(x_i)\n");
for(double i=-2;i<15;i=i+0.034){
	double fi=plot_data(i);
	printf("%10g %10g\n",i,fi);}
printf("\n\n");

gsl_vector_free(D);gsl_vector_free(u);gsl_vector_free(x);

	////plotning af effektiviteten af default som funktion af m for tre matricer med n=1000 og større og større indgange.////

printf("#index 3: Data for default effektivitet som funktion af m\n");
int L=1000;
gsl_vector* D1001=gsl_vector_alloc(L);
gsl_vector* u1001=gsl_vector_alloc(L);
gsl_vector* x1001=gsl_vector_alloc(L);
gsl_vector* D1005=gsl_vector_alloc(L);
gsl_vector* u1005=gsl_vector_alloc(L);
gsl_vector* x1005=gsl_vector_alloc(L);
gsl_vector* D1100=gsl_vector_alloc(L);
gsl_vector* u1100=gsl_vector_alloc(L);
gsl_vector* x1100=gsl_vector_alloc(L);
for(int i=0;i<L;i++){
	gsl_vector_set(D1001,i,RND);
	gsl_vector_set(u1001,i,RND);
	gsl_vector_set(D1005,i,RND*5);
	gsl_vector_set(u1005,i,RND*5);
	gsl_vector_set(D1100,i,RND*100);
	gsl_vector_set(u1100,i,RND*100);}
for(int i=0;i<=100;i=i+5){
	int eig1=0;int eig5=0;int eig100=0;
	secular_equation_default(D1001,u1001,25,x1001,0.0001,i);
	secular_equation_default(D1005,u1005,25,x1005,0.0001,i);
	secular_equation_default(D1100,u1100,25,x1100,0.0001,i);
	for(int j=0;j<L;j++){
		if(isnan(gsl_vector_get(x1001,j)))eig1++;
		if(isnan(gsl_vector_get(x1005,j)))eig5++;
		if(isnan(gsl_vector_get(x1100,j)))eig100++;}
	printf("%4i %4i %4i %4i\n",i,L-eig1,L-eig5,L-eig100);
	}


gsl_vector_free(D1001);gsl_vector_free(u1001);gsl_vector_free(x1001);
gsl_vector_free(D1005);gsl_vector_free(u1005);gsl_vector_free(x1005);
gsl_vector_free(D1100);gsl_vector_free(u1100);gsl_vector_free(x1100);

return 0;
}
