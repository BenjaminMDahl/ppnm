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

printf("#index 1:lambda_i     f(lambda_i)\n");
for(int i=0;i<x->size;i++){
	double xi=gsl_vector_get(x,i);
	double fi=plot_data(xi);
	printf("%10g %10g\n",xi,fi);}


printf("\n\n");
//Derefter til kurve
printf("#index 2:x_i    f(x_i)\n");
for(double i=-2;i<15;i=i+0.11){
	double fi=plot_data(i);
	printf("%10g %10g\n",i,fi);}

return 0;
}
