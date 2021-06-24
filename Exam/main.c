#include<stdio.h>
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

void make_rnd_vector(gsl_vector* v){
	int n=v->size;
	for(int i=0;i<n;i++)gsl_vector_set(v,i,RND);
}



void vector_print(char s[], gsl_vector* v);

void  matrix_print(char s[], gsl_matrix* A);

void secular_equation_default(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps,int m);

void secular_equation_GSL(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps);

void secular_equation_guess(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps);


int main(){


int n=19;
gsl_eigen_symm_workspace* w= gsl_eigen_symm_alloc(n);
gsl_vector* D=gsl_vector_alloc(n);
gsl_matrix* G=gsl_matrix_alloc(n,n);
gsl_vector* u=gsl_vector_alloc(n);
gsl_vector* x=gsl_vector_alloc(n);
gsl_vector* g=gsl_vector_alloc(n);
int p=2;
double tal[]={1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10};

for(int i=0;i<n;i++){
	gsl_vector_set(D,i,tal[i]);
	gsl_matrix_set(G,i,i,tal[i]);
	gsl_matrix_set(G,p,i,tal[i]);
	gsl_matrix_set(G,i,p,tal[i]);
	gsl_vector_set(u,i,tal[i]);
	gsl_vector_set(x,i,tal[i]+0.1);}
	gsl_matrix_set(G,p,p,tal[p]+2*tal[p]);
	matrix_print("Dette er den vi arbejder med",G);

	secular_equation_GSL(D,u,p,x,0.0000001);
	vector_print("res er",x);
//	secular_equation_guess(D,u,p,x,0.0001);
//	vector_print("med gæt",x);
	gsl_eigen_symm(G,g,w);
	gsl_sort_vector(g);
	vector_print("Fra GSL får vi",g);


//Tidstagning

clock_t start, stop;
double t1,t10;
n=200;

gsl_vector* D1=gsl_vector_alloc(n);
gsl_vector* u1=gsl_vector_alloc(n);
gsl_vector* x1=gsl_vector_alloc(n);
make_rnd_vector(D1);make_rnd_vector(u1);make_rnd_vector(x1);

start = clock();
secular_equation_default(D1,u1,5,x1,0.00001,0);
stop = clock();
t1 = ((double) (stop - start)) / CLOCKS_PER_SEC;


gsl_vector* D10=gsl_vector_alloc(n*10);
gsl_vector* u10=gsl_vector_alloc(n*10);
gsl_vector* x10=gsl_vector_alloc(n*10);
make_rnd_vector(D10);make_rnd_vector(u10);make_rnd_vector(x10);


start = clock();
secular_equation_default(D10,u10,5,x10,0.0001,0);
stop = clock();
t10 = ((double) (stop - start)) / CLOCKS_PER_SEC;

printf("Det tog for n=%i så lang tid: \n %g \n, og for n=%i så lang tid: \n %g \n de to delt med hinanden giver: \n %g \n Så en stigning på 10 i n gør det hele ca 100 gange langsommere, hvilket antyder at vi når i mål på O(n^2) operationer",n,t1,n*10,t10,t10/t1);

////////////////////////////////////

/* int L=1000;
gsl_vector* D1000=gsl_vector_alloc(L);
gsl_vector* u1000=gsl_vector_alloc(L);
gsl_vector* x1000=gsl_vector_alloc(L);

for(int i=0;i<L;i++){
	gsl_vector_set(D1000,i,RND);
	gsl_vector_set(u1000,i,RND);}

for(int i=500;i<3000;i=i+200){
	int eig=0;
	secular_equation_default(D1000,u1000,5,x1000,0.0001,i);
	for(int j=0;j<L;j++)if(isnan(gsl_vector_get(x1000,j)))eig++;
	printf("With %i guess' so many Eigenvalues were found: %i\n",i,L-eig);
	}*/


gsl_eigen_symm_free(w);
return 0;
}
