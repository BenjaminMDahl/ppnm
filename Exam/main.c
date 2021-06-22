#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_sort_vector.h>
#include<assert.h>
#define RND (double)rand()/RAND_MAX



void vector_print(char s[], gsl_vector* v);

void secular_equation_default(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps,int m);

int main(){

int n=19;
gsl_vector* D=gsl_vector_alloc(n);
gsl_vector* u=gsl_vector_alloc(n);
gsl_vector* x=gsl_vector_alloc(n);
int p=2;
double tal[]={1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10};

for(int i=0;i<n;i++){
	gsl_vector_set(D,i,tal[i]);
	gsl_vector_set(u,i,tal[i]);}



	secular_equation_default(D,u,p,x,0.0000001,10);
	vector_print("res er",x);

////////////////////////////////////

int L=1000;
gsl_vector* D1000=gsl_vector_alloc(L);
gsl_vector* u1000=gsl_vector_alloc(L);
gsl_vector* x1000=gsl_vector_alloc(L);

for(int i=0;i<L;i++){
	gsl_vector_set(D1000,i,RND);
	gsl_vector_set(u1000,i,RND);}

for(int i=0;i<21;i++){
	int eig=0;
	secular_equation_default(D1000,u1000,5,x1000,0.0001,i);
	for(int j=0;j<L;j++)if(isnan(gsl_vector_get(x1000,j)))eig++;
	printf("Eigenvalues found: %i\n",L-eig);
	}

return 0;
}
