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

void secular_equation_default(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps);

int main(){

clock_t start, end;
//double time_5,time_10,time_100,time_1000;
double time_5;

int n=5;
gsl_vector* D=gsl_vector_alloc(n);
gsl_vector* u=gsl_vector_alloc(n);
int p=0; //Overvej lige noget index her
double tal[]={1,1,3,4,5,6,7,8,9,10};

for(int i=0;i<n;i++){
	gsl_vector_set(D,i,tal[i]);
	gsl_vector_set(u,i,tal[i]);}
//gsl_vector_scale(u,1/sqrt(55));
//gsl_matrix_scale(D,1/sqrt(55));


gsl_vector* x=gsl_vector_alloc(n);

start = clock();

	secular_equation_default(D,u,p,x,0.0000001);
	vector_print("res er",x);

end = clock();
time_5= ((double) (end - start)) / CLOCKS_PER_SEC;
printf("det tog: \n %g \n",time_5);

// Prøver igen
/*
int m=10;
gsl_vector* d=gsl_vector_alloc(m);
gsl_vector* U=gsl_vector_alloc(m);
gsl_vector* xb=gsl_vector_alloc(m);
int P=1; //Overvej lige noget index her

for(int i=0;i<m;i++){
	gsl_vector_set(d,i,tal[i]);
	gsl_vector_set(U,i,tal[i]);}

//double b[]={1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5};

//start = clock();


secular_equation_solve(d,U,P,xb,0.0001);
vector_print("res er",xb);

//end = clock();
//time_10= ((double) (end - start)) / CLOCKS_PER_SEC;
*/

/*
// Prøver igen

int l=100;
gsl_matrix* D100=gsl_matrix_alloc(l,l);
gsl_vector* u100=gsl_vector_alloc(l);
int p100=25; //Overvej lige noget index her
double c[100];

for(int i=0;i<l;i++){
	gsl_matrix_set(D100,i,i,RND);
	gsl_vector_set(u100,i,RND);
	c[i]=RND;}

start = clock();
for(int i=0;i<l;i++){
	gsl_vector_set(x,0,c[i]);
	secular_equation_solve(D100,u100,p100,x,0.0001);
//	vector_print("res er",x);
	}
end = clock();

time_100= ((double) (end - start)) / CLOCKS_PER_SEC;

// prøver igen igen nu med 1000
*/
int L=1000;
gsl_vector* D1000=gsl_vector_alloc(L);
gsl_vector* u1000=gsl_vector_alloc(L);
gsl_vector* x1000=gsl_vector_alloc(L);

for(int i=0;i<L;i++){
	gsl_vector_set(D1000,i,RND);
	gsl_vector_set(u1000,i,RND);}

//start = clock();
	secular_equation_default(D1000,u1000,5,x1000,0.0001);
	vector_print("res er",x1000);

//end = clock();

//time_1000= ((double) (end - start)) / CLOCKS_PER_SEC;


// Tid
//printf("Tid1000/tid100 :\n %g \n",time_1000/time_100);

return 0;
}
