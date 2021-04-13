#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}

void matrix_print(char s[], gsl_matrix* A){
	int n=A->size1, m=A->size2;
	for(int i=0;i<n;i++){
		for(int j=0; j<m;j++){
			if(fabs(gsl_matrix_get(A,i,j))<10e-7)gsl_matrix_set(A,i,j,0);
		}
	}
	printf("%s\n",s);
	for(int i=0;i< n ;i++){							// Note til selv size1=vertical, size2=horisontal
		for(int j=0;j< m ;j++)printf("%10g ",gsl_matrix_get(A,i,j));
		printf("\n");}
	printf("\n");
}


void make_rand_sym_matrix(gsl_matrix* A){
	for(int i=0; i< A->size1; i++){
		double Aii=RND;
		gsl_matrix_set(A,i,i,Aii);
		for(int j=i+1; j<A->size2; j++){
			double Asym=RND;
			gsl_matrix_set(A,i,j,Asym);
			gsl_matrix_set(A,j,i,Asym);
			}
	}
}
