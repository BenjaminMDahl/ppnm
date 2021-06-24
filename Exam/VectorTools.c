#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#define RND (double)rand()/RAND_MAX

void make_rnd_vector(gsl_vector* v){
	int n=v->size;
	for(int i=0;i<n;i++)gsl_vector_set(v,i,RND);
}

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
	for(int i=0;i< n ;i++){
		for(int j=0;j< m ;j++)printf("%10g ",gsl_matrix_get(A,i,j));
		printf("\n");}
	printf("\n");
}

