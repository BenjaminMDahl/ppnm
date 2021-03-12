#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>


void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	}

int main(){
	int n=3;
	double Adata[] = {6.13,-2.90,5.86,8.08,-6.31,-3.89,-4.36,1.00,0.19};
	double bdata[] = {6.23,5.37,2.29};
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
	gsl_vector* b=gsl_vector_alloc(n);
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_calloc(n);	// callocs fylder ud med 0'er hvor alloc ikke gør som er hurtigere.

	for(int i=0; i< A->size1; i++) 		//size1=vertical size of A
		for(int j=0; j < A->size2; j++) //size2=horizontal size of A
		{
		double Aij=Adata[i*A->size1+j];
		gsl_matrix_set(A,i,j,Aij); 	// set definere hvad der skal ske i A på Aij indgangen.
		}
	gsl_matrix_memcpy(Acopy,A); 		//vi skal bruge et kopi af A fordi gsl_solve "ødelægger" A efter brug.
	for(int i=0; i< b->size; i++)	 	// vectors only have 1 size
		{
		double bi=bdata[i];
		gsl_vector_set(b,i,bi);
		}
	gsl_linalg_HH_solve(Acopy,b,x);	 	//Ax=b
	gsl_blas_dgemv(CblasNoTrans, 1 , A, x, 0 , y);
	vector_print("Our b vector:",b);
	vector_print("Check: A*x should be equal b, where x comes from HH_solve:",y);

gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(y);
return 0;
}
// hvis du glemmer at free et alloc space er det ikke katastrofe inden for en main function, fordi den gør det selv ved exit af main functionen
