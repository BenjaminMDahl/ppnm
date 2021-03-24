#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<assert.h>
#define RND (double)rand()/RAND_MAX

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));

}

void matrix_print(char s[], gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0;i< A->size1 ;i++){							// Note til selv size1=vertical, size2=horisontal
		for(int j=0;j< A->size2 ;j++)printf("%10g ",gsl_matrix_get(A,i,j));
		printf("\n");}
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



int main(){
	//Data laves //


	int n=5, m=3; assert(n>=m);
	gsl_matrix* A=gsl_matrix_alloc(n,m);
	gsl_matrix* A_pro=gsl_matrix_alloc(n,m);
	gsl_matrix* QtQ=gsl_matrix_alloc(m,m);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,m); //Jeg for brug for et copy da A skal blive til Q men jeg vil gerne kunne tjekke mit resultat
	gsl_matrix* R=gsl_matrix_alloc(m,m);
	gsl_vector* b=gsl_vector_alloc(n);
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_calloc(n);	// callocs fylder ud med 0'er hvor alloc ikke gør som er hurtigere.

	// Vi danner tilfældige A og b hvor indgangene alle er mellem 0 og 1 ligesom i forlæsningen omkring gsl matricer
	for(int i=0; i< A->size1; i++)
		for(int j=0; j<A->size2; j++)
		{
		double Aij=RND;
		gsl_matrix_set(A,i,j,Aij);
		}
	gsl_matrix_memcpy(Acopy,A);
	for(int i=0; i< b->size; i++)
		{
		double bi=RND;
		gsl_vector_set(b,i,bi);
		}




	//Opgave A del 1) Gram-Schmidt orthogonalization//

	GS_decomp(A,R);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1,A ,R,0, A_pro);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1,A ,A,0, QtQ);

	matrix_print("Vi printer A så vi har den til sammenligning",Acopy);
	printf("pause\n");
	matrix_print("Q som er den nye A",A);
	matrix_print("Det ses at R er upper triangular",R);
	matrix_print("Det ses at Q^(T)Q giver idenditet",QtQ);
	matrix_print("Det tjekkes om QR=A, som var den første matrix der blev printet",A_pro);




	//Opgave A del 2 ligning løsning//




// Vi rydder op//
gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(Acopy);
gsl_matrix_free(QtQ);
gsl_matrix_free(A_pro);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(y);
return 0;}
// hvis du glemmer at free et alloc space er det ikke katastrofe inden for en main function, fordi den gør det selv ved exit af main functionen
