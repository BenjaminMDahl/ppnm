#include<stdio.h>
#include<time.h>
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
	printf("\n");
}

void matrix_print(char s[], gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0;i< A->size1 ;i++){							// Note til selv size1=vertical, size2=horisontal
		for(int j=0;j< A->size2 ;j++)printf("%10g ",gsl_matrix_get(A,i,j));
		printf("\n");}
	printf("\n");
}

void make_random_matrix(gsl_matrix* A){
	for(int i=0; i< A->size1; i++)
			for(int j=0; j<A->size2; j++)
			{
			double Aij=RND;
			gsl_matrix_set(A,i,j,Aij);
			}
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


void backsub(gsl_matrix* R, gsl_vector* c){
	for(int i=c->size-1; i>=0; i--){
		double s=gsl_vector_get(c,i);
		for(int k=i+1; k< R->size1; k++)s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));
	}

}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x); //x=Q^(T)*b
	backsub(R,x);
}


void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* I){
	int n=Q->size1; assert(Q->size1==Q->size2); //Tjekker at Q og dermed A er en square matrix
	gsl_vector* bi=gsl_vector_calloc(n);
	gsl_vector* xi=gsl_vector_alloc(n);

	for(int i=0;i<n;i++){
		gsl_vector_set(bi,i,1);
		if(i>0)gsl_vector_set(bi,i-1,0);
		GS_solve(Q,R,bi,xi);
		for(int j=0;j<n;j++)gsl_matrix_set(I,j,i,gsl_vector_get(xi,j));
	}
	gsl_vector_free(bi);gsl_vector_free(xi);
}

int main(){
/* jeg prøvede at bruge dette til integratinos homework
	gsl_matrix* Q=gsl_matrix_alloc(3,4);
	double a1=(double) 1/6, a2=(double) 2/6, a3=(double) 4/6, a4=(double) 5/6;
	gsl_matrix_set(Q,0,0,1);
	gsl_matrix_set(Q,0,1,1);
	gsl_matrix_set(Q,0,2,1);
	gsl_matrix_set(Q,0,3,1);
	gsl_matrix_set(Q,1,0,a1);
	gsl_matrix_set(Q,1,1,a2);
	gsl_matrix_set(Q,1,2,a3);
	gsl_matrix_set(Q,1,3,a4);
	gsl_matrix_set(Q,2,0,a1*a1);
	gsl_matrix_set(Q,2,1,a2*a2);
	gsl_matrix_set(Q,2,2,a3*a3);
	gsl_matrix_set(Q,2,3,a4*a4);

	gsl_matrix* Rq=gsl_matrix_alloc(4,4);
	gsl_vector* bq=gsl_vector_alloc(3);
	double b2=(double) 1/2, b3=(double) 1/3, b4=(double) 1/4;
	gsl_vector_set(bq,0,1);
	gsl_vector_set(bq,1,b2);
	gsl_vector_set(bq,2,b3);

	gsl_vector* xq=gsl_vector_alloc(4);
	matrix_print("Tjekker om matricen er rigtig",Q);
	vector_print("ligeså for b",bq);

	GS_decomp(Q,Rq);
	GS_solve(Q,Rq,bq,xq);
	vector_print("Og så tjekker vi resultatet",xq);
*/
	//Opgave A del 1) Gram-Schmidt orthogonalization//
	printf("Opgave A del 1)\n\n");

	//Data laves
	int n=5, m=3; assert(n>=m);
	gsl_matrix* A=gsl_matrix_alloc(n,m);
	gsl_matrix* A_pro=gsl_matrix_alloc(n,m);
	gsl_matrix* QtQ=gsl_matrix_alloc(m,m);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,m); //Jeg for brug for et copy da A skal blive til Q men jeg vil gerne kunne tjekke mit resultat
	gsl_matrix* R=gsl_matrix_alloc(m,m);

	// Vi danner tilfældige A og b hvor indgangene alle er mellem 0 og 1 ligesom i forlæsningen omkring gsl matricer
	make_random_matrix(A);
	gsl_matrix_memcpy(Acopy,A);
	GS_decomp(A,R);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1,A ,R,0, A_pro);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1,A ,A,0, QtQ);

	matrix_print("Vi printer A så vi har den til sammenligning",Acopy);
	printf("pause\n");
	matrix_print("Q som er den nye A",A);
	matrix_print("Det ses at R er upper triangular",R);
	matrix_print("Det ses at Q^(T)Q giver idenditet",QtQ);
	matrix_print("Det tjekkes om QR=A, som var den første matrix der blev printet",A_pro);

// Vi rydder op
gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(Acopy);
gsl_matrix_free(QtQ);
gsl_matrix_free(A_pro);


	//Opgave A del 2 ligning løsning//
	printf("Opgave A del 2)\n\n");

	//Data laves
	int N=5;
	gsl_matrix* B=gsl_matrix_alloc(N,N);
	gsl_matrix* Bcopy=gsl_matrix_alloc(N,N);
	gsl_matrix* Rb=gsl_matrix_alloc(N,N);
	gsl_vector* b=gsl_vector_alloc(N);
	gsl_vector* x=gsl_vector_alloc(N);
	gsl_vector* y=gsl_vector_alloc(N);
	make_random_matrix(B);
	gsl_matrix_memcpy(Bcopy,B);
	for(int i=0; i< b->size; i++)//Makes random vektor
		{
		double bi=RND;
		gsl_vector_set(b,i,bi);
		}
	GS_decomp(B,Rb);
	GS_solve(B,Rb,b,x);
	gsl_blas_dgemv(CblasNoTrans,1,Bcopy,x,0,y);

	printf("Tjekker om b og Ax er ens\n");
	vector_print("b var fra starten af",b);
	vector_print("A*x giver",y);
	printf("Det ses at de to er ens\n");

// Vi rydder op
gsl_matrix_free(B);
gsl_matrix_free(Bcopy);
gsl_matrix_free(Rb);
gsl_vector_free(y);
gsl_vector_free(b);
gsl_vector_free(x);


	//Opgave B invers//
	printf("Opgave B)\n\n");

	//Laver data
	N=4;
	gsl_matrix* C=gsl_matrix_alloc(N,N);
	gsl_matrix* I=gsl_matrix_alloc(N,N);
	gsl_matrix* Ccopy=gsl_matrix_alloc(N,N);
	gsl_matrix* Rc=gsl_matrix_alloc(N,N);
	gsl_matrix* res=gsl_matrix_alloc(N,N);
	make_random_matrix(C);
	gsl_matrix_memcpy(Ccopy,C);


	GS_decomp(C,Rc);
	GS_inverse(C,Rc,I);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Ccopy ,I,0, res);

	matrix_print("Her er matricen vi kigger på i opgave B",Ccopy);
	matrix_print("Her er dens udregnede invers",I);
	matrix_print("Her ses at A*B giver identitets matricen",res);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,I,Ccopy,0, res);
	matrix_print("Det samme fås for B*A så B=A^(-1)",res);

// Vi rydder op
gsl_matrix_free(C);
gsl_matrix_free(Ccopy);
gsl_matrix_free(Rc);
gsl_matrix_free(I);
gsl_matrix_free(res);


	//Opgave C time
	printf("Opgave C)\n\n");


	int N_small=100, N_large=500;
	gsl_matrix* d=gsl_matrix_alloc(N_small,N_small);
	gsl_matrix* rd=gsl_matrix_alloc(N_small,N_small);
	gsl_matrix* D=gsl_matrix_alloc(N_large,N_large);
	gsl_matrix* Dcopy=gsl_matrix_alloc(N_large,N_large);
	gsl_matrix* rD=gsl_matrix_alloc(N_large,N_large);
	gsl_vector* T=gsl_vector_alloc(N_large);
	make_random_matrix(d);
	make_random_matrix(D);
	gsl_matrix_memcpy(Dcopy,D);
	clock_t start, end;
	double cpu_time_used_lillen, cpu_time_used_stortn, cpu_time_used_stortn_gsl;

	start = clock();
	GS_decomp(d,rd);
	end = clock();
	cpu_time_used_lillen = ((double) (end - start)) / CLOCKS_PER_SEC;

	printf("tiden det tog for at decompose en n=100 sqaure matrix=\n");
	printf("%10g\n",cpu_time_used_lillen);

	start = clock();
	GS_decomp(D,rD);
	end = clock();
	cpu_time_used_stortn = ((double) (end - start)) / CLOCKS_PER_SEC;

	printf("tiden det tog for at decompose en n=500 sqaure matrix=\n");
	printf("%10g\n",cpu_time_used_stortn);
	printf("Tidsforskellen for tid(n=500)/tid(n=100)=\n");
	printf("%10g\n",cpu_time_used_stortn/cpu_time_used_lillen);
	printf("Det ses at når n bliver 5 gange større er tiden blevet ca 5^3=125 gange større som forventet\n");

	start = clock();
	gsl_linalg_QR_decomp(Dcopy,T);
	end = clock();
	cpu_time_used_stortn_gsl = ((double) (end - start)) / CLOCKS_PER_SEC;

	printf("tiden det tog for at decompose en n=500 sqaure matrix med GSL=\n");
	printf("%10g\n",cpu_time_used_stortn_gsl);
	printf("Tidsforskellen for minkode/GSL\n");
	printf("%10g\n",cpu_time_used_stortn/cpu_time_used_stortn_gsl);
	printf("Det ses at GSL's metode er noget hurtigere end min");


	gsl_matrix_free(d);
	gsl_matrix_free(rd);
	gsl_matrix_free(D);
	gsl_matrix_free(Dcopy);
	gsl_matrix_free(rD);
	gsl_vector_free(T);
return 0;}
// hvis du glemmer at free et alloc space er det ikke katastrofe inden for en main function, fordi den gør det selv ved exit af main functionen
