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
			if(fabs(gsl_matrix_get(A,i,j))<10e-11)gsl_matrix_set(A,i,j,0);
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


// Fra opgaveformulering//
void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}


void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size2;i++){
		double new_api= c*gsl_matrix_get(A,p,i)+s*gsl_matrix_get(A,q,i);
		double new_aqi=-s*gsl_matrix_get(A,p,i)+c*gsl_matrix_get(A,q,i);
		gsl_matrix_set(A,p,i,new_api);
		gsl_matrix_set(A,q,i,new_aqi);
		}
}


void jacobi_diag(gsl_matrix* A , gsl_matrix* V){
		int changed;
		do{
			changed=0;
			for(int p=0;p<A->size1-1;p++){
			for(int q=p+1;q<A->size2;q++){
			double apq=gsl_matrix_get(A,p,q);
			double app=gsl_matrix_get(A,p,p);
			double aqq=gsl_matrix_get(A,q,q);
			double theta=0.5*atan2(2*apq,aqq-app);
			double c=cos(theta),s=sin(theta);
			double new_app=c*c*app-2*s*c*apq+s*s*aqq;
			double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
			if(new_app!=app || new_aqq!=aqq){         // Hvis C ikke kan se forskel på det nye og gamle diagonal element skal den ikke rotere igen, ellers skal den opdatere
				changed=1;
				timesJ(A,p,q, theta);
				Jtimes(A,p,q,-theta);
				timesJ(V,p,q,theta);		// Denne er til for også at få egenvektor matricen ud
		}
	}
		}

		}while(changed!=0);
}


			// ----- Her kommer Main ----- //

int main(){

// Data laves
int n=2;
gsl_matrix* A=gsl_matrix_alloc(n,n);
gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
gsl_matrix* V=gsl_matrix_alloc(n,n);
gsl_matrix* res1=gsl_matrix_alloc(n,n);
gsl_matrix* res12=gsl_matrix_alloc(n,n);
gsl_matrix* res2=gsl_matrix_alloc(n,n);
gsl_matrix* res22=gsl_matrix_alloc(n,n);
gsl_matrix* res3=gsl_matrix_alloc(n,n);
gsl_matrix_set_identity(V);
make_rand_sym_matrix(A);
gsl_matrix_memcpy(Acopy,A);


// Der laves matrix produkter og data printes
matrix_print("Min symmetriske tilfældige matrice",Acopy);
jacobi_diag(A,V);
matrix_print("Min A efter Jacobi algorithmen, hvilket skulles svare til D",A);
matrix_print("Min V matrice med egenvektorer",V);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,Acopy,0, res1);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,res1,V,0, res12);
matrix_print("Vi udregner V^(T)AV og ser at det give vores D igen:",res12);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,A,0, res2);
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,res2,V,0, res22);
matrix_print("Vi udregner VDV^(T) og ser at vi får A igen",res22);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,V,0,res3);
matrix_print("Vi udregner V^(T)V og ser at vi får en identitets matrice ud",res3);



gsl_matrix_free(A);gsl_matrix_free(Acopy);gsl_matrix_free(V);gsl_matrix_free(res1);gsl_matrix_free(res2);gsl_matrix_free(res3);

return 0;
}
