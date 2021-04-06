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
		gsl_matrix_set_identity(V);
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


//Optimeret udgave af overstående//

void timesJ_op(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<p+1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		}
	for(int i=0;i<q+1;i++){
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}


void Jtimes_op(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int j=p;j<A->size2;j++){
		double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		}
	for(int j=q;j<A->size2;j++){
		double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,q,j,new_aqj);
		}
}


void jacobi_diag_op(gsl_matrix* A , gsl_matrix* V){
		int changed;
		gsl_matrix_set_identity(V);
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
						timesJ_op(A,p,q, theta);
						Jtimes_op(A,p,q,-theta);
						timesJ(V,p,q,theta);
						matrix_print("A midtvejs",A);
					}
				}
			}
		}while(changed!=0);
}





			// ----- Her kommer Main ----- //

int main(){
/*
printf("#index0:stof der ikke skal plottes\n");
// Data laves

int n=4;
gsl_matrix* A=gsl_matrix_alloc(n,n);
gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
gsl_matrix* V=gsl_matrix_alloc(n,n);
gsl_matrix* res1=gsl_matrix_alloc(n,n);
gsl_matrix* res12=gsl_matrix_alloc(n,n);
gsl_matrix* res2=gsl_matrix_alloc(n,n);
gsl_matrix* res22=gsl_matrix_alloc(n,n);
gsl_matrix* res3=gsl_matrix_alloc(n,n);
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

// Opgave B
int N=20;
double s=1.0/(N+1);
gsl_matrix* H = gsl_matrix_alloc(N,N);
gsl_matrix* V_h = gsl_matrix_alloc(N,N);

for(int i=0;i<N-1;i++){
  gsl_matrix_set(H,i,i,-2);
  gsl_matrix_set(H,i,i+1,1);
  gsl_matrix_set(H,i+1,i,1);
  }
gsl_matrix_set(H,n-1,n-1,-2);
gsl_matrix_scale(H,-1/s/s);

matrix_print("Min H:",H);
jacobi_diag(H,V_h);
matrix_print("H efter jacobi",H);
matrix_print("Mit V som er eigenfunktioner:",V_h);
printf("k: calculated vs exact (hbar^2/2mL^2=1)\n");
for (int k=0; k < N/3; k++){
	double exact = M_PI*M_PI*(k+1)*(k+1);
	double calculated = gsl_matrix_get(H,k+1,k+1);
   	printf("%i %g %g\n",k,calculated,exact);}

	// Opgave C Comparison//
printf("\n Nu opgave C hvor vi undersøger tiden af vores operationer\n");
gsl_matrix* a_10=gsl_matrix_alloc(25,25);
gsl_matrix* v_10=gsl_matrix_alloc(25,25);
make_rand_sym_matrix(a_10);
gsl_matrix* a_100=gsl_matrix_alloc(250,250);
gsl_matrix* v_100=gsl_matrix_alloc(250,250);
gsl_matrix* a_100_copy=gsl_matrix_alloc(250,250);
gsl_matrix* v_100_copy=gsl_matrix_alloc(250,250);
make_rand_sym_matrix(a_100);
gsl_matrix* a_gsl=gsl_matrix_alloc(250,250);
gsl_matrix* v_gsl=gsl_matrix_alloc(250,250);
gsl_vector* vv=gsl_vector_alloc(250);
gsl_matrix_memcpy(a_gsl,a_100);
gsl_matrix_memcpy(a_100_copy,a_100);

clock_t start, end;
double cpu_time_used_10, cpu_time_used_100, cpu_time_used_gsl, cpu_time_used_copy;

	start = clock();
	jacobi_diag(a_10,v_10);
	end = clock();
	cpu_time_used_10 = ((double) (end - start)) / CLOCKS_PER_SEC;

	start = clock();
	jacobi_diag(a_100,v_100);
	end = clock();
	cpu_time_used_100 = ((double) (end - start)) / CLOCKS_PER_SEC;

	start = clock();
	gsl_linalg_SV_decomp_jacobi(a_gsl,v_gsl,vv);
	end = clock();
	cpu_time_used_gsl = ((double) (end - start)) / CLOCKS_PER_SEC;

printf("Vi undersøger tiden det tager for vores routine at køre, først for n=100 så n=1000\n");
printf("Det tog %5g i CPU tid at køre jacobi diag for n=25\n",cpu_time_used_10);
printf("Det tog %5g i CPU tid at køre jacobi diag for n=250\n",cpu_time_used_100);
printf("Tiden delt med hinanden giver %5g\n",cpu_time_used_100/cpu_time_used_10);
printf("Det ses at forskellen er næsten 1000, hvilket er at forvente, hvis der bruges N^3 operatrioner\n\n");

printf("Nu sammenligner vi vores rotine med GSL. For den samme n=250 matrice tog den %6g i CPU tid.\n",cpu_time_used_gsl);
printf("Tiden i forhold til hinanden hvor vi har MinMetode/GSL giver %5g\n",cpu_time_used_100/cpu_time_used_gsl);
printf("Det ses at GSL stadig er hurtigere, så vi prøver nu at optimere vores kode, ved kun at opdatere den øvre halvdel af matricen vi rotere\n");
*/

gsl_matrix* test=gsl_matrix_alloc(6,6);
gsl_matrix* vtest=gsl_matrix_alloc(6,6);
make_rand_sym_matrix(test);

//	start = clock()
	jacobi_diag_op(test,vtest);
//	end = clock();
//	cpu_time_used_copy = ((double) (end - start)) / CLOCKS_PER_SEC;

/*
matrix_print("For god ordens skyld starter vi med at tjekke at vores optimerede giver det rigtige, altså 0'er over diogonalen som har egenværdierne",a_100_copy);
matrix_print("Den ikke optimerede gave følgene matrice og det ses at egenværdierne er ens",a_100);
printf("Tiden det tog for den optimerede var %6g\n",cpu_time_used_copy);
printf("Forskellen var gammel/op=%6g\n",cpu_time_used_100/cpu_time_used_copy);


printf("\n\n");


printf("#index1: numerical vs analytical\n");
 	for(int i=0;i<N;i++){
	double k=(i+1.0)/(N+1);
	printf("%6g %6g %6g %6g %6g %6g %6g\n",k, gsl_matrix_get(V_h,i,0), gsl_matrix_get(V_h,i,1)+1.5, gsl_matrix_get(V_h,i,2)+3, sqrt(2/M_PI)*sin(M_PI*k/2),-sqrt(2/M_PI)*sin(M_PI*2*k/2)+1.5,sqrt(2/M_PI)*sin(M_PI*3*k/2)+3);
	}
*/




//gsl_matrix_free(A);gsl_matrix_free(Acopy);gsl_matrix_free(V);gsl_matrix_free(res1);gsl_matrix_free(res2);gsl_matrix_free(res3),gsl_matrix_free(H),gsl_matrix_free(V_h);
return 0;
}

