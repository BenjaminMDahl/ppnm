#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>

void matrix_print(char s[], gsl_matrix* A);

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



