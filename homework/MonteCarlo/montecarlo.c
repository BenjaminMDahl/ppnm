#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<complex.h>
#define RND (double)rand()/RAND_MAX


// plain monte carlo//

double complex plainMC(int dim,double f(int dim,double* x),double* a,double* b,int N){
	double sum=0,sum2=0,x[dim];
        double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
        for(int i=0;i<N;i++){
                for(int i=0;i<dim;i++)x[i]=a[i]+RND*(b[i]-a[i]); //Danner et pseduo tilfældigt punkt
                double fx=f(dim,x); sum+=fx; sum2+=fx*fx;	//udregner funktions værdien i punkten og summere den op, også for f(x)^2 for at kunne estimere fejlen.
                }
        double mean=sum/N, sigma=sqrt(sum2/N-mean*mean);
        double complex result=mean*V+I*sigma*V/sqrt(N);
        return result;
}

// Quasi random monte carlo//

// denne funktion giver dig det tilsvarende corput tal ud fra et hel tal n, i basen b.

double corput(int n, int base){
	double q=0;
	double bk=(double)1/base;
	while(n>0) { q += (n % base)*bk; n /= base; bk /= base; } // svarer til ligning 11 i monte carlo kapitlet
	return q;
}

//Har til formål at lave quasi tilfældige tal ved hjælp af corput
void halton1(int n, int dim, double *a, double *b, double *x){
	int base[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79};
	int dmax=sizeof(base)/sizeof(int); assert(dim <= dmax); //Vi sikre os at demensionen af vores integral ikke er større end vores base.
	for(int i=0;i<dim;i++)x[i]=a[i]+corput(n+1,base[i])*(b[i]-a[i]); //Her dannes der så quasi tilfældige tal ud fra basen
}

//Gør det samme som halton1 bare med en anden base, dette er nødvendigt for at kunne estimere fejlen af vores montecarlo integral
void halton2(int n, int dim, double *a, double *b, double *x){
	int base[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83};
	int dmax=sizeof(base)/sizeof(int); assert(dim <= dmax);
	for(int i=0;i<dim;i++)x[i]=a[i]+corput(n+1,base[i])*(b[i]-a[i]);
}

double complex quasiMC(int dim, double f(int,double*), double* a, double* b, int N){
	double sum1=0, sum2=0, x[dim];
	double vol=1; for(int i=0;i<dim;i++) vol*=b[i]-a[i];
	for(int i=0;i<N/2;i++){
		halton1(i,dim,a,b,x); sum1+=f(dim,x);
		}
	for(int i=0;i<N/2;i++){
		halton2(i,dim,a,b,x); sum2+=f(dim,x);
		}
	double integ=(sum1+sum2)/N*vol;
	double error=fabs(sum1-sum2)/N*vol;
	return integ+I*error;
	}
