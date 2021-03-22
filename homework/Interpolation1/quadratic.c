#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>


int binsearch(int n, double* x, double z){
	assert(n>1 && x[0]<=z && z<=x[n-1]); // Tjekker om Z faktisk er i intervallet og om vi har mere end 1 datapunkt, hvis ikke giver den fejl.
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	}


	// Qudratic spline //

typedef struct {int n ; double *x , *y , *b , *c;} qspline; 	// Der definieres en struct (qspline) hvor alt vores data skal opsamles i

qspline* qspline_alloc ( int n , double* x , double* y ){	//qspline bygges her
	qspline  *s = (qspline*) malloc (sizeof(qspline));		// Først laves der plads til alt data i hukommelsen
	s->b = (double*) malloc((n-1)*sizeof(double));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->x = (double*) malloc(n*sizeof(double)) ;
	s->y = (double*) malloc(n*sizeof(double)) ;
	s->n = n;
	int i; // Vi gidder ikke skrive int i alle vores forlykker
	for (i = 0; i<n; i++){s->x[i]=x[i]; s->y[i]=y[i];}		//Der henvises til vores data
	double a[n-1] , dx[n-1];
	for (i = 0; i<n-1; i++){
		dx[i]=x[i+1]-x[i];
		a[i]=(y[i+1]-y[i])/dx[i];}
	//Nu er alt sat op til at udregne en quadratic spline, som gøres på to måder, da starten ikke kendes
	s->c[0]=0;
	for( i=0; i<n-2; i++){
		s->c[i+1]=(a[i+1]-a[i]-s->c[i]*dx[i])/dx[i+1];}		//Sætning 13 fra kapitel 1
	s->c[n-2]/=2;
	for (i=n-3; i >=0;i--){
		s->c[i]=(a[i+1]-a[i]-s->c[i+1]*dx[i+1])/dx[i];}		//sætning 14 fra kapitel 1
	for (i =0; i<n-1; i++){
		s->b[i]=a[i]-s->c[i]*dx[i];}				//Sætning 15 fra kapitel 1
return s;
}


double qspline_eval(qspline *s, double z){			//Functionen som via qspline udregner s(z)
	int i=binsearch(s->n,s->x,z);
	double dx=z-s->x[i];
	double f_z=s->y[i]+dx*(s->b[i]+dx*s->c[i]);
return f_z;
}

double qua_integral(qspline *s, double z){
	int I=binsearch(s->n,s->x,z);
	double S=0;

	for(int i=0;i<I;i++){
	double dx=s->x[i+1]-s->x[i]; assert(dx>0);
		S+=dx*(s->y[i]+dx*(s->b[i]/2+dx*s->c[i]/3));
	}

	double dx=z-s->x[I]; assert(dx>0);
	double F_z=dx*(s->y[I]+dx*(s->b[I]/2+dx*s->c[I]/3));
	S+=F_z;
return S;
}

double qua_diff(qspline *s, double z){
	int i=binsearch(s->n,s->x,z);
	double dx=z-s->x[i];
	double df_dz=s->b[i]+2*dx*s->c[i];
return df_dz;
}

void qsplinefree(qspline *s){
	free(s->x);free(s->y);free(s->b);free(s->c);free(s);
}

int main(){
	//Rå data
	int n=5;
	double x[5]={1,2,3,4,5};
	double y1[5]={1,1,1,1,1};
	double yi[5]={1,2,3,4,5};
	double yi2[5]={1,4,9,16,25};


	//Linspace
	int i, N=79; double z[N-1];
	for(i=1; i<N;i++){
	z[i-1]=(double)(5*i+100)/(101);
	}

	// y=1
	double f_z_qua_1[N-1], F_z_qua_1[N-1], df_dz_qua_1[N-1];
	qspline* q_spline_y1=qspline_alloc(n,x,y1);
	printf("#index 0: quadratic spline data for y=1(x f(x) df_dx F(x))\n");
	for(i=0;i<N-1;i++){
		f_z_qua_1[i]=qspline_eval(q_spline_y1,z[i]);
		F_z_qua_1[i]=qua_integral(q_spline_y1,z[i]);
		df_dz_qua_1[i]=qua_diff(q_spline_y1,z[i]);
		printf("%g %g %g %g \n",z[i],f_z_qua_1[i],df_dz_qua_1[i],F_z_qua_1[i]);
	}
	printf("\n \n");
	qsplinefree(q_spline_y1);

	// y=i
	double f_z_qua_i[N-1], F_z_qua_i[N-1], df_dz_qua_i[N-1];
	qspline* q_spline_yi=qspline_alloc(n,x,yi);
	printf("#index 1: quadratic spline data for y=i(x f(x) df_dx F(x))\n");
	for(i=0;i<N-1;i++){
		f_z_qua_i[i]=qspline_eval(q_spline_yi,z[i]);
		F_z_qua_i[i]=qua_integral(q_spline_yi,z[i]);
		df_dz_qua_i[i]=qua_diff(q_spline_yi,z[i]);
		printf("%g %g %g %g \n",z[i],f_z_qua_i[i],df_dz_qua_i[i],F_z_qua_i[i]);
	}
	printf("\n \n");
	qsplinefree(q_spline_yi);


	// y=i^2
	double f_z_qua_i2[N-1], F_z_qua_i2[N-1], df_dz_qua_i2[N-1];
	qspline* q_spline_yi2=qspline_alloc(n,x,yi2);
	printf("#index 0: quadratic spline data for y=i^2(x f(x) df_dx F(x))\n");
	for(i=0;i<N-1;i++){
		f_z_qua_i2[i]=qspline_eval(q_spline_yi2,z[i]);
		F_z_qua_i2[i]=qua_integral(q_spline_yi2,z[i]);
		df_dz_qua_i2[i]=qua_diff(q_spline_yi2,z[i]);
		printf("%g %g %g %g \n",z[i],f_z_qua_i2[i],df_dz_qua_i2[i],F_z_qua_i2[i]);
	}
	printf("\n \n");
	qsplinefree(q_spline_yi2);

return 0;
}
