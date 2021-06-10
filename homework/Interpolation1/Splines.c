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


		// linear spline //

double linterp(int n, double x[], double y[], double z){
	int i=binsearch(n,x,z);

	double dy=y[i+1]-y[i];
	double dx=x[i+1]-x[i]; assert(dx>0);
	double f_z=y[i]+dy/dx*(z-x[i]);
return f_z;
}

double linterp_integral(int n, double x[], double y[], double z){
	int I=binsearch(n,x,z);
	double s=0;

	for(int i=0;i<I;i++){
	double dy=y[i+1]-y[i];
	double dx=x[i+1]-x[i]; assert(dx>0);
	double ai=dy/dx;
		s+=y[i]*(x[i+1]-x[i])+ai*pow((x[i+1]-x[i]),2)/2;
	}

	double dy=y[I+1]-y[I];
	double dx=x[I+1]-x[I]; assert(dx>0); //
	double az=dy/dx;
	double F_z=y[I]*(z-x[I])+az*pow((z-x[I]),2)/2;
	s+=F_z;
return s;
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


	// Cubic spline //

typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;					// Der skal igen dannes en ny strucktur med et parameter mere, men ellers gør vi meget det samme som ved qspline

cspline* cspline_alloc(int n, double *x, double *y){
	cspline* s = (cspline*)malloc(sizeof(cspline));
	s->x= (double*)malloc(n*sizeof(double));
	s->y= (double*)malloc(n*sizeof(double));
	s->b= (double*)malloc(n*sizeof(double));
	s->c= (double*)malloc((n-1)*sizeof(double));
	s->d= (double*)malloc((n-1)*sizeof(double));
	s->n= n;

	int i;
	for(i=0;i<n;i++){s->x[i]=x[i];s->y[i]=y[i];}

	double a[n-1] , dx[n-1];
	for (i = 0; i<n-1; i++){
		dx[i]=x[i+1]-x[i];
		a[i]=(y[i+1]-y[i])/dx[i];}
	// Her laves systemmet beskrevet ved sætning 21 fra kapitel 1
	double D[n], Q[n-1], B[n];
	D[0]=2; D[n-1]=2; Q[0]=1;
	for(i=0;i<n-2;i++)D[i+1]=2*dx[i]/dx[i+1]+2;
	for(i=0;i<n-2;i++)Q[i+1]=dx[i]/dx[i+1];
	for(i=0;i<n-2;i++)B[i+1]=3*(a[i]+a[i+1]*dx[i]/dx[i+1]);
	//solving the system:
	B[0]=3*a[0]; B[n-1]=3*a[n-2];
	for(i=1;i<n;i++){
		D[i]-=Q[i-1]/D[i-1];
		B[i]-=B[i-1]/D[i-1];}
	s->b[n-1]=B[n-1]/D[n-1];
	for(i=n-2;i>=0;i--)s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
	for(i=0;i<n-1;i++){
		s->c[i]=(-2*s->b[i]-s->b[i+1]+3*a[i])/dx[i];
		s->d[i]=(s->b[i]+s->b[i+1]-2*a[i])/dx[i]/dx[i];}
return s;
}

double cspline_eval(cspline *s,double z){
	int i=binsearch(s->n,s->x,z);
	double	dx=z-s->x[i];
	double f_z=s->y[i]+dx*(s->b[i]+dx*(s->c[i]+dx*s->d[i]));
return f_z;
}

double cub_integral(cspline *s, double z){
	int I=binsearch(s->n,s->x,z);
	double S=0;

	for(int i=0;i<I;i++){
	double dx=s->x[i+1]-s->x[i]; assert(dx>0);
		S+=dx*(s->y[i]+dx*(s->b[i]/2+dx*(s->c[i]/3+dx*s->d[i]/4)));
	}

	double dx=z-s->x[I]; assert(dx>0);
	double F_z=dx*(s->y[I]+dx*(s->b[I]/2+dx*(s->c[I]/3+dx*s->d[I]/4)));
	S+=F_z;
return S;
}

double cub_diff(cspline *s, double z){
	int i=binsearch(s->n,s->x,z);
	double dx=z-s->x[i];
	double df_dz=s->b[i]+dx*(2*s->c[i]+3*dx*s->d[i]);
return df_dz;
}


	// CLEAN UP //
	void qspline_free(qspline *s){
	free(s->x);free(s->y);free(s->b);free(s->c);free(s);
	}
	void cspline_free(cspline *s){
	free(s->x);free(s->y);free(s->b);free(s->c);free(s->d);free(s);
	}
