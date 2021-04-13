#include<math.h>
#include<assert.h>
#include<stdio.h>

//Vi tager udgangs punkt i et åbent interval og bruger 4 punkter, så vi kan bruge trapezium og rectangle rule fra bogen (vægtene fra sætning (49) og (50) i integrations kapitelet)


double quad(double f(double), double a, double b, double delta, double eps, double f2, double f3, int nrec)
{
assert(nrec<1000000);// tjekker at vi ikke får for mange recursions
double f1=f(a+(b-a)/6);
double f4=f(a+5*(b-a)/6);

double Q = (2*f1+f2+f3+2*f4)*(b-a)/6; // Baseret på (49) fra bogen kapitel Numerical integration
double q = (f1+f2+f3+f4)*(b-a)/4; // Baseret på (50) fra bogen kapitel numerical integration
double err=fabs(Q-q);
if (err < delta+eps*fabs(Q)) return Q;
else return quad(f,a,(a+b)/2,delta/sqrt(2),eps,f1,f2,nrec+1)+quad(f,(a+b)/2,b,delta/sqrt(2),eps,f3,f4,nrec+1);
}


double integrate(double f(double), double a, double b, double delta, double eps) // Dette første trin er for at kunne genbruge punkterne løbende
{
	int nrec=0;
	double f2=f(a+2*(b-a)/6);
	double f3=f(a+5*(b-a)/6);
	return quad(f,a,b,delta,eps,f2,f3,nrec);
}




