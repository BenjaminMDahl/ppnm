#include<stdio.h>
#include<math.h>
#include<complex.h>

int main(){
	double pi=M_PI;
	double g=tgamma(5);
	double b=y0(0.5);
	complex z1=csqrt(-2);
	complex z2= cexp (I * pi);
	complex z3= cexp (I);
	complex z4= pow(I,exp(1));
	complex z5= pow(I,I);


	printf("gamma(5)=%g\n",g);
	printf("bessel(0.5)=%g\n",b);
	printf("sqrt(-2)=%g+i%g\n",creal(z1),cimag(z1));
	printf("exp(i*pi)=%g+i%g\n",creal(z2),cimag(z2));
	printf("exp(i)=%g+i%g\n",creal(z3),cimag(z3));
	printf("i^e=%g+i%g\n",creal(z4),cimag(z4));
	printf("i^i=%g+i%g\n",creal(z5),cimag(z5));
return 0;
}


