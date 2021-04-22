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
	double complex z4= cpow(I,exp(1));
	double complex z5= cpow(I,I);

	printf("Via funktionerne fra math.h og complex.h har vi udregnet følgende størrelser:\n\n");

	printf("gamma(5)=%g\n",g);
	printf("bessel(0.5)=%g\n",b);
	printf("sqrt(-2)=%g+i%g\n",creal(z1),cimag(z1));
	printf("exp(i*pi)=%g+i%g\n",creal(z2),cimag(z2));
	printf("exp(i)=%g+i%g\n",creal(z3),cimag(z3));
	printf("i^e=%g+i%g\n",creal(z4),cimag(z4));
	printf("i^i=%g+i%g\n",creal(z5),cimag(z5));

	printf("\n Vi holder dette op imod, hvad man får hvis man bruger wolframalpha til at udregne det samme\n\n");

	//Det der bliver printet her er fra https://www.wolframalpha.com/:

	printf("gamma(5)=24\n");
	printf("bessel(0.5)=0.93847\n");
	printf("sqrt(-2)=i1.4142136\n");
	printf("exp(i*pi)=-1\n");
	printf("exp(i)=0.540302+i0.841471\n");
	printf("i^e=-0.42822-0,90367\n");
	printf("i^i=0.2078796\n");

	printf("\n Det ses at svarene er ret ens.");
return 0;
}


