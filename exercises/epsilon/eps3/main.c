#include<stdio.h>
#include<math.h>
#include<stdlib.h>


int equal(double a, double b, double tau, double epsilon);





int main(){
	double a=4.5;
	double b=3;
	double tau=5;
	double eps=10;

	printf("Vi prøver med a=%g b=%g tau=%g eps=%g. \n",a,b,tau,eps);
	printf("Vi forventer her at se et 0, da a+b=%g som ikke er lig hverken tau eller eps/2\n",a-b);
	int x=equal(a,b,tau,eps);
	printf("svaret er %i som forventet\n\n",x);

	double a2=3;
	double b2=2;
	double tau2=a2-b2;
	double eps2=3;
	printf("Vi prøver nu med a=%g b=%g tau=%g eps=%g. \n",a2,b2,tau2,eps2);
	printf("Her forventer vi et 1 da a-b=tau\n");
	int y=equal(a2,b2,tau2,eps2);
	printf("svaret er %i som forventet\n\n",y);

	double a3=1;
	double b3=0.5;
	double tau3=4;
	double eps3=(a3-b3)/(a3+b3)*2;
	printf("Vi prøver nu med a=%g b=%g tau=%g eps=%g. \n",a3,b3,tau3,eps3);
	printf("Her forventer vi et 1 da a-b=eps/2\n");
	int z=equal(a3,b3,tau3,eps3);
	printf("svaret er %i som forventet\n",z);


return 0;
}
