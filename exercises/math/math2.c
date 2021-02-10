#include<math.h>
#include<complex.h>
#include<stdio.h>

int main(){
float x_float = 1.f/9;
double x_double = 1./9;
long double x_long_double = 1.L/9;

	printf("float %25g\n",x_float);
	printf("double %251g\n",x_double);
	printf("long double %25Lg\n",x_long_double);
return 0;
}
