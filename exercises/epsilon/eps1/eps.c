
#include<limits.h>
#include<float.h>
#include<stdio.h>


int main(){
	int i=1; while(i+1>i) {i++;}
	printf("my max int = %i\n From while loop\n",i);

	int n=INT_MAX;
	int t=1;
	for(int k=1;k<n;k++)t++;
	printf("my max int = %i\n From for loop\n" ,t);

	int j=1;
	do
	{ 
	j++;
	}while (j+1>j);
	printf("my max int = %i\n From do while loop\n", j);

i=1; while(i>i-1) {i--;}
	printf("my min int = %i\n From while loop\n",i);

	n=INT_MIN;
	t=1;
	for(int k=1;k>n;k--)t--;
	printf("my min int = %i\n From for loop\n" ,t);

	j=1;
	do
	{ 
	j--;
	}while (j>j-1);
	printf("my min int = %i\n From do while loop\n", j);

float x=1; while(1+x!=1){x/=2;} x*=2;
printf("float difference=%g\n",x);
float X=FLT_EPSILON;
printf("value from float.h %g\n",X);
double y=1; while(1+y!=1){y/=2;} y*=2;
printf("double difference=%lg\n",y);
double Y=DBL_EPSILON;
printf("value from float.h %lg\n",Y);
long double z=1; while(1+z!=1){z/=2;} z*=2;
printf("long double difference%Lg\n",z);
long double Z=LDBL_EPSILON;
printf("value from float.h %Lg\n",Z);

float e; for(e=1; 1+e!=1; e/=2){} e*=2;
printf("from loop float differnce=%g\n", e);
double f; for(f=1; 1+f!=1; f/=2){} f*=2;
printf("from loop double differnce=%lg\n", f);
long double g; for(g=1; 1+g!=1; g/=2){} g*=2;
printf("from loop long double differnce=%Lg\n", g);

x=1;
do
{
x/=2;
}while (1+x!=1);
x*=2;
printf("from dowhile float difference=%g\n",x);


y=1;
do
{
y/=2;
}while (1+y!=1);
y*=2; printf("from dowhile double difference=%lg\n",y);


z=1;
do
{
z/=2;
}while (1+z!=1);
z*=2;
printf("from dowhile float difference=%Lg\n",z);

return 0;
}
