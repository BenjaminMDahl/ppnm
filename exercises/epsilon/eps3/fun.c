#include<math.h>
#include<stdio.h>


void change(int *y){
	*y=1;
}


// Da opgaven ikke beder om at få at vide om den retunere 1 pga absolut eller relativ value har vi udkommeteret nedstående, som ellers ville sige hvilken man var lig
int equal(double a, double b, double tau, double epsilon){

	int answer =0;
	double z=(fabs(a-b))/(fabs(a)+fabs(b));
	double d=fabs(a-b);
	if(d==tau){
		change(&answer);
//		printf("was equal to tau\n");
		}
	if(z==epsilon/2){
		change(&answer);
//		printf("was equal eps/2\n");
		}
	return answer;}


