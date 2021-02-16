#include<stdio.h>
#include<math.h>
#include<stdlib.h>


void change(int *y){
	*y=1;
}

int equal(double a, double b, double tau, double epsilon){

	int answer =0;
	double z=(fabs(a-b))/(fabs(a)+fabs(b));
	double d=fabs(a-b);
	if(d<tau){
		change(&answer);
		printf("from tau %lg\n",d);
}
	if(z<epsilon/2){
		change(&answer);
		printf("from eps%lg\n",z);

}

	return answer;
}




int main(){
	double a=4.5;
	double b=3;
	double tau=5;
	double eps=10;

	int x=equal(a,b,tau,eps);
	printf("answer is %i\n",x);
return 0;
}
