#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
	double x;
	int items;
	while(items!=EOF){
		items = scanf("%lg",&x); // from stdin into "x"
		printf("x=%3lg sin(x)=%10lg, cos(x)=%10lg \n",x,sin(x),cos(x));
	}
	printf("\n");
return 0;
}

