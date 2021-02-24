#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
	double x;
	int items;
	while(items!=EOF){
		items = scanf("%lg",&x); // from stdin into "x"
		printf("x=%lg sin(x)=%lg, cos(x)=%lg \n",x,sin(x),cos(x));
	}
return 0;
}

