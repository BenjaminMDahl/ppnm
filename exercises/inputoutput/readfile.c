#include<stdio.h>
#include<math.h>
int main(){
	double x;
	int items;
	FILE* out=fopen("opg3.out.txt","w");
	FILE* input=fopen("til3.txt","r");
	while (fscanf(input, "%lg", & x ) == 1 ){
		fprintf(out,"x= %3lg, sin(x)=%10lg, cos(x)=%10lg \n", x,sin(x),cos(x));
        }
	fclose(input);
	fclose(out);
	printf("\n");
return 0;
}
