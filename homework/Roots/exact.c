#include<stdio.h>
#include<math.h>

int main(){
int n=801;
double x[n],y[n];
for(int i=0;i<n;i++){
x[i]=i/100;
y[i]=x[i]*exp(-x[i]);
printf("%6g %6g\n",x[i],y[i]);}
}
