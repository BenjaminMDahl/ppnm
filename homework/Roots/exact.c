#include<stdio.h>
#include<math.h>

int main(){
int n=8001;
double x[n],y[n];
for(int i=0;i<n;i++){
x[i]=i/1000;
y[i]=x[i]*exp(-x[i]);
printf("%6g %6g\n",x[i],y[i]);}
}
