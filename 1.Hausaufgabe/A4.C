#include <math.h>
#include<tgmath.h>
#include<stdlib.h>
#include <stdio.h>

int sign(double x){
	if (x>=0){
		return 1;
	}
	if (x<0){
		return -1;
	}
	return 0;
}

double coshf(double x){
	return 0.5*(exp(x)+exp((-1)*x));
	}
double cosh_midpoint(double x, double delta_x){
	int N=(ceil(abs(x-0)/delta_x));
	double s=0;
	double m_i=0;
	for (int i=0;i<=N-1;i++){
		m_i=(i+i+1)*delta_x/2;
		s=s+cosh(m_i);
	}
	return sign(x)*s*delta_x;
}
double cosh_simpson(double x,double delta_x){
	double s=0;
	int N=(ceil(abs(x-0)/delta_x));
	double m_i=0;
	for (int i=0;i<=N-1;i++){
		m_i=(i+i+1)*delta_x/2;
		s=s+cosh(i*delta_x)+4*cosh(m_i)+cosh((i+1)*delta_x);
	}
	return sign(x)*s*delta_x/6;
}
double Integral(double x){
	return (0.5*(exp(x)-exp(-x))-0.5*(exp(0)-exp(-0)));
}
int main(){
	FILE*myFile = fopen("A4_values", "w");
	double dx[50];
	dx[0]=1e-6;
	for (int i=1;i<50;i++){
		dx[i]=1e-6*(i+1);
	}
	for (int i=0;i<50;i++){
		fprintf(myFile,"%g,%g,%g,%d\n",dx[i], cosh_simpson(1,dx[i]),cosh_midpoint(1,dx[i]), i);
	}	
	fclose(myFile); 
	
return 0;
}



