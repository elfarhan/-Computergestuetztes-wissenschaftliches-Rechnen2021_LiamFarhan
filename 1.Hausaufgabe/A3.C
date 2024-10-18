#include <math.h>
#include<tgmath.h>
#include<stdlib.h>
#include <stdio.h>
#define PI acos(-1.0)
int sign(double x){
	if (x>=0){
		return 1;
	}
	if (x<0){
		return -1;
	}
	return 0;
}
double erf_midpoint(double x, double delta_x){
	int N=(ceil(abs(x-0)/delta_x));
	double s=0;
	double m_i=0;
	for (int i=0;i<=N-1;i++){
		m_i=0+(i+i+1)*delta_x/2;
		s=s+exp((-1)*pow(m_i,2));
	}
	return sign(x)*s*delta_x*2/sqrt(PI);
}

double erf_simpson(double x, double delta_x){
	double s=0;
	int N=(ceil(abs(x-0)/delta_x));
	double m_i=0;
	for (int i=0;i<=N-1;i++){
		m_i=(i+i+1)*delta_x/2;
		s=s+exp((-1)*pow(m_i,2))*4 + exp((-1)*pow(i*delta_x,2)) + exp((-1)*pow((i+1)*delta_x,2));
	}

	return sign(x)*s*delta_x/6*2/sqrt(PI);
}

int main(){
	FILE*myFile = fopen("A3_values", "w");
	double x[100];
	for (int i=0;i<=100;i++){
		x[i]=-2+1e-2*i*4;		
	}
	for (int i=0;i<=100;i++){
		fprintf(myFile, "%g,%g,%g,%d\n", x[i], erf_simpson(x[i], 1e-4),erf_midpoint(x[i], 1e-4), i);
	}
	fclose(myFile);
	
	return 0;}

