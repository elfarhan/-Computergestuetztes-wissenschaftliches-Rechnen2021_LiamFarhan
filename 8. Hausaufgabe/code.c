#include <math.h>
#include<tgmath.h>
#include<stdlib.h>
#include <stdio.h>

// find char. polynomial
double p(int n, double a[n], double b[n], double lambda, int grade){
	double w=a[0];
	double t=b[0];
	if (n<0){
		printf("error!\n");
		return 0;
	}
	if (n==0){
		return 1;
	}
	if (n==1){
		return w-lambda;
	}
	return p(n-1, a, b, lambda, grade)*(w-lambda)- t*t* p(n-2,a, b, lambda, grade);
}

// find roots implementation

double find_root(double func(double), double x0, double delta, double rel_tol,int max_iter){
	if (max_iter<=0 || delta==0){
		return x0;
	}
	double letzter_wert;
	double func_ableitung;
	double x=x0;
	for (int n=1;n<=max_iter; n++){
		func_ableitung = ( func(x+delta) -func(x-delta) ) / (2*delta);
		letzter_wert=x;
		x = x - func(x)/func_ableitung;
		if ( abs((letzter_wert-x)/letzter_wert) < rel_tol){
			return x;
		}
	}
	return 0; 
}



// main program //
//
//
//
int main(){
// initilize Matrix
	int n=10;
	int grade=n;
	double a[n];
	double lambda[n];
	double b[n-1];
	double t=0.5;
	double w=1;
	for (int i=0;i<n;i++){
		a[i]=w;
	}
	for (int i=0;i<n-1;i++){
		b[i]=-t;
	}
// evaluate p for many different values of lambda and save it in a separate file	
	FILE*myFile = fopen("A3_values", "w");
	int k;
	double dt=0.001;
	double boundary_left=0;
	double boundary_right=2;
	k=(boundary_right-(-boundary_left))/dt;
	double x[k];
	for (int i=0;i<k;i++){
		x[i]=-boundary_left+i*dt;
		fprintf(myFile, "%g,%g\n",x[i], p(n, a, b, x[i], grade));
	}
	fclose(myFile);
// make a simple container for p to decrease number of arguments
	double f(double x){
		return p(n, a, b, x, grade);
	}
	

 	
// use find_root to find roots and output that:
	FILE*myFile2 = fopen("A4_values", "w");
	double lamda_guesses[10]={0.01, 0.15,0.30,0.55,0.8,1.1,1.4,1.6,1.8,1.97};
	double Null_stelle;
	fprintf(myFile2, "Die gefundene Nullstellen sind:\n");
	for (int i=0; i<10; i++){
		Null_stelle=find_root(f, lamda_guesses[i], 1e-9, 1e-5,1e9);
		lambda[i]=Null_stelle;
		fprintf(myFile2, "%lf\n",Null_stelle);
	
	}
	fclose(myFile2);
	return 0;
}
//
//
//////////

