#include<stdio.h>
#include<stdlib.h>

//
//ODE functions
//		



double f(double** y,int t){
	return y[0][t]*(1-y[0][t]);
	}
	
////////////////////////////////////////////
/////////////Euler algorithm////////////////
////////////////////////////////////////////
		
void derivs(int t, int Ordnung, double** y, double** F, double* params,double delta_t){
	F[Ordnung-2][t]=f(y,t);
//	for (int i=0; i<Ordnung-1;i++)	{
	//	F[i][t]=y[i+1][t];
	//	}
		
		return;
	}

void Eulerstep(double** y, double** F,int t,double mu, int order){
	y[0][t+1]=F[0][t]*4*mu;
	//for(int i=0; i<order-1;i++){
		//y[i][t+1]=F[i][t]*4*mu;// + y[i][t];
		//}
		//y[order-1][t+1]=f(y,t);
	return;
	}
	
void Eulersolver(double** y, double** F,int t, double delta_t, int order, double* params,int numberofsteps){
	F[0][0]=f(y,0);
	Eulerstep(y,F,0,delta_t,order);
	for(int t=0;t<numberofsteps-1;t++){
		derivs(t, order, y, F, params,delta_t);
		Eulerstep(y,F,t,delta_t,order);
		}
	return;
	}
		
		
		
		
		
		
/////////////'mainprogramm'#///////////
////////////////////////////////////////////		
				
int main(){
	FILE*myFile[] = {fopen("mu1_values", "w"),fopen("mu2_values", "w"),fopen("mu3_values", "w")};
//x_i between(0,...,100)#######
	int numberofsteps=101;
	int order=2;
	double r_0=1;
	double x_0=0.1;
	double params[]={r_0,x_0};
	// make dynamic 2d arrays
	double * y[order];
	double * F[order-1];
	for(int i=0;i<order;i++){
		y[i]=malloc(sizeof(double)*numberofsteps);
		F[i]=malloc(sizeof(double)*numberofsteps);
		if (y==NULL || F==NULL){return -1;}
	}
	y[0][0]=params[1];
	//Eulersolver(y, F, numberofsteps, 0.4, order, params,numberofsteps);
	double mu[]={0.4,0.74,0.77};
		
	for (int j=0;j<3;j++){
		Eulersolver(y, F, numberofsteps, mu[j], order, params,numberofsteps);
		for (int i=0;i<numberofsteps;i++){
			fprintf(myFile[j], "%g,%g\n",mu[j]*i,y[0][i]);
			}
		fclose(myFile[j]);
		}
//close allocated memory	
	for(int i=0;i<order;i++){
		free(y[i]);
		free(F[i]);
		}
	

//x_i between(0,...,1000)#######
	numberofsteps=1001;
	x_0=0.0000000001;
	// make dynamic 2d arrays
	for(int i=0;i<order;i++){
		y[i]=malloc(sizeof(double)*numberofsteps);
		F[i]=malloc(sizeof(double)*numberofsteps);
		if (y==NULL || F==NULL){return -1;}
	}
	y[0][0]=params[1];
	FILE* newFile=fopen("x1000vsmu", "w");
	double dmu=0.0001;
	for (int i=1;i<6001;i++){
		double muu=dmu*i+0.4;
		Eulersolver(y, F, numberofsteps, muu, order, params,numberofsteps);
		fprintf(newFile, "%g,%g\n",muu,y[0][1000]);
	}
	return 0;
	}
	
	
	
	
	 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

