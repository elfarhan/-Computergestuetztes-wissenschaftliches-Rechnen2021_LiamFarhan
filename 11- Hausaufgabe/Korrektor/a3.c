#include <math.h>
#include<tgmath.h>
#include<stdlib.h>
#include <stdio.h>
#include <lapacke.h>

double L=5;

int n(double delta_x, double L){ // Anzahl Gitterpunkte
	return 2*L/delta_x;

}

double V(double x){ //Potential
	return x*x/2;
}

double HH(double x,double g){ //Stör-Hamiltonoperator
	return g*x*x*x*x;
}

void Matrix_H(int N, double * H, double delta_x, double g){  //füllt Matrix H mit Einträgen
	double dummy;
	double x[N+1];
	for (int i=0;i<N+1;i++){
		x[i]=i*delta_x-L;;
	}
	// fülle Matrix mit 0
	for (int i=0;i<N+1;i++){
		for (int j=0;j<N+1;j++){
				H[i*(N+1) + j] = 0; // entspricht H[i][j]
		}
	}
	// erstelle Diagonale
	for (int i=0;i<N+1;i++){
		dummy = 1/(delta_x*delta_x) + V(x[i]) + HH(x[i],g);
		H[i*(N+1) + i] = dummy;
		//H[i][i]= H_D[i];
	}
	//erstelle Nebendiagonalen
	for (int i=0;i<N;i++){
		dummy = 1/(delta_x*delta_x*2);
		H[(i+1)*(N+1) + i] = dummy;
		H[i*(N+1) + (i+1)] = dummy;
		//H[i+1][i] = H_ND[i];
		//H[i][i+1] = H_ND[i];
	}	
	return ;
}


void k_th_Eigenzustand(double * H, double * psi_k,int N, int k){
		for(int i=0;i<N+1;i++){
			psi_k[i] = H[k*(N+1) + i];
		}

	return;
}


double skalarprodukt(double* a, double*b, int N, double delta_x){ // rechnet skalarprodukt zweier Funktionen
	double res=0;
	for (int i=0;i<N+1;i++){
		res=res+a[i]*b[i];
	}
	return res*delta_x;
}

void Matrix_vektor_produkt(double * H, double * u, double * res, int N){ //rechnet (H u) mit der komischen idex schreibweise und speichert das in res 
	double sum;
	for (int i = 0; i < N+1; i++){
		sum=0;
		for (int j = 0; j < N+1; j++){
			sum=sum + H[i*(N+1) + j] * u[j];
		}
		res[i] = sum;
	}
	return;
}



void Matrix_H_alpha_beta(double * H__, double * H_, double * psi ,int N, double delta_x,int M){ 
	 // füllt H__[a][b] mit <u_a|Hu_b>/ <u_a|u_b> aus
		double * u_a = malloc(sizeof(double)*(N+1)); // hier ist a-te Eigenzustand des ungestörten H Operators gespeichert
		double * u_b = malloc(sizeof(double)*(N+1));// hier ist b-te Eigenzustand des ungestörten H Operators gespeichert
		double * H_u_b= malloc(sizeof(double)*(N+1)); // in diesem Array wird (H u_b) temporär gespeichert
	
		for (int a = 0; a < M+1; a++){
			k_th_Eigenzustand(psi,  u_a, N, a);
			for (int b = 0; b < M+1; b++){
				k_th_Eigenzustand(psi,  u_b, N, b);
				Matrix_vektor_produkt(H_, u_b, H_u_b ,N); // H_ u_b ergebnis liegt in H_u_b. 
				H__[a*(M+1) + b] = skalarprodukt(u_a, H_u_b, N, 1);// / skalarprodukt(u_a, u_b, N, delta_x); 
			}
		}
	free(u_a);
	free(u_b);
	free(H_u_b);
	
	return;
}



void variations(FILE*someFile,  double * H1_, double * H_ ,int N, double delta_x,int M){
	// in H_ stehen die Eigenzustände
	// H_1 ist Hamiltonoperator mit Störung g=/=0 
	// H_a_b soll H mit Störung in der Eigenbasis H_ sein
	// (N+1) ist dimension von H_
	// (NN+1) soll Anzahl der verwendeten Eigenbasen sein
	// Die Funktion schreibt in einer Datei Anzahl der Eigenbasen sowie ersten vier Eigenwerte
	int info;
	double * H_a_b = malloc(sizeof(double)*(M+1)*(M+1)); // hier ist H [alpha] [beta]
	double *  w_a_b = malloc(sizeof(double)*(M+1)); // hier werden eigenwerte gespeichert
	Matrix_H_alpha_beta(H_a_b, H1_, H_, N, delta_x, M);
	info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'N', 'U', M+1, H_a_b, M+1, w_a_b);
	if (info!=0){
		printf("Fehler ist aufgetreten, Eigenwerte wurden nicht bestimmt. A3\n");
	}
	fprintf(someFile,"%d, ",M);
	for (int l=0;l<3;l++){
		fprintf(someFile,"%g, ",w_a_b[l]);
	}
	fprintf(someFile,"%g\n",w_a_b[3]);
	free(H_a_b);
	free(w_a_b);
	
	return;
}





int main(){
//Teilaufgabe 3
	int info;
	double delta_x = 0.01;
	int N = n(delta_x, L);
	double * H= malloc(sizeof(double)*(N+1)*(N+1));
	double * w = malloc(sizeof(double)*(N+1));
	double g=0;
	FILE*someFile = fopen("A3_Eigenvalues", "w");
	// finde ungestörten Operator
	Matrix_H(N, H, delta_x, g);
	//finde Eigenzustände nummerisch
	info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', N+1, H, N+1, w );
	if (info!=0){
		printf("Fehler ist aufgetreten, Eigenwerte wurden nicht bestimmt. Info ist %d, delta_x ist: %g \n",info, delta_x);
	}
	// find H gestört
	g = 0.5;
	double * H1 = malloc(sizeof(double)*(N+1)*(N+1));
	Matrix_H(N, H1, delta_x, g);
	// finde  H [alpha] [beta]
	int N_[] = {4,8,16,32,64,128};
	int M;
	for (int k=0;k<6;k++){
		M = N_[k];
		variations(someFile, H1, H, N, delta_x, M);
	}
	
	free(H);
	free(H1);
	free(w);
	fclose(someFile);
	return 0;


}



