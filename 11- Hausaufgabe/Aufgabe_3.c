#include <lapacke.h>
#include <math.h>
#include<tgmath.h>
#include<stdlib.h>
#include <stdio.h>


double V(double x){ //Potential
	return x*x/2;
}

double HH(double x,double g){ //Stör-Hamiltonoperator
	return g*x*x*x*x;
}

void Matrix_H(int N, double * H, double * x, double delta_x, double g){  //füllt Matrix H mit Einträgen
	double dummy;
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

int n(double delta_x, double L){ // Anzahl Gitterpunkte
	return 2*L/delta_x;

}

double skalarprodukt(double* a, double*b, int N, double delta_x){ // rechnet skalarprodukt zweier Funktionen
	double res=0;
	for (int i=0;i<N+1;i++){
		res=res+a[i]*b[i];
	}
	return res*delta_x;
}

void Matrix_vektor_produkt(double * H, double * psi, double * res, int N){ //rechnet (H psi) mit der komischen idex schreibweise und speichert das auf res 
	double sum;
	for (int i = 0; i < N+1; i++){
		sum=0;
		for (int j = 0; j < N+1; j++){
			sum=sum + H[i*(N+1) + j] * psi[j];
		}
		res[i] = sum;
	}
	return;
}

void Matrix_H_alpha_beta(double * H__, double * H_, double * psi ,int N, double delta_x,int NN){ 
	 // füllt H__[a][b] mit <u_a|Hu_b>/ <u_a|u_b> aus
		double * res=  malloc(sizeof(double)*(N+1)); // in diesem Array wird (H u_b) temporär gespeichert
		double * psii[N+1]; // in diesem werden die Eigenzustände gespeichert
		for (int k=0;k<N+1;k++){
			psii[k]= malloc(sizeof(double)*(N+1));;
		}
		for (int k=0;k<N+1;k++){
			for (int l=0;l<N+1;l++){
				psii[k][l]= psi[k*(N+1) + l];
			}
		}
		// psii[b] ist Eigenzustand b des ungestöreten Hamilton-Operators
		for (int a = 0; a < NN+1; a++){
			for (int b = 0; b < NN+1; b++){
				Matrix_vektor_produkt(H_, psii[b], res ,N); // H_ u_b ergebnis liegt in res. Es gilt also res=H_ u_b
				H__[a*(N+1) + b] = skalarprodukt(psii[a], res, N, delta_x) / skalarprodukt(psii[a],psii[b], N, delta_x); 
			}
		}
		for (int k=0;k<N+1;k++){
			free(psii[k]);
		}
		free(res);
	return;
}


void variations(FILE*someFile,  double * H1_, double * H_ ,int N, double delta_x,int NN){
	// in H_ stehen die Eigenzustände
	// H_1 ist Hamiltonoperator mit Störung g=/=0 
	// H_a_b soll H mit Störung in der Eigenbasis H_ sein
	// (N+1) ist dimension von H_
	// (NN+1) soll Anzahl der verwendeten Eigenbasen sein
	// Die Funktion schreibt in einer Datei Anzahl der Eigenbasen sowie ersten vier Eigenwerte
	int info;
	double * H_a_b; // hier ist H [alpha] [beta]
	double * w_a_b; // hier werden eigenwerte gespeichert
	double * zwischenspeicher = malloc(sizeof(double)*(N+1));
	H_a_b = malloc(sizeof(double)*(NN+1)*(NN+1));
	w_a_b = malloc(sizeof(double)*(NN+1)); // hier sind die Eigenwerte
	Matrix_H_alpha_beta(H_a_b, H1_, H_, N, delta_x, NN);
	info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'N', 'U', NN+1, H_a_b, NN+1, w_a_b);
	if (info!=0){
		printf("Fehler ist aufgetreten, Eigenwerte wurden nicht bestimmt. A3\n");
	}
	fprintf(someFile,"%d,",NN);
	for (int l=0;l<4;l++){
		fprintf(someFile,"%g,",w_a_b[l]);
	}
	fprintf(someFile,"\n");
	free(H_a_b);
	free(w_a_b);
	free(zwischenspeicher);
	return;
}




int main(){
//Teilaufgabe 3
	double delta_x;
	int info;
	double L=5;
	int N;
	double g;
	FILE*someFile = fopen("A3_Eigenvalues", "w");
	delta_x = 0.01;
	N = n(delta_x, L);
	double x_[N+1];
	double w_[N+1];
	double H_[(N+1)*(N+1)];
	// finde ungestörten Operator
	g = 0;
	for (int k=0;k<=N+1;k++){
		x_[k]=k*delta_x-L;
	}
	Matrix_H(N, H_,  x_, delta_x, g);
	//finde Eigenzustände nummerisch	
	info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', N+1, H_, N+1, w_ );
	if (info!=0){
		printf("Fehler ist aufgetreten, Eigenwerte wurden nicht bestimmt. Info ist %d, delta_x ist: %g \n",info, delta_x);
	}
	// find H gestört
	g = 0.5;
	double * H1_ = malloc(sizeof(double)*(N+1)*(N+1));
	Matrix_H(N, H1_,  x_, delta_x, g);
	// finde  H [alpha] [beta]
	int N_[] = {4,8,16,32,64,128};
	int NN;
	for (int k=0;k<6;k++){
		NN = N_[k];
		variations(someFile, H1_, H_, N, delta_x, NN);
	}
	
	free(H1_);		
	fclose(someFile);
	
	return 0;
}





