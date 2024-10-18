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


int main(){

//Teilaufgabe 1

	//diskretisiere x-Achse
	double L=5;
	int N=500;
	double delta_x=2*L/N;
	double x[N+1];
	for (int i=0; i<N+1;i++){
		x[i]=i*delta_x-L;
	}
	//Matrix H
	double H[(N+1)*(N+1)]; // Matrix erstellen
	double g=0;
	Matrix_H(N, H,  x, delta_x, g); // Einträge werden gefüllt
	//Diagonalisieren der Matrix
	int info;
	double w[N+1]; //hier werden die Eigenwerte gespeichert
	info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', N+1, H, N+1, w );
	if (info!=0){
		printf("Fehler ist aufgetreten, Eigenwerte wurden nicht bestimmt.\n");
	}
	//speichere 20 Eigenwerte in einer Datei
	FILE*myFile = fopen("A1_Eigenvalues", "w");
	for (int i=0;i<20;i++){ 
		fprintf(myFile,"%g \n",w[i]);
	}
	fclose(myFile);
	//speichere Komponenten der ersten 3 Eigenvektroren
	FILE*anotherFile = fopen("A1_Eigenvectors", "w");
	for (int j=0;j<3;j++){
		for (int i=0;i<N+1;i++){
			fprintf(anotherFile,"%g \n",H[j*(N+1) + i]);
		}
	}
	fclose(anotherFile);
	
//Teilaufgabe 2

	g=0.5;
	int itters= 1000;
	double * H_;
	double * x_;
	double * w_;
	FILE*yetanotherFile = fopen("A2_Eigenvalues", "w");
	for (int k=0;k<=itters;k++){
		delta_x= (1-0.01)/itters *k +0.01;
		N = n(delta_x, L);
		H_ = malloc(sizeof(double)*(N+1)*(N+1));
		x_ = malloc(sizeof(double)*(N+1));
		w_ = malloc(sizeof(double)*(N+1));
		for (int i=0; i<N+1;i++){
			x_[i]=i*delta_x-L;
			//printf("%g \n", x[_i]);
		}	
		Matrix_H(N, H_,  x_, delta_x, g);
		info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'N', 'U', N+1, H_, N+1, w_ );
		if (info!=0){
			printf("Fehler ist aufgetreten, Eigenwerte wurden nicht bestimmt. Info ist %d, delta_x ist: %g \n",info, delta_x);
		}
		fprintf(yetanotherFile,"%g, ",delta_x);
		for (int z=0; z<4; z++){
			fprintf(yetanotherFile,"%g,",w_[z]);
		}
		fprintf(yetanotherFile,"%g\n",w_[4]);
		free(H_);
		free(x_);
		free(w_);
	}
	

	return 0;
}



































