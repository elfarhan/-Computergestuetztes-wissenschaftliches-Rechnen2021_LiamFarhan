#include <math.h>
#include<tgmath.h>
#include<stdlib.h>
#include <stdio.h>
#include "code.c" 

int main(){
	printf("Bitte warten...\n");
	/* bereite zufallszahlerzeuger vor
	   viele Zufallszahlen werden generiert, sodass der Seed (last_number), welcher für die Erzeugung der Zufallszahlen verwendet wird, 
	zufällig aussieht
	*/
	long int last_number=42; 
	for (int k=0;k<4367222;k++){
		generate_r( &last_number);	
	}
//////		 //////	
////// Aufgabe 1 //////	
//////		 //////
	FILE*someFile = fopen("A1_values", "w"); // Ergebnisse von A1 werden hier gespeichert
	int length = 10; // Länge der Array, auf dem sich die Walker befinden werden
	double p = 0.2;
	double q = 0.5;
	int d[] = {1,3,2};
	int t = 500; // wie Lange dauert die Simulation
	int N = 1000; // wie viele Simulationen werden für die Berechnung der Mittelwerte verwendet
	double average = 0;
	double stdev = 0;
	double d_temp;
	for (int j=0;j<3;j++){
		d_temp=d[j];
		for (int k=0; k<=500; k++){
			p= (double) k/1000;
			average_density(p, q, &last_number, length, d_temp, t, N, &average, &stdev);	
			fprintf(someFile, "%g, %g, %g\n",p, average, stdev);
		
		}
	}
	fclose(someFile);
	printf("Aufgabe 1 fertig! \n");
	printf("Aufgabe 2 braucht viel mehr Zeit... \n");
//////		 //////	
////// Aufgabe 2 //////	
//////		 //////
	FILE*anotherFile = fopen("A2_values", "w"); // Ergebnisse von A2 werden hier gespeichert
	p = 0.05;
	d_temp = 3;
	for (int N=1; N<=3000; N++){
		average_density(p, q, &last_number, length, d_temp, t, N, &average, &stdev);
		fprintf(anotherFile, "%d, %g, %g\n", N, average, stdev);
		
		}
		fclose(anotherFile);
		printf("Aufgabe 2 fertig! \n");
	return 0;
}
