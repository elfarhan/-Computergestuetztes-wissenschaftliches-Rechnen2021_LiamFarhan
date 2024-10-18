#include <math.h>
#include<tgmath.h>
#include<stdlib.h>
#include <stdio.h>


int find_index(int i, int length){

/* Hilfsfunktion, welche ein Array zu einem Kreis macht 
   i: Index, der verarbeitet werden soll; length: Länge des Arrays, der zu einem Kreis gemacht werden soll
   Funktion gibt i verarbeitet zurück 
*/
	if (i>=0 && i<length){ // index is within boundaries -> do nothing
		return i;
	}
	if (i<0){ // index is out of boundaries from the left -> loop from the right
		while (i<0){
			i=length+i;
		}
		return i;
	}
	if (i>=length){// index is out of boundaries from the right -> loop from the left
		while (i>=length){
			i=i-length;
		}
		return i;
	}
	return i;
}


double generate_r( long int* last_number){

/* Funktion erzeugt eine pseudo-zufällige Zahl zwischen 0 und 1 und gibt sie zurück
   last_number : globales Seed der zufällig erstellten Zahlen, wird von der Funktion itteriert 
*/
	double r;
	int a =1664525;
	long int c = 1013904223;
	long int m = 4294967296; //2^32
	*(last_number) = ( *(last_number) * a +c) % m;
	r = (double) *last_number/ (m-1);
	return r;
}


int update_state(double p,double q, long int * last_number, int* L, int i, int length, int d){

/* Funktion macht eine zufällig ausgewählte Änderung zum ausgewählten Walker und gibt die totale Änderung der Anzahl der Walker zurück. 
   i: Ausgewählter Walker, L: Array wo die Walker sich befinden, length: Länge von L, d: wie weit soll der walker springen  
*/
	int change=0; // change in total number of walkers
	double r = generate_r(last_number);
	if (L[find_index( i, length)]==0){ // if chosen walker doesn't exisit
		printf("invalid index: There is no walker in %d\n",find_index( i, length));
		return -10;
	}
	if (r<0||r>1||1-p-q<0){//case zero: invalid input
		printf("Erorr: random number or probability overflow\n");
		change = -10;
	}
	
	if (r>=0 && r<p){ //case one: walker disappears
		i=find_index( i, length);
		L[i]=0;
		change = -1;
	}
	
	if (p<r && r< (1-q)){ //case two: walker is doplicated
		i=find_index( i, length);
		change = 2 - L[find_index( i-d, length)] - L[find_index( i+d, length)] -
		 L[i];
		L[i]=0;
		L[find_index( i-d, length)]=1;
		L[find_index( i+d, length)]=1;
	}
	
	if ( (1-q)<r && r< 1-q/2 ){ //case three: walker jumps to the right 
		i=find_index( i, length);
		L[i]=0;
		change = change - L[find_index( i+d, length)];
		L[find_index( i+d, length)]=1;
	}
	
	if ( 1-q/2<r && r< 1 ){ //case four: walker jumps to the left 
		i=find_index( i, length);
		L[i]=0;
		change = change - L[find_index( i-d, length)];
		L[find_index( i-d, length)]=1;
	}
	
	return change; 
}


int pick_a_walker(int number_of_walkers, long int* last_number, int* L, int length){

/* FunKtion wählt zufällig ein Walker aus und gibt seinen Index zurück. */

	int a =1664525;
	long int c = 1013904223;
	long int m = 4294967296; //2^32
	*last_number = ( (*last_number) * a +c) % m;
	int i =  *last_number % (number_of_walkers);
	i = i+1;
	int dummy=0;
	for (int k=0; k<length; k++){
		dummy=dummy+L[k];
		if (dummy==i){
			return k;
		}
	}
	return -1; //Falls irgendwas schief ausgelafen ist
}

void print_state(int * L, int length, int number_of_walkers){

/* Funktion existiert zur Debug-Zwecke. Sie stellten den aktuellen Zustand der Walker-Array dar. */

	for (int k=0;k<length-1;k++){
		printf("%d, ", L[k]);
	}
	printf("%d. Number of walkers: %d\n",L[length-1],number_of_walkers);
	
	return;
}


double density(double p,double q, long int * last_number, int length, int d, int t){

/* Funktion simuliert die Walker auf einem Array der Länge length für t Schritte und berechnet die Dichte der Walker
   nach ablauf der Simulation und gibt sie zurück.
   p,q,d: Simulationswariabeln (siehe dafür Aufgabenstellung), last_number: notwendig für Zufallszahlerstellung
*/
	// erzeuge L
	int number_of_walkers =1;
	int * L = malloc(sizeof(int)*(length));
	if (L==NULL){
		printf("Error: memory overflow\n");
		return -1;
	}
	L[0]=1;
	for (int k=1;k<length;k++){
		L[k]=0;
	}
	
	/* Durch mehrfaches Aufruf von pick_a_walker() wird last_number ständig geändert so, dass es zufällig aussieht
	for (int k=0;k<1000000;k++){ 
		pick_a_walker( number_of_walkers, last_number, L, length);
	}*/
	// durchlaufe die Simulation der Walker für t-Schritte
	int i=0;
	double tmp_change = 0; // hier wird die Änderung der Anzahl der Walker nach einem Schritt gespeichert
	for (int k=0;k<t;k++){
		tmp_change = update_state( p, q, last_number, L,  i, length, d);
		if (tmp_change ==-10){ // irgendwas ist in update_state schiefgelaufen
			return -1;
		}
		number_of_walkers= number_of_walkers + tmp_change;
		//print_state(L, length, number_of_walkers ); // Zur Debuging
		if (number_of_walkers==0){
			break;
		}
		i = pick_a_walker( number_of_walkers, last_number,  L, length);
	} 
	free(L);
	// berechne die Dichte
	double rho= (double) number_of_walkers / length;
	return rho;
}



double average_density(double p, double q, long int * last_number, int length, int d, int t, int N, double * average , double * stdev){

/* Funktion berechnet den Mittelwert oder die Standardabweichung der Dichte durch N-faches Aufrauf der funktion density()
  is_stdev: soll die Funktion die Standardabweichung berechnen? (ja=1, sonst wird der Mittelwert berechnet)
  N: über wie viele Aufrufe von density() soll die Dichte gemittelt werden, sonstige Variablen: siehe density()
*/	
	double average_sum = 0;
	double stdev_tmp = 0;
	double rho_tmp = 0;
	double rho_quadrat_sum = 0; // notwendig für stdev
	int valid_number_of_simulations = N;
	// Mittelwert und Standardabweichung der Dichte über N Simualtionen berechnen 
	for (int k=0;k<N;k++){
		rho_tmp = density(p , q, last_number,length, d, t); // speichere die Dcihte einer Simulation
		if (rho_tmp == -1){ // Überprüfe, ob diese falsch gerechnet wurde und nehme sie aus Datenset raus.
			rho_tmp = 0;
			valid_number_of_simulations = valid_number_of_simulations - 1;
		}
		average_sum = average_sum + rho_tmp; // summe über N smiulierten Dichten
		rho_quadrat_sum = rho_quadrat_sum + (double) pow(rho_tmp, 2); // summe über Quadrate von N simiulierten Dichten
	}
	average_sum = (double) average_sum/valid_number_of_simulations; //
	rho_quadrat_sum = (double) rho_quadrat_sum/valid_number_of_simulations;
	stdev_tmp = (rho_quadrat_sum - (double) pow(average_sum,2)) / valid_number_of_simulations;
	stdev_tmp = (double) pow(stdev_tmp, 0.5);
	// Ausgaben werden in den Variabeln gespeichert
	*average = average_sum;
	*stdev = stdev_tmp;
	return 1;
}

