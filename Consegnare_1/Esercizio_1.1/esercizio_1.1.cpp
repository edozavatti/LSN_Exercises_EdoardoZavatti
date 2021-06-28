#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "funzioni.h"

using namespace std;
 
	double error(double av, double av_2, int n) {
	
		if (n == 0) {
			return 0;
		}
		else {
			return sqrt((av_2 - pow(av, 2))/n);
		}
	}; 
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

	int N = 1000000;      //numero lanci Montecarlo
	int M = 100;         //numero blocchi	 
 	int L = int(N/M);    //numero lanci per blocco
 
	vector<double> media_medie; 
	vector<double> errore;
 
 	double sum_prog = 0.;
 	double sum_prog_2 = 0.;
 		
 	for (int i=0; i<M; i++) {
 	
 		double sum = 0.;
 		
 		for (int j=0; j<L; j++) {
 		
			double r = rnd.Rannyu();
			sum += r;
			
 		}

		double A_i = sum/L;
 		sum_prog += A_i;
 		sum_prog_2 += pow(A_i, 2);
 		
 		media_medie.push_back(sum_prog/(i+1) - 0.5);
 		double err = error(sum_prog/(i+1), sum_prog_2/(i+1), (i+1));
 		errore.push_back(err);
 				
 	}
 
 	Print("media.dat", media_medie);
 	Print("errore_media.dat", errore);
 	
 	vector<double> varianza_varianze; 
	vector<double> errore_varianza;
 
 	double sum_prog_var = 0.;
 	double sum_prog_2_var = 0.;
 		
 	for (int i=0; i<M; i++) {
 	
 		double sum = 0.;
 		
 		for (int j=0; j<L; j++) {
 		
			double r = rnd.Rannyu();
			sum += pow((r - 0.5), 2);
			
 		}

		double A_i = sum/L;
 		sum_prog_var += A_i;
 		sum_prog_2_var += pow(A_i, 2);
 		
 		varianza_varianze.push_back(sum_prog_var/(i+1) - 0.08333);
 		double err = error(sum_prog_var/(i+1), sum_prog_2_var/(i+1), (i+1));
 		errore_varianza.push_back(err);
 				
 	}
 
 	Print("varianza.dat", varianza_varianze);
 	Print("errore_varianza.dat", errore_varianza);

 	int intervalli = 100;  //numero intervalli in cui divido [0, 1]
 	int points = 10000;    //numero punti che genero 

 	vector<double> chi_2; 
 	
 
 	for (int j=0; j<intervalli; j++){
 	
		vector<int> n(intervalli, 0);   //vettore di contatori per contare 
		                                //quanti punti sono generati in ogni intervallo
		double sum = 0.;
 				
 		for (int i=0; i<points; i++){
 	
 			double n_i = rnd.Rannyu();
 			int a = n_i * intervalli;
 			n[a]++; 	                //riempio il contatore
 		}		
 	
 		for (int f=0; f<intervalli; f++){
 		
 			sum += pow(n[f]-(points/intervalli), 2)/(points/intervalli);
 				
 		}
 				
 		chi_2.push_back(sum);
 			
 	}

	Print("chi_2.dat", chi_2);
 
   rnd.SaveSeed();
   
   return 0;
}