/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "funzioni.h"
#include "random.h"

using namespace std;

	double dist(double x1, double x2){
	
		return sqrt(pow(x1, 2) + pow(x2, 2));
	
	}
	
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
   
   
   int N[4] = {1, 2, 10, 100}; //vettore di numeri di blocchi
   int M = 10000; //numero lanci Montecarlo
   
   vector<string> uniforme_filename = {"uniforme_1.dat", "uniforme_2.dat", "uniforme_10.dat", "uniforme_100.dat" };  //vettori per immagazzinare i
   vector<string> exp_filename = {"exp_1.dat", "exp_2.dat", "exp_10.dat", "exp_100.dat" };                           //nomi dei file di dati
   vector<string> lorenz_filename = {"lorenz_1.dat", "lorenz_2.dat", "lorenz_10.dat", "lorenz_100.dat" }; 
  	
    for(int k=0; k<M; k++){
        
        cout << k << endl;
   		for(int j=0; j<4; j++){
	
	 		double sum_unif = 0.;
	 		double sum_exp = 0.;
   	 		double sum_lorenz = 0.;
   	 
  			for (int i=0; i<N[j]; i++) {
   
  		 		double Lorenz = rnd.Lorenz(1., 0.);
   				double Exp = rnd.Exp(1.);
   				double unif = rnd.Rannyu();
   		
   				sum_lorenz += Lorenz;
   				sum_exp += Exp;		
   				sum_unif += unif;
   
   	 		}	
   	 		
   	 		double A_unif = sum_unif/N[j];
   	 		double A_exp = sum_exp/N[j];
   	 		double A_lorenz = sum_lorenz/N[j];
   	 
   	 		Print(uniforme_filename[j].c_str(), A_unif);
   	 		Print(exp_filename[j].c_str(), A_exp);
   	 		Print(lorenz_filename[j].c_str(), A_lorenz);
   	 		
   		}
   
   }
   
   
   
    rnd.SaveSeed();
   
   return 0;
}