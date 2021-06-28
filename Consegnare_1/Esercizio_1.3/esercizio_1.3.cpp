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

	double dist(double x1, double x2, double y1, double y2){
	
		return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
	
	};
	
	double error(double av, double av_2, int n) {
	
		if (n == 0) {
			return 0;
		}
		else {
			return sqrt((av_2 - pow(av, 2))/n);
		}
	}; 
	
	int int_part(double A){
	
		int B = (int)A;
		return A;
	
	}
	
	int sign(double A){
	
		if(A > 0){
		
			return 1;
	
		}else
		
			return -1;
	
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

	double l = 2.95; //lunghezza ago
	double d = 3.;  //distanza righe
//	double A = 3.;  //larghezza campo
	
	double N_throw, N_hit; 
	double PIE_AV, PIE_AV_2;
/*	
	int N = 100000;      //numero lanci Montecarlo
	int M = 100;         //numero blocchi	 
 	int L = int(N/M);    //numero lanci per blocco
 
	vector<double> media_medie; 
	vector<double> errore;
 
 	double sum_prog = 0.;
 	double sum_prog_2 = 0.;
 		
 	for (int i=0; i<M; i++) { //entro nel blocco i
		cout << i << endl;
		double P_i;
		double A_i = 0.;
		
		for (int k=0; k<L; k++){     //lancio nel blocco
		
			N_throw = 0.;
			N_hit = 0.;
		
			for (int j = 0; j<10000; j++){
		
				N_throw ++;
				
				double centro = rnd.Rannyu(0, d/2);            //genero centro e angolo del lancio
				double angolo = rnd.angolo();
	
				if (centro < l/2 * sin(angolo)) {        //condizione per colpire la riga
					
					N_hit ++;
				
				}
		
			}
			
			P_i = (N_hit/N_throw); //calcolo la probabilità nel blocco				
			A_i += 2. * l/(P_i * d); //calcolo il valore di pi nel blocco
		
		}
		
 		sum_prog += A_i/L; 		
 		sum_prog_2 += pow(A_i/L, 2);
 		
 		media_medie.push_back(sum_prog/(i+1) - M_PI);
 		double err = error(sum_prog/(i+1), sum_prog_2/(i+1), (i+1));
 		errore.push_back(err);
 				
 	}
 
 	Print("Media_pi.dat", media_medie);
 	Print("Errore_pi.dat", errore);
 	
*/

	int N = 100000000;      //numero lanci Montecarlo
	int M = 100;         //numero blocchi	 
 	int L = int(N/M);    //numero lanci per blocco
 	
 	vector<double> P; 
 	vector<double> Err;

	ofstream phi;
	
	/*for (int k=0; k<100 ; k++){
	
		cout << k<< endl;
		
		double sum = 0.;
		
		N = (k+1)*10;   //numero blocchi
		M = L*N;    //numero lanci totali*/
		
		double PIE ;
		double PIE_2 ;
		double errore_pi ;
		phi.open("angolo.dat", ios::app);
		for (int i=0; i<M; i++){   //entro nel blocco k_i
			cout << "Blocco " << i << endl;
			N_throw = 0;
			N_hit = 0;

			for (int j=0; j<L; j++){     //lancio nel blocco

				N_throw ++;
				
				double centro = rnd.Rannyu(0, d/2);            //genero centro e angolo del lancio
				double angolo = rnd.angolo();                  // genero angolo fra 0 e pi
				//phi << angolo << endl;

				if (centro < l/2 * sin(angolo)) {        //condizione per colpire la riga
				
					N_hit ++;
				
				}
				
			} //ho finito i lanci nel blocco
			
			phi.close();
			double P_i = (N_hit/N_throw); //calcolo la probabilità nel blocco
			//cout << P_i << endl;
			PIE += 2 * l/(P_i * d);       //calcolo il valore di pi nel blocco
			PIE_2 += pow(2 * l/(P_i * d), 2);
			
			PIE_AV = PIE/(i+1);                   //calcolo il valore medio di pi alla fine di tutti i blocchi
			PIE_AV_2 = PIE_2/(i+1);
			errore_pi = error(PIE_AV, PIE_AV_2, (i+1));  //errore statistico sulla media
		
			P.push_back(PIE_AV - M_PI);              //salvo la media di pi su N blocchi
			Err.push_back(errore_pi);         //salvo l'errore di pi su N blocchi			
		}
	
	//}
	
	Print("Media_pi.dat", P);
	Print("Errore_pi.dat", Err);
	
   rnd.SaveSeed();
	
	return 0;
	
}

