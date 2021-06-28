#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "funzioni.h"
#include "Posizione.h"

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


	int N = 10000;      //numero random walks
	int M = 100;         //numero blocchi	 
 	int L = int(N/M);    //numero random walks per blocco
	int passi = 101;
	double passo = 1.;
	
	vector<double> sum_r(passi);
	vector<double> ave_r(passi); //vettore di medie di r^2
	vector<double> ave_r_2(passi); 
	vector<double> err(passi);
/*	
	for(int k=0; k<passi; k++){
			
		ave_r[k] = 0; //inizializzo a zero il vettore di somme di somme di moduli quadri
				
	}
*/	
	for(int i=0; i<M; i++){ //ciclo sui blocchi
	
		for(int k=0; k<passi; k++){
			
			sum_r[k] = 0; //inizializzo a zero il vettore di somme di moduli quadri
				
		}
			
		for(int j=0; j<L; j++){ //entro nel blocco j
		
			vector<double> v;
			Posizione P; //genero punto nella posizione (0,0,0)
			
			for(int l=0; l<passi; l++){ //faccio random walk di 100 passi
			
				double r_2 = P.getR() * P.getR();
				v.push_back(r_2); //in questo modo il passo 0 è sempre 0

				double d = rnd.Rannyu();
			
				if(d<1./6.){P.setX(P.getX() + passo);
				}	
				else if(d>1./6. and d<2./6.){P.setX(P.getX() - passo);
				}
				else if(d>2./6. and d<3./6.){P.setY(P.getY() + passo);
				}
				else if(d>3./6. and d<4./6.){P.setY(P.getY() - passo);	
				}	
				else if(d>4./6. and d<5./6.){P.setZ(P.getZ() + passo);
				}
				else if(d>5./6.){P.setZ(P.getZ() - passo);	
				}
		
			}

			for(int k=0; k<passi; k++){
				
				sum_r[k] += v[k]; //salvo le somme dei moduli quadri di r per ogni passo k
				
			}
		
		}
			
		for(int k=0; k<passi; k++){
			
			double A_i = sum_r[k]/(L);
			ave_r[k] += A_i;
			ave_r_2[k] += pow(A_i, 2) ;
			err[k] = error(ave_r[k]/(M), ave_r_2[k]/(M), M);
		
		}
	
	}
	
	for (int i=0; i<passi; i++){
		
		Print("media_random_walk.dat", sqrt(ave_r[i]/M));
		Print("errore_random_walk.dat", err[i]/(2*sqrt(ave_r[i]/M)));
		
	}
		
	return 0;
	
}

