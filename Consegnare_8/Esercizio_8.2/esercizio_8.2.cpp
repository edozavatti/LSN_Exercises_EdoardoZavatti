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
#include <iomanip>
#include "random.h"
#include "FunzioneBase.h"
#include "Posizione.h"

using namespace std;

double error(double av, double av_2, int n) {
	
	if (n == 1) {
		return 0;
	}
	else {
		return sqrt((av_2 - pow(av, 2))/double(n-1));
	}
}
 
int main (int argc, char *argv[]){

   Random my_rand;

   int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
  		Primes >> p1 >> p2 ;
	} else std::cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
  		while ( !input.eof() ){
     		input >> property;
     		if( property == "RANDOMSEED" ){
        		input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        		my_rand.SetRandom(seed,p1,p2);
     		}
  		}
  		input.close();
	} else std::cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	my_rand.SaveSeed();

   ofstream E_ave, point;
   const int wd=12;

   double x, mu, sigma;
   double hbar, m;
   int Nstep = 1000000;
   int Nblocks = 100;
   int Nstep_block = int(Nstep/Nblocks);
   double passo = 3.;

   x=0.;
   mu=0.78;
   sigma=0.59;

   hbar = 1.;
   m = 1.;

   Psi_sigma_mu* f_1 = new Psi_sigma_mu();
   Potenziale_doppia_buca* V = new Potenziale_doppia_buca();
   H_Psi* H = new H_Psi(V, m, hbar);

   Posizione p(x);

   E_ave.open("Esercizio_8.2/output_Energia_precis.dat",ios::app);
   //point.open("Esercizio_8.2/output_sigma",ios::app);

   int accepted = 0;
   int attempted = 0;
   double E_best_new;

   while(mu<=0.82){
      cout << endl;
      cout << "mu = " << mu << endl;
      sigma = 0.59;
      while(sigma<=0.63){
         cout << "sigma = " << sigma << endl;
         double ave = 0.;
         double ave2 = 0.;
         for (int i = 1; i <= Nblocks; i++){
            double sum = 0.;
            //double sum2 = 0.;
            for(int step = 1; step <= Nstep_block; step++){
               //Metropolis algorithm
               double x_new = x + my_rand.Rannyu(-passo/2. , passo/2.);
               double z = f_1->Eval(x_new, mu, sigma) / (f_1->Eval(x, mu, sigma));
               double alpha = fmin(1., z);
               double ran = my_rand.Rannyu();
         
               if(ran <= alpha) { 
                  x = x_new;
                  accepted ++;
               }
               attempted ++;
               point << x << endl;         
               double E = H->Eval(x, mu, sigma); //energy mean value
               sum += E;          
            }
         ave += sum/double(Nstep_block);
         ave2 += pow(sum/double(Nstep_block), 2);

         }
         E_best_new = ave/double(Nblocks);
         double err = error(ave/double(Nblocks), ave2/double(Nblocks), Nblocks);
         E_ave << setw(wd) << mu << setw(wd) << sigma <<  setw(wd) << E_best_new << setw(wd) << err << endl;
      
         cout << "Energia = " << E_best_new << endl;
         cout << endl;
         cout << "----------------------------------------------"  << endl;
   
         sigma += 0.001;
      }

      mu += 0.001;
   }

   
   E_ave.close();
   point.close();

   cout << "Acceptance = " << double(accepted)/double(attempted) << endl;

   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
