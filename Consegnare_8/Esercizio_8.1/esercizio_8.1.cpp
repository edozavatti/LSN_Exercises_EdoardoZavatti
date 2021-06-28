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
   double passo = 6.0;

   x=1.;
   mu=2.;
   sigma=1.;

   hbar = 1.;
   m = 1.;
   
   double ave = 0.;
   double ave2 = 0.;

   Psi_sigma_mu* f_1 = new Psi_sigma_mu();
   Potenziale_doppia_buca* V = new Potenziale_doppia_buca();
   H_Psi* H = new H_Psi(V, m, hbar);

   Posizione p(x);

   E_ave.open("output_E_aveprova.dat",ios::app);
   point.open("output_x_new.dat",ios::app);

   int accepted = 0;
   int attempted = 0;

  for (int i = 1; i <= Nblocks; i++){
      double sum = 0.;
      for(int step = 1; step <= Nstep_block; step++){
         //Metropolis algorithm
         double x_new = x + my_rand.Rannyu(-passo/2. , passo/2.);
         double alpha = fmin(1., f_1->Eval(x_new, mu, sigma) / (f_1->Eval(x, mu, sigma)));
         if(my_rand.Rannyu() <= alpha) { 
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
   double err = error(ave/double(i), ave2/double(i), i);
   E_ave << setw(wd) << i <<  setw(wd) << sum/double(Nstep_block) << setw(wd) << ave/double(i) << setw(wd) << err << endl;
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
