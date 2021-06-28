#include "random.h"
#include "funzioni.h"
#include <cmath>
#include <fstream>
#include <vector>

double error(double av, double av_2, int n) {
	
		if (n == 0) {
			return 0;
		}
		else {
			return sqrt((av_2 - pow(av, 2))/n);
		}
	};

/*double N(double x) {

    return 0.5 * (1. + erf(x / sqrt(2.)));

};
*/
double Call(double r, double T, double S, double K){

    double a = S - K;
    //cout << exp( - r*T ) * max(0., a) << endl;
    return exp( - r*T ) * max(0., a);
    
};

double Put(double r, double T, double S, double K){

    double a = K - S;
    //cout << exp( - r*T ) * max(0., a) << endl;
    return exp( - r*T ) * max(0., a);

};

using namespace std;

int main() {

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

    double S0 = 100.; 
    double K = 100.;
    double T = 1.;
    double r = 0.1;
    double sigma = 0.25;
    
    int N = 100000; //numero lanci montecarlo
    int M = 1000; //numero blocchi
    int L = int(N/M); //numero lanci per blocco
    int passi = 100;

    double sum_prog_Call = 0., sum_prog_Call_2 = 0., sum_prog_Put = 0., sum_prog_Put_2 = 0., err_Call, err_Put;
    vector<double> media_Call, media_Put, errore_Call, errore_Put;
   

    //calcolo senza simulare il percorso

    for(int i = 0; i < M; i++){ //entro nel blocco

        double sumC = 0.;
        double sumP = 0.;

        for(int j = 0; j < L; j++){

            double t = T/passi;
            double S = S0;

            for (int k = 0; k < passi; k++){
                
                double Z = rnd.Gauss(0., 1.);
                S = S * exp((r - pow(sigma, 2)/2) * t + sigma * Z * sqrt(t));
           
            }

            sumC += exp(-r*T) * max(0. , S - K);
            sumP += exp(-r*T) * max(0. , K - S);

        }
        
        double A_i = sumC/L;
       // cout << A_i << endl;
        double B_i = sumP/L;
 		sum_prog_Call += A_i;
        sum_prog_Put += B_i;
 		sum_prog_Call_2 += pow(A_i, 2);
        sum_prog_Put_2 += pow(B_i, 2); 
 		
 		media_Call.push_back(sum_prog_Call/(i+1) - 14.975790778311286);
        media_Put.push_back(sum_prog_Put/(i+1) - 5.4595325819072364); 
 		err_Call = error(sum_prog_Call/(i+1), sum_prog_Call_2/(i+1), (i+1));
        err_Put = error(sum_prog_Put/(i+1), sum_prog_Put_2/(i+1), (i+1));
 		errore_Call.push_back(err_Call);
        errore_Put.push_back(err_Put);

    }

    Print("Call_bis.dat", media_Call);
 	Print("errore_Call_bis.dat", errore_Call);
    Print("Put_bis.dat", media_Put);
 	Print("errore_Put_bis.dat", errore_Put);

   rnd.SaveSeed();

    return 0;

}