#ifndef __integralMC_h__
#define __integralMC_h__

#include "FunzioneBase.h"
#include "random.h"
#include <iostream>
#include <fstream>

class integralMC {

	public:
	
	integralMC () {
	
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
            		m_my_rand.SetRandom(seed,p1,p2);
         		}
      		}
    		input.close();
   		} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	};
	
	~integralMC() {};
	
	double integralAVE(double xmin, double xmax, FunzioneBase* f, int punti);
	double integralAVE_par(double xmin, double xmax, FunzioneBase* f, int punti);

	private:
	
	Random m_my_rand;
	
};

#endif