#ifndef __Posizione_h__
#define __Posizione_h__

#include "random.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class Posizione {

public:

  // costruttori
  Posizione() {
	
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
            		m_my_rand.SetRandom(seed,p1,p2);
         		}
      		}
      		input.close();
   		} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
		m_x = 0.;
		m_y = 0.;
		m_z = 0.;
	
  };
  
  Posizione(double x, double y, double z); 
  // distruttore
  ~Posizione();
  // metodi
  double getX() const;       // Coordinate cartesiane
  double getY() const;
  double getZ() const;
  void setX(double) ;       // Coordinate cartesiane
  void setY(double) ;
  void setZ(double) ;
  double getR() const;       // Coordinate sferiche
  double getPhi() const;
  double getTheta() const;
  double getRho() const;     // raggio delle coordinate cilindriche
  double Distanza(const Posizione&) const; // distanza da un altro punto
  vector<double> Random_walk_discreto(int, double);

private:

  double m_x, m_y, m_z;   
  Random m_my_rand;

};

#endif // __posizione_h__


