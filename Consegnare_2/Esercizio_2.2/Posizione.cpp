#include "Posizione.h"
#include <cmath>
#include "random.h"
#include <iostream>

using namespace std;

	Posizione::Posizione(double x, double y, double z) {
	
	
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
   		} else std::cerr << "PROBLEM: Unable to open seed.in" << endl;
	
		m_x = x;
		m_y = y;
		m_z = z;
	}

	Posizione::~Posizione() {
	}

// Coordinate cartesiane
void Posizione::setX(double x) {
	m_x = x;
}
void Posizione::setY(double x) {
	m_y = x;
}
void Posizione::setZ(double x) {
	m_z = x;
}
double Posizione::getX() const{
	return m_x;
}
  double Posizione::getY() const{
	return m_y;
}
  double Posizione::getZ() const{
	return m_z;
}
  // Coordinate sferiche
double Posizione::getR() const {
	return sqrt(pow(m_x, 2)+pow(m_y, 2)+pow(m_z, 2));
}
  double Posizione::getPhi() const {
	return atan2               (m_y,m_x);
}
  double Posizione::getTheta() const {
	return acos(m_z/getR());
}
  double Posizione::getRho() const {     // raggio delle coordinate cilindriche
 	return sqrt(pow(m_x, 2)+pow(m_y, 2));
}
double Posizione::Distanza(const Posizione& b) const { // distanza da un altro punto
	return sqrt(pow(getX()-b.getX(), 2)+
		pow(getY()-b.getY(), 2)+
		pow(getZ()-b.getZ(), 2));
}

vector<double> Posizione::Random_walk_discreto(int n_passi, double passo){ //fa random walk discreto di N-passi e ritorna un vettore con il quadrato della distanza dall'origine a ogni passo

	vector<double> a;

	for(int l=0; l<n_passi; l++){ //faccio random walk di 100 passi
			
		double d = m_my_rand.Rannyu();

		ofstream outputFile("numeri.dat", std::ios_base::app);
	if (!outputFile) {
		cerr << "cannot create file" << endl;
	} else {
		outputFile << d << endl;
	}		
	outputFile.close();
   	

		if(d<1./6.){m_x = (m_x + passo);
		}	
		else if(d>1./6. and d<2./6.){m_x = m_x - passo;
		}
		else if(d>2./6. and d<3./6.){m_y = m_y + passo;
		}
		else if(d>3./6. and d<4./6.){m_y = m_y - passo;	
		}	
		else if(d>4./6. and d<5./6.){m_z = m_z + passo;
		}
		else if(d>5./6.){m_z = m_z - passo;	
		}
		
		double r_2 = getR() * getR(); //calcolo il quadrato della distanza dall'origine
		a.push_back(r_2); //salvo il valore nel vector
		
	}
	
	return a;

}