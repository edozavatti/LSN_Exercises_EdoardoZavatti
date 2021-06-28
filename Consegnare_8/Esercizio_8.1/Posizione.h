#ifndef __Posizione_h__
#define __Posizione_h__

#include "random.h"
#include "FunzioneBase.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class Posizione {

public:

  // costruttori
  Posizione();
  Posizione(double x, double y, double z);
  Posizione(double x); 
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
  double getAlpha() const;
  double getRho() const;     // raggio delle coordinate cilindriche
  double Distanza(const Posizione&) const; // distanza da un altro punto
  vector<double> Random_walk_discreto(int, double);
  void Metropolis(double, const FunzioneVettoriale&);
  void Metropolis(double, const FunzioneBase&);

private:

  double m_x, m_y, m_z, m_alpha;   
  Random m_my_rand;

};

#endif // __posizione_h__


