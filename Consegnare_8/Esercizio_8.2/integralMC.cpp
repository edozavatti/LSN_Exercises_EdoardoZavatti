#include "integralMC.h"
#include "random.h"
#include <iostream>

double integralMC :: integralAVE(double xmin, double xmax, FunzioneBase* f, int punti) {

	double x; 
	double sum = 0;	
   			
	for (int i=0; i<punti; i++) {
	
		x = m_my_rand.Rannyu(xmin, xmax);
		sum += f->Eval(x);
	
	}
	
	return (xmax-xmin)*sum/punti;

};

double integralMC :: integralAVE_par(double xmin, double xmax, FunzioneBase* f, int punti) {

	   	double x, y; 
		double sum = 0;	
		
		for (int i=0; i<punti; i++) {
	
			y = m_my_rand.Rannyu();
			x = 1 - sqrt(1-y);
			sum += f->Eval(x);
	
		}
	
		return (xmax-xmin)*sum/punti;

};
